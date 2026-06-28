"""Annotation plotting and GFF/GenBank generation."""

from __future__ import annotations

import logging
import os
import warnings
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.graph_objects as go
from Bio import BiopythonWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from BCBio import GFF
from dna_features_viewer import BiopythonTranslator

from pca.parsers import parse_hmmsearch, parse_trna_gff

__all__ = [
    "create_feature",
    "get_coordinates",
    "get_cordinates",
    "generate_plots_and_gff",
    "generate_plots_and_annotations",
]

logger = logging.getLogger(__name__)

DEFAULT_PLOT_FORMATS = ("pdf",)
SUPPORTED_PLOT_FORMATS = {"pdf", "png", "html"}


def get_coordinates(x: pd.Series) -> pd.Series:
    """Normalize tRNA coordinates and strand."""
    if x["begin"] > x["end"]:
        return pd.Series(
            [x["qname"], x["trna_no"], x["end"], x["begin"], -1, x["trna_type"], x["score"]],
            index=["contig", "trna_no", "begin", "end", "strand", "trna_type", "score"],
        )
    return pd.Series(
        [x["qname"], x["trna_no"], x["begin"], x["end"], 1, x["trna_type"], x["score"]],
        index=["contig", "trna_no", "begin", "end", "strand", "trna_type", "score"],
    )


# Backwards-compatible alias for the original typo.
get_cordinates = get_coordinates


def create_feature(x: pd.Series) -> SeqFeature:
    """Create a Biopython SeqFeature from a normalized tRNA row."""
    qualifiers = {
        "source": "tRNAscan-SE",
        "score": x["score"],
        "trna_type": x["trna_type"],
        "label": f"{x['trna_type']}_tRNA",
        "ID": f"trna_{x['trna_no']}",
    }
    return SeqFeature(
        FeatureLocation(int(x["begin"]), int(x["end"]), strand=int(x["strand"])),
        type="tRNA",
        id=str(x["trna_no"]),
        qualifiers=qualifiers,
    )


def _iter_gff_records(gff_path: str | os.PathLike[str]) -> Iterator[Any]:
    """Yield Biopython SeqRecord objects from a GFF file."""
    with open(gff_path) as handle:
        yield from GFF.parse(handle)


def _read_sequences(fasta_path: str | os.PathLike[str]) -> dict[str, Seq]:
    """Return a mapping of sequence ID to Biopython Seq from a FASTA file."""
    sequences: dict[str, Seq] = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences[record.id] = record.seq
    return sequences


def _attach_sequences(
    records: Iterator[Any],
    sequences: dict[str, Seq],
) -> Iterator[Any]:
    """Attach DNA sequences to GFF records and set GenBank annotations."""
    for record in records:
        if record.id in sequences:
            record.seq = sequences[record.id]
        else:
            logger.warning("no sequence found for contig %s", record.id)
            record.seq = Seq("")
        record.annotations.setdefault("molecule_type", "DNA")
        yield record


def _plotly_hover_text(feature: Any, record_id: str = "") -> str:
    """Build hover text for a feature in the interactive plot."""
    qualifiers = getattr(feature, "qualifiers", {})
    label = qualifiers.get("label") or record_id or "feature"
    if isinstance(label, list):
        label = label[0]
    lines = [f"<b>{label}</b>"]
    for key in ("phrog", "category", "eval", "score", "color", "trna_type"):
        value = qualifiers.get(key)
        if value is not None:
            if isinstance(value, list):
                value = value[0]
            lines.append(f"{key}: {value}")
    return "<br>".join(lines)


def _write_interactive_plot(
    record: Any,
    plots_dir: Path,
    contig_length: int,
) -> None:
    """Write an interactive Plotly HTML genome map for a contig."""
    traces: list[go.Scatter] = []

    for feature in record.features:
        start = int(feature.location.start)
        end = int(feature.location.end)
        strand = feature.location.strand
        color = feature.qualifiers.get("color", "#c9c9c9")
        if isinstance(color, list):
            color = color[0]

        y_base = 1 if strand == 1 else -1
        height = 0.3
        y = [y_base, y_base, y_base + height, y_base + height, y_base]
        x = [start, end, end, start, start]

        hover_text = _plotly_hover_text(feature, record.id)
        label = feature.qualifiers.get("label", "")
        if isinstance(label, list):
            label = label[0] if label else ""
        traces.append(
            go.Scatter(
                x=x,
                y=y,
                fill="toself",
                fillcolor=color,
                line={"color": color, "width": 1},
                mode="lines",
                name=str(label),
                text=hover_text,
                hoverinfo="text",
                showlegend=False,
            )
        )

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=record.id,
        xaxis_title="Position",
        yaxis={
            "tickvals": [-1, 1],
            "ticktext": ["Reverse strand", "Forward strand"],
            "range": [-1.8, 1.8],
            "zeroline": False,
        },
        xaxis={"range": [0, contig_length], "zeroline": False},
        hovermode="closest",
        plot_bgcolor="white",
        height=400,
        margin={"l": 80, "r": 40, "t": 60, "b": 60},
    )
    out_path = plots_dir / f"{record.id}.html"
    fig.write_html(str(out_path), include_plotlyjs="cdn")


def _write_static_plot(
    record: Any,
    plots_dir: Path,
    format: str,
) -> None:
    """Write a static matplotlib genome map for a contig."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    graphic_record = BiopythonTranslator().translate_record(record)
    for feat in graphic_record.features:
        if feat.label == "unknown function":
            feat.label = None

    fig, ax1 = plt.subplots(1, 1, figsize=(15, 4))
    fig.tight_layout(pad=2.5)
    ax, _ = graphic_record.plot(
        ax=ax1,
        strand_in_label_threshold=7,
        annotate_inline=False,
        figure_height=3,
    )
    ax.set_title(record.id)
    out_path = plots_dir / f"{record.id}.{format}"
    ax.figure.savefig(str(out_path), bbox_inches="tight", format=format)
    plt.clf()
    plt.close("all")


def _contig_length(record: Any) -> int:
    """Return the contig length from sequence or feature coordinates."""
    if record.seq:
        return len(record.seq)
    if record.features:
        return max(int(f.location.end) for f in record.features)
    return 0


def generate_plots_and_annotations(
    tmp_dir: str | os.PathLike[str],
    hmmsearch_dir: str | os.PathLike[str],
    trna_dir: str | os.PathLike[str],
    meta_dir: str | os.PathLike[str],
    gff_dir: str | os.PathLike[str],
    input_fasta: str | os.PathLike[str],
    plot_formats: Sequence[str] = DEFAULT_PLOT_FORMATS,
) -> None:
    """Generate annotated GFF/GenBank and per-contig plots from HMM results.

    Parameters
    ----------
    tmp_dir:
        Output directory for annotations, plots, and checkpoints.
    hmmsearch_dir:
        Path to the aggregated HMMER tblout file.
    trna_dir:
        Path to the tRNAscan-SE GFF file (may be empty).
    meta_dir:
        Path to the PHROG metadata TSV file.
    gff_dir:
        Path to the protein GFF file produced by gene calling.
    input_fasta:
        Path to the original nucleotide FASTA (used for GenBank and contig length).
    plot_formats:
        Iterable of plot formats to produce. Supported: ``pdf``, ``png``, ``html``.
        Default is ``("pdf",)``.
    """
    tmp_dir = Path(tmp_dir)
    plots_dir = tmp_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    requested_formats = set(plot_formats)
    unsupported = requested_formats - SUPPORTED_PLOT_FORMATS
    if unsupported:
        raise ValueError(f"Unsupported plot formats: {unsupported}")

    sequences = _read_sequences(input_fasta)

    logger.info("processing hmm results")
    search_results = pd.DataFrame(parse_hmmsearch(hmmsearch_dir))

    trna_path = Path(trna_dir)
    if trna_path.is_file() and trna_path.stat().st_size > 0:
        trna = pd.DataFrame(parse_trna_gff(trna_dir))
    else:
        trna = pd.DataFrame()

    if search_results.empty:
        raise ValueError("hmmsearch returned zero matches")

    if not trna.empty:
        trna = trna.apply(get_coordinates, axis=1).sort_values(by=["contig", "begin"])
    else:
        trna = None

    phrogs_anno = pd.read_table(meta_dir)
    phrogs_anno = phrogs_anno.fillna("unknown function")
    phrogs_anno["phrog"] = phrogs_anno["phrog"].apply(lambda x: f"phrog_{x}")

    results_with_annotate = search_results.merge(
        phrogs_anno, how="inner", left_on="tname", right_on="phrog"
    )
    results_with_annotate["position"] = results_with_annotate["qname"].apply(
        lambda x: int(x.split("_")[-1])
    )
    results_with_annotate["contig"] = results_with_annotate["qname"].apply(
        lambda x: x.rsplit("_", 1)[0]
    )

    results_filtered = results_with_annotate.loc[
        results_with_annotate.groupby("qname")["score"].idxmax()
    ].query("score > 50")

    logger.info("generating annotations and plots")
    gff_out_path = tmp_dir / "annotations.gff"
    gbk_out_path = tmp_dir / "annotations.gbk"

    annotated_records: list[Any] = []
    with open(gff_out_path, "w") as out_handle:
        for record in _attach_sequences(_iter_gff_records(gff_dir), sequences):
            tmp = results_filtered[results_filtered["contig"] == record.id]

            for pos, feature in enumerate(record.features, start=1):
                tmp_feature = tmp[tmp["position"] == pos][
                    ["category", "color", "annot", "phrog", "score", "eval"]
                ]
                if not tmp_feature.empty:
                    feature.qualifiers.update({"label": tmp_feature["annot"].values[0]})
                    feature.qualifiers.update({"category": tmp_feature["category"].values[0]})
                    feature.qualifiers.update({"eval": tmp_feature["eval"].values[0]})
                    feature.qualifiers.update({"score": tmp_feature["score"].values[0]})
                    feature.qualifiers.update({"color": tmp_feature["color"].values[0]})
                    feature.qualifiers.update({"phrog": tmp_feature["phrog"].values[0]})
                else:
                    feature.qualifiers.update({"label": "unknown function"})
                    feature.qualifiers.update({"color": "#c9c9c9"})

            if trna is not None:
                tmp_trna = trna[trna["contig"] == record.id]
                if not tmp_trna.empty:
                    tmp_trna = tmp_trna.reset_index(drop=True)
                    tmp_trna["feature"] = tmp_trna.apply(create_feature, axis=1)
                    record.features.extend(tmp_trna["feature"].to_list())

            annotated_records.append(record)
            GFF.write([record], out_handle)

            contig_length = _contig_length(record)

            if "html" in requested_formats:
                _write_interactive_plot(record, plots_dir, contig_length)

            static_formats = requested_formats & {"pdf", "png"}
            for fmt in static_formats:
                _write_static_plot(record, plots_dir, fmt)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        SeqIO.write(annotated_records, gbk_out_path, "genbank")
    logger.info("wrote %s and %s", gff_out_path, gbk_out_path)


# Backwards-compatible alias.
def generate_plots_and_gff(*args: Any, **kwargs: Any) -> None:
    """Deprecated alias for :func:`generate_plots_and_annotations`."""
    return generate_plots_and_annotations(*args, **kwargs)
