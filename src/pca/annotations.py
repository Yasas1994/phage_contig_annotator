"""Annotation plotting and GFF generation."""

from __future__ import annotations

import logging
import os
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pandas as pd
from Bio.SeqFeature import FeatureLocation, SeqFeature
from BCBio import GFF
from dna_features_viewer import BiopythonTranslator

__all__ = ["create_feature", "get_coordinates", "get_cordinates", "generate_plots_and_gff"]

logger = logging.getLogger(__name__)


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


def generate_plots_and_gff(
    tmp_dir: str | os.PathLike[str],
    hmmsearch_dir: str | os.PathLike[str],
    trna_dir: str | os.PathLike[str],
    meta_dir: str | os.PathLike[str],
    gff_dir: str | os.PathLike[str],
) -> None:
    """Generate annotated GFF and per-contig PNG maps from HMM search results."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from pca.parsers import parse_hmmsearch, parse_trna_gff

    tmp_dir = Path(tmp_dir)
    checkpoint_file = tmp_dir / "plotting_chkpt"
    if checkpoint_file.is_file():
        logger.info("plotting checkpoint found")
        return

    plots_dir = tmp_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

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

    logger.info("generating annotation plots")
    gff_out_dir = tmp_dir / "annotations.gff"
    with open(gff_out_dir, "w") as out_handle:
        for record in _iter_gff_records(gff_dir):
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

            graphic_record = BiopythonTranslator().translate_record(record)
            GFF.write([record], out_handle)

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
            out_name = plots_dir / f"{record.id}.png"
            ax.figure.savefig(out_name, bbox_inches="tight")
            plt.clf()
            plt.close("all")

    checkpoint_file.touch()
