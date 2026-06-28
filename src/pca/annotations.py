"""Annotation plotting and GFF/GenBank generation."""

from __future__ import annotations

import html
import json
import logging
import os
import warnings
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
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


def _qualifier_to_str(value: Any) -> str:
    """Return a qualifier value as a plain string."""
    if value is None:
        return ""
    if isinstance(value, list):
        return str(value[0]) if value else ""
    return str(value)


def _feature_to_dict(feature: Any) -> dict[str, Any]:
    """Return a JSON-serializable dict describing a SeqFeature."""
    return {
        "start": int(feature.location.start),
        "end": int(feature.location.end),
        "strand": feature.location.strand or 1,
        "label": _qualifier_to_str(feature.qualifiers.get("label")),
        "category": _qualifier_to_str(feature.qualifiers.get("category")),
        "color": _qualifier_to_str(feature.qualifiers.get("color", "#c9c9c9")),
        "phrog": _qualifier_to_str(feature.qualifiers.get("phrog")),
        "score": _qualifier_to_str(feature.qualifiers.get("score")),
        "eval": _qualifier_to_str(feature.qualifiers.get("eval")),
        "trna_type": _qualifier_to_str(feature.qualifiers.get("trna_type")),
    }


_D3_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>__TITLE__</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  body { font-family: Arial, Helvetica, sans-serif; margin: 20px; background: #fff; }
  #chart { width: 100%; overflow-x: auto; }
  .axis path, .axis line { fill: none; stroke: #333; shape-rendering: crispEdges; }
  .axis text { font-size: 12px; }
  .track-label { font-size: 12px; font-weight: bold; fill: #333; }
  .center-line { stroke: #999; stroke-width: 1; }
  .gene { stroke: #333; stroke-width: 1px; cursor: pointer; }
  .gene:hover { stroke: #000; stroke-width: 2px; }
  .label { font-size: 11px; fill: #000; pointer-events: none; }
  .tooltip {
    position: absolute;
    text-align: left;
    padding: 8px;
    font-size: 12px;
    background: rgba(255, 255, 255, 0.95);
    border: 1px solid #ccc;
    border-radius: 4px;
    pointer-events: none;
    box-shadow: 2px 2px 6px rgba(0,0,0,0.2);
  }
  .legend { font-size: 12px; }
  .legend-box { stroke: #333; stroke-width: 1px; }
  .legend-title { font-weight: bold; font-size: 13px; }
  .legend-label { font-size: 12px; line-height: 1.2; word-wrap: break-word; color: #000; }
</style>
</head>
<body>
<h3 style="text-align:center;">__TITLE__</h3>
<div id="chart"></div>
<script>
(function() {
  const data = __FEATURES_JSON__;
  const categoryColors = __COLORS_JSON__;
  const contigLength = __CONTIG_LENGTH__;

  const margin = {top: 60, right: 200, bottom: 70, left: 80};
  const plotWidth = Math.max(900, Math.min(1800, contigLength / 35));
  const width = plotWidth - margin.left - margin.right;
  const trackHeight = 55;
  const trackGap = 35;
  const forwardY = trackHeight / 2;
  const reverseY = trackHeight * 1.5 + trackGap;
  const trackHeightTotal = reverseY + trackHeight / 2 + 20;

  // Legend is placed to the right; ensure the SVG is tall enough.
  const legendRowHeight = 44;
  const legendHeight = legendRowHeight * Object.keys(categoryColors).length + 30;
  const height = Math.max(trackHeightTotal, legendHeight);

  const svg = d3.select("#chart")
    .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  const x = d3.scaleLinear()
    .domain([0, contigLength])
    .range([0, width]);

  // Axis
  svg.append("g")
    .attr("class", "axis")
    .attr("transform", "translate(0," + (height + 15) + ")")
    .call(d3.axisBottom(x).ticks(10).tickFormat(function(d) {
      if (d >= 1000000) return (d / 1000000).toFixed(1) + "M";
      if (d >= 1000) return (d / 1000).toFixed(d % 1000 === 0 ? 0 : 1) + "k";
      return d;
    }));

  svg.append("text")
    .attr("x", width / 2)
    .attr("y", height + 50)
    .style("text-anchor", "middle")
    .style("font-size", "14px")
    .text("Position (bp)");

  // Track labels
  svg.append("text")
    .attr("class", "track-label")
    .attr("x", -10)
    .attr("y", forwardY)
    .attr("text-anchor", "end")
    .attr("dominant-baseline", "middle")
    .text("Forward (+)");

  svg.append("text")
    .attr("class", "track-label")
    .attr("x", -10)
    .attr("y", reverseY)
    .attr("text-anchor", "end")
    .attr("dominant-baseline", "middle")
    .text("Reverse (-)");

  // Center lines
  svg.append("line")
    .attr("class", "center-line")
    .attr("x1", 0).attr("x2", width)
    .attr("y1", forwardY).attr("y2", forwardY);
  svg.append("line")
    .attr("class", "center-line")
    .attr("x1", 0).attr("x2", width)
    .attr("y1", reverseY).attr("y2", reverseY);

  const tooltip = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);

  const arrowHeadBp = Math.max(150, contigLength * 0.006);
  const featureHeight = 22;

  function genePath(d) {
    const start = x(d.start);
    const end = x(d.end);
    const width = Math.max(0, end - start);
    const headPixels = Math.min(x(arrowHeadBp) - x(0), width / 2);
    const cy = d.strand === 1 ? forwardY : reverseY;
    const y0 = cy - featureHeight / 2;
    const y1 = cy + featureHeight / 2;
    const mid = cy;

    if (d.strand === 1) {
      return "M" + start + "," + y0 +
             "H" + (end - headPixels) +
             "L" + end + "," + mid +
             "L" + (end - headPixels) + "," + y1 +
             "H" + start + "Z";
    } else {
      return "M" + (start + headPixels) + "," + y0 +
             "H" + end +
             "V" + y1 +
             "H" + (start + headPixels) +
             "L" + start + "," + mid + "Z";
    }
  }

  function fillColor(d) {
    return categoryColors[d.category] || d.color || "#c9c9c9";
  }

  // Draw genes
  svg.selectAll(".gene")
    .data(data)
    .enter()
    .append("path")
    .attr("class", "gene")
    .attr("d", genePath)
    .attr("fill", fillColor)
    .on("mouseover", function(event, d) {
      tooltip.transition().duration(150).style("opacity", 0.95);
      tooltip.html("<b>" + (d.label || "feature") + "</b>" +
        "<br>category: " + (d.category || "n/a") +
        "<br>phrog: " + (d.phrog || "n/a") +
        "<br>score: " + (d.score || "n/a") +
        "<br>e-value: " + (d.eval || "n/a"))
        .style("left", (event.pageX + 10) + "px")
        .style("top", (event.pageY - 28) + "px");
    })
    .on("mousemove", function(event) {
      tooltip.style("left", (event.pageX + 10) + "px")
             .style("top", (event.pageY - 28) + "px");
    })
    .on("mouseout", function() {
      tooltip.transition().duration(300).style("opacity", 0);
    });

  // Labels
  svg.selectAll(".label")
    .data(data.filter(function(d) {
      return d.label && d.label.toLowerCase() !== "unknown function";
    }))
    .enter()
    .append("text")
    .attr("class", "label")
    .attr("x", function(d) { return (x(d.start) + x(d.end)) / 2; })
    .attr("y", function(d) {
      const cy = d.strand === 1 ? forwardY : reverseY;
      return d.strand === 1 ? cy - featureHeight / 2 - 4 : cy + featureHeight / 2 + 12;
    })
    .attr("text-anchor", "middle")
    .text(function(d) { return d.label; });

  // Legend
  const legendItems = Object.entries(categoryColors);
  if (legendItems.length > 0) {
    const legend = svg.append("g")
      .attr("class", "legend")
      .attr("transform", "translate(" + (width + 20) + ",20)");

    legend.append("text")
      .attr("class", "legend-title")
      .attr("x", 0)
      .attr("y", 0)
      .text("Category");
    const rows = legend.selectAll(".legend-item")
      .data(legendItems)
      .enter()
      .append("g")
      .attr("class", "legend-item")
      .attr("transform", function(d, i) { return "translate(0," + ((i + 1) * legendRowHeight) + ")"; });

    rows.append("rect")
      .attr("width", 12)
      .attr("height", 12)
      .attr("fill", function(d) { return d[1]; })
      .attr("class", "legend-box");

    rows.append("foreignObject")
      .attr("x", 18)
      .attr("y", -2)
      .attr("width", 170)
      .attr("height", legendRowHeight)
      .append("xhtml:div")
      .attr("class", "legend-label")
      .text(function(d) { return d[0]; });
  }
})();
</script>
</body>
</html>
"""


def _write_interactive_plot(
    record: Any,
    plots_dir: Path,
    contig_length: int,
    category_colors: dict[str, str],
) -> None:
    """Write an interactive D3.js HTML genome map for a contig.

    Genes are drawn on two tracks: forward strand above and reverse strand
    below. A category color legend and hover tooltips are included.
    """
    features = [_feature_to_dict(f) for f in record.features]
    features.sort(key=lambda f: (f["start"], f["end"]))

    present_categories = {f["category"] for f in features if f["category"]}
    ordered_colors = {
        cat: color for cat, color in category_colors.items() if cat in present_categories
    }

    features_json = json.dumps(features)
    colors_json = json.dumps(ordered_colors)

    safe_title = html.escape(record.id)
    html_content = (
        _D3_HTML_TEMPLATE
        .replace("__TITLE__", safe_title)
        .replace("__FEATURES_JSON__", features_json)
        .replace("__COLORS_JSON__", colors_json)
        .replace("__CONTIG_LENGTH__", str(contig_length))
    )

    out_path = plots_dir / f"{record.id}.html"
    out_path.write_text(html_content, encoding="utf-8")


def _write_static_plot(
    record: Any,
    plots_dir: Path,
    format: str,
    category_colors: dict[str, str],
) -> None:
    """Write a static matplotlib genome map for a contig with a category legend."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

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

    present_categories = {
        _qualifier_to_str(f.qualifiers.get("category"))
        for f in record.features
        if _qualifier_to_str(f.qualifiers.get("category"))
    }
    if present_categories:
        legend_handles = [
            Patch(facecolor=color, edgecolor="#333333", label=category)
            for category, color in category_colors.items()
            if category in present_categories
        ]
        if legend_handles:
            ax.legend(
                handles=legend_handles,
                loc="upper left",
                bbox_to_anchor=(1.01, 1.0),
                title="Category",
                frameon=True,
            )
            fig.subplots_adjust(right=0.82)

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
    annotations_dir = tmp_dir / "annotations"
    annotations_dir.mkdir(parents=True, exist_ok=True)

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
    category_colors = dict(zip(phrogs_anno["category"], phrogs_anno["color"]))

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
    gff_out_path = annotations_dir / "annotations.gff"
    gbk_out_path = annotations_dir / "annotations.gbk"

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
                    feature.qualifiers.update({"category": "unknown"})
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
                _write_interactive_plot(record, plots_dir, contig_length, category_colors)

            static_formats = requested_formats & {"pdf", "png"}
            for fmt in static_formats:
                _write_static_plot(record, plots_dir, fmt, category_colors)

    # Write per-format sentinel files so workflow engines can depend on
    # specific formats rather than the whole plots directory.
    for fmt in requested_formats:
        (plots_dir / f".{fmt}_done").touch()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", BiopythonWarning)
        SeqIO.write(annotated_records, gbk_out_path, "genbank")
    logger.info("wrote %s and %s", gff_out_path, gbk_out_path)


# Backwards-compatible alias.
def generate_plots_and_gff(*args: Any, **kwargs: Any) -> None:
    """Deprecated alias for :func:`generate_plots_and_annotations`."""
    return generate_plots_and_annotations(*args, **kwargs)
