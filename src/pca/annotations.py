"""Annotation plotting and GFF/GenBank generation."""

from __future__ import annotations

import html
import json
import logging
import os
import re
import warnings
from collections.abc import Iterator, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import BiopythonWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqUtils import gc_fraction
from BCBio import GFF
from dna_features_viewer import BiopythonTranslator
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter, MaxNLocator
from matplotlib.transforms import Bbox

from pca.parsers import parse_hmmsearch, parse_trf_dat, parse_trna_gff

__all__ = [
    "convert_to_html",
    "create_feature",
    "get_coordinates",
    "get_cordinates",
    "generate_plots_and_gff",
    "generate_plots_and_annotations",
]

logger = logging.getLogger(__name__)

DEFAULT_PLOT_FORMATS = ("pdf",)
SUPPORTED_PLOT_FORMATS = {"pdf", "png", "html"}

# Display label/category used when a CDS has no PHROG HMM hit.
_NO_PHROG_MATCH_CATEGORY = "no phrogs match"

# Display category used when a CDS has a PHROG hit but the PHROG metadata
# does not assign it a functional category.
_UNKNOWN_FUNCTION_CATEGORY = "unknown function"

# Fallback colors used for standalone HTML conversion when no PHROG metadata
# CSV is supplied. The palette matches the default PHROG category colors.
_DEFAULT_CATEGORY_COLORS: dict[str, str] = {
    "head and packaging": "#3e83f6",
    "tail": "#07e9a2",
    "lysis": "#f35f49",
    "DNA, RNA and nucleotide metabolism": "#ffdf59",
    "integration and excision": "#fea328",
    "transcription regulation": "#a861e3",
    "connector": "#35d7ff",
    "moron, auxiliary metabolic gene and host takeover": "#ff59f5",
    "other": "#808000",
    _UNKNOWN_FUNCTION_CATEGORY: "#c9c9c9",
    _NO_PHROG_MATCH_CATEGORY: "#a0a0a0",
    "tRNA": "#e377c2",
    "tandem repeat": "#17becf",
    # Optional extra database categories
    "antimicrobial resistance": "#d62728",
    "virulence factor": "#e377c2",
    "anti-phage defense": "#9467bd",
    "anti-CRISPR": "#bcbd22",
    "toxin-antitoxin": "#8c564b",
    "diversity-generating retroelement": "#17becf",
}

# Map optional extra-database names to their display category.
_EXTRA_DB_CATEGORIES: dict[str, str] = {
    "card": "antimicrobial resistance",
    "vfdb": "virulence factor",
    "defensefinder": "anti-phage defense",
    "anticrisprdb": "anti-CRISPR",
    "netflax": "toxin-antitoxin",
    "dgr": "diversity-generating retroelement",
}

# Static plot layout tunables.
_STATIC_FIGURE_WIDTH = 24  # inches
_STATIC_MAX_BPS_PER_LINE = 50_000
_STATIC_PNG_DPI = 300
_STATIC_LABEL_FONTSIZE = 9
_STATIC_LABEL_MAX_LEN = 35


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
        "category": "tRNA",
        "color": "#e377c2",
        "ID": f"trna_{x['trna_no']}",
    }
    return SeqFeature(
        FeatureLocation(int(x["begin"]), int(x["end"]), strand=int(x["strand"])),
        type="tRNA",
        id=str(x["trna_no"]),
        qualifiers=qualifiers,
    )


def _create_trf_feature(x: pd.Series) -> SeqFeature:
    """Create a Biopython SeqFeature from a normalized TRF repeat row."""
    qualifiers = {
        "source": "TRF",
        "score": x["score"],
        "label": f"{x['period']} bp tandem repeat",
        "category": "tandem repeat",
        "color": _DEFAULT_CATEGORY_COLORS["tandem repeat"],
        "ID": f"trf_{x['trf_no']}",
        "period": x["period"],
        "copies": x["copies"],
        "consensus": x["consensus"],
        "entropy": x["entropy"],
    }
    return SeqFeature(
        FeatureLocation(int(x["begin"]), int(x["end"]), strand=0),
        type="tandem_repeat",
        id=str(x["trf_no"]),
        qualifiers=qualifiers,
    )


def _parse_translation_tables(gff_path: str | os.PathLike[str]) -> dict[str, int]:
    """Extract per-contig translation tables from a prodigal-gv GFF header.

    Prodigal-gv writes the selected ``transl_table`` in the ``# Model Data:``
    line, even when running in automatic metagenomic mode. The preceding
    ``# Sequence Data:`` line gives the contig header via ``seqhdr=...``.
    """
    tables: dict[str, int] = {}
    current_header: str | None = None
    with open(gff_path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith("# Sequence Data:"):
                match = re.search(r'seqhdr="([^"]+)"', line)
                current_header = match.group(1) if match else None
            elif line.startswith("# Model Data:"):
                match = re.search(r'transl_table=(\d+)', line)
                if match and current_header is not None:
                    tables[current_header] = int(match.group(1))
    return tables


def _genbank_translation_table(record: Any) -> int | None:
    """Return the NCBI translation table used in a GenBank record.

    The table is first looked up in the record-level annotations, then in the
    first CDS feature qualifier. ``None`` is returned when no table is found.
    """
    table = record.annotations.get("transl_table")
    if table is not None:
        try:
            return int(table)
        except ValueError:
            pass
    for feature in record.features:
        if feature.type == "CDS":
            value = _qualifier_to_str(feature.qualifiers.get("transl_table"))
            if value:
                try:
                    return int(value)
                except ValueError:
                    continue
    return None


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
    result: dict[str, Any] = {
        "start": int(feature.location.start),
        "end": int(feature.location.end),
        "strand": feature.location.strand or 1,
        "label": _qualifier_to_str(feature.qualifiers.get("label")),
        "locus_tag": _qualifier_to_str(feature.qualifiers.get("locus_tag")),
        "product": _qualifier_to_str(feature.qualifiers.get("product")),
        "category": _qualifier_to_str(feature.qualifiers.get("category")),
        "color": _qualifier_to_str(feature.qualifiers.get("color", "#c9c9c9")),
        "phrog": _qualifier_to_str(feature.qualifiers.get("phrog")),
        "score": _qualifier_to_str(feature.qualifiers.get("score")),
        "eval": _qualifier_to_str(feature.qualifiers.get("eval")),
        "hostDomain": _qualifier_to_str(feature.qualifiers.get("hostDomain")),
        "Avg_#AA": _qualifier_to_str(feature.qualifiers.get("Avg_#AA")),
        "#ColMSA": _qualifier_to_str(feature.qualifiers.get("#ColMSA")),
        "#prot": _qualifier_to_str(feature.qualifiers.get("#prot")),
        "RefSeq": _qualifier_to_str(feature.qualifiers.get("RefSeq")),
        "Pfam_hit": _qualifier_to_str(feature.qualifiers.get("Pfam_hit")),
        "GO_hit": _qualifier_to_str(feature.qualifiers.get("GO_hit")),
        "KO_hit": _qualifier_to_str(feature.qualifiers.get("KO_hit")),
        "trna_type": _qualifier_to_str(feature.qualifiers.get("trna_type")),
        "period": _qualifier_to_str(feature.qualifiers.get("period")),
        "copies": _qualifier_to_str(feature.qualifiers.get("copies")),
        "consensus": _qualifier_to_str(feature.qualifiers.get("consensus")),
        "entropy": _qualifier_to_str(feature.qualifiers.get("entropy")),
    }
    for db_name in _EXTRA_DB_CATEGORIES:
        result[f"{db_name}_hit"] = _qualifier_to_str(feature.qualifiers.get(f"{db_name}_hit"))
        result[f"{db_name}_evalue"] = _qualifier_to_str(feature.qualifiers.get(f"{db_name}_evalue"))
        result[f"{db_name}_score"] = _qualifier_to_str(feature.qualifiers.get(f"{db_name}_score"))
    result["defensefinder_hit"] = _qualifier_to_str(feature.qualifiers.get("defensefinder_hit"))
    result["defensefinder_system"] = _qualifier_to_str(feature.qualifiers.get("defensefinder_system"))
    result["defensefinder_type"] = _qualifier_to_str(feature.qualifiers.get("defensefinder_type"))
    return result


_TABLE_COLUMNS: list[tuple[str, str]] = [
    ("phrog", "PHROG"),
    ("label", "Label"),
    ("locus_tag", "Locus tag"),
    ("product", "Product"),
    ("category", "Category"),
    ("hostDomain", "hostDomain"),
    ("Avg_#AA", "Avg #AA"),
    ("#ColMSA", "#ColMSA"),
    ("#prot", "# Prot"),
    ("RefSeq", "RefSeq"),
    ("Pfam_hit", "Pfam hit"),
    ("GO_hit", "GO hit"),
    ("KO_hit", "KO hit"),
    ("period", "Period"),
    ("copies", "Copies"),
    ("consensus", "Consensus"),
    ("entropy", "Entropy"),
    *[(f"{db}_hit", db.upper()) for db in _EXTRA_DB_CATEGORIES],
]
_TABLE_COLUMNS += [
    ("defensefinder_hit", "DefenseFinder"),
    ("defensefinder_system", "Defense system"),
    ("defensefinder_type", "Defense type"),
    ("start", "Start"),
    ("end", "End"),
    ("strand", "Strand"),
    ("score", "Score"),
    ("eval", "E-value"),
]

# Map CSV column names (after rename) to Biopython SeqFeature qualifier keys.
_METADATA_QUALIFIER_MAP: dict[str, str] = {
    "annot": "label",
    "category": "category",
    "color": "color",
    "phrog": "phrog",
    "score": "score",
    "eval": "eval",
    "hostDomain": "hostDomain",
    "Avg_#AA": "Avg_#AA",
    "#ColMSA": "#ColMSA",
    "#prot": "#prot",
    "RefSeq": "RefSeq",
    "Pfam_hit": "Pfam_hit",
    "GO_hit": "GO_hit",
    "KO_hit": "KO_hit",
}


_D3_HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<script src="https://d3js.org/d3.v7.min.js"></script>
<style>
  :root {
    --bg: #ffffff;
    --text: #1a1a18;
    --text-muted: #666;
    --link: #185fa5;
    --border: #d3d1c7;
    --header-bg: #ebebea;
    --header-text: #888780;
    --body-bg: #ffffff;
    --metric-bg: #f5f5f3;
    --chip-bg: #e6f1fb;
    --chip-text: #185fa5;
    --chip-border: #b5d4f4;
    --btn-bg: #fff;
    --btn-text: #1a1a18;
    --btn-border: #c0bfb8;
    --btn-hover: #ebebea;
    --overview-bg: #fafafa;
    --overview-border: #ddd;
    --chart-bg: #ffffff;
    --chart-border: #eee;
    --chart-axis: #333;
    --chart-axis-text: #555;
    --chart-text: #333;
    --gene-stroke: #333;
    --gene-stroke-strong: #000;
    --center-line: #999;
    --gc-chart-bg: #fafaf9;
    --gc-chart-border: #e3e1db;
    --gc-axis: #e0ded8;
    --gc-axis-text: #8c8a82;
    --gc-zero: #c9c7c1;
    --gc-crosshair: #8c8a82;
    --gc-title: #6b6961;
    --tooltip-bg: rgba(255, 255, 255, 0.95);
    --tooltip-border: #ccc;
    --tooltip-text: #1a1a18;
    --tooltip-shadow: rgba(0,0,0,0.2);
    --table-bg: #fff;
    --table-header-bg: #f0efea;
    --table-header-text: #666;
    --table-row-even: #fafaf9;
    --table-row-hover: rgba(0,0,0,0.03);
    --highlight-overlay: rgba(0,0,0,0.35);
    --legend-box-stroke: #333;
    --overview-viewport-fill: rgba(70,130,180,0.15);
    --overview-viewport-stroke: #4682b4;
  }
  body.dark {
    --bg: #1a1a1a;
    --text: #e8e6e1;
    --text-muted: #a09c90;
    --link: #7ab8ff;
    --border: #3f3e3a;
    --header-bg: #2a2a28;
    --header-text: #9c9a92;
    --body-bg: #1e1e1e;
    --metric-bg: #2a2a28;
    --chip-bg: #1f3a55;
    --chip-text: #7ab8ff;
    --chip-border: #365c82;
    --btn-bg: #2a2a28;
    --btn-text: #e8e6e1;
    --btn-border: #5f5e5a;
    --btn-hover: #3f3e3a;
    --overview-bg: #222220;
    --overview-border: #3f3e3a;
    --chart-bg: #1e1e1e;
    --chart-border: #3f3e3a;
    --chart-axis: #8c8a82;
    --chart-axis-text: #a09c90;
    --chart-text: #c9c7c1;
    --gene-stroke: #8c8a82;
    --gene-stroke-strong: #e8e6e1;
    --center-line: #5f5e5a;
    --gc-chart-bg: #222220;
    --gc-chart-border: #3f3e3a;
    --gc-axis: #3f3e3a;
    --gc-axis-text: #9c9a92;
    --gc-zero: #5f5e5a;
    --gc-crosshair: #a09c90;
    --gc-title: #9c9a92;
    --tooltip-bg: rgba(40, 40, 40, 0.95);
    --tooltip-border: #5f5e5a;
    --tooltip-text: #e8e6e1;
    --tooltip-shadow: rgba(0,0,0,0.5);
    --table-bg: #222220;
    --table-header-bg: #2e2e2c;
    --table-header-text: #a09c90;
    --table-row-even: #2a2a28;
    --table-row-hover: rgba(255,255,255,0.05);
    --highlight-overlay: rgba(255,255,255,0.2);
    --legend-box-stroke: #8c8a82;
    --overview-viewport-fill: rgba(100,160,220,0.25);
    --overview-viewport-stroke: #6aa9e6;
  }
  body {
    font-family: Arial, Helvetica, sans-serif;
    margin: 20px;
    background: var(--bg);
    color: var(--text);
    transition: background 0.2s, color 0.2s;
  }
  #theme-toggle {
    background: transparent;
    border: 0.5px solid var(--border);
    border-radius: 6px;
    padding: 4px 10px;
    cursor: pointer;
    font-size: 12px;
    color: var(--header-text);
    line-height: 1;
  }
  #theme-toggle:hover { background: var(--btn-hover); }
  #controls { margin-bottom: 10px; display: flex; flex-wrap: wrap; gap: 6px; align-items: center; }
  #controls button, #controls label {
    padding: 4px 8px;
    background: var(--btn-bg);
    color: var(--btn-text);
    border: 0.5px solid var(--btn-border);
    border-radius: 6px;
    font-size: 12px;
  }
  #controls button { cursor: pointer; }
  #controls button:hover { background: var(--btn-hover); }
  #controls input[type="checkbox"] { accent-color: var(--link); }
  #overview { width: 100%; height: 55px; margin-bottom: 6px; border: 1px solid var(--overview-border); background: var(--overview-bg); }
  .overview-axis text { font-size: 10px; fill: var(--chart-axis-text); }
  #chart { width: 100%; overflow: hidden; border: 1px solid var(--chart-border); background: var(--chart-bg); cursor: grab; touch-action: none; }
  #chart:active { cursor: grabbing; }
  .axis path, .axis line { fill: none; stroke: var(--chart-axis); shape-rendering: crispEdges; }
  .axis text { font-size: 12px; fill: var(--chart-axis-text); }
  .axis-label { font-size: 14px; fill: var(--chart-axis-text); }
  .track-label { font-size: 12px; font-weight: bold; fill: var(--chart-text); }
  .center-line { stroke: var(--center-line); stroke-width: 1; }
  .gene { stroke: var(--gene-stroke); stroke-width: 1px; cursor: pointer; }
  .gene:hover { stroke: var(--gene-stroke-strong); stroke-width: 2px; }
  .gene.selected { stroke: var(--gene-stroke-strong); stroke-width: 3px; }
  .label { font-size: 11px; fill: var(--chart-text); pointer-events: none; }
  .label.hidden { display: none; }
  .tooltip {
    position: absolute;
    text-align: left;
    padding: 8px;
    font-size: 12px;
    background: var(--tooltip-bg);
    border: 1px solid var(--tooltip-border);
    border-radius: 4px;
    pointer-events: none;
    box-shadow: 2px 2px 6px var(--tooltip-shadow);
    color: var(--tooltip-text);
  }
  .legend { font-size: 12px; }
  .legend-box { stroke: var(--legend-box-stroke); stroke-width: 1px; }
  .legend-title { font-weight: bold; font-size: 13px; fill: var(--text); }
  .legend-label { font-size: 12px; line-height: 1.2; word-wrap: break-word; color: var(--text); }
  #legend { margin-top: 18px; }
  .legend-grid-title {
    font-weight: bold;
    font-size: 13px;
    color: var(--text);
    margin-bottom: 8px;
  }
  .legend-grid {
    display: grid;
    grid-template-columns: repeat(4, minmax(0, 1fr));
    gap: 6px 14px;
  }
  .legend-item {
    display: flex;
    align-items: center;
    gap: 6px;
    font-size: 12px;
    color: var(--text);
    min-width: 0;
  }
  .legend-swatch {
    width: 12px;
    height: 12px;
    border: 0.5px solid var(--legend-box-stroke);
    flex-shrink: 0;
  }
  .legend-text {
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
  }
  #annotation-table { margin-top: 16px; }
  #annotation-table h4 {
    margin: 0 0 10px 0;
    font-size: 13px;
    font-weight: 600;
    color: var(--text);
    text-transform: uppercase;
    letter-spacing: 0.06em;
  }
  .table-wrap {
    border: 0.5px solid var(--border);
    border-radius: 10px;
    overflow-x: auto;
    overflow-y: hidden;
    background: var(--table-bg);
  }
  #annotation-table table {
    border-collapse: collapse;
    width: 100%;
    font-size: 12px;
    background: transparent;
  }
  #annotation-table th {
    background: var(--table-header-bg);
    font-weight: 600;
    font-size: 11px;
    color: var(--table-header-text);
    text-transform: uppercase;
    letter-spacing: 0.05em;
    padding: 8px 10px;
    text-align: left;
    border-bottom: 0.5px solid var(--border);
    white-space: nowrap;
  }
  #annotation-table td {
    padding: 7px 10px;
    text-align: left;
    border-bottom: 0.5px solid var(--border);
    max-width: 180px;
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
    vertical-align: top;
    cursor: pointer;
    color: var(--text);
  }
  #annotation-table td.expander-col {
    width: 80px;
    min-width: 80px;
    max-width: 80px;
    text-align: center;
    vertical-align: middle;
    cursor: default;
  }
  .row-expander {
    padding: 2px 6px;
    font-size: 11px;
    border: 0.5px solid var(--btn-border);
    border-radius: 6px;
    background: var(--btn-bg);
    cursor: pointer;
    color: var(--link);
    white-space: nowrap;
  }
  .row-expander:hover { background: var(--btn-hover); }
  #annotation-table td[data-col="hostDomain"],
  #annotation-table td[data-col="RefSeq"],
  #annotation-table td[data-col="Pfam_hit"],
  #annotation-table td[data-col="GO_hit"],
  #annotation-table td[data-col="KO_hit"],
  #annotation-table td[data-col="label"],
  #annotation-table td[data-col="category"] {
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    min-width: 120px;
    max-width: 220px;
  }
  #annotation-table tr.expanded td[data-col="hostDomain"],
  #annotation-table tr.expanded td[data-col="RefSeq"],
  #annotation-table tr.expanded td[data-col="Pfam_hit"],
  #annotation-table tr.expanded td[data-col="GO_hit"],
  #annotation-table tr.expanded td[data-col="KO_hit"],
  #annotation-table tr.expanded td[data-col="label"],
  #annotation-table tr.expanded td[data-col="category"] {
    white-space: normal;
    word-break: break-word;
    overflow: visible;
    text-overflow: unset;
    max-width: 600px;
  }
  #annotation-table tr:last-child td { border-bottom: none; }
  #annotation-table tr { background-color: var(--row-bg, transparent); transition: background 0.1s; }
  #annotation-table tr:nth-child(even) { background-color: var(--row-bg-even, var(--table-bg)); }
  body.dark #annotation-table tr,
  body.dark #annotation-table tr:nth-child(even) { background-color: var(--table-bg) !important; }
  #annotation-table tr:not(.annotation-row-highlight):hover td { background: var(--table-row-hover); }
  #annotation-table tr.annotation-row-highlight td {
    box-shadow: inset 0 0 0 1000px var(--highlight-overlay);
    font-weight: 600;
  }
  .pagination {
    margin-top: 10px;
    display: flex;
    align-items: center;
    gap: 6px;
    font-size: 12px;
    color: var(--text-muted);
  }
  .pagination button {
    padding: 4px 10px;
    font-size: 12px;
    border: 0.5px solid var(--btn-border);
    border-radius: 6px;
    background: var(--btn-bg);
    cursor: pointer;
    color: var(--btn-text);
  }
  .pagination button:hover:not(:disabled) { background: var(--btn-hover); }
  .pagination button:disabled { opacity: 0.4; cursor: default; }
  .pagination .page-info { margin: 0 4px; }
  .zoom-hint { font-size: 11px; color: var(--text-muted); margin-left: auto; }
  /* ── Report header ────────────────────────────────────────── */
  .ph-header {
    font-family: "SFMono-Regular", "Consolas", "Liberation Mono", monospace;
    border-radius: 12px;
    overflow: hidden;
    margin-bottom: 20px;
    background: var(--body-bg);
    border: 0.5px solid var(--border);
  }
  .ph-topbar {
    display: flex;
    align-items: center;
    justify-content: space-between;
    padding: 10px 18px;
    background: var(--header-bg);
    gap: 10px;
  }
  .ph-tool-badge {
    display: flex;
    align-items: center;
    gap: 8px;
    font-size: 11px;
    color: var(--header-text);
    letter-spacing: 0.04em;
    text-transform: uppercase;
  }
  .ph-dot {
    width: 7px;
    height: 7px;
    border-radius: 50%;
    background: var(--link);
    flex-shrink: 0;
  }
  .ph-timestamp {
    font-size: 11px;
    color: var(--header-text);
    letter-spacing: 0.03em;
  }
  .ph-body {
    padding: 16px 20px 18px;
    background: var(--body-bg);
  }
  .ph-contig-id {
    font-size: 17px;
    font-weight: 500;
    color: var(--text);
    letter-spacing: -0.01em;
    margin: 0 0 14px 0;
    display: flex;
    align-items: center;
    gap: 10px;
    word-break: break-all;
  }
  .ph-chip {
    font-size: 10px;
    font-weight: 500;
    background: var(--chip-bg);
    color: var(--chip-text);
    border: 0.5px solid var(--chip-border);
    border-radius: 4px;
    padding: 2px 7px;
    white-space: nowrap;
    letter-spacing: 0.05em;
    text-transform: uppercase;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  }
  .ph-metrics {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(120px, 1fr));
    gap: 10px;
  }
  .ph-metric {
    background: var(--metric-bg);
    border: 0.5px solid var(--border);
    border-radius: 8px;
    padding: 9px 12px;
  }
  .ph-metric-label {
    font-size: 10px;
    color: var(--text-muted);
    text-transform: uppercase;
    letter-spacing: 0.07em;
    margin-bottom: 4px;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
  }
  .ph-metric-value {
    font-size: 15px;
    font-weight: 500;
    color: var(--text);
    letter-spacing: -0.01em;
    white-space: nowrap;
  }
  .ph-metric-value .unit {
    font-size: 11px;
    color: var(--text-muted);
    font-weight: 400;
    margin-left: 2px;
  }
  #gc-plots {
    margin: 8px 0;
  }
  .gc-chart {
    width: 100%;
    height: 55px;
    margin-bottom: 3px;
    border: 0.5px solid var(--gc-chart-border);
    border-radius: 4px;
    background: var(--gc-chart-bg);
    position: relative;
    overflow: hidden;
  }
  .gc-chart:last-child { margin-bottom: 0; }
  .gc-chart svg { display: block; background: transparent; }
  .gc-chart path.line {
    fill: none;
    stroke-width: 1.25px;
    stroke-linejoin: round;
    stroke-linecap: round;
  }
  .gc-chart .axis text { font-size: 9px; fill: var(--gc-axis-text); }
  .gc-chart .axis path, .gc-chart .axis line { stroke: var(--gc-axis); }
  .gc-chart .zero-line { stroke: var(--gc-zero); stroke-dasharray: 2,2; pointer-events: none; }
  .gc-chart .crosshair {
    stroke: var(--gc-crosshair);
    stroke-dasharray: 2,2;
    stroke-width: 1px;
    pointer-events: none;
  }
  .gc-chart .gc-title {
    font-size: 9px;
    font-weight: 500;
    fill: var(--gc-title);
    pointer-events: none;
  }
</style>
</head>
<body class="__THEME_CLASS__" data-theme="__THEME__">

<div class="ph-header">
  <div class="ph-topbar">
    <div class="ph-tool-badge">
      <div class="ph-dot"></div>
      phage_contig_annotator · PHROG annotation pipeline
    </div>
    <button id="theme-toggle" title="Toggle dark mode" aria-label="Toggle dark mode">🌙</button>
    <div class="ph-timestamp">__HEADER_TIMESTAMP__</div>
  </div>
  <div class="ph-body">
    <div class="ph-contig-id">
      __CONTIG_ID__
      <span class="ph-chip">phage contig</span>
    </div>
    <div class="ph-metrics">
      <div class="ph-metric">
        <div class="ph-metric-label">Contig length</div>
        <div class="ph-metric-value">__CONTIG_LENGTH_FORMATTED__ <span class="unit">bp</span></div>
      </div>
      <div class="ph-metric">
        <div class="ph-metric-label">GC content</div>
        <div class="ph-metric-value">__GC_CONTENT__</div>
      </div>
      <div class="ph-metric">
        <div class="ph-metric-label">Strand bias</div>
        <div class="ph-metric-value">__STRAND_BIAS__</div>
      </div>
      <div class="ph-metric">
        <div class="ph-metric-label">__FEATURE_SUMMARY_LABEL__</div>
        <div class="ph-metric-value">__FEATURE_SUMMARY_VALUE__</div>
      </div>
      <div class="ph-metric">
        <div class="ph-metric-label">Translation table</div>
        <div class="ph-metric-value">__TRANSLATION_TABLE__</div>
      </div>
    </div>
  </div>
</div>
<div id="controls">
  <button id="pan-left">&larr; Pan</button>
  <button id="pan-right">Pan &rarr;</button>
  <button id="zoom-in">Zoom in</button>
  <button id="zoom-out">Zoom out</button>
  <button id="zoom-reset">Reset</button>
  <label><input type="checkbox" id="show-labels" checked> Show labels</label>
  <span class="zoom-hint">Wheel = zoom &middot; Shift+wheel / horizontal swipe = pan</span>
</div>
<div id="overview"></div>
<div id="gc-plots">
  <div class="gc-chart" id="gc-content-chart"></div>
  <div class="gc-chart" id="gc-skew-chart"></div>
  <div class="gc-chart" id="gc-cumskew-chart"></div>
</div>
<div id="chart"></div>
<div id="legend"></div>
<div id="annotation-table"></div>

<script>
(function() {
  const data = __FEATURES_JSON__;
  const categoryColors = __COLORS_JSON__;
  const tableColumns = __COLUMNS_JSON__;
  const contigLength = __CONTIG_LENGTH__;
  const gcPositions = __GC_POSITIONS_JSON__;
  const gcContent = __GC_CONTENT_JSON__;
  const gcSkew = __GC_SKEW_JSON__;
  const gcCumSkew = __GC_CUM_SKEW_JSON__;

  data.forEach(function(d, i) { d._idx = i; });

  // Theme toggle (persisted in localStorage)
  const docBody = document.body;
  const themeToggle = document.getElementById("theme-toggle");
  function updateThemeIcon() {
    themeToggle.textContent = docBody.classList.contains("dark") ? "\u2600" : "\u263E";
  }
  // Respect the theme chosen when the file was generated. Only fall back
  // to a saved browser preference when no explicit theme was baked in.
  const initialTheme = docBody.dataset.theme;
  if (initialTheme === "dark") {
    docBody.classList.add("dark");
  } else if (initialTheme === "light") {
    docBody.classList.remove("dark");
  } else {
    try {
      const savedTheme = localStorage.getItem("pca-theme");
      if (savedTheme === "dark") docBody.classList.add("dark");
      else if (savedTheme === "light") docBody.classList.remove("dark");
    } catch (e) {}
  }
  updateThemeIcon();
  themeToggle.addEventListener("click", function() {
    docBody.classList.toggle("dark");
    try {
      localStorage.setItem("pca-theme", docBody.classList.contains("dark") ? "dark" : "light");
    } catch (e) {}
    updateThemeIcon();
  });

  const margin = {top: 55, right: 20, bottom: 55, left: 80};
  const trackHeight = 55;
  const trackGap = 35;
  const featureHeight = 22;
  const laneGap = 6;
  const laneStep = featureHeight + laneGap;

  // Stack overlapping features on the same strand into separate lanes so that
  // consecutive overlapping genes do not occlude each other.
  function assignLanes(features) {
    const byStrand = {1: [], "-1": []};
    features.forEach(function(d) {
      byStrand[d.strand === 1 ? 1 : -1].push(d);
    });
    [1, -1].forEach(function(strand) {
      const feats = byStrand[strand].sort(function(a, b) {
        return a.start - b.start || a.end - b.end;
      });
      const laneEnds = [];
      feats.forEach(function(d) {
        let lane = -1;
        for (let i = 0; i < laneEnds.length; i++) {
          if (d.start >= laneEnds[i]) {
            lane = i;
            break;
          }
        }
        if (lane === -1) {
          lane = laneEnds.length;
          laneEnds.push(0);
        }
        d.lane = lane;
        laneEnds[lane] = Math.max(laneEnds[lane] || 0, d.end);
      });
    });
  }
  assignLanes(data);

  function featureCenterY(d) {
    const base = d.strand === 1 ? forwardY : reverseY;
    const offset = (d.lane || 0) * laneStep;
    return d.strand === 1 ? base - offset : base + offset;
  }

  const maxForwardLane = d3.max(data.filter(function(d) { return d.strand === 1; }), function(d) { return d.lane || 0; }) || 0;
  const maxReverseLane = d3.max(data.filter(function(d) { return d.strand === -1; }), function(d) { return d.lane || 0; }) || 0;
  const maxLaneOffset = Math.max(maxForwardLane, maxReverseLane) * laneStep;
  const lanePadding = maxLaneOffset + 10;

  const forwardY = lanePadding + trackHeight / 2;
  const reverseY = forwardY + trackHeight + trackGap;
  const trackHeightTotal = reverseY + trackHeight / 2 + lanePadding;
  const plotHeight = trackHeightTotal;
  const overviewHeight = 55;

  const chartDiv = document.getElementById("chart");
  let width = Math.max(600, chartDiv.clientWidth) - margin.left - margin.right;

  const x = d3.scaleLinear().domain([0, contigLength]).range([0, width]);
  let currentXScale = x.copy();
  let currentTransform = d3.zoomIdentity;

  let showLabels = true;
  const minZoom = 1;
  const maxZoom = 50;
  const zoomStep = 1.3;
  const panStep = 0.2;

  // --- Main SVG ---
  const svg = d3.select("#chart")
    .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", plotHeight + margin.top + margin.bottom)
    .style("display", "block");

  const plotWrapper = svg.append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  plotWrapper.append("defs").append("clipPath")
      .attr("id", "plot-clip")
    .append("rect")
      .attr("width", width)
      .attr("height", plotHeight);

  const plotArea = plotWrapper.append("g")
      .attr("class", "plot-area")
      .attr("clip-path", "url(#plot-clip)");

  const staticGroup = plotWrapper.append("g").attr("class", "static-group");

  // --- Zoom behavior ---
  const zoom = d3.zoom()
    .scaleExtent([minZoom, maxZoom])
    .extent([[0, 0], [width, plotHeight]])
    .translateExtent([[0, 0], [width, plotHeight]])
    .on("zoom", function(event) {
      currentTransform = event.transform;
      currentXScale = event.transform.rescaleX(x);
      updateView();
    });

  // The root SVG captures pan/zoom events; genes and labels are drawn
  // on top so their own click/hover handlers are reached first.
  const overlay = svg;
  svg.call(zoom);

  const geneGroup = plotArea.append("g").attr("class", "gene-group");
  const labelGroup = plotArea.append("g").attr("class", "label-group");
  const axisGroup = plotArea.append("g").attr("class", "axis");

  staticGroup.append("text")
    .attr("class", "track-label")
    .attr("x", -10)
    .attr("y", forwardY)
    .attr("text-anchor", "end")
    .attr("dominant-baseline", "middle")
    .text("Forward (+)");

  staticGroup.append("text")
    .attr("class", "track-label")
    .attr("x", -10)
    .attr("y", reverseY)
    .attr("text-anchor", "end")
    .attr("dominant-baseline", "middle")
    .text("Reverse (-)");

  const axisLabel = plotArea.append("text")
    .attr("class", "axis-label")
    .attr("y", plotHeight + 40)
    .style("text-anchor", "middle")
    .text("Position (bp)");

  const tooltip = d3.select("body").append("div")
    .attr("class", "tooltip")
    .style("opacity", 0);

  const arrowHeadPx = 7;  // target arrow-head width in screen pixels
  const targetBp = x.invert(arrowHeadPx) - x.invert(0);

  // --- Overview bar ---
  const overviewDiv = d3.select("#overview");
  const overviewWidth = overviewDiv.node().clientWidth;
  const overviewSvg = overviewDiv.append("svg")
    .attr("width", overviewWidth)
    .attr("height", overviewHeight);
  const overviewX = d3.scaleLinear().domain([0, contigLength]).range([0, overviewWidth]);
  const overviewGroup = overviewSvg.append("g");

  const overviewAxis = d3.axisTop(overviewX)
    .ticks(Math.max(2, Math.floor(overviewWidth / 80)))
    .tickSize(4)
    .tickFormat(d3.format("~s"));
  overviewGroup.append("g")
    .attr("class", "overview-axis")
    .attr("transform", "translate(0,22)")
    .call(overviewAxis)
    .selectAll("text")
    .attr("dy", "-2px");

  const overviewGeneHeight = 8;
  const overviewForwardY = 26;
  const overviewReverseY = 38;
  overviewGroup.selectAll(".overview-gene")
    .data(data.filter(function(d) { return d.phrog || d.label; }))
    .enter()
    .append("rect")
    .attr("class", "overview-gene")
    .attr("x", function(d) { return overviewX(d.start); })
    .attr("y", function(d) { return d.strand === 1 ? overviewForwardY : overviewReverseY; })
    .attr("width", function(d) { return Math.max(1, overviewX(d.end) - overviewX(d.start)); })
    .attr("height", overviewGeneHeight)
    .attr("fill", function(d) { return categoryColors[d.category] || d.color || "#c9c9c9"; })
    .attr("stroke", "none");

  const viewportRect = overviewGroup.append("rect")
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", overviewWidth)
    .attr("height", overviewHeight)
    .attr("fill", "var(--overview-viewport-fill)")
    .attr("stroke", "var(--overview-viewport-stroke)")
    .attr("stroke-width", 2)
    .style("pointer-events", "none");

  function updateOverview() {
    const visibleMin = Math.max(0, currentXScale.invert(0));
    const visibleMax = Math.min(contigLength, currentXScale.invert(width));
    viewportRect
      .attr("x", overviewX(visibleMin))
      .attr("width", Math.max(2, overviewX(visibleMax) - overviewX(visibleMin)));
  }

  overviewSvg.on("click", function(event) {
    const [mx] = d3.pointer(event, overviewGroup.node());
    const bp = Math.max(0, Math.min(contigLength, overviewX.invert(mx)));
    const targetK = Math.max(1, currentTransform.k);
    const tx = width / 2 - targetK * x(bp);
    const newTransform = d3.zoomIdentity.translate(tx, 0).scale(targetK);
    overlay.call(zoom.transform, clampTransform(newTransform));
  });

  // --- GC metrics line charts ---
  function drawGcLineChart(containerId, values, color, title, yFormat, showXAxis) {
    if (!values || values.length === 0) return;
    const w = overviewWidth;
    const chartMargin = {top: 12, right: 0, bottom: showXAxis ? 14 : 2, left: 0};
    const h = Math.max(0, overviewHeight - chartMargin.top - chartMargin.bottom);
    const svg = d3.select("#" + containerId).append("svg")
      .attr("width", w)
      .attr("height", overviewHeight);
    const g = svg.append("g")
      .attr("transform", "translate(0," + chartMargin.top + ")");

    const xScale = d3.scaleLinear().domain([0, contigLength]).range([0, w]);
    const yMin = d3.min(values);
    const yMax = d3.max(values);
    let yDomain;
    if (title === "GC content") {
      yDomain = [0, 100];
    } else {
      const yPad = (yMax - yMin) * 0.05 || Math.abs(yMax) * 0.05 || 1;
      yDomain = [yMin - yPad, yMax + yPad];
    }
    const yScale = d3.scaleLinear()
      .domain(yDomain)
      .range([h, 0])
      .nice();

    const series = values.map(function(v, i) {
      return {pos: gcPositions[i], value: v};
    });

    if (showXAxis) {
      g.append("g")
        .attr("class", "axis")
        .attr("transform", "translate(0," + h + ")")
        .call(d3.axisBottom(xScale)
          .ticks(Math.max(2, Math.floor(w / 80)))
          .tickSize(3)
          .tickFormat(d3.format("~s")));
    } else {
      g.append("line")
        .attr("x1", 0).attr("x2", w)
        .attr("y1", h).attr("y2", h)
        .attr("stroke", "var(--gc-axis)")
        .attr("stroke-width", 0.5);
    }

    if (title !== "GC content" && yMin <= 0 && yMax >= 0) {
      g.append("line")
        .attr("class", "zero-line")
        .attr("x1", 0).attr("x2", w)
        .attr("y1", yScale(0)).attr("y2", yScale(0));
    }

    const line = d3.line()
      .x(function(d) { return xScale(d.pos); })
      .y(function(d) { return yScale(d.value); });

    g.append("path")
      .datum(series)
      .attr("class", "line")
      .attr("d", line)
      .attr("stroke", color);

    svg.append("text")
      .attr("class", "gc-title")
      .attr("x", 5)
      .attr("y", 10)
      .text(title);

    const crosshair = g.append("line")
      .attr("class", "crosshair")
      .attr("y1", 0).attr("y2", h)
      .style("opacity", 0);

    const bisect = d3.bisector(function(d) { return d.pos; }).left;
    const overlay = g.append("rect")
      .attr("width", w).attr("height", h)
      .style("fill", "none")
      .style("pointer-events", "all");

    overlay.on("mousemove", function(event) {
      const mx = d3.pointer(event, g.node())[0];
      const bp = Math.max(0, Math.min(contigLength, xScale.invert(mx)));
      const idx = bisect(series, bp, 1);
      const d0 = series[idx - 1];
      const d1 = series[idx];
      const d = (d0 && d1) ? (bp - d0.pos < d1.pos - bp ? d0 : d1) : (d0 || d1);
      if (!d) return;
      crosshair.attr("x1", xScale(d.pos)).attr("x2", xScale(d.pos)).style("opacity", 1);
      tooltip.transition().duration(50).style("opacity", 0.95);
      tooltip.html("<b>" + title + "</b><br>position: " + d3.format(",")(Math.round(d.pos)) + " bp<br>value: " + (yFormat || d3.format(".3g"))(d.value))
        .style("left", (event.pageX + 10) + "px")
        .style("top", (event.pageY - 28) + "px");
    }).on("mouseout", function() {
      crosshair.style("opacity", 0);
      tooltip.transition().duration(300).style("opacity", 0);
    });
  }

  if (gcPositions.length > 0) {
    drawGcLineChart("gc-content-chart", gcContent, "#5b9bd5", "GC content", d3.format(".1f"), false);
    drawGcLineChart("gc-skew-chart", gcSkew, "#e29642", "GC skew", d3.format(".1f"), false);
    drawGcLineChart("gc-cumskew-chart", gcCumSkew, "#6cbf6b", "Cum. GC skew", d3.format(".1f"), true);
  }

  function genePath(d) {
    const start = currentXScale(d.start);
    const end = currentXScale(d.end);
    const w = Math.max(0, end - start);
    const headPixels = Math.min(x(targetBp), w / 2);
    const cy = featureCenterY(d);
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

  function drawGenes() {
    const genes = geneGroup.selectAll(".gene").data(data);
    const genesEnter = genes.enter().append("path").attr("class", "gene");
    genesEnter.merge(genes)
      .attr("d", genePath)
      .attr("fill", fillColor)
      .attr("data-index", function(d) { return d._idx; })
      .on("click", function(event, d) {
        d3.selectAll(".gene").classed("selected", false);
        d3.select(this).classed("selected", true);
        highlightTableRow(d._idx);
      })
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
    genes.exit().remove();
  }

  function hideOverlappingLabels() {
    if (!showLabels) return;
    const labels = labelGroup.selectAll(".label").nodes();
    const padding = 4;
    const items = labels.map(function(node) {
      const bbox = node.getBBox();
      return {node: node, left: bbox.x, right: bbox.x + bbox.width};
    }).sort(function(a, b) { return a.left - b.left; });

    let lastRight = -Infinity;
    items.forEach(function(item) {
      if (item.left < lastRight + padding) {
        d3.select(item.node).classed("hidden", true);
      } else {
        d3.select(item.node).classed("hidden", false);
        lastRight = item.right;
      }
    });
  }

  function drawLabels() {
    if (!showLabels) {
      labelGroup.selectAll(".label").remove();
      return;
    }
    const labels = labelGroup.selectAll(".label").data(data.filter(function(d) {
      return d.label && d.label.toLowerCase() !== "unknown function"
             && d.label.toLowerCase() !== "no phrogs match";
    }));
    const labelsEnter = labels.enter().append("text").attr("class", "label");
    labelsEnter.merge(labels)
      .attr("x", function(d) { return (currentXScale(d.start) + currentXScale(d.end)) / 2; })
      .attr("y", function(d) {
        const cy = featureCenterY(d);
        return d.strand === 1 ? cy - featureHeight / 2 - 4 : cy + featureHeight / 2 + 12;
      })
      .attr("text-anchor", "middle")
      .classed("hidden", false)
      .text(function(d) { return d.label; });
    labels.exit().remove();

    if (labelsEnter.size() > 0 || labels.size() > 0) {
      requestAnimationFrame(hideOverlappingLabels);
    }
  }

  function drawLegend() {
    const legendItems = Object.entries(categoryColors);
    if (legendItems.length === 0) return;

    const legend = d3.select("#legend");
    legend.selectAll("*").remove();
    legend.append("div")
      .attr("class", "legend-grid-title")
      .text("Category");

    const grid = legend.append("div").attr("class", "legend-grid");
    const items = grid.selectAll(".legend-item")
      .data(legendItems)
      .enter()
      .append("div")
      .attr("class", "legend-item");

    items.append("div")
      .attr("class", "legend-swatch")
      .style("background-color", function(d) { return d[1]; });

    items.append("div")
      .attr("class", "legend-text")
      .attr("title", function(d) { return d[0]; })
      .text(function(d) { return d[0]; });
  }

  function updateView() {
    axisGroup.attr("transform", "translate(0," + (plotHeight + 15) + ")")
      .call(d3.axisBottom(currentXScale).ticks(10).tickFormat(function(d) {
        if (d >= 1000000) return (d / 1000000).toFixed(1) + "M";
        if (d >= 1000) return (d / 1000).toFixed(d % 1000 === 0 ? 0 : 1) + "k";
        return d;
      }));

    axisLabel.attr("x", width / 2);

    geneGroup.selectAll(".gene").attr("d", genePath);
    drawLabels();
    updateOverview();
  }

  function setZoom(level, centerPixel) {
    const k = Math.max(minZoom, Math.min(maxZoom, level));
    const cx = centerPixel === undefined ? width / 2 : centerPixel;
    // New transform: scale around cx, keep current translation as much as possible.
    const tx = cx - k * (cx - currentTransform.x) / currentTransform.k;
    const newTransform = d3.zoomIdentity.translate(tx, 0).scale(k);
    overlay.call(zoom.transform, clampTransform(newTransform));
  }

  function clampTransform(t) {
    const tx = Math.min(0, Math.max(width * (1 - t.k), t.x));
    return d3.zoomIdentity.translate(tx, 0).scale(t.k);
  }

  function panBy(fraction) {
    const visibleSpan = currentXScale.invert(width) - currentXScale.invert(0);
    const shiftBp = fraction * visibleSpan;
    const tx = currentTransform.x - currentTransform.k * x(shiftBp);
    const newTransform = d3.zoomIdentity.translate(tx, 0).scale(currentTransform.k);
    overlay.call(zoom.transform, clampTransform(newTransform));
  }

  d3.select("#zoom-in").on("click", function() {
    setZoom(currentTransform.k * zoomStep, width / 2);
  });
  d3.select("#zoom-out").on("click", function() {
    setZoom(currentTransform.k / zoomStep, width / 2);
  });
  d3.select("#zoom-reset").on("click", function() {
    overlay.call(zoom.transform, d3.zoomIdentity);
  });
  d3.select("#pan-left").on("click", function() { panBy(-panStep); });
  d3.select("#pan-right").on("click", function() { panBy(panStep); });
  d3.select("#show-labels").on("change", function() {
    showLabels = this.checked;
    drawLabels();
  });
  // Wheel: vertical = zoom, horizontal/shift+wheel = pan.
  overlay.on("wheel", function(event) {
    event.preventDefault();
    const isPan = event.shiftKey || Math.abs(event.deltaX) > Math.abs(event.deltaY);
    const point = d3.pointer(event, overlay.node());
    if (isPan) {
      const dx = event.deltaX || event.deltaY;
      const tx = currentTransform.x - dx;
      const newTransform = d3.zoomIdentity.translate(tx, 0).scale(currentTransform.k);
      overlay.call(zoom.transform, clampTransform(newTransform));
    } else {
      const factor = event.deltaY < 0 ? zoomStep : 1 / zoomStep;
      const k = Math.max(minZoom, Math.min(maxZoom, currentTransform.k * factor));
      const tx = point[0] - k * (point[0] - currentTransform.x) / currentTransform.k;
      const newTransform = d3.zoomIdentity.translate(tx, 0).scale(k);
      overlay.call(zoom.transform, clampTransform(newTransform));
    }
  });

  function pastelColor(hex) {
    const rgb = d3.rgb(hex);
    return d3.rgb(
      Math.round(255 - (255 - rgb.r) * 0.25),
      Math.round(255 - (255 - rgb.g) * 0.25),
      Math.round(255 - (255 - rgb.b) * 0.25)
    ).formatHex();
  }

  const tableContainer = d3.select("#annotation-table");
  tableContainer.append("h4")
    .text("Annotations");
  const tableWrap = tableContainer.append("div").attr("class", "table-wrap");
  const table = tableWrap.append("table");
  const theadRow = table.append("thead").append("tr");
  theadRow.append("th").attr("class", "expander-col").text("");
  tableColumns.forEach(function(col) {
    theadRow.append("th").text(col.header);
  });
  const tbody = table.append("tbody");
  const paginationDiv = tableContainer.append("div").attr("class", "pagination");
  paginationDiv.append("button")
    .attr("id", "page-prev")
    .text("Previous")
    .on("click", function() {
      if (currentPage > 0) {
        currentPage--;
        renderTable();
        renderPagination();
      }
    });
  paginationDiv.append("span").attr("class", "page-info");
  paginationDiv.append("button")
    .attr("id", "page-next")
    .text("Next")
    .on("click", function() {
      if ((currentPage + 1) * rowsPerPage < tableData().length) {
        currentPage++;
        renderTable();
        renderPagination();
      }
    });

  const rowsPerPage = 10;
  let currentPage = 0;

  function tableData() {
    return data;
  }

  // Click any table cell to highlight the corresponding feature on the map.
  tbody.on("click", function(event) {
    if (event.target.closest(".row-expander")) return;
    const cell = event.target.closest("td");
    if (!cell) return;
    const row = cell.closest("tr");
    const d = row && row.__data__;
    if (d) highlightFeature(d._idx);
  });

  function renderTable() {
    const tData = tableData();
    const pageData = tData.slice(currentPage * rowsPerPage, (currentPage + 1) * rowsPerPage);
    const rows = tbody.selectAll("tr").data(pageData, function(d) { return d._idx; });
    rows.exit().remove();
    const rowsEnter = rows.enter().append("tr").style("cursor", "pointer");
    const allRows = rowsEnter.merge(rows);
    function rowColor(d) {
      return pastelColor(categoryColors[d.category] || d.color || "#ffffff");
    }
    allRows
      .style("--row-bg", rowColor)
      .style("--row-bg-even", rowColor)
      .each(function(d) {
        const row = d3.select(this);
        row.selectAll("td").remove();
        row.append("td")
          .attr("class", "expander-col")
          .append("button")
          .attr("class", "row-expander")
          .text(d._expanded ? "see less" : "see more")
          .on("click", function(event) {
            event.stopPropagation();
            d._expanded = !d._expanded;
            d3.select(this).text(d._expanded ? "see less" : "see more");
            row.classed("expanded", d._expanded);
          });
        tableColumns.forEach(function(col) {
          let val = d[col.key];
          if (val === undefined || val === null || val === "") val = "-";
          row.append("td")
            .attr("data-col", col.key)
            .text(String(val));
        });
      });
  }

  function renderPagination() {
    const tData = tableData();
    const totalPages = Math.ceil(tData.length / rowsPerPage) || 1;
    if (currentPage >= totalPages) currentPage = totalPages - 1;
    paginationDiv.select(".page-info")
      .text("Page " + (currentPage + 1) + " of " + totalPages + " (" + tData.length + " features)");
    paginationDiv.select("#page-prev").property("disabled", currentPage === 0);
    paginationDiv.select("#page-next").property("disabled", (currentPage + 1) >= totalPages);
  }

  function highlightFeature(index) {
    const feature = data[index];
    const centerBp = (feature.start + feature.end) / 2;
    const visibleSpan = currentXScale.invert(width) - currentXScale.invert(0);
    // If feature is outside the visible window, center on it.
    if (centerBp < currentXScale.invert(0) || centerBp > currentXScale.invert(width) || visibleSpan > contigLength * 0.9) {
      const targetK = currentTransform.k > 2 ? currentTransform.k : 3;
      const tx = width / 2 - targetK * x(centerBp);
      const newTransform = d3.zoomIdentity.translate(tx, 0).scale(targetK);
      overlay.call(zoom.transform, clampTransform(newTransform));
    }
    d3.selectAll(".gene").classed("selected", false);
    d3.select('.gene[data-index="' + index + '"]').classed("selected", true);
    highlightTableRow(index);
  }

  function highlightTableRow(index) {
    const tData = tableData();
    const localIndex = tData.findIndex(function(d) { return d._idx === index; });
    if (localIndex >= 0) {
      const page = Math.floor(localIndex / rowsPerPage);
      if (page !== currentPage) {
        currentPage = page;
        renderTable();
        renderPagination();
      }
      tbody.selectAll("tr").classed("annotation-row-highlight", function(d) {
        return d._idx === index;
      });
      const rowNode = tbody.selectAll("tr").filter(function(d) { return d._idx === index; }).node();
      if (rowNode) rowNode.scrollIntoView({ behavior: "smooth", block: "center" });
    }
  }

  drawGenes();
  drawLegend();
  updateView();
  renderTable();
  renderPagination();
})();
</script>
</body>
</html>
"""


def _compute_gc_metrics(
    seq: str | Seq, num_windows: int = 200
) -> tuple[list[int], list[float], list[float], list[float]]:
    """Return window centers and GC content, GC skew, cumulative GC skew.

    Metrics are computed over contiguous, non-overlapping windows.  GC skew is
    ``(G - C) / (G + C) * 100`` and cumulative GC skew is the running sum of
    the per-window GC skew values.
    """
    seq_str = str(seq).upper()
    length = len(seq_str)
    if length == 0:
        return [], [], [], []

    window_size = max(1, length // num_windows)
    positions: list[int] = []
    gc_values: list[float] = []
    skew_values: list[float] = []

    for start in range(0, length, window_size):
        end = min(start + window_size, length)
        window = seq_str[start:end]
        a = window.count("A")
        t = window.count("T")
        g = window.count("G")
        c = window.count("C")
        valid = a + t + g + c
        positions.append((start + end) // 2)
        gc_values.append((g + c) / valid * 100 if valid else 0.0)
        denom = g + c
        skew_values.append((g - c) / denom * 100 if denom else 0.0)

    cumulative_skew = list(pd.Series(skew_values).cumsum())
    return positions, gc_values, skew_values, cumulative_skew


def _write_interactive_plot(
    record: Any,
    plots_dir: Path,
    contig_length: int,
    category_colors: dict[str, str],
    translation_table: int | None = None,
    theme: str = "light",
    out_path: Path | None = None,
    feature_summary: tuple[str, str] | None = None,
) -> None:
    """Write an interactive D3.js HTML genome map for a contig.

    Genes are drawn on two tracks: forward strand above and reverse strand
    below. A category color legend and hover tooltips are included. A header
    card with contig metrics is rendered above the map.

    Parameters
    ----------
    record:
        Biopython SeqRecord to plot.
    plots_dir:
        Directory used when ``out_path`` is not provided.
    contig_length:
        Length of the contig in base pairs.
    category_colors:
        Mapping of category names to hex colors.
    translation_table:
        NCBI translation table to display in the report header.
    theme:
        ``"light"`` or ``"dark"``.
    out_path:
        Optional explicit output path. When omitted, the file is written to
        ``plots_dir / f"{record.id}.html"``.
    feature_summary:
        Optional ``(label, value_html)`` tuple for the feature summary metric
        in the header. Defaults to ``"PHROG hits"`` / ``"N/M annotated"``.
    """
    features = [_feature_to_dict(f) for f in record.features]
    features.sort(key=lambda f: (f["start"], f["end"]))

    present_categories = {f["category"] for f in features if f["category"]}
    ordered_colors = {
        cat: color for cat, color in category_colors.items() if cat in present_categories
    }

    features_json = json.dumps(features)
    colors_json = json.dumps(ordered_colors)

    # Only include optional metadata columns for which at least one feature has
    # a non-empty value, keeping the predefined column order.
    present_meta_cols = {
        key
        for key, _header in _TABLE_COLUMNS
        if any(feature.get(key) for feature in features)
    }
    # PHROG/score/e-value are pipeline-specific; only show them when present.
    always_cols = {"label", "category", "start", "end", "strand"}
    pipeline_cols = {"phrog", "score", "eval"} & present_meta_cols
    table_columns = [
        {"key": key, "header": header}
        for key, header in _TABLE_COLUMNS
        if key in always_cols | pipeline_cols
        or key in present_meta_cols
    ]
    columns_json = json.dumps(table_columns)

    # Header metrics
    timestamp = datetime.now(timezone.utc).astimezone().strftime(
        "%Y-%m-%d · %H:%M:%S UTC%z"
    )

    if record.seq:
        gc_value = f"{gc_fraction(record.seq) * 100:.1f}<span class=\"unit\">%</span>"
    else:
        gc_value = "— <span class=\"unit\" style=\"font-size:10px\">awaited</span>"

    forward = sum(1 for f in record.features if f.location.strand == 1)
    reverse = sum(1 for f in record.features if f.location.strand == -1)
    if forward > reverse:
        strand_bias = 'Forward <span class="unit">dominant</span>'
    elif reverse > forward:
        strand_bias = 'Reverse <span class="unit">dominant</span>'
    else:
        strand_bias = "Balanced"

    total_genes = sum(1 for f in record.features if f.type == "CDS")
    phrog_hits = sum(
        1
        for f in record.features
        if f.type == "CDS"
        and _qualifier_to_str(f.qualifiers.get("category"))
        not in {_UNKNOWN_FUNCTION_CATEGORY, _NO_PHROG_MATCH_CATEGORY}
    )

    if feature_summary is None:
        feature_summary_label = "PHROG hits"
        feature_summary_value = (
            f'{phrog_hits}/{total_genes} <span class="unit">annotated</span>'
        )
    else:
        feature_summary_label, feature_summary_value = feature_summary

    if translation_table is None:
        translation_value = "Auto <span class=\"unit\">pyrodigal-gv</span>"
    else:
        table_labels = {
            1: "standard",
            4: "mycoplasma",
            11: "bacterial",
            15: "alternative yeasts",
        }
        label = table_labels.get(translation_table, "NCBI table")
        translation_value = (
            f"{translation_table} <span class=\"unit\">{label}</span>"
        )

    if record.seq:
        positions, gc_values, skew_values, cumulative_skew = _compute_gc_metrics(
            record.seq
        )
    else:
        positions, gc_values, skew_values, cumulative_skew = [], [], [], []

    safe_title = html.escape(record.id)
    theme_class = "dark" if theme.lower() == "dark" else ""
    theme_value = "dark" if theme.lower() == "dark" else "light"
    html_content = (
        _D3_HTML_TEMPLATE
        .replace("__TITLE__", safe_title)
        .replace("__FEATURES_JSON__", features_json)
        .replace("__COLORS_JSON__", colors_json)
        .replace("__COLUMNS_JSON__", columns_json)
        .replace("__CONTIG_LENGTH__", str(contig_length))
        .replace("__HEADER_TIMESTAMP__", timestamp)
        .replace("__CONTIG_ID__", safe_title)
        .replace("__CONTIG_LENGTH_FORMATTED__", f"{contig_length:,}")
        .replace("__GC_CONTENT__", gc_value)
        .replace("__STRAND_BIAS__", strand_bias)
        .replace("__FEATURE_SUMMARY_LABEL__", feature_summary_label)
        .replace("__FEATURE_SUMMARY_VALUE__", feature_summary_value)
        .replace("__TRANSLATION_TABLE__", translation_value)
        .replace("__GC_POSITIONS_JSON__", json.dumps(positions))
        .replace("__GC_CONTENT_JSON__", json.dumps(gc_values))
        .replace("__GC_SKEW_JSON__", json.dumps(skew_values))
        .replace("__GC_CUM_SKEW_JSON__", json.dumps(cumulative_skew))
        .replace("__THEME_CLASS__", theme_class)
        .replace("__THEME__", theme_value)
    )

    if out_path is None:
        out_path = plots_dir / f"{record.id}.html"
    out_path.write_text(html_content, encoding="utf-8")


def _load_category_colors(meta_path: str | os.PathLike[str] | None) -> dict[str, str]:
    """Return a category-to-color map from a PHROG metadata CSV or defaults."""
    if meta_path is None:
        return _DEFAULT_CATEGORY_COLORS

    df = pd.read_csv(meta_path, low_memory=False)
    df = df.rename(
        columns={"#phrog": "phrog", "Annotation": "annot", "Category": "category"}
    )
    if "category" not in df.columns or "color" not in df.columns:
        return _DEFAULT_CATEGORY_COLORS
    colors = dict(zip(df["category"].astype(str), df["color"].astype(str)))
    return colors if colors else _DEFAULT_CATEGORY_COLORS


_IGNORED_FEATURE_TYPES = frozenset(
    {"source", "region", "assembly_gap", "gap", "sequence_feature"}
)


def _feature_location_key(feature: Any) -> tuple[int, int, int | None]:
    """Return a hashable key for a feature location."""
    return (int(feature.location.start), int(feature.location.end), feature.location.strand)


def _normalize_features(record: Any) -> None:
    """Ensure SeqFeature qualifiers needed by the HTML report are present.

    Features whose type is ``source``, ``region``, ``gap``, etc. are dropped.
    When both ``gene`` and ``CDS``/RNA features share the same location, the
    ``gene`` feature is dropped so each gene is represented by a single row.
    For the remaining features a ``label``, ``category`` and ``color`` qualifier
    is added when missing.
    """
    record.features = [
        f for f in record.features if f.type not in _IGNORED_FEATURE_TYPES
    ]

    non_gene_locs = {
        _feature_location_key(f)
        for f in record.features
        if f.type != "gene"
    }
    record.features = [
        f
        for f in record.features
        if not (f.type == "gene" and _feature_location_key(f) in non_gene_locs)
    ]

    for feature in record.features:
        quals = feature.qualifiers

        label = _qualifier_to_str(quals.get("label"))
        if not label:
            for candidate in ("gene", "product", "Name", "locus_tag", "ID"):
                label = _qualifier_to_str(quals.get(candidate))
                if label:
                    break
            if not label:
                label = feature.type or "feature"
            quals["label"] = [label]

        category = _qualifier_to_str(quals.get("category"))
        if not category:
            if feature.type == "tRNA":
                category = "tRNA"
            else:
                category = _NO_PHROG_MATCH_CATEGORY
            quals["category"] = [category]

        color = _qualifier_to_str(quals.get("color"))
        if not color:
            quals["color"] = [_DEFAULT_CATEGORY_COLORS.get(category, "#c9c9c9")]


def convert_to_html(
    input_path: str | os.PathLike[str],
    output_path: str | os.PathLike[str],
    *,
    fasta_path: str | os.PathLike[str] | None = None,
    meta_path: str | os.PathLike[str] | None = None,
    theme: str = "light",
) -> list[Path]:
    """Convert a GenBank or GFF file into an interactive HTML genome report.

    Parameters
    ----------
    input_path:
        Path to a ``.gb``, ``.gbk``, ``.genbank``, ``.gff`` or ``.gff3`` file.
    output_path:
        Target ``.html`` file (for single-record inputs) or output directory
        (for multi-record inputs).
    fasta_path:
        Optional nucleotide FASTA file used to attach sequences to GFF records.
        Required for GFF inputs that do not contain an embedded ``##FASTA`` block.
    meta_path:
        Optional PHROG metadata CSV. When omitted, built-in category colors are
        used.
    theme:
        ``"light"`` or ``"dark"``.

    Returns
    -------
    List of generated HTML file paths.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    suffix = input_path.suffix.lower()

    if suffix in {".gb", ".gbk", ".genbank"}:
        records = list(SeqIO.parse(input_path, "genbank"))
    elif suffix in {".gff", ".gff3"}:
        sequences = _read_sequences(fasta_path) if fasta_path else {}
        records = list(_attach_sequences(_iter_gff_records(input_path), sequences))
    else:
        raise ValueError(
            f"Unsupported input format '{suffix}'. "
            "Expected .gb, .gbk, .genbank, .gff or .gff3."
        )

    if not records:
        raise ValueError(f"No records found in {input_path}")

    category_colors = _load_category_colors(meta_path)

    if suffix in {".gff", ".gff3"}:
        translation_tables = _parse_translation_tables(input_path)
    else:
        translation_tables = {}

    if output_path.suffix == ".html":
        if len(records) > 1:
            raise ValueError(
                f"Input contains {len(records)} records; "
                "specify an output directory instead of a single .html file."
            )
        output_path.parent.mkdir(parents=True, exist_ok=True)
        out_paths = [output_path]
    else:
        output_path.mkdir(parents=True, exist_ok=True)
        out_paths = [output_path / f"{record.id}.html" for record in records]

    for record, out_path in zip(records, out_paths):
        _normalize_features(record)
        contig_length = _contig_length(record)
        if suffix in {".gff", ".gff3"}:
            translation_table = translation_tables.get(record.id)
        else:
            translation_table = _genbank_translation_table(record)
        feature_count = len(record.features)
        _write_interactive_plot(
            record,
            plots_dir=out_path.parent,
            contig_length=contig_length,
            category_colors=category_colors,
            translation_table=translation_table,
            theme=theme,
            out_path=out_path,
            feature_summary=(
                "FEATURES",
                f'{feature_count} <span class="unit">features</span>',
            ),
        )
        logger.info("wrote interactive HTML report to %s", out_path)

    return out_paths


def _assign_static_feature_levels(features: Iterable[Any]) -> tuple[int, int]:
    """Assign non-overlapping y-levels to graphic features by strand.

    Same-strand overlapping features are stacked into separate lanes; forward
    features are placed above the axis and reverse features below it.  The
    returned ``(min_level, max_level)`` can be used to set a symmetric y-axis
    range.
    """
    by_strand: dict[int | None, list[Any]] = {1: [], -1: [], None: []}
    for feature in features:
        by_strand.setdefault(feature.strand, []).append(feature)

    def _assign(strand_features: list[Any], sign: int) -> None:
        lane_ends: list[int] = []
        for feature in sorted(strand_features, key=lambda f: (f.start, f.end)):
            assigned = -1
            for i, end in enumerate(lane_ends):
                if feature.start >= end:
                    assigned = i
                    lane_ends[i] = feature.end
                    break
            if assigned == -1:
                assigned = len(lane_ends)
                lane_ends.append(feature.end)
            feature.data["fixed_level"] = sign * (assigned + 1)

    _assign(by_strand.get(1, []), 1)
    _assign(by_strand.get(-1, []), -1)
    # Strand-less features (e.g. unstranded annotations) are treated like forward.
    _assign(by_strand.get(None, []), 1)

    levels = [f.data.get("fixed_level", 0) for f in features]
    return (min(levels), max(levels)) if levels else (0, 0)


def _format_genome_coord(x: float) -> str:
    """Format a genome coordinate with K/M suffixes."""
    if x >= 1_000_000:
        value = x / 1_000_000
        formatted = f"{value:.1f}".rstrip("0").rstrip(".")
        return f"{formatted}M"
    if x >= 1_000:
        value = x / 1_000
        formatted = f"{value:.1f}".rstrip("0").rstrip(".")
        return f"{formatted}K"
    return str(int(x))


def _truncate_label(label: str | None, max_len: int = _STATIC_LABEL_MAX_LEN) -> str | None:
    """Return ``label`` truncated to ``max_len`` characters with an ellipsis."""
    if not label:
        return label
    if len(label) > max_len:
        return label[:max_len].rstrip() + "..."
    return label


def _hide_label_leader_lines(
    ax: Any, positions: Sequence[tuple[float, float]], tol: float = 0.5
) -> None:
    """Hide leader lines attached to label positions that have been hidden."""
    for line in ax.lines:
        xdata, ydata = line.get_data()
        if len(xdata) != 2:
            continue
        for tx, ty in positions:
            if (
                (abs(xdata[0] - tx) < tol and abs(ydata[0] - ty) < tol)
                or (abs(xdata[1] - tx) < tol and abs(ydata[1] - ty) < tol)
            ):
                line.set_visible(False)
                break


def _drop_overlapping_labels(ax: Any, padding_pt: float = 3.0) -> list[tuple[float, float]]:
    """Hide labels whose bounding boxes overlap, keeping lower-left labels first.

    Returns the positions of hidden labels so their leader lines can be removed.
    """
    renderer = ax.figure.canvas.get_renderer()
    hidden_positions: list[tuple[float, float]] = []
    kept: list[tuple[Any, Any]] = []
    for text in sorted(
        ax.texts, key=lambda t: (t.get_position()[1], t.get_position()[0])
    ):
        if not text.get_visible():
            continue
        bbox = text.get_window_extent(renderer)
        padded = bbox.padded(padding_pt)
        if any(padded.overlaps(other[1]) for other in kept):
            text.set_visible(False)
            hidden_positions.append(text.get_position())
        else:
            kept.append((text, bbox))
    return hidden_positions


def _align_labels_to_feature_side(
    ax: Any, features: Sequence[Any], feature_level_height: float
) -> None:
    """Make labels horizontal, move them to the strand side, and keep them off arrows."""
    by_center = {round(f.x_center, 6): f for f in features}
    for text in list(ax.texts):
        if not text.get_visible():
            continue
        tx, ty = text.get_position()
        matched = None
        matched_line = None
        text_idx = -1
        for line in ax.lines:
            xdata, ydata = line.get_data()
            if len(xdata) != 2:
                continue
            if abs(xdata[0] - tx) < 0.5 and abs(ydata[0] - ty) < 0.5:
                matched_line = line
                text_idx = 0
            elif abs(xdata[1] - tx) < 0.5 and abs(ydata[1] - ty) < 0.5:
                matched_line = line
                text_idx = 1
            if matched_line is not None:
                other_x = xdata[1] if text_idx == 0 else xdata[0]
                matched = by_center.get(round(other_x, 6))
                if matched is None:
                    for center, feat in by_center.items():
                        if abs(center - other_x) < 0.5:
                            matched = feat
                            break
                break
        if matched is None:
            continue
        sign = -1 if matched.strand == -1 else 1
        level = matched.data.get("fixed_level", 0) or 0
        min_label_y = sign * (abs(level) + 1.5) * feature_level_height
        if sign * ty < sign * min_label_y:
            ty = min_label_y
        text.set_position((tx, ty))
        if matched_line is not None:
            xdata, ydata = matched_line.get_data()
            new_xdata = list(xdata)
            new_ydata = list(ydata)
            new_xdata[text_idx] = tx
            new_ydata[text_idx] = ty
            matched_line.set_data(new_xdata, new_ydata)
        text.set_rotation(0)
        text.set_fontsize(_STATIC_LABEL_FONTSIZE)
        text.set_bbox(None)
        text.set_horizontalalignment("center")
        if matched.strand == -1:
            text.set_verticalalignment("top")
        else:
            text.set_verticalalignment("bottom")


def _drop_dense_labels(ax: Any, min_gap_pt: float = 4.0) -> list[tuple[float, float]]:
    """Hide labels that are horizontally too close to an already-kept label."""
    renderer = ax.figure.canvas.get_renderer()
    hidden_positions: list[tuple[float, float]] = []
    kept: list[tuple[Any, Any]] = []
    for text in sorted(ax.texts, key=lambda t: t.get_position()[0]):
        if not text.get_visible():
            continue
        bbox = text.get_window_extent(renderer)
        x0, _, x1, _ = bbox.extents
        overlapping = False
        for other_text, other_bbox in kept:
            ox0, _, ox1, _ = other_bbox.extents
            if x1 < ox0 or x0 > ox1:
                continue
            gap = min(abs(x1 - ox0), abs(ox1 - x0))
            if gap < min_gap_pt:
                overlapping = True
                break
        if overlapping:
            text.set_visible(False)
            hidden_positions.append(text.get_position())
        else:
            kept.append((text, bbox))
    return hidden_positions


def _drop_labels_overlapping_patches(
    ax: Any, features: Sequence[Any]
) -> list[tuple[float, float]]:
    """Hide labels whose bounding boxes overlap visible feature patches."""
    renderer = ax.figure.canvas.get_renderer()
    hidden_positions: list[tuple[float, float]] = []
    patch_bboxes = [
        patch.get_window_extent(renderer)
        for patch in ax.patches
        if patch.get_visible()
    ]
    for text in ax.texts:
        if not text.get_visible():
            continue
        tbbox = text.get_window_extent(renderer)
        for pbbox in patch_bboxes:
            if tbbox.overlaps(pbbox):
                text.set_visible(False)
                hidden_positions.append(text.get_position())
                break
    return hidden_positions


def _write_static_plot(
    record: Any,
    plots_dir: Path,
    format: str,
    category_colors: dict[str, str],
    theme: str = "light",
) -> None:
    """Write a static matplotlib genome map for a contig with a category legend.

    Labels are kept horizontal, placed on the same side as their feature's
    strand, and overlapping/dense labels are hidden. Long contigs are split
    into multiple horizontal rows. PNG output is rendered at 300 DPI.
    """
    is_dark = theme.lower() == "dark"
    style = "dark_background" if is_dark else "default"

    def _feature_properties(_feature: Any) -> dict[str, Any]:
        if is_dark:
            return {"linecolor": "white", "label_link_color": "white"}
        return {}

    class _TallTranslator(BiopythonTranslator):
        graphic_record_parameters = {"feature_level_height": 0.4, "labels_spacing": 25}

        def compute_feature_label(self, feature: Any) -> str | None:
            category = _qualifier_to_str(feature.qualifiers.get("category"))
            if category in {
                _UNKNOWN_FUNCTION_CATEGORY,
                _NO_PHROG_MATCH_CATEGORY,
                "tRNA",
            }:
                return None
            return super().compute_feature_label(feature)

    with plt.style.context(style):
        graphic_record = _TallTranslator(
            features_properties=_feature_properties
        ).translate_record(record)

        for feat in graphic_record.features:
            feat.label = _truncate_label(feat.label)

        # Compute global y range so every row has the same height and the
        # bottom margin is large enough for tick labels below reverse features.
        min_level, max_level = _assign_static_feature_levels(graphic_record.features)
        y_min = (min_level - 2.5) * graphic_record.feature_level_height
        y_max = (max_level + 2.5) * graphic_record.feature_level_height

        seq_length = _contig_length(record)
        formatter = FuncFormatter(lambda x, pos: _format_genome_coord(x))

        axes: list[Any] = []
        line_features: list[list[Any]] = []
        if seq_length > _STATIC_MAX_BPS_PER_LINE:
            n_lines = (seq_length + _STATIC_MAX_BPS_PER_LINE - 1) // _STATIC_MAX_BPS_PER_LINE
            nucl_per_line = (seq_length + n_lines - 1) // n_lines
            row_height = max(2.5, (max_level - min_level + 5.5) * graphic_record.feature_level_height)
            fig, raw_axes = plt.subplots(
                n_lines,
                1,
                figsize=(_STATIC_FIGURE_WIDTH, n_lines * row_height),
                gridspec_kw={"hspace": 0.1},
            )
            if n_lines == 1:
                axes = [raw_axes]
            else:
                axes = list(raw_axes)
            for idx, ax in enumerate(axes):
                start = graphic_record.first_index + idx * nucl_per_line
                end = min(graphic_record.last_index, start + nucl_per_line)
                line_record = graphic_record.crop((start, end))
                line_features.append(line_record.features)
                line_record.plot(
                    ax=ax,
                    x_lim=(start, end),
                    with_ruler=True,
                    draw_line=True,
                    annotate_inline=False,
                    elevate_outline_annotations=False,
                    strand_in_label_threshold=None,
                )
        else:
            fig, ax = plt.subplots(1, 1, figsize=(_STATIC_FIGURE_WIDTH, 4))
            graphic_record.plot(
                ax=ax,
                with_ruler=True,
                draw_line=True,
                annotate_inline=False,
                elevate_outline_annotations=False,
                strand_in_label_threshold=None,
            )
            axes = [ax]
            line_features = [graphic_record.features]

        axis_color = "white" if is_dark else "black"
        for ax, feats in zip(axes, line_features):
            ax.xaxis.set_major_formatter(formatter)
            ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
            # Place the x-axis ruler at y=0 (the strand separator) and draw
            # tick marks/labels just below it, matching the reference layout.
            ax.spines["bottom"].set_position(("data", 0))
            ax.spines["bottom"].set_color(axis_color)
            ax.spines["top"].set_visible(False)
            ax.spines["left"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.xaxis.set_ticks_position("bottom")
            ax.tick_params(
                axis="x",
                direction="out",
                length=5,
                pad=4,
                labelsize=_STATIC_LABEL_FONTSIZE,
                colors=axis_color,
                bottom=True,
                top=False,
                labelbottom=True,
                labeltop=False,
            )
            for tick_label in ax.get_xticklabels():
                tick_label.set_rotation(0)
                tick_label.set_horizontalalignment("center")
                tick_label.set_verticalalignment("top")

            _align_labels_to_feature_side(
                ax, feats, graphic_record.feature_level_height
            )
            hidden_positions: list[tuple[float, float]] = []
            hidden_positions.extend(_drop_overlapping_labels(ax))
            hidden_positions.extend(_drop_dense_labels(ax))
            hidden_positions.extend(_drop_labels_overlapping_patches(ax, feats))
            _hide_label_leader_lines(ax, hidden_positions)

            ax.set_ylim(y_min, y_max)
            ax.set_yticks([])

        fig.suptitle(record.id, fontsize=14)

        present_categories = {
            _qualifier_to_str(f.qualifiers.get("category"))
            for f in record.features
            if _qualifier_to_str(f.qualifiers.get("category"))
        }
        if present_categories:
            legend_edge = "#dddddd" if is_dark else "#333333"
            legend_handles = [
                Patch(facecolor=color, edgecolor=legend_edge, label=category)
                for category, color in category_colors.items()
                if category in present_categories
            ]
            if legend_handles:
                fig.legend(
                    handles=legend_handles,
                    loc="lower center",
                    bbox_to_anchor=(0.5, -0.02),
                    ncol=min(4, len(legend_handles)),
                    frameon=True,
                    fontsize="small",
                    title="Category",
                )
                fig.subplots_adjust(bottom=0.12)

        out_path = plots_dir / f"{record.id}.{format}"
        savefig_kwargs = {"bbox_inches": "tight", "format": format, "pad_inches": 0.3}
        if format == "png":
            savefig_kwargs["dpi"] = _STATIC_PNG_DPI
        fig.savefig(str(out_path), **savefig_kwargs)
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
    translation_table: int | None = None,
    theme: str = "light",
    trf_path: str | os.PathLike[str] | None = None,
    extra_hits: dict[str, str | os.PathLike[str]] | None = None,
    defensefinder_tsv: str | os.PathLike[str] | None = None,
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
        Path to the PHROG metadata CSV file.
    gff_dir:
        Path to the protein GFF file produced by gene calling.
    input_fasta:
        Path to the original nucleotide FASTA (used for GenBank and contig length).
    plot_formats:
        Iterable of plot formats to produce. Supported: ``pdf``, ``png``, ``html``.
        Default is ``("pdf",)``.
    translation_table:
        NCBI translation table used for gene calling. Displayed in the HTML
        report header. When the protein GFF contains prodigal-gv ``Model Data``
        headers, the actual per-contig table is extracted and this value is
        used only as a fallback. ``None`` means pyrodigal-gv's metagenomic
        auto-selection.
    trf_path:
        Optional path to a Tandem Repeats Finder ``.dat`` file.  When provided,
        repeat regions are added to the GFF/GenBank records and plots.
    extra_hits:
        Optional mapping of extra database name (e.g. ``"card"``) to a CSV of
        HMMER/phmmer hits produced by :func:`pca.pipeline.search_extra_db`.
    defensefinder_tsv:
        Optional path to a DefenseFinder ``defense_finder_genes.tsv`` file.
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
    translation_tables = _parse_translation_tables(gff_dir)
    if translation_tables:
        logger.info("parsed translation tables from GFF: %s", translation_tables)

    logger.info("processing hmm results")
    search_results = pd.DataFrame(parse_hmmsearch(hmmsearch_dir))

    trna_path = Path(trna_dir)
    if trna_path.is_file() and trna_path.stat().st_size > 0:
        trna = pd.DataFrame(parse_trna_gff(trna_dir))
    else:
        trna = pd.DataFrame()

    trf_dat = Path(trf_path) if trf_path else None
    if trf_dat is not None and trf_dat.is_file() and trf_dat.stat().st_size > 0:
        repeats = pd.DataFrame(parse_trf_dat(trf_dat))
    else:
        repeats = pd.DataFrame()

    if search_results.empty:
        raise ValueError("hmmsearch returned zero matches")

    if not trna.empty:
        trna = trna.apply(get_coordinates, axis=1).sort_values(by=["contig", "begin"])
    else:
        trna = None

    if not repeats.empty:
        repeats = repeats.sort_values(by=["contig", "begin"]).reset_index(drop=True)
        repeats["trf_no"] = repeats.index + 1
    else:
        repeats = None

    extra_hits = extra_hits or {}
    extra_results: dict[str, pd.DataFrame] = {}
    for db_name, csv_path in extra_hits.items():
        csv_path = Path(csv_path)
        if not csv_path.is_file() or csv_path.stat().st_size == 0:
            extra_results[db_name] = pd.DataFrame()
            continue
        try:
            df = pd.read_csv(csv_path)
        except pd.errors.EmptyDataError:
            extra_results[db_name] = pd.DataFrame()
            continue
        if not df.empty:
            df = df.sort_values(
                ["qname", "score", "eval"], ascending=[True, False, True]
            ).drop_duplicates("qname")
        extra_results[db_name] = df

    defensefinder_results = pd.DataFrame()
    if defensefinder_tsv:
        df_path = Path(defensefinder_tsv)
        if df_path.is_file() and df_path.stat().st_size > 0:
            defensefinder_results = pd.DataFrame(
                parse_defensefinder_genes(defensefinder_tsv)
            )
            if not defensefinder_results.empty:
                defensefinder_results = defensefinder_results.sort_values(
                    ["qname"]
                ).drop_duplicates("qname")

    phrogs_anno = pd.read_csv(meta_dir, low_memory=False)
    phrogs_anno = phrogs_anno.rename(
        columns={"#phrog": "phrog", "Annotation": "annot", "Category": "category"}
    )
    phrogs_anno = phrogs_anno.fillna("unknown function")
    category_colors = {
        **_DEFAULT_CATEGORY_COLORS,
        **{
            cat: color
            for cat, color in zip(phrogs_anno["category"], phrogs_anno["color"])
            if pd.notna(color)
        },
    }
    if trna is not None:
        category_colors["tRNA"] = "#e377c2"
    if repeats is not None:
        category_colors["tandem repeat"] = _DEFAULT_CATEGORY_COLORS["tandem repeat"]
    category_colors["other"] = _DEFAULT_CATEGORY_COLORS["other"]

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
                qname = f"{record.id}_{pos}"
                available_meta_cols = [
                    col for col in _METADATA_QUALIFIER_MAP if col in tmp.columns
                ]
                tmp_feature = tmp[tmp["position"] == pos][available_meta_cols]
                has_phrog = not tmp_feature.empty
                if has_phrog:
                    for col in available_meta_cols:
                        qualifier_key = _METADATA_QUALIFIER_MAP[col]
                        feature.qualifiers.update(
                            {qualifier_key: tmp_feature[col].values[0]}
                        )
                else:
                    feature.qualifiers.update({"label": _NO_PHROG_MATCH_CATEGORY})
                    feature.qualifiers.update({"category": _NO_PHROG_MATCH_CATEGORY})
                    feature.qualifiers.update(
                        {"color": _DEFAULT_CATEGORY_COLORS[_NO_PHROG_MATCH_CATEGORY]}
                    )

                for db_name, db_df in extra_results.items():
                    if db_df.empty or "qname" not in db_df.columns:
                        continue
                    db_hit = db_df[db_df["qname"] == qname]
                    if not db_hit.empty:
                        hit = db_hit.iloc[0]
                        feature.qualifiers[f"{db_name}_hit"] = str(hit["tname"])
                        feature.qualifiers[f"{db_name}_evalue"] = f"{float(hit['eval']):.2g}"
                        feature.qualifiers[f"{db_name}_score"] = f"{float(hit['score']):.1f}"
                        if not has_phrog:
                            category = _EXTRA_DB_CATEGORIES.get(
                                db_name, _UNKNOWN_FUNCTION_CATEGORY
                            )
                            feature.qualifiers["category"] = category
                            feature.qualifiers["color"] = _DEFAULT_CATEGORY_COLORS.get(
                                category, "#c9c9c9"
                            )
                            feature.qualifiers["label"] = str(hit["tname"])

                if not defensefinder_results.empty:
                    df_hit = defensefinder_results[defensefinder_results["qname"] == qname]
                    if not df_hit.empty:
                        hit = df_hit.iloc[0]
                        feature.qualifiers["defensefinder_hit"] = str(hit["gene_name"])
                        feature.qualifiers["defensefinder_system"] = str(hit["system"])
                        feature.qualifiers["defensefinder_type"] = str(hit["type"])
                        if not has_phrog and not any(
                            feature.qualifiers.get(f"{db_name}_hit")
                            for db_name in _EXTRA_DB_CATEGORIES
                        ):
                            category = _EXTRA_DB_CATEGORIES["defensefinder"]
                            feature.qualifiers["category"] = category
                            feature.qualifiers["color"] = _DEFAULT_CATEGORY_COLORS.get(
                                category, "#c9c9c9"
                            )
                            feature.qualifiers["label"] = str(hit["gene_name"])

                # Ensure the plotted/written color matches the chosen category
                # palette (overrides any stale color from the metadata CSV).
                category = _qualifier_to_str(feature.qualifiers.get("category"))
                if category:
                    feature.qualifiers["color"] = category_colors.get(
                        category, _DEFAULT_CATEGORY_COLORS.get(category, "#c9c9c9")
                    )

            if trna is not None:
                tmp_trna = trna[trna["contig"] == record.id]
                if not tmp_trna.empty:
                    tmp_trna = tmp_trna.reset_index(drop=True)
                    tmp_trna["feature"] = tmp_trna.apply(create_feature, axis=1)
                    record.features.extend(tmp_trna["feature"].to_list())

            if repeats is not None:
                tmp_repeats = repeats[repeats["contig"] == record.id]
                if not tmp_repeats.empty:
                    tmp_repeats = tmp_repeats.reset_index(drop=True)
                    tmp_repeats["feature"] = tmp_repeats.apply(_create_trf_feature, axis=1)
                    record.features.extend(tmp_repeats["feature"].to_list())

            annotated_records.append(record)
            GFF.write([record], out_handle)

            contig_length = _contig_length(record)

            if "html" in requested_formats:
                record_translation_table = translation_tables.get(
                    record.id, translation_table
                )
                _write_interactive_plot(
                    record, plots_dir, contig_length, category_colors, record_translation_table, theme
                )

            static_formats = requested_formats & {"pdf", "png"}
            for fmt in static_formats:
                _write_static_plot(record, plots_dir, fmt, category_colors, theme)

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
