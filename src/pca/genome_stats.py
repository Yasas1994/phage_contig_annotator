"""Genomic statistics calculator for phage contigs.

This module computes per-contig metrics such as strand-switch rates and
tandem-repeat density. It is designed to work on the gene-calling GFF
produced by the pipeline together with the input FASTA and the TRF ``.dat``
output.
"""

from __future__ import annotations

import logging
from collections.abc import Sequence
from itertools import groupby
from pathlib import Path
from typing import Any

import pandas as pd
from BCBio import GFF

from pca.io import get_compressed_file_handle, read_fasta
from pca.multicopy import analyze
from pca.parsers import parse_crispr_cas_finder_gff, parse_minced_gff, parse_trf_dat

logger = logging.getLogger(__name__)

__all__ = [
    "strand_switch_frequency",
    "strand_run_stats",
    "weighted_strand_switch_frequency",
    "strand_bias_index",
    "tandem_repeat_stats",
    "contig_stats_from_gff",
    "compute_genome_stats",
]


def _normalize_strand(strand: int | str | None) -> str:
    """Return ``+`` or ``-`` for a Biopython or GFF strand value."""
    if strand in (1, "+", "plus"):
        return "+"
    if strand in (-1, "-", "minus"):
        return "-"
    return "+"


def strand_switch_frequency(strands: Sequence[str]) -> float:
    """Return the strand switch frequency (SSF) for an ordered strand list.

    SSF = number of strand switches / (N - 1). Returns ``0.0`` for fewer
    than two genes.
    """
    if len(strands) < 2:
        return 0.0
    switches = sum(1 for i in range(len(strands) - 1) if strands[i] != strands[i + 1])
    return switches / (len(strands) - 1)


def strand_run_stats(strands: Sequence[str]) -> dict[str, Any]:
    """Return run-length statistics for an ordered strand list.

    Metrics returned:

    - ``n_runs``: number of consecutive same-strand runs
    - ``mean_run_length``: average run length
    - ``max_run_length``: longest run length
    """
    if not strands:
        return {
            "n_runs": 0,
            "mean_run_length": 0.0,
            "max_run_length": 0,
        }
    runs = [list(g) for _, g in groupby(strands)]
    run_lengths = [len(r) for r in runs]
    n = len(strands)
    return {
        "n_runs": len(runs),
        "mean_run_length": sum(run_lengths) / len(run_lengths),
        "max_run_length": max(run_lengths),
    }


def strand_bias_index(strands: Sequence[str]) -> float:
    """Return the strand bias index (SBI): (N+ - N-) / (N+ + N-).

    Values are in ``[-1, 1]``. ``0`` means balanced strands; values near
    ``+1`` or ``-1`` mean strong strand bias.
    """
    n = len(strands)
    if n == 0:
        return 0.0
    n_plus = sum(1 for s in strands if s == "+")
    return (n_plus - (n - n_plus)) / n


def weighted_strand_switch_frequency(genes: Sequence[dict[str, Any]]) -> float:
    """Return strand switch frequency weighted by inter-gene span.

    ``genes`` must be dicts with integer ``start`` and ``strand`` keys, sorted
    by start position. Longer inter-gene distances contribute more to the
    metric than short ones.
    """
    if len(genes) < 2:
        return 0.0
    total_span = 0.0
    switched_span = 0.0
    for i in range(len(genes) - 1):
        span = genes[i + 1]["start"] - genes[i]["start"]
        if span <= 0:
            continue
        total_span += span
        if genes[i]["strand"] != genes[i + 1]["strand"]:
            switched_span += span
    return switched_span / total_span if total_span > 0 else 0.0


def tandem_repeat_stats(
    trf_path: str | Path,
    contig_id: str,
    contig_length: int,
) -> dict[str, Any]:
    """Return tandem-repeat statistics for ``contig_id`` from a TRF ``.dat`` file."""
    count = 0
    total_length = 0
    for record in parse_trf_dat(trf_path):
        if record.get("contig") != contig_id:
            continue
        count += 1
        total_length += int(record["end"]) - int(record["begin"]) + 1
    density = total_length / contig_length if contig_length > 0 else 0.0
    return {
        "n_tandem_repeats": count,
        "total_tandem_repeat_length": total_length,
        "tandem_repeat_density": density,
    }


def crispr_stats(
    crispr_gff_path: str | Path,
    contig_id: str,
) -> dict[str, Any]:
    """Return CRISPR array statistics for ``contig_id`` from a GFF file.

    Accepts both CRISPRCasFinder GFFs (feature type ``CRISPR``) and MinCED
    GFFs (feature type ``repeat_region`` with ``rpt_family=CRISPR``).
    """
    count = 0
    total_length = 0
    for parser in (parse_crispr_cas_finder_gff, parse_minced_gff):
        for record in parser(crispr_gff_path):
            if record.get("qname") != contig_id:
                continue
            count += 1
            total_length += record["end"] - record["begin"] + 1
    return {
        "n_crispr_arrays": count,
        "total_crispr_array_length": total_length,
    }


def _is_gene_feature(feature_type: str) -> bool:
    """Return True if ``feature_type`` is a gene-level feature."""
    return feature_type.lower() in {"cds", "gene", "exon", "mrna"}


def _genes_from_gff(gff_path: str | Path) -> dict[str, list[dict[str, Any]]]:
    """Return a mapping contig -> sorted gene records from a GFF file.

    Each gene record is a dict with ``start``, ``end``, and ``strand``.
    """
    genes: dict[str, list[dict[str, Any]]] = {}
    with get_compressed_file_handle(gff_path) as handle:
        for record in GFF.parse(handle):
            contig_genes: list[dict[str, Any]] = []
            for feature in record.features:
                if not _is_gene_feature(feature.type):
                    continue
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = _normalize_strand(feature.location.strand)
                contig_genes.append({"start": start, "end": end, "strand": strand})
            if contig_genes:
                contig_genes.sort(key=lambda g: (g["start"], g["end"]))
                genes[record.id] = contig_genes
    return genes


def contig_stats_from_gff(
    gff_path: str | Path,
    contig_id: str | None = None,
) -> dict[str, dict[str, Any]]:
    """Compute strand-based statistics from a GFF file.

    If ``contig_id`` is provided, only that contig is returned. Otherwise all
    contigs present in the GFF are returned.
    """
    genes_by_contig = _genes_from_gff(gff_path)
    result: dict[str, dict[str, Any]] = {}
    for cid, genes in genes_by_contig.items():
        if contig_id is not None and cid != contig_id:
            continue
        strands = [g["strand"] for g in genes]
        stats: dict[str, Any] = {
            "n_genes": len(genes),
            "strand_switch_frequency": strand_switch_frequency(strands),
            "weighted_strand_switch_frequency": weighted_strand_switch_frequency(genes),
            "strand_bias_index": strand_bias_index(strands),
        }
        stats.update(strand_run_stats(strands))
        result[cid] = stats
    return result


def compute_genome_stats(
    fasta_path: str | Path,
    gff_path: str | Path,
    trf_path: str | Path | None = None,
    crispr_path: str | Path | None = None,
    minced_path: str | Path | None = None,
    multicopy_path: str | Path | None = None,
) -> pd.DataFrame:
    """Compute a per-contig statistics table.

    Parameters
    ----------
    fasta_path:
        Path to a nucleotide FASTA file.
    gff_path:
        Path to a gene-calling GFF file. Used to count genes and compute
        strand-switch metrics.
    trf_path:
        Optional Tandem Repeats Finder ``.dat`` file. If provided, tandem repeat
        statistics are included.
    crispr_path:
        Optional CRISPRCasFinder GFF file. If provided, CRISPR array counts
        are included.
    minced_path:
        Optional MinCED GFF file. If provided, CRISPR array counts are
        included and combined with any CRISPRCasFinder results.
    multicopy_path:
        Optional multi-copy detection TSV. If provided, copy-number
        estimates and confidence are merged into the table.

    Returns
    -------
    A ``pandas.DataFrame`` with one row per contig and columns for length,
    gene count, strand metrics, and (when optional inputs are given)
    tandem-repeat and CRISPR metrics.
    """
    lengths = {seq_id: len(seq) for seq_id, seq in read_fasta(fasta_path)}
    gff_stats = contig_stats_from_gff(gff_path)

    rows: list[dict[str, Any]] = []
    for contig_id in sorted(gff_stats):
        length = lengths.get(contig_id, 0)
        stats = gff_stats[contig_id]
        stats["contig_id"] = contig_id
        stats["length"] = length
        if trf_path is not None and length > 0:
            stats.update(tandem_repeat_stats(trf_path, contig_id, length))
        if crispr_path is not None:
            crispr_stats_ = crispr_stats(crispr_path, contig_id)
            stats["n_crispr_arrays"] = stats.get("n_crispr_arrays", 0) + crispr_stats_["n_crispr_arrays"]
            stats["total_crispr_array_length"] = stats.get("total_crispr_array_length", 0) + crispr_stats_["total_crispr_array_length"]
        if minced_path is not None:
            minced_stats_ = crispr_stats(minced_path, contig_id)
            stats["n_crispr_arrays"] = stats.get("n_crispr_arrays", 0) + minced_stats_["n_crispr_arrays"]
            stats["total_crispr_array_length"] = stats.get("total_crispr_array_length", 0) + minced_stats_["total_crispr_array_length"]
        rows.append(stats)

    df = pd.DataFrame(rows)

    multicopy_df = pd.DataFrame()
    if multicopy_path is not None:
        multicopy_path_obj = Path(multicopy_path)
        if multicopy_path_obj.is_file() and multicopy_path_obj.stat().st_size > 0:
            multicopy_df = pd.read_csv(multicopy_path_obj, sep="\t", low_memory=False)
            if not multicopy_df.empty:
                keep_cols = [
                    "contig_id",
                    "mean_kmer_freq",
                    "copies_kmer",
                    "validator",
                    "validator_score",
                    "validator_snr",
                    "validator_ok",
                    "copies_final",
                    "confidence",
                    "flag",
                    "note",
                ]
                multicopy_df = multicopy_df[[c for c in keep_cols if c in multicopy_df.columns]]
                multicopy_df = multicopy_df.rename(
                    columns={"confidence": "multicopy_confidence", "flag": "multicopy_flag"}
                )
                df = df.merge(multicopy_df, on="contig_id", how="left")

    if not df.empty:
        column_order = [
            "contig_id",
            "length",
            "n_genes",
            "strand_switch_frequency",
            "weighted_strand_switch_frequency",
            "strand_bias_index",
            "n_runs",
            "mean_run_length",
            "max_run_length",
            "n_tandem_repeats",
            "total_tandem_repeat_length",
            "tandem_repeat_density",
            "n_crispr_arrays",
            "total_crispr_array_length",
            "mean_kmer_freq",
            "copies_kmer",
            "validator",
            "validator_score",
            "validator_snr",
            "validator_ok",
            "copies_final",
            "multicopy_confidence",
            "multicopy_flag",
            "note",
        ]
        # Only include columns that actually exist (optional inputs may be missing).
        df = df[[c for c in column_order if c in df.columns]]
    return df
