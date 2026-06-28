"""Parsers for external-tool tabular outputs."""

from __future__ import annotations

import logging
from collections.abc import Iterator
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "parse_blastp",
    "parse_hmmsearch",
    "parse_trf_dat",
    "parse_trna",
    "parse_trna_gff",
]


def parse_blastp(path: str | Path) -> Iterator[dict[str, Any]]:
    """Parse BLAST tabular format 6 output."""
    names = [
        "qname", "tname", "pid", "aln", "mis", "gap",
        "qstart", "qstop", "tstart", "tstop", "eval", "score",
    ]
    formats = [str, str, float, int, int, int, int, int, int, int, float, float]
    with open(path) as f:
        for line in f:
            values = line.split()
            if len(values) < 12:
                continue
            yield {
                names[i]: formats[i](values[i])
                for i in range(12)
            }


def parse_hmmsearch(path: str | Path) -> Iterator[dict[str, Any]]:
    """Parse HMMER tblout output."""
    names = [
        "qname", "qacc", "tname", "tacc", "eval",
        "score", "bias", "beval", "bscore", "bbias",
    ]
    formats = [str, str, str, str, float, float, float, float, float, float]
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            values = line.split()
            if len(values) < 10:
                logger.debug("skipping malformed hmmsearch line: %s", line.rstrip())
                continue
            try:
                yield {names[i]: formats[i](values[i]) for i in range(10)}
            except (ValueError, IndexError):
                logger.debug("skipping erroneous hmmsearch line: %s", line.rstrip())


def parse_trf_dat(path: str | Path) -> Iterator[dict[str, Any]]:
    """Parse Tandem Repeats Finder (TRF) ``.dat`` output.

    TRF writes one block per input sequence.  Each block starts with a
    ``Sequence:`` header, followed by a ``Parameters:`` line and then one line
    per tandem repeat with the 15 space-delimited fields described in the TRF
    documentation.
    """
    names = [
        "begin",
        "end",
        "period",
        "copies",
        "consensus_size",
        "matches",
        "indels",
        "score",
        "a",
        "c",
        "g",
        "t",
        "entropy",
        "consensus",
        "repeat_seq",
    ]
    formats = [
        int,
        int,
        int,
        float,
        int,
        int,
        int,
        int,
        float,
        float,
        float,
        float,
        float,
        str,
        str,
    ]

    contig: str | None = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("Sequence:"):
                contig = line.split(":", 1)[1].strip()
                continue
            if contig is None or not line or line.startswith("Parameters"):
                continue
            values = line.split()
            if len(values) < 15:
                continue
            try:
                record = {names[i]: formats[i](values[i]) for i in range(15)}
            except (ValueError, IndexError) as exc:
                logger.debug("skipping erroneous TRF line: %s (%s)", line, exc)
                continue
            record["contig"] = contig
            record["strand"] = 0
            yield record


def parse_trna(path: str | Path) -> Iterator[dict[str, Any]]:
    """Parse tRNAscan-SE tabular output."""
    names = [
        "qname", "trna_no", "begin", "end", "trna_type",
        "anticodon", "intron_begin", "intron_end", "score",
    ]
    formats = [str, int, int, int, str, str, int, int, float]
    with open(path) as f:
        for line in f:
            if line.startswith(("Sequence", "Name", "-----")):
                continue
            values = line.split()
            if len(values) < 9:
                continue
            try:
                yield {names[i]: formats[i](values[i]) for i in range(9)}
            except (ValueError, IndexError) as exc:
                logger.debug("skipping erroneous tRNA line: %s (%s)", line.rstrip(), exc)


def parse_trna_gff(path: str | Path) -> Iterator[dict[str, Any]]:
    """Parse tRNAscan-SE GFF output."""
    names = [
        "qname", "begin", "end", "score", "strand",
        "trna_no", "trna_type", "anticodon",
    ]
    indices = [0, 3, 4, 5, 6]
    formats = [str, int, int, float, str, int, str, str]
    with open(path) as f:
        for line in f:
            if line.startswith("#") or "exon" in line:
                continue
            values = line.split()
            if len(values) < 9:
                continue
            v = [values[i] for i in indices]
            meta = values[-1].split(";")
            for field in meta:
                if field.startswith("ID"):
                    v.append(int(field.split("trna")[-1]))
                elif field.startswith("isotype"):
                    v.append(field.split("=")[-1])
                elif field.startswith("anticodon"):
                    v.append(field.split("=")[-1])
            if len(v) != 8:
                logger.debug("skipping malformed tRNA GFF line: %s", line.rstrip())
                continue
            try:
                yield {names[i]: formats[i](v[i]) for i in range(8)}
            except (ValueError, IndexError) as exc:
                logger.debug("skipping erroneous tRNA GFF line: %s (%s)", line.rstrip(), exc)
