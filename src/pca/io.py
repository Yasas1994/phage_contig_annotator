"""Input/output helpers for FASTA and compressed files."""

from __future__ import annotations

import bz2
import gzip
import logging
import lzma
import os
import re
import shutil
from collections import defaultdict
from collections.abc import Iterable, Iterator
from enum import Enum, auto
from pathlib import Path
from typing import Any, TextIO

import pandas as pd
from Bio import SeqIO
from BCBio import GFF

logger = logging.getLogger(__name__)

__all__ = [
    "Compression",
    "check_fasta",
    "convert_to_fasta",
    "detect_sequence_format",
    "get_compressed_file_handle",
    "is_compressed",
    "read_fasta",
    "split_csv_by_contig",
    "split_fasta",
    "split_fasta_by_contig",
    "split_genbank_by_contig",
    "split_gff_by_contig",
    "split_hmmsearch_by_contig",
    "split_trf_dat_by_contig",
    "split_trna_tsv_by_contig",
    "split_tsv_by_contig",
]


class Compression(Enum):
    """Supported input compression formats."""

    gzip = auto()
    bzip2 = auto()
    xz = auto()
    noncompressed = auto()


def is_compressed(filepath: str | os.PathLike[str]) -> Compression:
    """Check whether a file is gzip, bzip2, xz, or uncompressed."""
    with open(filepath, "rb") as fin:
        signature = fin.peek(8)[:8]

    if tuple(signature[:2]) == (0x1F, 0x8B):
        return Compression.gzip
    if tuple(signature[:3]) == (0x42, 0x5A, 0x68):
        return Compression.bzip2
    if tuple(signature[:7]) == (0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00, 0x00):
        return Compression.xz
    return Compression.noncompressed


def get_compressed_file_handle(path: str | os.PathLike[str]) -> TextIO:
    """Return a text-mode file handle for a possibly-compressed file."""
    compression = is_compressed(path)
    if compression == Compression.gzip:
        return gzip.open(path, "rt")
    if compression == Compression.bzip2:
        return bz2.open(path, "rt")
    if compression == Compression.xz:
        return lzma.open(path, "rt")
    return open(path, "r")


def _first_non_blank_line(path: str | os.PathLike[str]) -> str | None:
    """Return the first non-empty line of a possibly-compressed text file."""
    with get_compressed_file_handle(path) as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                return stripped
    return None


def detect_sequence_format(path: str | os.PathLike[str]) -> str:
    """Detect whether ``path`` is FASTA, GenBank, or EMBL.

    Returns one of ``"fasta"``, ``"genbank"``, ``"embl"``, or ``"unknown"``.
    Detection is based on the first non-empty line of the (possibly
    compressed) file.
    """
    first_line = _first_non_blank_line(path)
    if first_line is None:
        return "unknown"
    if first_line.startswith(">"):
        return "fasta"
    if first_line.startswith("LOCUS"):
        return "genbank"
    if first_line.startswith("ID   "):
        return "embl"
    return "unknown"


def convert_to_fasta(
    input_path: str | os.PathLike[str],
    output_path: str | os.PathLike[str],
) -> None:
    """Convert a GenBank or EMBL file to a multi-FASTA file.

    Raises ``ValueError`` if the input format cannot be determined or if no
    records are found.
    """
    fmt = detect_sequence_format(input_path)
    if fmt == "unknown":
        raise ValueError(
            f"Unable to determine sequence format for {input_path}. "
            "Expected FASTA, GenBank, or EMBL."
        )
    if fmt == "fasta":
        shutil.copy2(input_path, output_path)
        return

    seq_format = "genbank" if fmt == "genbank" else "embl"
    count = SeqIO.write(
        SeqIO.parse(get_compressed_file_handle(input_path), seq_format),
        output_path,
        "fasta",
    )
    if count == 0:
        raise ValueError(f"No sequences found in {input_path}.")


def read_fasta(path: str | os.PathLike[str]) -> Iterator[tuple[str, str]]:
    """Yield (sequence id, sequence) tuples from a possibly compressed FASTA file."""
    with get_compressed_file_handle(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            name = record.id
            seq = str(record.seq).upper()
            if name and seq:
                yield name, seq


def split_fasta(
    input_path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    n_chunks: int,
) -> list[Path]:
    """Split a multi-FASTA into ``n_chunks`` chunk files.

    Chunks are written to ``output_dir`` as ``chunk_0.fasta``,
    ``chunk_1.fasta``, and so on. If ``n_chunks`` exceeds the number of
    records, it is capped so that no empty chunk is produced. Records are kept
    in their original order within each chunk.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    records = list(SeqIO.parse(get_compressed_file_handle(input_path), "fasta"))
    n_chunks = max(1, min(n_chunks, len(records))) if records else 1
    if n_chunks == 1:
        single = output_dir / "chunk_0.fasta"
        SeqIO.write(records, single, "fasta")
        return [single]

    chunk_size = (len(records) + n_chunks - 1) // n_chunks
    chunks: list[list[Any]] = [[] for _ in range(n_chunks)]
    for idx, record in enumerate(records):
        chunk_idx = idx // chunk_size
        chunks[chunk_idx].append(record)

    paths: list[Path] = []
    for i, chunk_records in enumerate(chunks):
        chunk_path = output_dir / f"chunk_{i}.fasta"
        SeqIO.write(chunk_records, chunk_path, "fasta")
        paths.append(chunk_path)

    return paths


def _derive_contig_id(protein_id: str) -> str:
    """Return the contig ID embedded in a ``{contig}_{position}`` protein ID.

    The position suffix is the trailing underscore-separated numeric component.
    Contig IDs that contain underscores are handled correctly because the
    position is always the final component.
    """
    parts = protein_id.rsplit("_", 1)
    if len(parts) == 2 and parts[1].isdigit():
        return parts[0]
    return protein_id


def split_fasta_by_contig(
    input_path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    filename: str = "proteins.faa",
    contig_ids: Iterable[str] | None = None,
) -> dict[str, Path]:
    """Split a protein FASTA into per-contig files.

    Protein record IDs are expected to be ``{contig_id}_{position}``. Each
    contig gets a subdirectory under ``output_dir`` containing ``filename``.
    If ``contig_ids`` is supplied, only those contigs are emitted; otherwise
    contigs are discovered from the protein IDs.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    allowed = set(contig_ids) if contig_ids is not None else None
    handles: dict[str, Any] = {}
    paths: dict[str, Path] = {}

    try:
        with get_compressed_file_handle(input_path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                contig = _derive_contig_id(record.id)
                if allowed is not None and contig not in allowed:
                    continue
                if contig not in handles:
                    contig_dir = output_dir / contig
                    contig_dir.mkdir(parents=True, exist_ok=True)
                    out_path = contig_dir / filename
                    handles[contig] = open(out_path, "w")
                    paths[contig] = out_path
                SeqIO.write([record], handles[contig], "fasta")
    finally:
        for fh in handles.values():
            fh.close()

    return paths


def split_gff_by_contig(
    input_path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    filename: str = "proteins.gff",
) -> dict[str, Path]:
    """Split a GFF file into per-contig files using the record ID."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: dict[str, Path] = {}
    with open(input_path) as handle:
        for record in GFF.parse(handle):
            contig_dir = output_dir / record.id
            contig_dir.mkdir(parents=True, exist_ok=True)
            out_path = contig_dir / filename
            with open(out_path, "w") as out_handle:
                GFF.write([record], out_handle)
            paths[record.id] = out_path

    return paths


def split_hmmsearch_by_contig(
    input_path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    filename: str = "hmmsearch.txt",
) -> dict[str, Path]:
    """Split an HMMER tblout file into per-contig files.

    The query name (first column) is expected to be a ``{contig}_{position}``
    protein ID. Comment lines are omitted from the per-contig outputs.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    handles: dict[str, Any] = {}
    paths: dict[str, Path] = {}

    try:
        with open(input_path) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                qname = line.split("\t", 1)[0] if "\t" in line else line.split(None, 1)[0]
                contig = _derive_contig_id(qname)
                if contig not in handles:
                    contig_dir = output_dir / contig
                    contig_dir.mkdir(parents=True, exist_ok=True)
                    out_path = contig_dir / filename
                    handles[contig] = open(out_path, "w")
                    paths[contig] = out_path
                handles[contig].write(line)
    finally:
        for fh in handles.values():
            fh.close()

    return paths


def split_csv_by_contig(
    input_path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    filename: str = "hmmsearch.csv",
    contig_column: str = "contig",
    sep: str = ",",
) -> dict[str, Path]:
    """Split a CSV/TSV into per-contig files by a column value.

    For protein-derived tables (e.g. HMMER CSV) the column may be ``qname``;
    the contig is then derived from ``{contig}_{position}`` IDs. For tables
    that already contain a contig name, pass that column name.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: dict[str, Path] = {}
    try:
        df = pd.read_csv(input_path, sep=sep)
    except pd.errors.EmptyDataError:
        return paths

    if df.empty:
        return paths

    if contig_column not in df.columns:
        logger.warning("split_csv_by_contig: column %r not found in %s", contig_column, input_path)
        return paths

    if contig_column == "qname":
        df["contig"] = df[contig_column].apply(_derive_contig_id)
        contig_column = "contig"

    for contig, group in df.groupby(contig_column, sort=False):
        contig = str(contig)
        contig_dir = output_dir / contig
        contig_dir.mkdir(parents=True, exist_ok=True)
        out_path = contig_dir / filename
        group.to_csv(out_path, sep=sep, index=False)
        paths[contig] = out_path

    return paths


def split_tsv_by_contig(
    input_path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    filename: str = "trna.tsv",
    contig_column: str = "qname",
) -> dict[str, Path]:
    """Split a TSV into per-contig files by a column value."""
    return split_csv_by_contig(input_path, output_dir, filename, contig_column, sep="\t")


def split_trna_tsv_by_contig(
    input_path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    filename: str = "trna.tsv",
) -> dict[str, Path]:
    """Split a tRNAscan-SE tabular output into per-contig files.

    The tRNAscan-SE TSV has a multi-line header block ending with a dashed
    separator line. The header is copied into each per-contig file and data
    rows are routed based on the first column (sequence name).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    header_lines: list[str] = []
    header_done = False
    handles: dict[str, Any] = {}
    paths: dict[str, Path] = {}

    try:
        with open(input_path) as handle:
            for line in handle:
                if not line.strip():
                    continue
                if not header_done:
                    header_lines.append(line)
                    if line.startswith("-") or line.lstrip().startswith("-"):
                        header_done = True
                    continue

                values = line.split("\t")
                if not values:
                    continue
                contig = values[0].strip()
                if contig not in handles:
                    contig_dir = output_dir / contig
                    contig_dir.mkdir(parents=True, exist_ok=True)
                    out_path = contig_dir / filename
                    fh = open(out_path, "w")
                    fh.writelines(header_lines)
                    handles[contig] = fh
                    paths[contig] = out_path
                handles[contig].write(line)
    finally:
        for fh in handles.values():
            fh.close()

    return paths


def split_trf_dat_by_contig(
    input_path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    filename: str = "trf.dat",
) -> dict[str, Path]:
    """Split TRF ``.dat`` output into per-contig files by ``Sequence:`` header."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Preserve the file header (everything before the first Sequence: line).
    header_lines: list[str] = []
    header_done = False
    current_contig: str | None = None
    handles: dict[str, Any] = {}
    paths: dict[str, Path] = {}

    try:
        with open(input_path) as handle:
            for line in handle:
                if line.startswith("Sequence:"):
                    header_done = True
                    current_contig = line.split()[1].strip()
                    if current_contig not in handles:
                        contig_dir = output_dir / current_contig
                        contig_dir.mkdir(parents=True, exist_ok=True)
                        out_path = contig_dir / filename
                        fh = open(out_path, "w")
                        fh.writelines(header_lines)
                        handles[current_contig] = fh
                        paths[current_contig] = out_path
                    handles[current_contig].write(line)
                elif header_done and current_contig is not None:
                    handles[current_contig].write(line)
                else:
                    header_lines.append(line)
    finally:
        for fh in handles.values():
            fh.close()

    return paths


def split_genbank_by_contig(
    input_path: str | os.PathLike[str],
    output_dir: str | os.PathLike[str],
    filename: str = "annotations.gbk",
) -> dict[str, Path]:
    """Split a GenBank file into per-contig files using the record ID."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: dict[str, Path] = {}
    for record in SeqIO.parse(input_path, "genbank"):
        contig_dir = output_dir / record.id
        contig_dir.mkdir(parents=True, exist_ok=True)
        out_path = contig_dir / filename
        SeqIO.write([record], out_path, "genbank")
        paths[record.id] = out_path

    return paths


def check_fasta(path: str | os.PathLike[str], tmp_dir: str | os.PathLike[str]) -> None:
    """Validate an input FASTA file and check for duplicate sequence IDs."""
    checkpoint_file = Path(tmp_dir) / "input_validation_checkpoint"
    if checkpoint_file.is_file():
        return

    seen: set[str] = set()
    duplicates: set[str] = set()
    with get_compressed_file_handle(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq_id = record.id
            if seq_id in seen:
                duplicates.add(seq_id)
            else:
                seen.add(seq_id)

    if not seen:
        raise ValueError("Input FASTA file is empty or not properly formatted.")

    if duplicates:
        raise ValueError(
            "Please remove duplicated sequence IDs from the input FASTA file: "
            f"{', '.join(sorted(duplicates))}"
        )

    checkpoint_file.touch()
