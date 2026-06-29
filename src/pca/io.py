"""Input/output helpers for FASTA and compressed files."""

from __future__ import annotations

import bz2
import gzip
import logging
import lzma
import os
import shutil
from collections.abc import Iterator
from enum import Enum, auto
from pathlib import Path
from typing import TextIO

from Bio import SeqIO

logger = logging.getLogger(__name__)

__all__ = [
    "Compression",
    "check_fasta",
    "convert_to_fasta",
    "detect_sequence_format",
    "get_compressed_file_handle",
    "is_compressed",
    "read_fasta",
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
