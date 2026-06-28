"""Input/output helpers for FASTA and compressed files."""

from __future__ import annotations

import bz2
import gzip
import logging
import lzma
import os
from collections.abc import Iterator
from enum import Enum, auto
from pathlib import Path
from typing import TextIO

from Bio import SeqIO

logger = logging.getLogger(__name__)

__all__ = ["Compression", "check_fasta", "get_compressed_file_handle", "is_compressed", "read_fasta"]


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
