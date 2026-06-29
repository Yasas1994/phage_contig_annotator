"""Backward-compatible public API facade.

New code should import from the focused submodules (``io``, ``validation``,
``pipeline``, ``annotations``, ``parsers``, ``externals``, ``parallel``,
``logutils``) instead of this module.
"""

from __future__ import annotations

from pca.annotations import (
    convert_to_html,
    create_feature,
    generate_plots_and_annotations,
    generate_plots_and_gff,
    get_coordinates,
    get_cordinates,
)
from pca.externals import run_trnascan
from pca.io import (
    Compression,
    check_fasta,
    convert_to_fasta,
    detect_sequence_format,
    get_compressed_file_handle,
    is_compressed,
    read_fasta,
)
from pca.logutils import get_logger
from pca.parallel import async_parallel
from pca.parsers import (
    parse_blastp,
    parse_defensefinder_genes,
    parse_hmmsearch,
    parse_trna,
    parse_trna_gff,
)
from pca.pipeline import (
    call_genes_phanotate,
    call_genes_pyrodigal,
    search_extra_db,
    search_hmms_pyhmmer,
    search_phmmer_pyhmmer,
)
from pca.validation import (
    check_executables,
    dbname,
    is_valid_dir,
    is_valid_file_path,
)

__all__ = [
    # io
    "Compression",
    "check_fasta",
    "convert_to_fasta",
    "detect_sequence_format",
    "get_compressed_file_handle",
    "is_compressed",
    "read_fasta",
    # validation
    "check_executables",
    "dbname",
    "is_valid_dir",
    "is_valid_file_path",
    # logging
    "get_logger",
    # parallel
    "async_parallel",
    # externals
    "run_trnascan",
    # parsers
    "parse_blastp",
    "parse_defensefinder_genes",
    "parse_hmmsearch",
    "parse_trna",
    "parse_trna_gff",
    # annotations
    "convert_to_html",
    "create_feature",
    "generate_plots_and_annotations",
    "generate_plots_and_gff",
    "get_coordinates",
    "get_cordinates",
    # pipeline
    "call_genes_phanotate",
    "call_genes_pyrodigal",
    "search_extra_db",
    "search_hmms_pyhmmer",
    "search_phmmer_pyhmmer",
]
