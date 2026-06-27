"""High-level pipeline steps (gene calling, HMM search, result aggregation)."""

from __future__ import annotations

import logging
import os
from collections.abc import Iterator
from pathlib import Path

import pandas as pd
import pyhmmer
import pyrodigal_gv
from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.plan7 import HMM, HMMFile

from pca.io import read_fasta
from pca.parsers import parse_hmmsearch

__all__ = ["call_genes_pyrodigal", "search_hmms_pyhmmer"]

logger = logging.getLogger(__name__)


def _iter_hmms(db_dir: Path) -> Iterator[HMM]:
    """Yield HMMs from every ``.hmm`` file in ``db_dir``."""
    for path in sorted(db_dir.glob("*.hmm")):
        with HMMFile(path) as hmm_file:
            yield from hmm_file


def search_hmms_pyhmmer(
    in_faa: str | os.PathLike[str],
    db_dir: str | os.PathLike[str],
    out_txt: str | os.PathLike[str],
    out_csv: str | os.PathLike[str],
    threads: int = 1,
    evalue: float = 10.0,
) -> None:
    """Search proteins against all HMMs in ``db_dir`` using pyhmmer.

    Writes both an HMMER-style tblout file and a CSV summary.
    """
    db_dir = Path(db_dir)
    out_txt = Path(out_txt)
    out_csv = Path(out_csv)

    logger.info("loading HMMs from %s", db_dir)
    hmms = list(_iter_hmms(db_dir))
    if not hmms:
        raise RuntimeError(f"No .hmm files found in {db_dir}")

    logger.info("loading protein sequences from %s", in_faa)
    alphabet = Alphabet.amino()
    with SequenceFile(in_faa, digital=True, alphabet=alphabet) as seq_file:
        sequences = seq_file.read_block()

    logger.info("running hmmsearch against %d HMMs", len(hmms))
    cpus = threads if threads > 0 else None
    all_hits = pyhmmer.hmmsearch(hmms, sequences, cpus=cpus)

    out_txt.parent.mkdir(parents=True, exist_ok=True)
    with open(out_txt, "w") as fh:
        for top_hits in all_hits:
            query = top_hits.query
            query_name = query.name.decode() if isinstance(query.name, bytes) else query.name
            query_acc = (
                query.accession.decode()
                if isinstance(query.accession, bytes)
                else query.accession or "-"
            )
            for hit in top_hits:
                if not hit.domains:
                    continue
                domain = hit.domains[0]
                hit_name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                hit_acc = (
                    hit.accession.decode()
                    if isinstance(hit.accession, bytes)
                    else hit.accession or "-"
                )
                fh.write(
                    f"{hit_name}\t{hit_acc}\t{query_name}\t{query_acc}\t"
                    f"{hit.evalue:.2g}\t{hit.score:.1f}\t{hit.bias:.1f}\t"
                    f"{domain.c_evalue:.2g}\t{domain.score:.1f}\t{domain.bias:.1f}\n"
                )

    logger.info("writing hmmsearch CSV summary")
    search_results = pd.DataFrame(parse_hmmsearch(out_txt))
    search_results.to_csv(out_csv, index=False)


def call_genes_pyrodigal(
    in_fna: str | os.PathLike[str],
    out_faa: str | os.PathLike[str],
    out_gff: str | os.PathLike[str],
    threads: int = 1,
) -> None:
    """Predict genes with pyrodigal-gv and write proteins + GFF."""
    in_fna = Path(in_fna)
    out_faa = Path(out_faa)
    out_gff = Path(out_gff)
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    out_gff.parent.mkdir(parents=True, exist_ok=True)

    logger.info("gene calling with pyrodigal-gv: %s", in_fna)
    orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True)

    with open(out_faa, "w") as faa_fh, open(out_gff, "w") as gff_fh:
        for seq_id, seq in read_fasta(in_fna):
            pred = orf_finder.find_genes(seq.encode())
            pred.write_translations(faa_fh, sequence_id=seq_id)
            pred.write_gff(gff_fh, sequence_id=seq_id, header=False)

    logger.info("wrote %s and %s", out_faa, out_gff)
