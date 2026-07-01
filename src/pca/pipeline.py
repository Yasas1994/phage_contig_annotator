"""High-level pipeline steps (gene calling, HMM search, result aggregation)."""

from __future__ import annotations

import logging
import os
import tempfile
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import pandas as pd
import pyhmmer
import pyrodigal
import pyrodigal_gv
from BCBio import GFF
from Bio.Seq import Seq
from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.plan7 import HMM, HMMFile, TopHits

from pca.externals import run_phanotate, run_phanotate_rs
from pca.io import read_fasta
from pca.parsers import parse_hmmsearch

__all__ = [
    "call_genes_phanotate",
    "call_genes_phanotate_rs",
    "call_genes_pyrodigal",
    "search_extra_db",
    "search_hmms_pyhmmer",
    "search_phmmer_pyhmmer",
]

logger = logging.getLogger(__name__)


def _iter_hmms(db_path: Path) -> Iterator[HMM]:
    """Yield HMMs from a single ``.hmm`` file or every ``.hmm`` file under ``db_path``."""
    if db_path.is_file():
        with HMMFile(db_path) as hmm_file:
            yield from hmm_file
        return
    for path in sorted(db_path.rglob("*.hmm")):
        with HMMFile(path) as hmm_file:
            yield from hmm_file


def _write_hmmer_tblout(
    top_hits_iter: Iterator[TopHits],
    out_txt: Path,
) -> None:
    """Write HMMER tblout-style lines from an iterable of ``TopHits``.

    For ``hmmsearch`` the query is an HMM and the hit is a protein; we emit
    the protein (``hit``) first so that ``qname`` in downstream parsing is the
    input protein ID.  For ``phmmer`` the query is a protein and the hit is a
    database protein; we emit the query first for the same reason.
    """
    out_txt.parent.mkdir(parents=True, exist_ok=True)
    with open(out_txt, "w") as fh:
        for top_hits in top_hits_iter:
            query = top_hits.query
            query_name = query.name.decode() if isinstance(query.name, bytes) else query.name
            query_acc = (
                query.accession.decode()
                if isinstance(query.accession, bytes)
                else query.accession or "-"
            )
            is_hmm_query = isinstance(query, HMM)
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
                if is_hmm_query:
                    name1, acc1 = hit_name, hit_acc
                    name2, acc2 = query_name, query_acc
                else:
                    name1, acc1 = query_name, query_acc
                    name2, acc2 = hit_name, hit_acc
                fh.write(
                    f"{name1}\t{acc1}\t{name2}\t{acc2}\t"
                    f"{hit.evalue:.2g}\t{hit.score:.1f}\t{hit.bias:.1f}\t"
                    f"{domain.c_evalue:.2g}\t{domain.score:.1f}\t{domain.bias:.1f}\n"
                )


def search_hmms_pyhmmer(
    in_faa: str | os.PathLike[str],
    db_path: str | os.PathLike[str],
    out_txt: str | os.PathLike[str],
    out_csv: str | os.PathLike[str],
    threads: int = 1,
    evalue: float = 10.0,
) -> None:
    """Search proteins against all HMMs in ``db_path`` using pyhmmer.

    ``db_path`` may be a directory containing ``.hmm`` files, a single ``.hmm``
    file, or a pressed HMM database. Writes both an HMMER-style tblout file
    and a CSV summary.
    """
    db_path = Path(db_path)
    out_txt = Path(out_txt)
    out_csv = Path(out_csv)

    logger.info("loading HMMs from %s", db_path)
    hmms = list(_iter_hmms(db_path))
    if not hmms:
        raise RuntimeError(f"No .hmm files found in {db_path}")

    logger.info("loading protein sequences from %s", in_faa)
    alphabet = Alphabet.amino()
    with SequenceFile(in_faa, digital=True, alphabet=alphabet) as seq_file:
        sequences = seq_file.read_block()

    logger.info("running hmmsearch against %d HMMs", len(hmms))
    cpus = threads if threads > 0 else None
    all_hits = pyhmmer.hmmsearch(hmms, sequences, cpus=cpus, E=evalue)

    _write_hmmer_tblout(all_hits, out_txt)

    logger.info("writing hmmsearch CSV summary")
    search_results = pd.DataFrame(parse_hmmsearch(out_txt))
    search_results.to_csv(out_csv, index=False)


def search_phmmer_pyhmmer(
    in_faa: str | os.PathLike[str],
    db_fasta: str | os.PathLike[str],
    out_txt: str | os.PathLike[str],
    out_csv: str | os.PathLike[str],
    threads: int = 1,
    evalue: float = 1e-5,
) -> None:
    """Search proteins against a protein sequence database using pyhmmer/phmmer.

    ``db_fasta`` should be a multi-FASTA file of protein reference sequences
    (e.g. CARD, VFDB, Anti-CRISPRdb). Writes both an HMMER-style tblout file
    and a CSV summary.
    """
    db_fasta = Path(db_fasta)
    out_txt = Path(out_txt)
    out_csv = Path(out_csv)

    if not db_fasta.is_file():
        raise RuntimeError(f"Protein database FASTA not found: {db_fasta}")

    logger.info("loading query proteins from %s", in_faa)
    alphabet = Alphabet.amino()
    with SequenceFile(in_faa, digital=True, alphabet=alphabet) as seq_file:
        queries = seq_file.read_block()

    logger.info("running phmmer against %s", db_fasta)
    cpus = threads if threads > 0 else None
    with SequenceFile(db_fasta, digital=True, alphabet=alphabet) as targets:
        all_hits = pyhmmer.phmmer(queries, targets, cpus=cpus, E=evalue)
        _write_hmmer_tblout(all_hits, out_txt)

    logger.info("writing phmmer CSV summary")
    search_results = pd.DataFrame(parse_hmmsearch(out_txt))
    search_results.to_csv(out_csv, index=False)


def search_extra_db(
    in_faa: str | os.PathLike[str],
    db_path: str | os.PathLike[str],
    out_txt: str | os.PathLike[str],
    out_csv: str | os.PathLike[str],
    threads: int = 1,
    evalue: float = 1e-5,
) -> None:
    """Dispatch an extra database search to ``hmmsearch`` or ``phmmer``.

    If ``db_path`` is a directory, ``.hmm`` files take precedence over protein
    FASTA files. If ``db_path`` is a file, its extension selects the program.
    """
    db_path = Path(db_path)
    fasta_suffixes = {".faa", ".fasta", ".fa"}

    if db_path.is_dir():
        hmm_files = list(db_path.glob("*.hmm"))
        fasta_files = [
            p for p in db_path.iterdir() if p.suffix.lower() in fasta_suffixes
        ]
        if hmm_files:
            search_hmms_pyhmmer(in_faa, db_path, out_txt, out_csv, threads, evalue)
        elif fasta_files:
            search_phmmer_pyhmmer(in_faa, fasta_files[0], out_txt, out_csv, threads, evalue)
        else:
            raise RuntimeError(
                f"Extra database directory {db_path} contains no .hmm or protein FASTA files"
            )
    elif db_path.suffix.lower() == ".hmm":
        search_hmms_pyhmmer(in_faa, db_path, out_txt, out_csv, threads, evalue)
    elif db_path.suffix.lower() in fasta_suffixes:
        search_phmmer_pyhmmer(in_faa, db_path, out_txt, out_csv, threads, evalue)
    else:
        raise RuntimeError(
            f"Unsupported extra database path {db_path}: expected a .hmm file or protein FASTA"
        )


def call_genes_pyrodigal(
    in_fna: str | os.PathLike[str],
    out_faa: str | os.PathLike[str],
    out_gff: str | os.PathLike[str],
    threads: int = 1,
    translation_table: int | None = None,
) -> None:
    """Predict genes with pyrodigal-gv and write proteins + GFF.

    When ``translation_table`` is given, run in single-genome mode with that
    table. Otherwise use pyrodigal-gv's default metagenomic mode, which
    automatically selects among viral models.
    """
    in_fna = Path(in_fna)
    out_faa = Path(out_faa)
    out_gff = Path(out_gff)
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    out_gff.parent.mkdir(parents=True, exist_ok=True)

    logger.info("gene calling with pyrodigal-gv: %s", in_fna)
    if translation_table is None:
        orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True)
    else:
        training_info = pyrodigal.TrainingInfo(
            gc=0.5,
            translation_table=translation_table,
        )
        orf_finder = pyrodigal_gv.ViralGeneFinder(
            training_info=training_info,
            meta=False,
        )

    with open(out_faa, "w") as faa_fh, open(out_gff, "w") as gff_fh:
        for seq_id, seq in read_fasta(in_fna):
            pred = orf_finder.find_genes(seq.encode())
            pred.write_translations(faa_fh, sequence_id=seq_id)
            pred.write_gff(gff_fh, sequence_id=seq_id, header=False)

    logger.info("wrote %s and %s", out_faa, out_gff)


def _normalize_phanotate_gff_coordinates(
    gff_path: str | os.PathLike[str], out_path: str | os.PathLike[str]
) -> None:
    """Normalize PHANOTATE/PHANOTATE-rs GFF coordinates for BCBio.GFF.

    PHANOTATE-rs (and some PHANOTATE versions) emit negative-strand CDS
    features with the start column greater than the end column. GFF requires
    start <= end, with the strand sign carrying the orientation. This helper
    swaps those coordinates in place before parsing.
    """
    with open(gff_path) as in_fh, open(out_path, "w") as out_fh:
        for line in in_fh:
            if line.startswith("#") or not line.strip():
                out_fh.write(line)
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5:
                try:
                    start = int(parts[3])
                    end = int(parts[4])
                except ValueError:
                    pass
                else:
                    if start > end:
                        parts[3], parts[4] = str(end), str(start)
            out_fh.write("\t".join(parts) + "\n")


def _parse_phanotate_gff(
    gff_path: str | os.PathLike[str],
    sequences: dict[str, Seq],
    translation_table: int = 11,
) -> list[Any]:
    """Parse a PHANOTATE GFF3 and return normalized SeqRecords.

    CDS features are enumerated in file order and assigned deterministic
    IDs (``<contig>_<position>``) so that the protein FASTA headers match
    the GFF features consumed downstream.
    """
    records: list[Any] = []
    with open(gff_path) as handle:
        for record in GFF.parse(handle):
            seq = sequences.get(record.id)
            if seq is not None:
                record.seq = seq

            cds_count = 0
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                cds_count += 1
                qname = f"{record.id}_{cds_count}"
                old_id = feature.qualifiers.get("ID", [None])[0]
                feature.qualifiers["ID"] = [qname]
                if old_id and "Name" not in feature.qualifiers:
                    feature.qualifiers["Name"] = [old_id]

                if seq is not None:
                    try:
                        nuc = feature.extract(record.seq)
                        # Trim incomplete trailing codons before translating.
                        trim = len(nuc) % 3
                        if trim:
                            nuc = nuc[:-trim]
                        protein = str(nuc.translate(table=translation_table, to_stop=True))
                        logger.warning("DEBUG translated %s len=%d protein=%r", qname, len(nuc), protein)
                    except Exception as exc:
                        logger.warning("failed to translate %s: %s", qname, exc)
                        protein = ""
                    feature.qualifiers["translation"] = [protein]

            records.append(record)

    return records


def call_genes_phanotate_rs(
    in_fna: str | os.PathLike[str],
    out_faa: str | os.PathLike[str],
    out_gff: str | os.PathLike[str],
    threads: int = 1,
    translation_table: int | None = None,
) -> None:
    """Predict genes with PHANOTATE-rs and write proteins + GFF.

    PHANOTATE-rs is a Rust reimplementation of PHANOTATE. The wrapper
    requests GFF3 output, parses it, and translates the proteins using
    the selected NCBI translation table.
    """
    in_fna = Path(in_fna)
    out_faa = Path(out_faa)
    out_gff = Path(out_gff)
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    out_gff.parent.mkdir(parents=True, exist_ok=True)

    table = translation_table if translation_table is not None else 11

    logger.info("gene calling with PHANOTATE-rs: %s", in_fna)
    with tempfile.TemporaryDirectory(prefix="phanotate_rs_") as tmp:
        tmp_path = Path(tmp)
        raw_gff = tmp_path / "phanotate_rs.gff"
        normalized_gff = tmp_path / "phanotate_rs_normalized.gff"
        success = run_phanotate_rs(
            str(in_fna), str(raw_gff), threads=threads, translation_table=table
        )
        if not success or not raw_gff.is_file():
            raise RuntimeError(f"PHANOTATE-rs gene calling failed for {in_fna}")

        _normalize_phanotate_gff_coordinates(raw_gff, normalized_gff)
        sequences = {seq_id: Seq(seq) for seq_id, seq in read_fasta(in_fna)}
        records = _parse_phanotate_gff(normalized_gff, sequences, translation_table=table)

    with open(out_faa, "w") as faa_fh:
        for record in records:
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                qname = feature.qualifiers.get("ID", ["unknown"])[0]
                protein = feature.qualifiers.get("translation", [""])[0]
                faa_fh.write(f">{qname}\n{protein}\n")

    with open(out_gff, "w") as gff_fh:
        GFF.write(records, gff_fh)

    logger.info("wrote %s and %s", out_faa, out_gff)


def call_genes_phanotate(
    in_fna: str | os.PathLike[str],
    out_faa: str | os.PathLike[str],
    out_gff: str | os.PathLike[str],
    threads: int = 1,
    translation_table: int | None = None,
) -> None:
    """Predict genes with PHANOTATE and write proteins + GFF.

    ``translation_table`` is accepted for API compatibility but defaults to
    NCBI table 11, which PHANOTATE itself does not allow changing.
    """
    in_fna = Path(in_fna)
    out_faa = Path(out_faa)
    out_gff = Path(out_gff)
    out_faa.parent.mkdir(parents=True, exist_ok=True)
    out_gff.parent.mkdir(parents=True, exist_ok=True)

    table = translation_table if translation_table is not None else 11

    logger.info("gene calling with PHANOTATE: %s", in_fna)
    with tempfile.TemporaryDirectory(prefix="phanotate_") as tmp:
        tmp_path = Path(tmp)
        raw_gff = tmp_path / "phanotate.gff"
        normalized_gff = tmp_path / "phanotate_normalized.gff"
        success = run_phanotate(str(in_fna), str(raw_gff), threads=threads)
        if not success or not raw_gff.is_file():
            raise RuntimeError(f"PHANOTATE gene calling failed for {in_fna}")

        _normalize_phanotate_gff_coordinates(raw_gff, normalized_gff)
        sequences = {seq_id: Seq(seq) for seq_id, seq in read_fasta(in_fna)}
        records = _parse_phanotate_gff(normalized_gff, sequences, translation_table=table)

    with open(out_faa, "w") as faa_fh:
        for record in records:
            for feature in record.features:
                if feature.type != "CDS":
                    continue
                qname = feature.qualifiers.get("ID", ["unknown"])[0]
                protein = feature.qualifiers.get("translation", [""])[0]
                faa_fh.write(f">{qname}\n{protein}\n")

    with open(out_gff, "w") as gff_fh:
        GFF.write(records, gff_fh)

    logger.info("wrote %s and %s", out_faa, out_gff)
