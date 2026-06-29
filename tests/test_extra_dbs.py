"""Unit tests for optional extra database searches (CARD, VFDB, etc.)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import Background, Builder

from pca.annotations import generate_plots_and_annotations
from pca.parsers import parse_defensefinder_genes, parse_hmmsearch
from pca.pipeline import search_extra_db


_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def _digital_sequence(name: bytes, sequence: str) -> TextSequence:
    alphabet = Alphabet.amino()
    return TextSequence(name=name, sequence=sequence.encode()).digitize(alphabet)


def _build_hmm(sequence: str, name: bytes = b"test_hmm") -> Path:
    """Return the path to a tiny HMM built from a single sequence."""
    alphabet = Alphabet.amino()
    seq = _digital_sequence(name, sequence)
    builder = Builder(alphabet)
    bg = Background(alphabet)
    hmm, _, _ = builder.build(seq, bg)
    hmm.name = name
    return hmm


def _write_hmm(path: Path, sequence: str, name: bytes = b"test_hmm") -> None:
    hmm = _build_hmm(sequence, name)
    with open(path, "wb") as fh:
        hmm.write(fh)


def _write_fasta(path: Path, entries: dict[bytes, str]) -> None:
    with open(path, "w") as fh:
        for name, seq in entries.items():
            fh.write(f">{name.decode()}\n{seq}\n")


def test_search_extra_db_with_hmm_profile(tmp_path: Path) -> None:
    """``search_extra_db`` dispatches to hmmsearch when given ``.hmm`` files."""
    seq = "MKTLLILFRKAV" + _AMINO_ACIDS
    query_faa = tmp_path / "queries.faa"
    _write_fasta(query_faa, {b"contig1_1": seq, b"contig1_2": _AMINO_ACIDS})

    db_dir = tmp_path / "db"
    db_dir.mkdir()
    _write_hmm(db_dir / "test.hmm", seq)

    out_txt = tmp_path / "out.txt"
    out_csv = tmp_path / "out.csv"
    search_extra_db(query_faa, db_dir, out_txt, out_csv, threads=1, evalue=10.0)

    assert out_txt.exists()
    assert out_csv.exists()
    rows = list(parse_hmmsearch(out_txt))
    assert any(row["qname"] == "contig1_1" and row["tname"] == "test_hmm" for row in rows)
    df = pd.read_csv(out_csv)
    assert not df.empty


def test_search_extra_db_with_protein_fasta(tmp_path: Path) -> None:
    """``search_extra_db`` dispatches to phmmer when given a protein FASTA."""
    seq = "MKTLLILFRKAV" + _AMINO_ACIDS
    query_faa = tmp_path / "queries.faa"
    _write_fasta(query_faa, {b"contig1_1": seq})

    db_dir = tmp_path / "db"
    db_dir.mkdir()
    _write_fasta(db_dir / "ref.faa", {b"card_target": seq})

    out_txt = tmp_path / "out.txt"
    out_csv = tmp_path / "out.csv"
    search_extra_db(query_faa, db_dir, out_txt, out_csv, threads=1, evalue=10.0)

    rows = list(parse_hmmsearch(out_txt))
    assert any(row["qname"] == "contig1_1" and row["tname"] == "card_target" for row in rows)


def test_parse_defensefinder_genes(tmp_path: Path) -> None:
    """Parse DefenseFinder gene output into simple dict records."""
    tsv = tmp_path / "defense_finder_genes.tsv"
    tsv.write_text(
        "replicon\thit_id\tgene_name\tsys_id\tsubtype\n"
        "contig1\tcontig1_1\tCAS_Class1-Subtype-I-E\tsys_1\tCAS\n"
        "contig1\tcontig1_2\tCAS_Class1-Subtype-I-E\tsys_1\tCAS\n"
    )
    rows = list(parse_defensefinder_genes(tsv))
    assert len(rows) == 2
    assert rows[0] == {
        "qname": "contig1_1",
        "gene_name": "CAS_Class1-Subtype-I-E",
        "system": "sys_1",
        "type": "CAS",
    }


def test_generate_plots_and_annotations_with_extra_hits(tmp_path: Path) -> None:
    """Extra database hits are attached as CDS qualifiers and written to outputs."""
    tmp_path = Path(tmp_path)
    sequences = {"contig1": Seq("A" * 300)}
    fasta = tmp_path / "input.fna"
    with open(fasta, "w") as fh:
        for seq_id, seq in sequences.items():
            fh.write(f">{seq_id}\n{str(seq)}\n")

    gff = tmp_path / "genes.gff"
    gff.write_text(
        "##gff-version 3\n"
        "contig1\tprodigal-gv\tCDS\t1\t300\t.\t+\t0\tID=1_1\n"
    )

    meta = tmp_path / "meta.csv"
    meta.write_text("#phrog,Annotation,Category,color\nphrog_1,test protein,other,#bcbd22\n")

    hmmsearch_txt = tmp_path / "hmmsearch.txt"
    hmmsearch_txt.write_text(
        "contig1_1\t-\tphrog_1\t-\t1e-50\t200.0\t0.0\t1e-50\t200.0\t0.0\n"
    )

    extra_csv = tmp_path / "card.csv"
    extra_csv.write_text(
        "qname,qacc,tname,tacc,eval,score,bias,beval,bscore,bbias\n"
        "contig1_1,-,TetM,-,1e-20,150.0,0.0,1e-20,150.0,0.0\n"
    )

    generate_plots_and_annotations(
        tmp_dir=tmp_path,
        hmmsearch_dir=hmmsearch_txt,
        trna_dir=tmp_path / "empty_trna.gff",
        meta_dir=meta,
        gff_dir=gff,
        input_fasta=fasta,
        plot_formats=[],
        extra_hits={"card": extra_csv},
    )

    gff_out = tmp_path / "annotations" / "annotations.gff"
    gbk_out = tmp_path / "annotations" / "annotations.gbk"
    assert gff_out.exists()
    assert gbk_out.exists()

    gff_text = gff_out.read_text()
    assert "card_hit=TetM" in gff_text or "card_hit" in gff_text

    record = next(SeqIO.parse(gbk_out, "genbank"))
    cds = [f for f in record.features if f.type == "CDS"][0]
    assert cds.qualifiers.get("card_hit") == ["TetM"]
    assert float(cds.qualifiers["card_score"][0]) == 150.0
    assert cds.qualifiers.get("category") == ["other"]


def test_generate_plots_and_annotations_sets_category_from_extra_hit(tmp_path: Path) -> None:
    """A CDS with no PHROG hit but an extra DB hit gets the extra category."""
    fasta = tmp_path / "input.fna"
    fasta.write_text(">contig1\n" + "A" * 300 + "\n")
    gff = tmp_path / "genes.gff"
    gff.write_text(
        "##gff-version 3\n"
        "contig1\tprodigal-gv\tCDS\t1\t300\t.\t+\t0\tID=1_1\n"
    )
    meta = tmp_path / "meta.csv"
    meta.write_text("#phrog,Annotation,Category,color\nphrog_1,test protein,other,#bcbd22\n")

    hmmsearch_txt = tmp_path / "hmmsearch.txt"
    hmmsearch_txt.write_text("")
    # Function raises on empty PHROG hits, so we give a dummy hit with score <=50
    # so it is filtered out and the CDS remains unannotated by PHROG.
    hmmsearch_txt.write_text(
        "contig1_1\t-\tphrog_1\t-\t1e-2\t10.0\t0.0\t1e-2\t10.0\t0.0\n"
    )

    extra_csv = tmp_path / "vfdb.csv"
    extra_csv.write_text(
        "qname,qacc,tname,tacc,eval,score,bias,beval,bscore,bbias\n"
        "contig1_1,-,IcaA,-,1e-10,80.0,0.0,1e-10,80.0,0.0\n"
    )

    generate_plots_and_annotations(
        tmp_dir=tmp_path,
        hmmsearch_dir=hmmsearch_txt,
        trna_dir=tmp_path / "empty_trna.gff",
        meta_dir=meta,
        gff_dir=gff,
        input_fasta=fasta,
        plot_formats=[],
        extra_hits={"vfdb": extra_csv},
    )

    gbk_out = tmp_path / "annotations" / "annotations.gbk"
    record = next(SeqIO.parse(gbk_out, "genbank"))
    cds = [f for f in record.features if f.type == "CDS"][0]
    assert cds.qualifiers.get("vfdb_hit") == ["IcaA"]
    assert cds.qualifiers.get("category") == ["virulence factor"]
