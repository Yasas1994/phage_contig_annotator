"""Unit tests for pca.io splitting helpers."""

from __future__ import annotations

from pathlib import Path

from BCBio import GFF
from Bio import SeqIO

from pca import io


class TestSplitFastaByContig:
    def test_split_protein_fasta_by_contig(self, tmp_path: Path) -> None:
        fasta = tmp_path / "proteins.faa"
        fasta.write_text(
            ">contig1_1\nMKT\n>contig1_2\nMLR\n>contig2_1\nMST\n>contig2_2\nMVS\n"
        )
        paths = io.split_fasta_by_contig(fasta, tmp_path / "out", filename="proteins.faa")

        assert set(paths.keys()) == {"contig1", "contig2"}
        assert (tmp_path / "out" / "contig1" / "proteins.faa").read_text() == (
            ">contig1_1\nMKT\n>contig1_2\nMLR\n"
        )
        assert (tmp_path / "out" / "contig2" / "proteins.faa").read_text() == (
            ">contig2_1\nMST\n>contig2_2\nMVS\n"
        )

    def test_split_with_underscore_contig_id(self, tmp_path: Path) -> None:
        fasta = tmp_path / "proteins.faa"
        fasta.write_text(
            ">my_contig_1_1\nMKT\n>my_contig_1_2\nMLR\n>my_contig_2_1\nMST\n"
        )
        paths = io.split_fasta_by_contig(fasta, tmp_path / "out", filename="proteins.faa")

        assert set(paths.keys()) == {"my_contig_1", "my_contig_2"}
        assert ">my_contig_1_1" in (
            tmp_path / "out" / "my_contig_1" / "proteins.faa"
        ).read_text()
        assert ">my_contig_2_1" in (
            tmp_path / "out" / "my_contig_2" / "proteins.faa"
        ).read_text()


class TestSplitGffByContig:
    def test_split_gff_by_contig(self, tmp_path: Path) -> None:
        gff = tmp_path / "proteins.gff"
        gff.write_text(
            "##gff-version 3\n"
            "##sequence-region contig1 1 100\n"
            "contig1\tphage\tCDS\t1\t30\t.\t+\t0\tID=1\n"
            "##sequence-region contig2 1 100\n"
            "contig2\tphage\tCDS\t1\t30\t.\t+\t0\tID=2\n"
        )
        paths = io.split_gff_by_contig(gff, tmp_path / "out", filename="proteins.gff")

        assert set(paths.keys()) == {"contig1", "contig2"}
        for contig in ("contig1", "contig2"):
            out_path = tmp_path / "out" / contig / "proteins.gff"
            assert out_path.is_file()
            records = list(GFF.parse(out_path))
            assert len(records) == 1
            assert records[0].id == contig


class TestSplitHmmsearchByContig:
    def test_split_hmmsearch_by_contig(self, tmp_path: Path) -> None:
        txt = tmp_path / "hmmsearch.txt"
        txt.write_text(
            "# hmmsearch :: search for hidden markov models\n"
            "contig1_1\t-\tPHROG1\t-\t1e-50\t100.0\t0.0\t1e-50\t100.0\t0.0\n"
            "contig1_2\t-\tPHROG2\t-\t2e-40\t90.0\t0.0\t2e-40\t90.0\t0.0\n"
            "contig2_1\t-\tPHROG3\t-\t3e-30\t80.0\t0.0\t3e-30\t80.0\t0.0\n"
        )
        paths = io.split_hmmsearch_by_contig(txt, tmp_path / "out", filename="hmmsearch.txt")

        assert set(paths.keys()) == {"contig1", "contig2"}
        c1 = (tmp_path / "out" / "contig1" / "hmmsearch.txt").read_text()
        c2 = (tmp_path / "out" / "contig2" / "hmmsearch.txt").read_text()
        assert "contig1_1" in c1
        assert "contig1_2" in c1
        assert "contig2_1" not in c1
        assert "contig2_1" in c2
        assert "#" not in c1


class TestSplitCsvByContig:
    def test_split_csv_by_contig_column(self, tmp_path: Path) -> None:
        csv = tmp_path / "table.csv"
        csv.write_text("contig,score\ncontig1,10\ncontig1,20\ncontig2,30\n")
        paths = io.split_csv_by_contig(csv, tmp_path / "out", filename="table.csv")

        assert set(paths.keys()) == {"contig1", "contig2"}
        c1 = (tmp_path / "out" / "contig1" / "table.csv").read_text()
        c2 = (tmp_path / "out" / "contig2" / "table.csv").read_text()
        assert c1.count("contig1") == 2
        assert "contig2" not in c1
        assert "contig2" in c2

    def test_split_csv_by_qname_column(self, tmp_path: Path) -> None:
        csv = tmp_path / "table.csv"
        csv.write_text("qname,tname,score\ncontig1_1,PHROG1,10\ncontig1_2,PHROG2,20\ncontig2_1,PHROG3,30\n")
        paths = io.split_csv_by_contig(
            csv, tmp_path / "out", filename="table.csv", contig_column="qname"
        )

        assert set(paths.keys()) == {"contig1", "contig2"}
        c1 = (tmp_path / "out" / "contig1" / "table.csv").read_text()
        assert "contig1_1" in c1
        assert "contig1_2" in c1
        assert "contig2_1" not in c1


class TestSplitTsvByContig:
    def test_split_tsv_by_contig(self, tmp_path: Path) -> None:
        tsv = tmp_path / "table.tsv"
        tsv.write_text("qname\tscore\ncontig1_1\t10\ncontig2_1\t20\n")
        paths = io.split_tsv_by_contig(tsv, tmp_path / "out", filename="table.tsv")

        assert set(paths.keys()) == {"contig1", "contig2"}
        c1 = (tmp_path / "out" / "contig1" / "table.tsv").read_text()
        assert "contig1_1" in c1
        assert "contig2_1" not in c1


class TestSplitTrnaTsvByContig:
    def test_split_trna_tsv_by_contig(self, tmp_path: Path) -> None:
        tsv = tmp_path / "trna.tsv"
        tsv.write_text(
            "Sequence\t\ttRNA\tBounds\ttRNA\tAnti\tIntron Bounds\tInf\n"
            "Name\t\ttRNA #\tBegin\tEnd\tType\tCodon\tBegin\tEnd\tScore\tNote\n"
            "--------\t\t------\t-----\t------\t----\t-----\t-----\t----\t------\t------\n"
            "contig1\t1\t10\t80\tGln\tTTG\t0\t0\t62.4\n"
            "contig1\t2\t90\t160\tMet\tCAT\t0\t0\t68.6\n"
            "contig2\t1\t10\t80\tArg\tACG\t0\t0\t55.0\n"
        )
        paths = io.split_trna_tsv_by_contig(tsv, tmp_path / "out", filename="trna.tsv")

        assert set(paths.keys()) == {"contig1", "contig2"}
        c1 = (tmp_path / "out" / "contig1" / "trna.tsv").read_text()
        c2 = (tmp_path / "out" / "contig2" / "trna.tsv").read_text()
        assert "contig1" in c1
        assert "contig2" not in c1
        assert "contig2" in c2
        assert "--------" in c1
        assert "--------" in c2


class TestSplitTrfDatByContig:
    def test_split_trf_dat_by_contig(self, tmp_path: Path) -> None:
        dat = tmp_path / "trf.dat"
        dat.write_text(
            "Tandem Repeats Finder Program\n\n"
            "Sequence: contig1\n"
            "Parameters: 2 7 7 80 10 50 2000\n"
            "1 10 4 2.0 4 10 0 20 5 5 5 5 2.00 ATGC ATGC\n"
            "Sequence: contig2\n"
            "Parameters: 2 7 7 80 10 50 2000\n"
            "11 20 4 2.0 4 10 0 20 5 5 5 5 2.00 CGTA CGTA\n"
        )
        paths = io.split_trf_dat_by_contig(dat, tmp_path / "out", filename="trf.dat")

        assert set(paths.keys()) == {"contig1", "contig2"}
        c1 = (tmp_path / "out" / "contig1" / "trf.dat").read_text()
        c2 = (tmp_path / "out" / "contig2" / "trf.dat").read_text()
        assert "Sequence: contig1" in c1
        assert "1 10" in c1
        assert "Sequence: contig2" not in c1
        assert "Sequence: contig2" in c2
        assert "11 20" in c2


class TestSplitGenbankByContig:
    def test_split_genbank_by_contig(self, tmp_path: Path) -> None:
        gbk = tmp_path / "annotations.gbk"
        gbk.write_text(
            "LOCUS       contig1        4 bp    DNA     linear       01-JAN-2025\n"
            "ORIGIN\n"
            "        1 acgt\n"
            "//\n"
            "LOCUS       contig2        4 bp    DNA     linear       01-JAN-2025\n"
            "ORIGIN\n"
            "        1 tgca\n"
            "//\n"
        )
        paths = io.split_genbank_by_contig(gbk, tmp_path / "out", filename="annotations.gbk")

        assert set(paths.keys()) == {"contig1", "contig2"}
        for contig in ("contig1", "contig2"):
            out_path = tmp_path / "out" / contig / "annotations.gbk"
            assert out_path.is_file()
            records = list(SeqIO.parse(out_path, "genbank"))
            assert len(records) == 1
            assert records[0].id == contig
