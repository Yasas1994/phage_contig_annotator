"""Unit tests for pca.utils and submodules."""

import argparse
import gzip
import logging
from pathlib import Path

import pytest

from pca import (
    annotations,
    externals,
    io,
    logutils,
    parsers,
    pipeline,
    utils,
    validation,
)


class TestModuleFacades:
    def test_io_symbols_re_exported(self) -> None:
        assert utils.Compression is io.Compression
        assert utils.read_fasta is io.read_fasta
        assert utils.is_compressed is io.is_compressed

    def test_validation_symbols_re_exported(self) -> None:
        assert utils.is_valid_dir is validation.is_valid_dir
        assert utils.check_executables is validation.check_executables

    def test_pipeline_symbols_re_exported(self) -> None:
        assert utils.call_genes_pyrodigal is pipeline.call_genes_pyrodigal
        assert utils.search_hmms_pyhmmer is pipeline.search_hmms_pyhmmer

    def test_annotations_symbols_re_exported(self) -> None:
        assert utils.get_coordinates is annotations.get_coordinates
        assert utils.get_cordinates is annotations.get_cordinates

    def test_externals_symbols_re_exported(self) -> None:
        assert utils.run_trnascan is externals.run_trnascan

    def test_parsers_symbols_re_exported(self) -> None:
        assert utils.parse_hmmsearch is parsers.parse_hmmsearch
        assert utils.parse_trna_gff is parsers.parse_trna_gff

    def test_logger_symbol_re_exported(self) -> None:
        assert utils.get_logger is logutils.get_logger


class TestLogger:
    def test_get_logger_colors_level_name_when_tty(self, monkeypatch) -> None:
        monkeypatch.setattr("sys.stderr.isatty", lambda: True)
        log = logutils.get_logger()
        handler = log.handlers[0]
        record = logging.LogRecord(
            name="pca",
            level=logging.ERROR,
            pathname="",
            lineno=0,
            msg="test",
            args=(),
            exc_info=None,
        )
        formatted = handler.format(record)
        assert "\033[31mERROR\033[0m" in formatted

    def test_get_logger_does_not_color_when_not_tty(self, monkeypatch) -> None:
        monkeypatch.setattr("sys.stderr.isatty", lambda: False)
        log = logutils.get_logger()
        handler = log.handlers[0]
        record = logging.LogRecord(
            name="pca",
            level=logging.ERROR,
            pathname="",
            lineno=0,
            msg="test",
            args=(),
            exc_info=None,
        )
        formatted = handler.format(record)
        assert "ERROR" in formatted
        assert "\033[" not in formatted


class TestCompression:
    def test_is_compressed_plain_text(self, tmp_path: Path) -> None:
        plain = tmp_path / "plain.fna"
        plain.write_text(">seq1\nACGT\n")
        assert utils.is_compressed(plain) == utils.Compression.noncompressed

    def test_is_compressed_gzip(self, tmp_path: Path) -> None:
        gz = tmp_path / "plain.fna.gz"
        with gzip.open(gz, "wt") as fh:
            fh.write(">seq1\nACGT\n")
        assert utils.is_compressed(gz) == utils.Compression.gzip


class TestReadFasta:
    def test_read_plain_fasta(self) -> None:
        records = list(utils.read_fasta("test/bin.460.fna"))
        assert len(records) == 9
        header, seq = records[0]
        assert header.startswith("NODE_94_length_44776_cov_27.159388")
        assert set(seq).issubset({"A", "C", "G", "T", "N"})


class TestValidationHelpers:
    def test_is_valid_dir_creates_missing_directory(self, tmp_path: Path) -> None:
        target = tmp_path / "new_dir"
        result = utils.is_valid_dir(str(target))
        assert result == str(target)
        assert target.is_dir()

    def test_is_valid_dir_rejects_existing_file(self, tmp_path: Path) -> None:
        existing = tmp_path / "existing"
        existing.write_text("not a dir")
        with pytest.raises(argparse.ArgumentTypeError):
            utils.is_valid_dir(str(existing))

    def test_is_valid_file_path_accepts_file(self, tmp_path: Path) -> None:
        existing = tmp_path / "file.fna"
        existing.write_text(">seq\nACGT")
        assert utils.is_valid_file_path(str(existing)) == str(existing)

    def test_is_valid_file_path_rejects_missing(self, tmp_path: Path) -> None:
        with pytest.raises(argparse.ArgumentTypeError):
            utils.is_valid_file_path(str(tmp_path / "missing.fna"))


class TestDbName:
    def test_dbname_extracts_prefix(self, tmp_path: Path) -> None:
        (tmp_path / "phrogs.dmnd").write_text("dummy")
        assert utils.dbname(str(tmp_path)) == "phrogs"


class TestParsers:
    def test_parse_hmmsearch_skips_comments(self, tmp_path: Path) -> None:
        hmmout = tmp_path / "out.tbl"
        hmmout.write_text(
            "# header line\n"
            "q1 - phrog_1 - 1e-10 50.0 0.1 1e-5 30.0 0.0\n"
        )
        rows = list(utils.parse_hmmsearch(str(hmmout)))
        assert len(rows) == 1
        assert rows[0]["qname"] == "q1"
        assert rows[0]["tname"] == "phrog_1"
        assert rows[0]["score"] == 50.0

    def test_parse_trna_gff(self, tmp_path: Path) -> None:
        gff = tmp_path / "trna.gff"
        gff.write_text(
            "##gff-version 3\n"
            "contig\ttRNAscan-SE\ttRNA\t100\t120\t50.0\t+\t.\t"
            "ID=trna1;isotype=Met;anticodon=CAT\n"
        )
        rows = list(utils.parse_trna_gff(str(gff)))
        assert len(rows) == 1
        assert rows[0]["qname"] == "contig"
        assert rows[0]["trna_type"] == "Met"
        assert rows[0]["anticodon"] == "CAT"

    def test_parse_blastp(self, tmp_path: Path) -> None:
        blast = tmp_path / "out.blast"
        blast.write_text(
            "q1\tt1\t95.0\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200.0\n"
        )
        rows = list(utils.parse_blastp(str(blast)))
        assert len(rows) == 1
        assert rows[0]["qname"] == "q1"
        assert rows[0]["score"] == 200.0


class TestCoordinates:
    def test_get_cordinates_forward(self) -> None:
        import pandas as pd

        row = pd.Series(
            {"qname": "c1", "trna_no": 1, "begin": 100, "end": 120, "trna_type": "Met", "score": 50.0}
        )
        result = utils.get_cordinates(row)
        assert result["begin"] == 100
        assert result["end"] == 120
        assert result["strand"] == 1

    def test_get_cordinates_reverse(self) -> None:
        import pandas as pd

        row = pd.Series(
            {"qname": "c1", "trna_no": 1, "begin": 120, "end": 100, "trna_type": "Met", "score": 50.0}
        )
        result = utils.get_cordinates(row)
        assert result["begin"] == 100
        assert result["end"] == 120
        assert result["strand"] == -1

    def test_get_coordinates_alias(self) -> None:
        import pandas as pd

        row = pd.Series(
            {"qname": "c1", "trna_no": 1, "begin": 100, "end": 120, "trna_type": "Met", "score": 50.0}
        )
        assert utils.get_coordinates is utils.get_cordinates
        result = utils.get_coordinates(row)
        assert result["strand"] == 1


class TestPyrodigal:
    def test_call_genes_pyrodigal(self, tmp_path: Path) -> None:
        fna = tmp_path / "input.fna"
        fna.write_text(">contig1\n" + "ATG" * 200 + "TAA\n")
        out_faa = tmp_path / "proteins.faa"
        out_gff = tmp_path / "proteins.gff"

        pipeline.call_genes_pyrodigal(fna, out_faa, out_gff)

        assert out_faa.is_file()
        assert out_gff.is_file()
        faa_text = out_faa.read_text()
        assert ">contig1_1" in faa_text

    def test_call_genes_pyrodigal_translation_table(self, tmp_path: Path) -> None:
        fna = tmp_path / "input.fna"
        fna.write_text(">contig1\n" + "ATG" * 200 + "TAA\n")
        out_faa = tmp_path / "proteins.faa"
        out_gff = tmp_path / "proteins.gff"

        pipeline.call_genes_pyrodigal(
            fna, out_faa, out_gff, translation_table=4
        )

        assert out_faa.is_file()
        assert out_gff.is_file()
        faa_text = out_faa.read_text()
        assert ">contig1_1" in faa_text


class TestPyhmmer:
    def test_search_hmms_pyhmmer(self, tmp_path: Path) -> None:
        from pyhmmer.easel import Alphabet, TextSequence
        from pyhmmer.plan7 import Builder, Background

        alphabet = Alphabet.amino()
        background = Background(alphabet)
        builder = Builder(alphabet)

        seq = TextSequence(
            name=b"query1",
            sequence="MKTLLILTLGVMMMVSAKSSDKDLSQNLQINLKEQLNQLRQKKMQQNKENIEKQLTQLEASLQQAQQNREQEDRLRTEIEALLAKNGNPDEITAVAEAAMKIAQSATDVTVEIGMGKSGGIAAKILAEAGFDVTAISVDNGDNVLQAGESKVTAQLTEQGVKLINDNPQGEIPLAAGSGGFLGGLLAPAKK",
        ).digitize(alphabet)
        hmm, _, _ = builder.build(seq, background)
        hmm.name = b"test_hmm"

        db_dir = tmp_path / "hmmdb"
        db_dir.mkdir()
        with open(db_dir / "test.hmm", "wb") as f:
            hmm.write(f)

        faa = tmp_path / "proteins.faa"
        seq2 = TextSequence(
            name=b"seq1",
            description=b"seq1",
            sequence="MKTLLILTLGVMMMVSAKSSDKDLSQNLQINLKEQLNQLRQKKMQQNKENIEKQLTQLEASLQQAQQNREQEDRLRTEIEALLAKNGNPDEITAVAEAAMKIAQSATDVTVEIGMGKSGGIAAKILAEAGFDVTAISVDNGDNVLQAGESKVTAQLTEQGVKLINDNPQGEIPLAAGSGGFLGGLLAPAKK",
        )
        with open(faa, "w") as fh:
            seq_name = seq2.name.decode() if isinstance(seq2.name, bytes) else seq2.name
            fh.write(f">{seq_name}\n{seq2.sequence}\n")

        out_txt = tmp_path / "hmmsearch.txt"
        out_csv = tmp_path / "hmmsearch.csv"

        pipeline.search_hmms_pyhmmer(faa, db_dir, out_txt, out_csv, threads=1)

        assert out_txt.is_file()
        assert out_csv.is_file()
        rows = list(utils.parse_hmmsearch(out_txt))
        assert len(rows) >= 1
        assert rows[0]["qname"] == "seq1"
        assert rows[0]["tname"] == "test_hmm"


class TestAnnotations:
    def test_generate_plots_and_annotations_default_pdf_and_genbank(
        self, tmp_path: Path
    ) -> None:
        from Bio import SeqIO

        from pca.annotations import generate_plots_and_annotations

        # Create a minimal protein GFF and matching FASTA
        fna = tmp_path / "input.fna"
        fna.write_text(">contig1\n" + "ATG" * 300 + "TAA\n")

        proteins_gff = tmp_path / "proteins.gff"
        proteins_gff.write_text(
            "##gff-version 3\n"
            "##sequence-region contig1 1 903\n"
            "contig1\tpyrodigal\tCDS\t1\t903\t100.0\t+\t0\tID=contig1_1\n"
        )

        hmmsearch_txt = tmp_path / "hmmsearch.txt"
        hmmsearch_txt.write_text(
            "contig1_1\t-\tphrog_1\t-\t1e-10\t100.0\t0.0\t1e-5\t50.0\t0.0\n"
        )

        meta = tmp_path / "meta.csv"
        meta.write_text("#phrog,Annotation,Category,color\nphrog_1,test,unknown,#c9c9c9\n")

        trna_gff = tmp_path / "trna.gff"
        trna_gff.write_text("##gff-version 3\n")

        out_dir = tmp_path / "out"
        out_dir.mkdir()

        generate_plots_and_annotations(
            tmp_dir=out_dir,
            hmmsearch_dir=hmmsearch_txt,
            trna_dir=trna_gff,
            meta_dir=meta,
            gff_dir=proteins_gff,
            input_fasta=fna,
        )

        assert (out_dir / "annotations" / "annotations.gff").is_file()
        assert (out_dir / "annotations" / "annotations.gbk").is_file()
        assert (out_dir / "plots" / "contig1.pdf").is_file()

        gbk_records = list(SeqIO.parse(out_dir / "annotations" / "annotations.gbk", "genbank"))
        assert len(gbk_records) == 1
        assert gbk_records[0].id == "contig1"
        assert len(gbk_records[0].seq) == 903

    def test_generate_plots_and_annotations_html_output(self, tmp_path: Path) -> None:
        from pca.annotations import generate_plots_and_annotations

        fna = tmp_path / "input.fna"
        fna.write_text(">contig1\n" + "ATG" * 300 + "TAA\n")

        proteins_gff = tmp_path / "proteins.gff"
        proteins_gff.write_text(
            "##gff-version 3\n"
            "##sequence-region contig1 1 903\n"
            "contig1\tpyrodigal\tCDS\t1\t903\t100.0\t+\t0\tID=contig1_1\n"
        )

        hmmsearch_txt = tmp_path / "hmmsearch.txt"
        hmmsearch_txt.write_text(
            "contig1_1\t-\tphrog_1\t-\t1e-10\t100.0\t0.0\t1e-5\t50.0\t0.0\n"
        )

        meta = tmp_path / "meta.csv"
        meta.write_text("#phrog,Annotation,Category,color\nphrog_1,test,unknown,#c9c9c9\n")

        trna_gff = tmp_path / "trna.gff"
        trna_gff.write_text("##gff-version 3\n")

        out_dir = tmp_path / "out"
        out_dir.mkdir()

        generate_plots_and_annotations(
            tmp_dir=out_dir,
            hmmsearch_dir=hmmsearch_txt,
            trna_dir=trna_gff,
            meta_dir=meta,
            gff_dir=proteins_gff,
            input_fasta=fna,
            plot_formats=["html"],
        )

        assert (out_dir / "annotations" / "annotations.gff").is_file()
        assert (out_dir / "annotations" / "annotations.gbk").is_file()
        assert (out_dir / "plots" / "contig1.html").is_file()

    def test_convert_to_html_from_genbank(self, tmp_path: Path) -> None:
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, FeatureLocation
        from Bio.SeqRecord import SeqRecord

        gbk = tmp_path / "genome.gbk"
        record = SeqRecord(
            Seq("ATG" * 200 + "TAA"),
            id="contig1",
            annotations={"molecule_type": "DNA"},
        )
        record.features.append(
            SeqFeature(
                FeatureLocation(0, 300, strand=1),
                type="CDS",
                qualifiers={"locus_tag": ["gene_1"], "product": ["tail protein"]},
            )
        )
        SeqIO.write(record, gbk, "genbank")

        out_html = tmp_path / "report.html"
        paths = annotations.convert_to_html(gbk, out_html)

        assert len(paths) == 1
        assert paths[0] == out_html
        assert out_html.is_file()
        content = out_html.read_text()
        assert "contig1" in content
        assert "tail protein" in content

    def test_convert_to_html_from_gff_with_fasta(self, tmp_path: Path) -> None:
        from pca.annotations import convert_to_html

        fasta = tmp_path / "genome.fna"
        fasta.write_text(">contig1\n" + "ATG" * 100 + "TAA\n")

        gff = tmp_path / "genome.gff"
        gff.write_text(
            "##gff-version 3\n"
            "##sequence-region contig1 1 303\n"
            "contig1\tprodigal\tCDS\t1\t300\t.\t+\t0\tID=gene_1;Name=head protein\n"
        )

        out_dir = tmp_path / "html_out"
        paths = convert_to_html(gff, out_dir, fasta_path=fasta)

        assert len(paths) == 1
        assert paths[0].name == "contig1.html"
        assert paths[0].is_file()
        assert "head protein" in paths[0].read_text()

