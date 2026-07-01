"""Unit tests for pca.parsers."""

from pathlib import Path

from pca.parsers import parse_crispr_cas_finder_gff, parse_minced_gff, parse_trf_dat


def test_parse_crispr_cas_finder_gff_extracts_crispr_features(tmp_path: Path) -> None:
    gff = tmp_path / "crispr.gff"
    gff.write_text(
        "##gff-version 3\n"
        "contig1\tCRISPRCasFinder\tCRISPR\t10\t30\t.\t+\t.\tID=CRISPR1;Name=CRISPR1\n"
        "contig1\tCRISPRCasFinder\tCas\t50\t70\t.\t+\t.\tName=cas9\n"
        "contig2\tCRISPRCasFinder\tCRISPR\t5\t25\t23.5\t-\t.\tID=CRISPR2\n"
    )
    records = list(parse_crispr_cas_finder_gff(gff))
    assert len(records) == 2
    assert records[0]["qname"] == "contig1"
    assert records[0]["begin"] == 10
    assert records[0]["end"] == 30
    assert records[0]["strand"] == "+"
    assert records[0]["score"] == 0.0
    assert records[1]["qname"] == "contig2"
    assert records[1]["score"] == 23.5
    assert records[1]["strand"] == "-"


def test_parse_crispr_cas_finder_gff_skips_header_and_blank_lines(
    tmp_path: Path,
) -> None:
    gff = tmp_path / "crispr.gff"
    gff.write_text(
        "# CRISPRCasFinder output\n"
        "\n"
        "contig1\tCRISPRCasFinder\tCRISPR\t10\t30\t.\t+\t.\tID=CRISPR1\n"
    )
    records = list(parse_crispr_cas_finder_gff(gff))
    assert len(records) == 1


def test_parse_minced_gff_extracts_crispr_repeat_regions(tmp_path: Path) -> None:
    gff = tmp_path / "minced.gff"
    gff.write_text(
        "##gff-version 3\n"
        "contig1\tminced:0.4.2\trepeat_region\t10\t30\t5\t.\t.\tID=CRISPR1;rpt_type=direct;rpt_family=CRISPR;rpt_unit_seq=GTTCC\n"
        "contig1\tminced:0.4.2\trepeat_region\t50\t70\t4\t+\t.\tID=CRISPR2;rpt_family=CRISPR\n"
        "contig1\tminced:0.4.2\trepeat_region\t90\t110\t3\t.\t.\tID=other;rpt_family=other\n"
        "contig2\tminced:0.4.2\tCDS\t1\t100\t.\t+\t.\tID=gene1\n"
    )
    records = list(parse_minced_gff(gff))
    assert len(records) == 2
    assert records[0]["qname"] == "contig1"
    assert records[0]["begin"] == 10
    assert records[0]["end"] == 30
    assert records[0]["score"] == 5.0
    assert records[0]["strand"] == "."
    assert records[1]["begin"] == 50
    assert records[1]["score"] == 4.0


def test_parse_minced_gff_skips_header_and_non_crispr_features(
    tmp_path: Path,
) -> None:
    gff = tmp_path / "minced.gff"
    gff.write_text(
        "# MinCED output\n"
        "\n"
        "contig1\tminced:0.4.2\trepeat_region\t10\t30\t.\t.\t.\tID=CRISPR1;rpt_family=CRISPR\n"
    )
    records = list(parse_minced_gff(gff))
    assert len(records) == 1

def test_parse_trf_dat_extracts_repeats_per_sequence(tmp_path: Path) -> None:
    dat = tmp_path / "repeats.dat"
    dat.write_text(
        "Tandem Repeats Finder Program written by:\n\n"
        "Gary Benson\n"
        "Program in Bioinformatics\n"
        "Boston University\n"
        "Version 4.07b\n\n"
        "Sequence: contig_1\n\n"
        "Parameters: 2 7 7 80 10 50 2000\n\n"
        "3404 3433 14 2.1 14 93 0 51 30 16 30 23 1.96 "
        "AGCAGCACTTGAGT AGCAGTACTTGAGTAGCAGCACTTGAGTAG\n"
        "6344 6383 18 2.2 18 95 0 71 45 2 10 42 1.51 "
        "ATTATTATAAGGAATTAC ATTATTATAAGGAATTATATTATTATAAGGAATTACATTA\n\n"
        "Sequence: contig_2\n\n"
        "Parameters: 2 7 7 80 10 50 2000\n\n"
        "15030 15056 14 1.9 14 100 0 54 66 0 0 33 0.92 "
        "TTAAATAAATAAAT TTAAATAAATAAATTTAAATAAATAAA\n"
    )
    records = list(parse_trf_dat(dat))
    assert len(records) == 3

    first, second, third = records
    assert first["contig"] == "contig_1"
    assert first["begin"] == 3404
    assert first["end"] == 3433
    assert first["period"] == 14
    assert first["copies"] == 2.1
    assert first["consensus_size"] == 14
    assert first["matches"] == 93
    assert first["indels"] == 0
    assert first["score"] == 51
    assert first["entropy"] == 1.96
    assert first["consensus"] == "AGCAGCACTTGAGT"
    assert first["strand"] == 0

    assert second["contig"] == "contig_1"
    assert third["contig"] == "contig_2"
    assert third["begin"] == 15030


def test_parse_trf_dat_skips_malformed_lines(tmp_path: Path) -> None:
    dat = tmp_path / "repeats.dat"
    dat.write_text(
        "Sequence: contig_1\n\n"
        "Parameters: 2 7 7 80 10 50 2000\n\n"
        "100 130 6 4.0 6 80 0 30 20 20 20 20 1.4 ACGTAC ACGTACACGTAC\n"
        "this is not a valid data line\n"
    )
    records = list(parse_trf_dat(dat))
    assert len(records) == 1
    assert records[0]["period"] == 6
