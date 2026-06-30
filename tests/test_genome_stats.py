"""Unit tests for pca.genome_stats."""

from __future__ import annotations

import pytest
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from pca import genome_stats


class TestStrandSwitchFrequency:
    def test_empty_returns_zero(self) -> None:
        assert genome_stats.strand_switch_frequency([]) == 0.0

    def test_single_gene_returns_zero(self) -> None:
        assert genome_stats.strand_switch_frequency(["+"]) == 0.0

    def test_all_same_strand(self) -> None:
        assert genome_stats.strand_switch_frequency(["+", "+", "+"]) == 0.0

    def test_alternating_strands(self) -> None:
        assert genome_stats.strand_switch_frequency(["+", "-", "+", "-"]) == 1.0

    def test_half_switch_rate(self) -> None:
        assert genome_stats.strand_switch_frequency(["+", "+", "-", "-"]) == 1 / 3


class TestStrandRunStats:
    def test_empty(self) -> None:
        assert genome_stats.strand_run_stats([]) == {
            "n_runs": 0,
            "mean_run_length": 0.0,
            "max_run_length": 0,
            "switch_rate": 0.0,
        }

    def test_single_run(self) -> None:
        assert genome_stats.strand_run_stats(["+", "+", "+"]) == {
            "n_runs": 1,
            "mean_run_length": 3.0,
            "max_run_length": 3,
            "switch_rate": 0.0,
        }

    def test_multiple_runs(self) -> None:
        stats = genome_stats.strand_run_stats(["+", "+", "-", "-", "+"])
        assert stats["n_runs"] == 3
        assert stats["mean_run_length"] == pytest.approx(5 / 3)
        assert stats["max_run_length"] == 2
        assert stats["switch_rate"] == pytest.approx(2 / 4)


class TestStrandBiasIndex:
    def test_empty(self) -> None:
        assert genome_stats.strand_bias_index([]) == 0.0

    def test_balanced(self) -> None:
        assert genome_stats.strand_bias_index(["+", "-"]) == 0.0

    def test_all_plus(self) -> None:
        assert genome_stats.strand_bias_index(["+", "+", "+"]) == 1.0

    def test_all_minus(self) -> None:
        assert genome_stats.strand_bias_index(["-", "-"]) == -1.0

    def test_unbalanced(self) -> None:
        assert genome_stats.strand_bias_index(["+", "+", "-"]) == pytest.approx(1 / 3)


class TestWeightedStrandSwitchFrequency:
    def test_empty(self) -> None:
        assert genome_stats.weighted_strand_switch_frequency([]) == 0.0

    def test_single_gene(self) -> None:
        assert genome_stats.weighted_strand_switch_frequency(
            [{"start": 1, "strand": "+"}]
        ) == 0.0

    def test_no_switch(self) -> None:
        genes = [
            {"start": 100, "strand": "+"},
            {"start": 1000, "strand": "+"},
        ]
        assert genome_stats.weighted_strand_switch_frequency(genes) == 0.0

    def test_switch_with_span(self) -> None:
        genes = [
            {"start": 100, "strand": "+"},
            {"start": 1000, "strand": "-"},
        ]
        assert genome_stats.weighted_strand_switch_frequency(genes) == 1.0

    def test_weighted_average(self) -> None:
        genes = [
            {"start": 100, "strand": "+"},
            {"start": 200, "strand": "-"},
            {"start": 1000, "strand": "+"},
        ]
        # spans: 100 and 800; switched: 100 and 800 -> 900 / 900 = 1.0
        assert genome_stats.weighted_strand_switch_frequency(genes) == 1.0


class TestTandemRepeatStats:
    def test_no_repeats(self, tmp_path: Path) -> None:
        dat = tmp_path / "empty.dat"
        dat.write_text("")
        stats = genome_stats.tandem_repeat_stats(dat, "contig1", 1000)
        assert stats["n_tandem_repeats"] == 0
        assert stats["total_tandem_repeat_length"] == 0
        assert stats["tandem_repeat_density"] == 0.0

    def test_counts_only_matching_contig(self, tmp_path: Path) -> None:
        dat = tmp_path / "repeats.dat"
        dat.write_text(
            "Sequence: contig1\n"
            "Parameters: 2 7 7 80 10 50 2000\n"
            "10 30 5 4.0 5 20 0 50 15.0 20.0 10.0 5.0 1.5 AAAAA AAAAA\n"
            "100 130 10 3.0 10 25 1 60 10.0 10.0 10.0 10.0 1.5 CCCCC CCCCC\n"
            "Sequence: contig2\n"
            "Parameters: 2 7 7 80 10 50 2000\n"
            "1 20 4 5.0 4 18 0 40 12.0 8.0 10.0 10.0 1.4 GGGGG GGGGG\n"
        )
        stats = genome_stats.tandem_repeat_stats(dat, "contig1", 1000)
        # lengths: 21 and 31 -> total 52
        assert stats["n_tandem_repeats"] == 2
        assert stats["total_tandem_repeat_length"] == 52
        assert stats["tandem_repeat_density"] == pytest.approx(52 / 1000)


class TestContigStatsFromGff:
    def _make_gff(self, tmp_path: Path, features: list[SeqFeature]) -> Path:
        record = SeqRecord(seq=Seq("ACGTACGTACGT"), id="contig1")
        record.features = features
        gff_path = tmp_path / "test.gff"
        with open(gff_path, "w") as fh:
            from BCBio import GFF

            GFF.write([record], fh)
        return gff_path

    def test_single_contig(self, tmp_path: Path) -> None:
        features = [
            SeqFeature(FeatureLocation(10, 20, strand=1), type="CDS"),
            SeqFeature(FeatureLocation(25, 35, strand=-1), type="CDS"),
            SeqFeature(FeatureLocation(40, 50, strand=-1), type="CDS"),
            SeqFeature(FeatureLocation(55, 65, strand=1), type="CDS"),
        ]
        gff_path = self._make_gff(tmp_path, features)
        stats = genome_stats.contig_stats_from_gff(gff_path)
        assert set(stats.keys()) == {"contig1"}
        contig = stats["contig1"]
        assert contig["n_genes"] == 4
        # strands: + - - + -> switches at positions 0, 2 -> 2/3
        assert contig["strand_switch_frequency"] == pytest.approx(2 / 3)
        assert contig["strand_bias_index"] == 0.0
        assert contig["n_runs"] == 3

    def test_filter_by_contig_id(self, tmp_path: Path) -> None:
        features = [
            SeqFeature(FeatureLocation(10, 20, strand=1), type="CDS"),
        ]
        gff_path = self._make_gff(tmp_path, features)
        stats = genome_stats.contig_stats_from_gff(gff_path, contig_id="missing")
        assert stats == {}


    def test_crispr_stats(self, tmp_path: Path) -> None:
        gff = tmp_path / "crispr.gff"
        gff.write_text(
            "##gff-version 3\n"
            "contig1\tCRISPRCasFinder\tCRISPR\t10\t30\t.\t+\t.\tID=CRISPR1\n"
            "contig1\tCRISPRCasFinder\tCRISPR\t50\t70\t.\t+\t.\tID=CRISPR2\n"
            "contig2\tCRISPRCasFinder\tCRISPR\t5\t25\t.\t+\t.\tID=CRISPR3\n"
        )
        stats = genome_stats.crispr_stats(gff, "contig1")
        assert stats["n_crispr_arrays"] == 2
        assert stats["total_crispr_array_length"] == 42

    def test_crispr_stats_no_match(self, tmp_path: Path) -> None:
        gff = tmp_path / "crispr.gff"
        gff.write_text(
            "##gff-version 3\n"
            "contig1\tCRISPRCasFinder\tCRISPR\t10\t30\t.\t+\t.\tID=CRISPR1\n"
        )
        stats = genome_stats.crispr_stats(gff, "missing")
        assert stats["n_crispr_arrays"] == 0
        assert stats["total_crispr_array_length"] == 0


class TestComputeGenomeStats:
    def test_basic(self, tmp_path: Path) -> None:
        fasta = tmp_path / "input.fna"
        fasta.write_text(
            ">contig1\nACGTACGTACGTACGTACGT\n>contig2\nAAAA\n"
        )
        gff = tmp_path / "genes.gff"
        record = SeqRecord(seq=Seq("ACGTACGTACGTACGTACGT"), id="contig1")
        record.features = [
            SeqFeature(FeatureLocation(10, 20, strand=1), type="CDS"),
            SeqFeature(FeatureLocation(25, 35, strand=-1), type="CDS"),
        ]
        with open(gff, "w") as fh:
            from BCBio import GFF

            GFF.write([record], fh)

        df = genome_stats.compute_genome_stats(fasta, gff)
        assert len(df) == 1
        assert df.iloc[0]["contig_id"] == "contig1"
        assert df.iloc[0]["length"] == 20
        assert df.iloc[0]["n_genes"] == 2


class TestComputeGenomeStatsWithCrispr:
    def test_basic(self, tmp_path: Path) -> None:
        fasta = tmp_path / "input.fna"
        fasta.write_text(
            ">contig1\nACGTACGTACGTACGTACGT\n>contig2\nAAAA\n"
        )
        gff = tmp_path / "genes.gff"
        record = SeqRecord(seq=Seq("ACGTACGTACGTACGTACGT"), id="contig1")
        record.features = [
            SeqFeature(FeatureLocation(10, 20, strand=1), type="CDS"),
            SeqFeature(FeatureLocation(25, 35, strand=-1), type="CDS"),
        ]
        with open(gff, "w") as fh:
            from BCBio import GFF

            GFF.write([record], fh)

        crispr_gff = tmp_path / "crispr.gff"
        crispr_gff.write_text(
            "contig1\tCRISPRCasFinder\tCRISPR\t1\t5\t.\t+\t.\tID=CRISPR1\n"
        )

        df = genome_stats.compute_genome_stats(fasta, gff, crispr_path=crispr_gff)
        assert len(df) == 1
        assert df.iloc[0]["n_crispr_arrays"] == 1
        assert df.iloc[0]["total_crispr_array_length"] == 5
