"""Unit tests for pca.annotations."""

from pathlib import Path

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from pca.annotations import (
    _compute_gc_metrics,
    _create_trf_feature,
    _parse_translation_tables,
    _write_static_plot,
)


def test_compute_gc_metrics_for_simple_sequence() -> None:
    seq = "GCGCGCATAT" * 10  # 60% GC, equal G and C -> skew 0
    positions, gc_values, skew_values, cumulative = _compute_gc_metrics(seq, num_windows=10)
    assert len(positions) == len(gc_values) == len(skew_values) == len(cumulative) == 10
    assert all(59 < g < 61 for g in gc_values)
    assert all(abs(s) < 1e-9 for s in skew_values)
    assert all(abs(c) < 1e-9 for c in cumulative)


def test_parse_translation_tables_extracts_per_contig_tables(tmp_path: Path) -> None:
    gff = tmp_path / "genes.gff"
    gff.write_text(
        '# Sequence Data: seqnum=1;seqlen=4421;seqhdr="NC_003438.1"\n'
        '# Model Data: version=pyrodigal.v3.7.1;run_type=Single;model="Ab initio";'
        'gc_cont=50.00;transl_table=4;uses_sd=1\n'
        '# Sequence Data: seqnum=2;seqlen=6000;seqhdr="NC_999999.1"\n'
        '# Model Data: version=pyrodigal.v3.7.1;run_type=Single;model="Ab initio";'
        'gc_cont=40.00;transl_table=11;uses_sd=1\n'
        'NC_003438.1\tprodigal-gv\tCDS\t1\t100\t.\t+\t0\tID=1_1\n'
    )
    tables = _parse_translation_tables(gff)
    assert tables == {"NC_003438.1": 4, "NC_999999.1": 11}


def test_parse_translation_tables_returns_empty_for_missing_headers(tmp_path: Path) -> None:
    gff = tmp_path / "genes.gff"
    gff.write_text("NC_003438.1\tprodigal\tCDS\t1\t100\t.\t+\t0\tID=1_1\n")
    assert _parse_translation_tables(gff) == {}


def test_write_static_plot_creates_pdf_and_png(tmp_path: Path) -> None:
    plots_dir = tmp_path / "plots"
    plots_dir.mkdir()
    record = SeqRecord(Seq("A" * 5000), id="test_contig")
    record.features.append(
        SeqFeature(
            FeatureLocation(100, 900, strand=1),
            type="CDS",
            qualifiers={
                "label": "test protein",
                "category": "other",
                "color": "#9467bd",
            },
        )
    )
    category_colors = {"other": "#9467bd", "unknown": "#c9c9c9"}
    _write_static_plot(record, plots_dir, "pdf", category_colors)
    _write_static_plot(record, plots_dir, "png", category_colors)
    pdf_path = plots_dir / "test_contig.pdf"
    png_path = plots_dir / "test_contig.png"
    assert pdf_path.exists() and pdf_path.stat().st_size > 0
    assert png_path.exists() and png_path.stat().st_size > 0


def test_create_trf_feature_builds_repeat_feature() -> None:
    row = pd.Series(
        {
            "begin": 100,
            "end": 200,
            "period": 12,
            "copies": 8.3,
            "score": 123,
            "entropy": 1.4,
            "consensus": "ATATATATATAT",
            "trf_no": 1,
        }
    )
    feature = _create_trf_feature(row)
    assert feature.type == "tandem_repeat"
    assert int(feature.location.start) == 100
    assert int(feature.location.end) == 200
    assert feature.location.strand == 0
    assert feature.qualifiers["category"] == "tandem repeat"
    assert feature.qualifiers["period"] == 12
    assert feature.qualifiers["copies"] == 8.3
    assert feature.qualifiers["consensus"] == "ATATATATATAT"
    assert feature.qualifiers["ID"] == "trf_1"
