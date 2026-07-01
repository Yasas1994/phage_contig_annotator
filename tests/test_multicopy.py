"""Unit tests for pca.multicopy."""
from pathlib import Path

from pca.multicopy import (
    _canonical_kmer,
    _reverse_complement,
    analyze_contig,
    detect_multicopy,
    fft_dominant_period,
    mean_kmer_frequency,
    reconcile_copy_number,
    sliding_gc_content,
)


def test_reverse_complement() -> None:
    assert _reverse_complement("ATCG") == "CGAT"
    assert _reverse_complement("AAAA") == "TTTT"


def test_canonical_kmer_prefers_lexicographically_smaller() -> None:
    assert _canonical_kmer("ATCG") == "ATCG"
    # Reverse complement of ATCG is CGAT, which is larger.
    assert _canonical_kmer("CGAT") == "ATCG"


def test_mean_kmer_frequency_for_single_copy() -> None:
    # A random-looking 100 bp sequence should have most canonical 21-mers unique.
    seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTC"
    freq = mean_kmer_frequency(seq, k=21)
    assert freq is not None
    assert 1.0 <= freq <= 2.0


def test_mean_kmer_frequency_for_repeated_sequence() -> None:
    # Repeating a 60 bp motif many times produces multi-copy k-mers.
    motif = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTC"
    seq = motif * 10  # 600 bp
    freq = mean_kmer_frequency(seq, k=21)
    assert freq is not None
    assert freq >= 8.0


def test_mean_kmer_frequency_returns_none_for_short_sequence() -> None:
    assert mean_kmer_frequency("ATCG", k=21) is None


def test_sliding_gc_computes_windowed_fractions() -> None:
    seq = "GC" * 50 + "AT" * 50  # 50% GC then 0% GC, total 200 bp
    gc = sliding_gc_content(seq, window=100, step=100)
    assert len(gc) == 2
    assert abs(gc[0] - 1.0) < 1e-9
    assert abs(gc[1] - 0.0) < 1e-9


def test_fft_dominant_period_finds_periodicity() -> None:
    # Synthetic 3000 bp sequence: 1000 bp blocks of different GC content.
    block = "GC" * 250 + "AT" * 250  # 1000 bp, 50% GC then 0% GC
    seq = block * 3
    gc = sliding_gc_content(seq, window=100, step=100)
    period, snr = fft_dominant_period(gc, step=100, contig_len=3000, min_period_bp=500)
    assert period is not None
    assert snr is not None
    # Period should be close to 1000 bp.
    assert 900 <= period <= 1100


def test_fft_dominant_period_returns_none_for_short_signal() -> None:
    assert fft_dominant_period(
        sliding_gc_content("ATCG" * 10, window=10, step=10),
        step=10,
        contig_len=40,
    ) == (None, None)


def test_reconcile_single_copy() -> None:
    result = reconcile_copy_number(
        contig_len=5000,
        mean_kmer_freq=1.1,
        fft_period_bp=None,
        fft_snr=None,
    )
    assert result["copies_final"] == 1
    assert result["confidence"] == "high"
    assert result["flag"] == "single_copy"


def test_reconcile_agreement() -> None:
    result = reconcile_copy_number(
        contig_len=6000,
        mean_kmer_freq=3.0,
        fft_period_bp=2000.0,
        fft_snr=5.0,
    )
    assert result["copies_final"] == 3
    assert result["confidence"] == "high"
    assert result["flag"] == "agreement"
    assert result["genome_unit_bp"] == 2000


def test_reconcile_disagreement() -> None:
    result = reconcile_copy_number(
        contig_len=6000,
        mean_kmer_freq=3.0,
        fft_period_bp=1500.0,
        fft_snr=5.0,
    )
    assert result["copies_final"] == 3
    assert result["confidence"] == "low"
    assert "disagreement" in result["flag"]


def test_analyze_contig_returns_expected_keys() -> None:
    seq = "ATCG" * 1000
    result = analyze_contig(seq, k=21, window=100, step=100)
    expected_keys = {
        "contig_len",
        "mean_kmer_freq",
        "fft_period_bp",
        "fft_snr",
        "gc_windows",
        "copies_kmer",
        "copies_fft",
        "copies_final",
        "genome_unit_bp",
        "confidence",
        "flag",
    }
    assert set(result.keys()) == expected_keys


def test_detect_multicopy_filters_short_contigs(tmp_path: Path) -> None:
    fasta = tmp_path / "input.fna"
    fasta.write_text(
        ">short\n" + "ATCG" * 100 + "\n"  # 400 bp, below default 2000
        ">long\n" + "ATCG" * 600 + "\n"    # 2400 bp, passes filter
    )

    df = detect_multicopy(fasta, k=21, window=100, step=100, min_len=2000)
    assert len(df) == 1
    assert df.iloc[0]["contig_id"] == "long"


def test_detect_multicopy_empty_result() -> None:
    fasta = Path("/tmp/pca_multicopy_empty.fna")
    fasta.write_text(">short\n" + "ATCG" * 100 + "\n")
    df = detect_multicopy(fasta, min_len=2000)
    assert df.empty
    assert "contig_id" in df.columns
    assert "copies_final" in df.columns
