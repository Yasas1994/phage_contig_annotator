"""Unit tests for pca.multicopy."""
import pytest
from pathlib import Path
from unittest.mock import patch

from pca.multicopy import (
    _canonical,
    _count_kmers_python,
    _rc,
    _safe_float,
    analyze,
    detect_multicopy,
    mean_kmer_freq,
    minimap2_validator,
)


def test_rc() -> None:
    assert _rc("ATCG") == "CGAT"
    assert _rc("AAAA") == "TTTT"
    assert _rc("ACGTacgt") == "acgtACGT"


def test_canonical_prefers_lexicographically_smaller() -> None:
    assert _canonical("ATCG") == "ATCG"
    # Reverse complement of ATCG is CGAT, which is larger.
    assert _canonical("CGAT") == "ATCG"


def test_count_kmers_python_counts_canonical_kmers() -> None:
    seq = "ATCGATCG"  # 8 bp, k=5 gives 4 canonical 5-mers
    counts = _count_kmers_python(seq, k=5)
    assert sum(counts.values()) == 4
    # Each key is a canonical k-mer string of the requested length.
    assert all(isinstance(k, str) and len(k) == 5 for k in counts)


def test_mean_kmer_freq_for_single_copy() -> None:
    # A random-looking 100 bp sequence should have most canonical 21-mers unique.
    seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTC"
    freq = mean_kmer_freq(seq, k=21)
    assert freq is not None
    assert 1.0 <= freq <= 2.0


def test_mean_kmer_freq_for_repeated_sequence() -> None:
    # Repeating a 60 bp motif many times produces multi-copy k-mers.
    motif = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTC"
    seq = motif * 10  # 600 bp
    freq = mean_kmer_freq(seq, k=21)
    assert freq is not None
    assert freq >= 8.0


def test_mean_kmer_freq_returns_none_for_short_sequence() -> None:
    assert mean_kmer_freq("ATCG", k=21) is None


def test_safe_float() -> None:
    assert _safe_float("1.5") == 1.5
    assert _safe_float(2.0) == 2.0
    assert _safe_float(None) is None
    assert _safe_float("not a number") is None
    assert _safe_float(float("nan")) is None
    assert _safe_float(float("inf")) is None


def test_analyze_single_copy() -> None:
    seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTC"
    result = analyze(seq, k=21)
    assert result["copies_final"] == 1
    assert result["confidence"] == "high"
    assert result["flag"] == "PASS"
    assert result["validator"] == "none"
    assert result["validator_ok"] is False


def test_analyze_repeated_sequence_with_mocked_minimap2() -> None:
    motif = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTC"
    seq = motif * 10  # ~600 bp, but k=21 may not give enough. Use k=11.
    with patch("pca.multicopy.minimap2_validator") as mock_validator:
        mock_validator.return_value = {
            "ok": True,
            "mean_identity": 0.99,
            "boundaries": [60, 120, 180, 240, 300, 360, 420, 480, 540],
            "repeat_fraction": 0.95,
            "inferred_copies": 10,
            "n_alignments": 50,
            "note": "minimap2_confirmed",
        }
        result = analyze(seq, k=11)
    assert result["copies_final"] == 10
    assert result["validator_ok"] is True
    assert result["validator"] == "minimap2"
    assert result["confidence"] == "high"
    assert result["flag"] == "FAIL"


def test_minimap2_validator_no_alignments() -> None:
    with patch("pca.multicopy._run_minimap2_self", return_value=[]) as _:
        result = minimap2_validator("ATCG" * 100, 400, 2)
    assert result["ok"] is False
    assert result["note"] == "no_alignments"


def test_detect_multicopy_returns_dataframe_for_fasta(tmp_path: Path) -> None:
    fasta = tmp_path / "input.fna"
    fasta.write_text(
        ">short\n" + "ATCG" * 10 + "\n"  # 40 bp
        ">long\n" + "ATCG" * 600 + "\n"  # 2400 bp
    )
    df = detect_multicopy(fasta, k=21)
    assert len(df) == 2
    assert set(df.columns) >= {
        "contig_id", "contig_len", "mean_kmer_freq", "copies_kmer",
        "validator", "validator_score", "validator_snr", "validator_ok",
        "copies_final", "confidence", "flag", "note",
    }
    long_row = df[df["contig_id"] == "long"].iloc[0]
    assert long_row["contig_len"] == 2400
