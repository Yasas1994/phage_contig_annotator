"""Multi-copy viral contig detection.

Detects contigs that likely contain multiple copies of a genome unit using two
complementary methods:

1. Mean k-mer frequency (CheckV-style) – estimates copy number from k-mer
   redundancy. Fast and robust for uniform-GC genomes.
2. FFT on sliding-window GC content – detects periodicity in sequence
   composition and estimates the genome unit length. Works on chimeric or
   partial contigs.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from Bio import SeqIO

logger = None  # logging is optional; imported lazily where needed

# Optional kcounter (CheckV uses this). Falls back to pure Python if absent.
try:
    import kcounter

    HAS_KCOUNTER = True
except ImportError:  # pragma: no cover
    HAS_KCOUNTER = False


_RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]


def _canonical_kmer(kmer: str) -> str:
    """Return the lexicographically smaller of a k-mer and its reverse complement."""
    rc = _reverse_complement(kmer)
    return kmer if kmer <= rc else rc


def _count_kmers_python(seq: str, k: int, canonical_kmers: bool = True) -> Counter:
    """Pure-Python k-mer counter (fallback when kcounter is unavailable)."""
    seq = seq.upper()
    counts: Counter = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" in kmer:
            continue
        counts[_canonical_kmer(kmer) if canonical_kmers else kmer] += 1
    return counts


def mean_kmer_frequency(seq: str, k: int = 21) -> Optional[float]:
    """
    CheckV-style mean k-mer frequency.

    For a single-copy genome most k-mers appear once, so the mean is ≈ 1.0.
    For N tandem copies the same k-mers appear roughly N times, so the mean
    is ≈ N.

    Returns ``None`` if the sequence is too short to count k-mers.
    """
    if len(seq) < k:
        return None

    if HAS_KCOUNTER:
        counts = list(kcounter.count_kmers(seq, k, canonical_kmers=True).values())
    else:
        counts = list(_count_kmers_python(seq, k, canonical_kmers=True).values())

    if not counts:
        return None
    return float(np.mean(counts))


def sliding_gc_content(seq: str, window: int, step: int) -> np.ndarray:
    """Return GC fraction in sliding windows across ``seq``."""
    seq = seq.upper()
    gc_vals = []
    for i in range(0, len(seq) - window + 1, step):
        w = seq[i : i + window]
        gc = (w.count("G") + w.count("C")) / len(w)
        gc_vals.append(gc)
    return np.array(gc_vals)


def fft_dominant_period(
    signal: np.ndarray,
    step: int,
    contig_len: int,
    min_period_bp: int = 1000,
) -> tuple[Optional[float], Optional[float]]:
    """
    Run FFT on a 1D signal and return the dominant period and its SNR.

    Returns ``(None, None)`` if the signal has fewer than 4 points or no
    plausible period is found.
    """
    n = len(signal)
    if n < 4:
        return None, None

    detrended = signal - signal.mean()
    fft_mag = np.abs(np.fft.rfft(detrended))
    freqs = np.fft.rfftfreq(n, d=step)  # cycles per bp

    # Exclude DC and convert to periods in bp.
    fft_mag_no_dc = fft_mag[1:]
    freqs_no_dc = freqs[1:]

    max_period_bp = contig_len / 1.5
    valid = (1.0 / freqs_no_dc >= min_period_bp) & (1.0 / freqs_no_dc <= max_period_bp)

    if not valid.any():
        return None, None

    valid_mag = fft_mag_no_dc.copy()
    valid_mag[~valid] = 0.0

    peak_idx = int(np.argmax(valid_mag))
    peak_period_bp = 1.0 / freqs_no_dc[peak_idx]
    peak_power_ratio = float(
        valid_mag[peak_idx] / (np.mean(fft_mag_no_dc[valid]) + 1e-9)
    )

    return peak_period_bp, peak_power_ratio


def reconcile_copy_number(
    contig_len: int,
    mean_kmer_freq: Optional[float],
    fft_period_bp: Optional[float],
    fft_snr: Optional[float],
    kmer_threshold: float = 1.3,
    fft_snr_threshold: float = 3.0,
) -> dict:
    """
    Combine k-mer and FFT signals into a final copy-number call.

    Returns ``copies_kmer``, ``copies_fft``, ``copies_final``,
    ``genome_unit_bp``, ``confidence`` and ``flag``.
    """
    result = {
        "copies_kmer": None,
        "copies_fft": None,
        "copies_final": 1,
        "genome_unit_bp": None,
        "confidence": "low",
        "flag": "insufficient_data",
    }

    kmer_copies = None
    if mean_kmer_freq is not None:
        kmer_copies = max(1, round(mean_kmer_freq))
        result["copies_kmer"] = kmer_copies

    fft_copies = None
    if fft_period_bp is not None and fft_snr is not None:
        if fft_snr >= fft_snr_threshold:
            fft_copies = max(1, round(contig_len / fft_period_bp))
            result["copies_fft"] = fft_copies
            result["genome_unit_bp"] = int(round(fft_period_bp))

    kmer_positive = (kmer_copies is not None) and (mean_kmer_freq >= kmer_threshold)
    fft_positive = fft_copies is not None and fft_copies >= 2

    if not kmer_positive and not fft_positive:
        result.update(copies_final=1, confidence="high", flag="single_copy")
    elif kmer_positive and fft_positive:
        if kmer_copies == fft_copies:
            result.update(
                copies_final=kmer_copies,
                genome_unit_bp=result["genome_unit_bp"],
                confidence="high",
                flag="agreement",
            )
        else:
            result.update(
                copies_final=kmer_copies,
                genome_unit_bp=result["genome_unit_bp"],
                confidence="low",
                flag=f"disagreement_kmer={kmer_copies}_fft={fft_copies}",
            )
    elif kmer_positive and not fft_positive:
        result.update(
            copies_final=kmer_copies,
            confidence="medium",
            flag="kmer_only_uniform_gc_suspected",
        )
    elif not kmer_positive and fft_positive:
        result.update(
            copies_final=fft_copies,
            genome_unit_bp=result["genome_unit_bp"],
            confidence="medium",
            flag="fft_only_possible_chimera",
        )

    return result


def analyze_contig(
    seq: str,
    k: int = 21,
    window: int = 500,
    step: int = 100,
    min_period_bp: int = 1000,
) -> dict:
    """Run both multi-copy detection methods and return a combined result."""
    seq = seq.upper().replace(" ", "")
    contig_len = len(seq)

    mean_kf = mean_kmer_frequency(seq, k)
    gc = sliding_gc_content(seq, window, step)
    fft_period, fft_snr = fft_dominant_period(gc, step, contig_len, min_period_bp)
    rec = reconcile_copy_number(contig_len, mean_kf, fft_period, fft_snr)

    return {
        "contig_len": contig_len,
        "mean_kmer_freq": round(mean_kf, 4) if mean_kf is not None else None,
        "fft_period_bp": int(round(fft_period)) if fft_period is not None else None,
        "fft_snr": round(fft_snr, 2) if fft_snr is not None else None,
        "gc_windows": len(gc),
        **rec,
    }


def detect_multicopy(
    fasta_path: str | Path,
    k: int = 21,
    window: int = 500,
    step: int = 100,
    min_period_bp: int = 1000,
    min_len: int = 2000,
) -> pd.DataFrame:
    """
    Run multi-copy detection on every contig in ``fasta_path``.

    Returns a ``pandas.DataFrame`` with one row per contig and columns for
    copy-number estimates, confidence, and flag.
    """
    rows = []
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        seq = str(rec.seq)
        if len(seq) < min_len:
            continue
        row = analyze_contig(seq, k=k, window=window, step=step, min_period_bp=min_period_bp)
        row["contig_id"] = rec.id
        rows.append(row)

    if not rows:
        return pd.DataFrame(
            columns=[
                "contig_id",
                "contig_len",
                "mean_kmer_freq",
                "copies_kmer",
                "fft_period_bp",
                "fft_snr",
                "copies_fft",
                "copies_final",
                "genome_unit_bp",
                "confidence",
                "flag",
            ]
        )

    df = pd.DataFrame(rows)
    column_order = [
        "contig_id",
        "contig_len",
        "mean_kmer_freq",
        "copies_kmer",
        "fft_period_bp",
        "fft_snr",
        "copies_fft",
        "copies_final",
        "genome_unit_bp",
        "confidence",
        "flag",
    ]
    return df[[c for c in column_order if c in df.columns]]
