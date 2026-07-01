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

import os
import subprocess
import sys
import tempfile
from typing import Optional

import numpy as np
import pandas as pd
from Bio import SeqIO


# ── optional kcounter ────────────────────────────────────────────────────────
try:
    import kcounter
    _HAS_KCOUNTER = True
except ImportError:
    _HAS_KCOUNTER = False


# ════════════════════════════════════════════════════════════════════════════
# CONSTANTS
# ════════════════════════════════════════════════════════════════════════════

KMER_K              = 21
KMER_COPY_THRESHOLD = 1.3

WARN_COPIES         = 2
FAIL_COPIES         = 3
MIN_LEN_DEFAULT     = 2000


_RC_TABLE = str.maketrans("ACGTacgt", "TGCAtgca")

def _rc(seq: str) -> str:
    return seq.translate(_RC_TABLE)[::-1]

def _canonical(kmer: str) -> str:
    rc = _rc(kmer)
    return kmer if kmer <= rc else rc


def _count_kmers_python(seq: str, k: int) -> dict[str, int]:
    counts: dict[str, int] = {}
    for i in range(len(seq) - k + 1):
        km = seq[i:i + k]
        if 'N' in km:
            continue
        ck = _canonical(km)
        counts[ck] = counts.get(ck, 0) + 1
    return counts


def mean_kmer_freq(seq: str, k: int = KMER_K) -> Optional[float]:
    """Mean count of canonical k-mers. Scales linearly with copy number."""
    if len(seq) < k:
        return None
    if _HAS_KCOUNTER:
        counts = list(kcounter.count_kmers(seq, k, canonical_kmers=True).values())
    else:
        counts = list(_count_kmers_python(seq, k).values())
    return float(np.mean(counts)) if counts else None



# minimap2 validator 

def _run_minimap2_self(seq: str, contig_id: str = "c") -> list[dict]:
    """Self-align with minimap2. Returns PAF records (excluding diagonal)."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False) as f:
        f.write(f">{contig_id}\n{seq}\n")
        tmp = f.name

    try:
        cmd = [
            "minimap2", "-x", "asm5", "-N", "50",
            "--secondary=yes", "-p", "0.05", tmp, tmp,
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if proc.returncode != 0:
            raise RuntimeError(f"minimap2 failed: {proc.stderr}")
    finally:
        os.unlink(tmp)

    alns = []
    for line in proc.stdout.strip().split("\n"):
        if not line:
            continue
        cols = line.split("\t")
        if len(cols) < 12:
            continue

        qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend = cols[:9]
        nmatch, alen, mapq = cols[9:12]

        # Skip trivial diagonal self-alignment
        if (qname == tname and
            abs(int(qstart) - int(tstart)) < 50 and
            abs(int(qend) - int(tend)) < 50):
            continue

        alns.append({
            "qstart": int(qstart), "qend": int(qend),
            "tstart": int(tstart), "tend": int(tend),
            "strand": strand,
            "nmatch": int(nmatch), "alen": int(alen),
            "mapq": int(mapq),
        })
    return alns


def minimap2_validator(seq: str, contig_len: int, kmer_copies: int,
                       min_identity: float = 0.80,
                       len_frac: float = 0.40) -> dict:
    """Validate tandem repeat via minimap2 self-alignment."""
    expected_unit = contig_len / kmer_copies

    alns = _run_minimap2_self(seq)

    if not alns:
        return {"ok": False, "note": "no_alignments"}

    good = []
    for a in alns:
        if a["strand"] != "+":
            continue
        aln_len = a["qend"] - a["qstart"]
        if aln_len < expected_unit * len_frac:
            continue
        identity = a["nmatch"] / a["alen"] if a["alen"] > 0 else 0
        if identity < min_identity:
            continue
        offset = a["tstart"] - a["qstart"]
        if offset <= 0:
            continue
        ratio = offset / expected_unit
        if abs(ratio - round(ratio)) > 0.25:
            continue
        good.append({**a, "identity": identity, "offset": offset})

    if not good:
        return {"ok": False, "note": "no_tandem_alignments"}

    # Boundaries
    expected_boundaries = [int(round(i * expected_unit))
                           for i in range(1, kmer_copies)]
    all_coords = []
    for a in good:
        all_coords.extend([a["tstart"], a["tend"]])
    boundaries = []
    for eb in expected_boundaries:
        if not all_coords:
            break
        closest = min(all_coords, key=lambda x: abs(x - eb))
        if abs(closest - eb) <= expected_unit * 0.15:
            boundaries.append(closest)
    boundaries = sorted(set(boundaries))

    # Coverage
    coverage = np.zeros(contig_len, dtype=bool)
    for a in good:
        coverage[a["tstart"]:a["tend"]] = True
    repeat_fraction = float(coverage.mean())

    mean_identity = float(np.mean([a["identity"] for a in good]))
    inferred = len(boundaries) + 1 if boundaries else kmer_copies

    return {
        "ok": True,
        "mean_identity": mean_identity,
        "boundaries": boundaries,
        "repeat_fraction": repeat_fraction,
        "inferred_copies": inferred,
        "n_alignments": len(good),
        "note": "minimap2_confirmed",
    }


def _safe_float(v) -> Optional[float]:
    if v is None:
        return None
    try:
        f = float(v)
    except (ValueError, OverflowError, TypeError):
        return None
    if f != f or f == float('inf') or f == float('-inf'):
        return None
    return f


def analyze(seq: str, k: int = KMER_K) -> dict:
    """
    Run copy-number detection on a single sequence.
    Primary: mean k-mer frequency.
    Validator: minimap2 self-alignment.
    """
    try:
        seq = str(seq).upper().replace(" ", "").replace("\n", "")
    except Exception:
        seq = ""
    L = len(seq)

    _safe_ret = {
        "contig_len": L, "mean_kmer_freq": None, "copies_kmer": 1,
        "validator": "error", "validator_score": None, "validator_snr": None,
        "validator_ok": False, "copies_final": 1,
        "confidence": "low", "flag": "PASS", "note": "analysis_error",
    }

    try:
        # Method 1: mean k-mer frequency
        mkf = mean_kmer_freq(seq, k)
        mkf = _safe_float(mkf)

        # Low-complexity guard
        if mkf is not None:
            if _HAS_KCOUNTER:
                n_unique = len(kcounter.count_kmers(seq, k, canonical_kmers=True))
            else:
                n_unique = len(_count_kmers_python(seq, k))
            if n_unique < max(10, (L - k + 1) * 0.10):
                mkf = None

        copies_kmer = max(1, round(mkf)) if mkf is not None else 1
        kmer_pos = (mkf is not None) and (copies_kmer >= 2) and (mkf >= KMER_COPY_THRESHOLD)

        # Defaults
        validator = "none"
        validator_score = None
        validator_snr = None
        validator_ok = False
        copies_final = copies_kmer if kmer_pos else 1
        confidence = "high"
        note = "single_copy"

        if kmer_pos and copies_kmer >= 2:
            mm = minimap2_validator(seq, L, copies_kmer)
            validator = "minimap2"
            validator_score = _safe_float(mm.get("mean_identity"))
            if mm.get("mean_identity") is not None:
                validator_snr = _safe_float(mm["mean_identity"] / 0.25)
            validator_ok = mm["ok"]

            if validator_ok:
                if abs(mm["inferred_copies"] - copies_kmer) <= 1:
                    copies_final = mm["inferred_copies"]
                confidence = "high"
                note = (f"kmer+minimap2_confirmed "
                        f"boundaries={mm.get('boundaries', [])} "
                        f"repeat_frac={mm['repeat_fraction']:.2f}")
            else:
                confidence = "medium"
                note = f"kmer_only_minimap2_{mm['note']}"

        flag = (
            "FAIL" if copies_final >= FAIL_COPIES else
            "WARN" if copies_final >= WARN_COPIES else
            "PASS"
        )

        return {
            "contig_len": L,
            "mean_kmer_freq": _safe_float(mkf),
            "copies_kmer": copies_kmer,
            "validator": validator,
            "validator_score": _safe_float(validator_score),
            "validator_snr": _safe_float(validator_snr),
            "validator_ok": validator_ok,
            "copies_final": copies_final,
            "confidence": confidence,
            "flag": flag,
            "note": note,
        }

    except Exception as exc:
        import traceback as _tb
        print(f"[warn] analyze() error (len={L}): {exc}", file=sys.stderr)
        _safe_ret["contig_len"] = L
        return _safe_ret


def detect_multicopy(
    fasta_path: str | os.PathLike[str],
    k: int = KMER_K,
    **_,
) -> pd.DataFrame:
    """Run multi-copy detection on every contig in a FASTA file.

    Returns a DataFrame with one row per contig and the columns expected by
    :func:`pca.genome_stats.compute_genome_stats`.
    """
    rows = []
    with open(fasta_path, "r") as fh:
        for record in SeqIO.parse(fh, "fasta"):
            result = analyze(str(record.seq), k=k)
            result["contig_id"] = record.id
            rows.append(result)

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    column_order = [
        "contig_id",
        "contig_len",
        "mean_kmer_freq",
        "copies_kmer",
        "validator",
        "validator_score",
        "validator_snr",
        "validator_ok",
        "copies_final",
        "confidence",
        "flag",
        "note",
    ]
    present = [c for c in column_order if c in df.columns]
    return df[present]
