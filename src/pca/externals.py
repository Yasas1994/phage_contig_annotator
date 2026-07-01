"""Wrappers for external bioinformatics binaries.

``prodigal-gv`` and ``hmmsearch`` have been replaced by Python libraries
(``pyrodigal-gv`` and ``pyhmmer``); ``tRNAscan-SE``, Tandem Repeats Finder
(``trf``), PHANOTATE, DefenseFinder, CRISPRCasFinder and MinCED are the
remaining external dependencies.
"""

from __future__ import annotations

import logging
import os
import shlex
import shutil
import subprocess as sp
import tempfile
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "run_crispr_cas_finder",
    "run_defensefinder",
    "run_minced",
    "run_phanotate",
    "run_trf",
    "run_trnascan",
]


def _phanotate_executable() -> str:
    """Return the PHANOTATE command available on PATH, preferring ``phanotate.py``."""
    for cmd in ("phanotate.py", "phanotate"):
        if shutil.which(cmd):
            return cmd
    return "phanotate.py"


def _crispr_cas_finder_executable() -> str | None:
    """Return the CRISPRCasFinder executable on PATH, or None if missing."""
    for cmd in ("CRISPRCasFinder.pl", "CRISPRCasFinder", "crisprcasfinder"):
        if shutil.which(cmd):
            return cmd
    return None


def run_minced(out: str, in_: str) -> bool:
    """Run MinCED on a nucleotide FASTA file.

    MinCED detects CRISPR arrays and writes a table file to ``{out}.crisprs``
    and a GFF3 file to ``{out}.gff``. The ``out`` argument is treated as a
    prefix.
    """
    if shutil.which("minced") is None:
        logger.error("minced not found on PATH; skipping CRISPR array detection")
        return False

    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    table_out = f"{out}.crisprs"
    gff_out = f"{out}.gff"
    cmd = ["minced", in_, table_out, gff_out]
    return _run_tool(cmd, f"{out}.log", f"{out}.cmd")


def _run_tool(
    cmd: list[str],
    log_path: str | os.PathLike[str],
    cmd_path: str | os.PathLike[str] | None = None,
    out_path: str | os.PathLike[str] | None = None,
    **popen_kwargs: Any,
) -> bool:
    """Run an external command, logging stdout/stderr to ``log_path``.

    If ``out_path`` is given, command stdout is written there; otherwise it is
    merged into the log.
    """
    if cmd_path is not None:
        with open(cmd_path, "w") as fh:
            fh.write(" ".join(shlex.quote(part) for part in cmd) + "\n")

    try:
        if out_path is not None:
            with open(out_path, "w") as stdout_target, open(log_path, "w") as stderr_target:
                proc = sp.Popen(
                    cmd, stdout=stdout_target, stderr=stderr_target, **popen_kwargs
                )
                return proc.wait() == 0

        with open(log_path, "w") as log_fh:
            proc = sp.Popen(cmd, stdout=log_fh, stderr=sp.STDOUT, **popen_kwargs)
            return proc.wait() == 0
    except FileNotFoundError as exc:
        with open(log_path, "w") as log_fh:
            log_fh.write(f"Command not found: {exc}\n")
        return False


def run_crispr_cas_finder(out: str, in_: str, threads: int = 1) -> bool:
    """Run CRISPRCasFinder on a nucleotide FASTA file.

    CRISPRCasFinder detects CRISPR arrays and Cas genes. It writes a results
    directory under ``out`` containing a GFF, a summary CSV and additional
    auxiliary files. This wrapper runs the tool with ``-so`` 1 (output
    orientation) and ``-plotNoTOC`` to suppress the HTML table of contents, and
    copies the resulting GFF to ``{out}.gff`` for easier downstream access.
    """
    executable = _crispr_cas_finder_executable()
    if executable is None:
        logger.error("CRISPRCasFinder not found on PATH; skipping CRISPR-Cas detection")
        return False

    out_path = Path(out)
    out_path.mkdir(parents=True, exist_ok=True)
    cmd = [
        executable,
        "-i", in_,
        "-out", str(out_path),
        "-so", "1",
        "-plotNoTOC",
    ]
    if threads > 1:
        cmd.extend(["-cpu", str(threads)])

    success = _run_tool(cmd, f"{out}.log", f"{out}.cmd")
    if not success:
        return False

    # Copy the GFF next to the output directory for stable downstream access.
    gff_candidates = list(out_path.glob("*.gff")) + list(out_path.glob("*CRISPR*gff"))
    target_gff = out_path / "crisprcasfinder.gff"
    if gff_candidates and not target_gff.exists():
        shutil.copy2(gff_candidates[0], target_gff)
    return True


def run_defensefinder(out: str, in_: str, threads: int = 1) -> bool:
    """Run DefenseFinder on a protein FASTA file.

    DefenseFinder detects anti-phage defense systems using profile HMMs and
    system-specific rules. It writes files named ``<input>_defense_finder_*.tsv``
    under ``out``; this wrapper renames them to stable ``defense_finder_*.tsv``
    names so the rest of the pipeline can depend on a fixed path.
    """
    out_path = Path(out)
    out_path.mkdir(parents=True, exist_ok=True)
    in_path = Path(in_)
    cmd = [
        "defense-finder",
        "run",
        "-w", str(threads),
        "--out-dir", str(out_path),
        in_,
    ]
    success = _run_tool(cmd, f"{out}.log", f"{out}.cmd")
    if not success:
        return False

    # DefenseFinder names outputs after the input basename (e.g.
    # proteins_defense_finder_genes.tsv). Rename to stable names.
    expected = {
        f"{in_path.stem}_defense_finder_genes.tsv": "defense_finder_genes.tsv",
        f"{in_path.stem}_defense_finder_systems.tsv": "defense_finder_systems.tsv",
    }
    for source_name, target_name in expected.items():
        source = out_path / source_name
        target = out_path / target_name
        if source.exists() and not target.exists():
            source.rename(target)
        elif source.exists() and target.exists():
            # Should not happen, but prefer the just-produced file.
            target.unlink()
            source.rename(target)

    return (out_path / "defense_finder_genes.tsv").exists()


def _run_trf_single(out: str, in_: str, flags: list[str] | None = None) -> Path:
    """Run TRF on a single FASTA file and return the path to the ``.dat`` output."""
    in_path = Path(in_).resolve()
    out_path = Path(f"{out}.dat")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    params = ["2", "7", "7", "80", "10", "50", "2000"]
    dot_params = ".".join(params)
    trf_flags = ["-d", "-h"] if flags is None else flags
    cmd = ["trf", "{input}", *params, *trf_flags]

    with tempfile.TemporaryDirectory(dir=str(out_path.parent)) as run_dir:
        local_fasta = Path(run_dir) / in_path.name
        shutil.copy2(in_path, local_fasta)
        cmd[1] = str(local_fasta)
        success = _run_tool(
            cmd, f"{out}.log", f"{out}.cmd", cwd=str(run_dir)
        )
        expected = Path(run_dir) / f"{in_path.name}.{dot_params}.dat"
        if success and expected.exists():
            expected.rename(out_path)
        elif success and not out_path.exists():
            out_path.touch()

    return out_path


def run_trf(out: str, in_: str, threads: int = 1, flags: list[str] | None = None) -> bool:
    """Run Tandem Repeats Finder (TRF) on a nucleotide FASTA file.

    ``out`` is treated as a prefix; the final TRF table is written to
    ``{out}.dat``. TRF writes its output next to the input file, so the input
    is copied into the output directory, the command is run there, and the
    resulting ``.dat`` file is moved to ``{out}.dat``. If TRF produces no
    repeats and therefore no ``.dat`` file, an empty ``{out}.dat`` is created
    so downstream steps have a stable input.

    When ``threads`` is greater than 1 and the input contains multiple
    sequences, the input is split into at most ``threads`` chunks and TRF is
    run on each chunk in parallel. The resulting ``.dat`` files are concatenated
    to produce the final output.

    ``flags`` controls the TRF output format. The default is ``["-d", "-h"]``
    (data file, no HTML). Pass ``["-f", "-d", "-m"]`` to produce flanking
    sequence and masked sequence output in addition to the ``.dat`` table.
    """
    from pca.io import read_fasta, split_fasta
    from pca.parallel import async_parallel

    out_path = Path(f"{out}.dat")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    n_records = sum(1 for _ in read_fasta(in_))
    n_chunks = max(1, min(threads, n_records))
    if n_chunks == 1:
        return _run_trf_single(out, in_).exists()

    chunk_dir = out_path.parent / f"{out_path.stem}_chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)
    chunk_paths = split_fasta(in_, chunk_dir, n_chunks)

    chunk_outs = [str(chunk_dir / f"chunk_{i}") for i in range(len(chunk_paths))]
    args = [
        [chunk_out, str(chunk_path), flags]
        for chunk_out, chunk_path in zip(chunk_outs, chunk_paths)
    ]
    async_parallel(_run_trf_single, args, threads=n_chunks)

    with open(out_path, "w") as out_fh:
        for chunk_out in chunk_outs:
            chunk_dat = Path(f"{chunk_out}.dat")
            if chunk_dat.exists() and chunk_dat.stat().st_size > 0:
                content = chunk_dat.read_text()
                out_fh.write(content)
                if not content.endswith("\n"):
                    out_fh.write("\n")

    shutil.rmtree(chunk_dir, ignore_errors=True)
    return out_path.exists()


def run_phanotate(in_: str, out: str, threads: int = 1) -> bool:
    """Run PHANOTATE on a nucleotide FASTA file and write a GFF3 file.

    PHANOTATE is a gene caller specialized for phage genomes. The wrapper
    invokes ``phanotate.py`` (or ``phanotate`` as a fallback) and requests
    GFF3 output.
    """
    out_path = Path(out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        _phanotate_executable(),
        in_,
        "-o", str(out_path),
        "-f", "gff3",
    ]
    return _run_tool(cmd, f"{out}.log", f"{out}.cmd")


def run_trnascan(out: str, in_: str, threads: int = 1) -> bool:
    """Run tRNAscan-SE on a nucleotide FASTA file."""
    cmd = [
        "tRNAscan-SE",
        "-G",
        "-p", "meta",
        "-o", f"{out}.tsv",
        "-j", f"{out}.gff",
        "-i", in_,
        "--thread", str(threads),
    ]
    return _run_tool(cmd, f"{out}.log", f"{out}.cmd")
