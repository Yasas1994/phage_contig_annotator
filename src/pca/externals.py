"""Wrappers for external bioinformatics binaries.

``prodigal-gv`` and ``hmmsearch`` have been replaced by Python libraries
(``pyrodigal-gv`` and ``pyhmmer``); only ``tRNAscan-SE`` remains as an
external dependency.
"""

from __future__ import annotations

import os
import shlex
import subprocess as sp
from typing import Any

__all__ = ["run_trnascan"]


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
