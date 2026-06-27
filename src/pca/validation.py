"""CLI argument validation and executable availability checks."""

from __future__ import annotations

import argparse
import glob
import os
import shutil
from collections.abc import Iterable
from pathlib import Path

__all__ = ["check_executables", "dbname", "is_valid_dir", "is_valid_file_path"]


def is_valid_dir(path: str | os.PathLike[str]) -> str:
    """Validate an output directory path, creating it if necessary."""
    path_obj = Path(path)
    if path_obj.is_dir():
        return str(path_obj)
    if path_obj.is_file():
        raise argparse.ArgumentTypeError(
            f"ERROR: a file named '{path_obj}' exists; cannot use it as output directory"
        )
    path_obj.mkdir(parents=True, exist_ok=True)
    return str(path_obj)


def is_valid_file_path(path: str | os.PathLike[str]) -> str:
    """Validate that a file exists."""
    path_obj = Path(path)
    if path_obj.is_file():
        return str(path_obj)
    raise argparse.ArgumentTypeError(f"ERROR: {path_obj} is not a valid file")


def check_executables(requirements: Iterable[str]) -> None:
    """Verify that required executables are available on PATH."""
    missing = [program for program in requirements if shutil.which(program) is None]
    if missing:
        raise RuntimeError(
            "Error: required program(s) not executable or not found on $PATH:\n"
            + "\n".join(f"  - {p}" for p in missing)
        )


def dbname(path: str | os.PathLike[str]) -> str:
    """Infer a database name from the first file inside ``path``."""
    files = glob.glob(str(Path(path) / "*"))
    if not files:
        return ""
    return Path(files[0]).stem
