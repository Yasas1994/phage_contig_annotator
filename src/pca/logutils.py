"""Logging helpers for the phage contig annotator pipeline."""

from __future__ import annotations

import logging
import sys

__all__ = ["get_logger"]


def get_logger(quiet: bool = False) -> logging.Logger:
    """Configure and return the module logger."""
    log = logging.getLogger("pca")
    log.setLevel(logging.WARNING if quiet else logging.INFO)
    formatter = logging.Formatter(
        fmt="%(asctime)s : %(levelname)s : %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
    )
    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)
    log.handlers.clear()
    log.addHandler(stream_handler)
    return log
