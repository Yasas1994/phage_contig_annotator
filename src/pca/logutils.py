"""Logging helpers for the phage contig annotator pipeline."""

from __future__ import annotations

import logging
import sys

__all__ = ["get_logger"]


_LEVEL_COLORS: dict[int, str] = {
    logging.DEBUG: "\033[36m",       # cyan
    logging.INFO: "\033[32m",        # green
    logging.WARNING: "\033[33m",     # yellow
    logging.ERROR: "\033[31m",       # red
    logging.CRITICAL: "\033[1;31m",  # bold red
}
_RESET = "\033[0m"


class _ColoredFormatter(logging.Formatter):
    """Formatter that colorizes the level name when writing to a TTY."""

    def __init__(self, fmt: str, datefmt: str | None = None) -> None:
        super().__init__(fmt=fmt, datefmt=datefmt)
        self._use_color = sys.stderr.isatty()

    def format(self, record: logging.LogRecord) -> str:
        levelno = record.levelno
        if self._use_color and levelno in _LEVEL_COLORS:
            record.levelname = (
                f"{_LEVEL_COLORS[levelno]}{record.levelname}{_RESET}"
            )
        return super().format(record)


def get_logger(quiet: bool = False) -> logging.Logger:
    """Configure and return the module logger."""
    log = logging.getLogger("pca")
    log.setLevel(logging.WARNING if quiet else logging.INFO)
    formatter = _ColoredFormatter(
        fmt="%(asctime)s : %(levelname)s : %(message)s",
        datefmt="%d-%b-%y %H:%M:%S",
    )
    stream_handler = logging.StreamHandler(sys.stderr)
    stream_handler.setFormatter(formatter)
    log.handlers.clear()
    log.addHandler(stream_handler)
    return log
