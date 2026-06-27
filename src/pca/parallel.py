"""Process-pool parallel execution helpers."""

from __future__ import annotations

import multiprocessing as mp
import os
import signal
import time
from collections.abc import Callable
from typing import Any

__all__ = ["async_parallel", "init_worker", "terminate_tree"]


def init_worker() -> None:
    """Ignore SIGINT inside pool workers; the parent handles cleanup."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)


def terminate_tree(pid: int, including_parent: bool = True) -> None:
    """Terminate a process and its descendants."""
    import psutil

    parent = psutil.Process(pid)
    for child in parent.children(recursive=True):
        child.terminate()
    if including_parent:
        parent.terminate()


def async_parallel(
    function: Callable[..., Any],
    argument_list: list[list[Any]],
    threads: int,
) -> list[Any]:
    """Run ``function`` over ``argument_list`` using a process pool."""
    pool = mp.Pool(threads, init_worker)
    try:
        results = [pool.apply_async(function, args=args) for args in argument_list]
        pool.close()
        while True:
            if all(r.ready() for r in results):
                return [r.get() for r in results]
            time.sleep(1)
    except KeyboardInterrupt:
        pid = os.getpid()
        terminate_tree(pid)
        raise
