"""Run independent, self-contained analysis steps concurrently.

Used by the Data Analysis and Heatmap upload paths so the post-load work on a
merged dataframe — protein-info extraction, per-protein abundance totals,
per-function peptide counts — happens simultaneously instead of one step after
another. On large files that sequential chain is the bulk of the upload wait.

Threads, not processes, on purpose:

* The heavy steps are pandas / numpy vectorised operations (``explode``,
  ``groupby``, ``to_numeric``, ``value_counts``) whose inner loops release the
  GIL, so real wall-clock overlap is available without separate interpreters.
* Every task reads the same in-memory dataframe. Handing that to worker
  *processes* would pickle a multi-hundred-MB frame once per task — reliably
  more expensive than the parallelism saves for a single request.
* Django serves these requests on a worker thread already; a bounded local pool
  keeps the extra concurrency predictable.
"""
from __future__ import annotations

import concurrent.futures
import os

# Small pool: there are only ever a handful of tasks and the box also has to
# serve other requests. Override with PEPTILINE_ANALYSIS_WORKERS if needed.
try:
    _MAX_WORKERS = max(1, int(os.environ.get("PEPTILINE_ANALYSIS_WORKERS", "4")))
except ValueError:
    _MAX_WORKERS = 4


def run_tasks(tasks: dict, max_workers: int | None = None) -> dict:
    """Execute ``{key: zero-arg callable}`` concurrently; return ``{key: result}``.

    - 0 or 1 tasks run inline (no pool spun up).
    - Exceptions propagate: the first task to raise wins and its exception is
      re-raised here, so callers behave exactly as with a sequential call.
    """
    if not tasks:
        return {}
    if len(tasks) == 1:
        (key, fn), = tasks.items()
        return {key: fn()}

    workers = min(max_workers or _MAX_WORKERS, len(tasks))
    results: dict = {}
    with concurrent.futures.ThreadPoolExecutor(
        max_workers=workers, thread_name_prefix="peptiline-analysis"
    ) as executor:
        future_to_key = {executor.submit(fn): key for key, fn in tasks.items()}
        for future in concurrent.futures.as_completed(future_to_key):
            # .result() re-raises inside the caller's thread; the context
            # manager still joins the rest before we leave.
            results[future_to_key[future]] = future.result()
    return results
