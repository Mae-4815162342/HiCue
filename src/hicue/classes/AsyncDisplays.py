"""
AsyncDisplays.py
================
Thread-based asynchronous display workers for HiCue.

These classes wrap async matplotlib display functions inside ``threading.Thread``
workers so that figure rendering can be parallelised with the rest of the
extraction pipeline without blocking the main thread.

.. warning::
    ``matplotlib.use('Agg')`` **must** be called before importing any other
    matplotlib module (it is set in ``imports.py``).  The non-interactive Agg
    backend is mandatory here because these threads do not have access to a
    display server and rendering must be purely file-based.

Classes
-------
Display
    Generic single-item display worker.
DisplayBatch
    Batching display worker that accumulates items before rendering them
    together in a single figure.

Bug-fix – matplotlib non-thread-safe state
------------------------------------------
``asyncio.run()`` creates and destroys a new event loop on every call.
Inside a long-lived thread this causes a teardown/creation race with
``plt.close()`` called at the end of the previous coroutine, leading to:

* blank or partially-drawn output figures,
* wrong data rendered on a figure (corrupted "current figure" singleton),
* occasional ``RuntimeError`` from the Agg renderer.

Fix: each Display/DisplayBatch worker creates **one** event loop
(``self._loop = asyncio.new_event_loop()``) that persists for the lifetime of
the thread.  ``loop.run_until_complete(coro)`` is used instead of
``asyncio.run(coro)``.  The loop is closed in an overridden ``join()`` so
resources are always released.
"""

from hicue.imports import *


def initialize_globals():
    """Initialise module-level globals (no-op placeholder for the worker loader)."""
    pass


class Display(threading.Thread):
    """Generic asynchronous display worker.

    Reads parameter dicts from *input_queue*, calls *function* with each set
    of parameters (merged with *static_params*), and forwards every raw value
    it receives to each queue in *output_queues* **before** rendering so that
    downstream pipeline stages are not stalled by slow matplotlib rendering.

    The sentinel value ``"DONE"`` causes the worker to stop.

    Parameters
    ----------
    input_queue : queue.Queue
        Source of ``dict`` parameter objects.  ``"DONE"`` signals end-of-stream.
    output_queues : list[queue.Queue]
        Queues to which every value (including ``"DONE"``) is forwarded before
        processing.
    function : coroutine function
        An ``async def`` function that produces a matplotlib figure and saves
        it to disk.  It receives the merged parameters from *input_queue* and
        *static_params*.
    **static_params
        Additional keyword arguments forwarded to *function* on every call
        (e.g. ``positions``, ``chromsizes``).

    Notes
    -----
    The thread starts immediately on construction.
    """

    def __init__(self, input_queue, output_queues, function, **static_params):
        super(Display, self).__init__()
        self._input_queue = input_queue
        self._output_queues = output_queues
        self._func = function
        self._statics = static_params
        # Bug-fix: one persistent event loop per worker thread instead of
        # creating/destroying a loop on every asyncio.run() call.
        self._loop = asyncio.new_event_loop()
        self.start()

    def call_display(self, **params):
        """Execute the display coroutine synchronously on this thread's event loop.

        Parameters
        ----------
        **params
            Per-item parameters merged with ``self._statics`` and forwarded to
            the wrapped coroutine.
        """
        # Bug-fix: run_until_complete reuses the existing loop; asyncio.run()
        # would tear it down after each call and race with plt.close().
        self._loop.run_until_complete(self._func(**params, **self._statics))

    def run(self):
        """Main loop: consume queue items and render each one.

        Terminates when the ``"DONE"`` sentinel is received or the queue raises
        ``Empty``.
        """
        while True:
            try:
                value = self._input_queue.get()
                for queue in self._output_queues:
                    queue.put(value)
            except Empty:
                break
            if value == "DONE":
                break
            params = value
            self.call_display(**params)

    def join(self, timeout=None):
        """Join the thread and close the event loop.

        Parameters
        ----------
        timeout : float or None
            Passed through to ``threading.Thread.join``.
        """
        super().join(timeout)
        self._loop.close()


class DisplayBatch(threading.Thread):
    """Batching asynchronous display worker.

    Accumulates items with the same ``(window, sep_id)`` key into batches of
    size *batch_size* and renders each full batch as a single multi-panel
    figure.  Partial batches remaining at end-of-stream are flushed during
    the cleanup loop so no submatrix is silently dropped.

    Parameters
    ----------
    input_queue : queue.Queue
        Source of parameter dicts.  Each dict must contain at least ``window``
        and ``pair["Sep_id"]`` keys.  ``"DONE"`` signals end-of-stream.
    output_queues : list[queue.Queue]
        Queues to which every raw value is forwarded before processing.
    function : coroutine function
        An ``async def`` function that renders a batch of submatrices.
    batch_size : int
        Number of submatrices to accumulate before triggering a render.
    params_to_batch : list[str], optional
        Parameter names that should be promoted to *static_params* on first
        encounter (e.g. ``["track_unit"]``), avoiding redundancy in every dict.
    **static_params
        Additional keyword arguments forwarded to *function* on every batch
        call.

    Notes
    -----
    The thread starts immediately on construction.
    """

    def __init__(
        self,
        input_queue,
        output_queues,
        function,
        batch_size,
        params_to_batch=[],
        **static_params,
    ):
        super(DisplayBatch, self).__init__()
        self._input_queue = input_queue
        self._output_queues = output_queues
        self._func = function
        self._batch_size = batch_size
        self._current_batch = {}
        self._statics = static_params
        self._params_to_batch = params_to_batch
        self._outfolders = {}
        self._keys = []
        self._nb_batch = {}
        # Bug-fix: same persistent-loop strategy as Display.
        self._loop = asyncio.new_event_loop()
        self.start()

    # ------------------------------------------------------------------
    # Batch management helpers
    # ------------------------------------------------------------------

    def add_to_batch(self, size_metric, sep_id, **params):
        """Register a new submatrix in the batch for ``(window, sep_id)``.

        If a key listed in *params_to_batch* is found in *params* and is not
        yet present in ``self._statics``, it is promoted so the rendering
        function can access it.

        Parameters
        ----------
        size_metric : int
            Extraction expected or window size in pixels or base-pairs respectively.
        sep_id : str
            Separation group identifier.
        **params
            Per-item display parameters, must include ``outfolder``.
        """
        for param in self._params_to_batch:
            if param not in self._statics and param in params:
                self._statics[param] = params[param]
        key = f"{size_metric}.{sep_id}"
        if key not in self._current_batch:
            self._current_batch[key] = []
            self._nb_batch[key] = 1
            self._outfolders[key] = params["outfolder"]
            self._keys.append((size_metric, sep_id))
        self._current_batch[key].append(params)

    def empty_batch(self, size_metric, sep_id):
        """Clear accumulated items for ``(window, sep_id)`` after rendering.

        Parameters
        ----------
        size_metric : int
            Extraction expected or window size in pixels or base-pairs respectively.
        sep_id : str
            Separation group identifier.
        """
        self._current_batch[f"{size_metric}.{sep_id}"] = []

    def get_current_batch_size(self, size_metric, sep_id):
        """Return the number of items currently queued for ``(window, sep_id)``.

        Parameters
        ----------
        size_metric : int
            Extraction expected or window size in pixels or base-pairs respectively.
        sep_id : str
            Separation group identifier.

        Returns
        -------
        int
        """
        return len(self._current_batch[f"{size_metric}.{sep_id}"])

    def is_empty_batch(self, size_metric, sep_id):
        """Return ``True`` if there are no pending items for ``(window, sep_id)``.

        Parameters
        ----------
        size_metric : int
            Extraction expected or window size in pixels or base-pairs respectively.
        sep_id : str
            Separation group identifier.

        Returns
        -------
        bool
        """
        return self.get_current_batch_size(size_metric, sep_id) == 0

    # ------------------------------------------------------------------
    # Rendering
    # ------------------------------------------------------------------

    def call_batch_display(self, size_metric, sep_id):
        """Render the current batch for ``(size_metric, sep_id)`` and save it.

        Parameters
        ----------
        size_metric : int
            Extraction expected or window size in pixels or base-pairs respectively.
        sep_id : str
            Separation group identifier.
        """
        key = f"{size_metric}.{sep_id}"
        # Bug-fix: run_until_complete on the persistent loop avoids the
        # teardown/creation overhead and the plt.close() race condition.
        self._loop.run_until_complete(
            self._func(
                self._current_batch[key],
                batch_size=self._batch_size,
                outfolder=self._outfolders[key],
                title=f"batch#{self._nb_batch[key]}",
                size_metric=size_metric,
                **self._statics,
            )
        )

    # ------------------------------------------------------------------
    # Thread main loop
    # ------------------------------------------------------------------

    def run(self):
        """Consume queue items, accumulating them into batches by ``(window, sep_id)`` or ``(expected_size, sep_id)``.

        When a batch reaches *batch_size* it is rendered immediately.  After
        the ``"DONE"`` sentinel, any remaining non-empty partial batches are
        flushed so that all submatrices produce output.
        """
        while True:
            try:
                value = self._input_queue.get()
                for queue in self._output_queues:
                    queue.put(value)
            except Empty:
                break
            if value == "DONE":
                break

            size_metric = value["size_metric"]

            sep_id = value["pair"]["Sep_id"]
            self.add_to_batch(**value, sep_id=sep_id)
            if self.get_current_batch_size(size_metric, sep_id) == self._batch_size:
                self.call_batch_display(size_metric, sep_id)
                self.empty_batch(size_metric, sep_id)
                self._nb_batch[f"{size_metric}.{sep_id}"] += 1

        # Flush partial batches that never reached batch_size.
        for size_metric, sep_id in self._keys:
            if not self.is_empty_batch(size_metric, sep_id):
                self.call_batch_display(size_metric, sep_id)

    def join(self, timeout=None):
        """Join the thread and release the event loop.

        Parameters
        ----------
        timeout : float or None
            Passed through to ``threading.Thread.join``.
        """
        super().join(timeout)
        self._loop.close()