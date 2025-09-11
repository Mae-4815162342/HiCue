from .utils import *

def initialize_globals():
    pass

class AggregatorScheduler(threading.Thread):
    
    def __init__(self, input_queue, output_queues, **subfargs):
        super(AggregatorScheduler, self).__init__()

        self._input_queue = input_queue
        self._subfargs = subfargs

        self.start()

    def run(self):
        aggregator = Aggregator(**self._subfargs)
        while True:
            try:
                val = self._input_queue.get()
            except Empty:
                break
            if val == 'DONE':
                break

            _, window, sep_id, submatrix = val
            aggregator.add(sep_id, window, submatrix)

class Aggregator():
    """Class for aggregating submatrices in a Pileup object."""
    def __init__(self, pileups):
        self._pileups = pileups

    def add(self, sep_id, window,  matrix):
        self._pileups[sep_id].add_submatrix(window, matrix)