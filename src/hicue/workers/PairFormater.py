from .imports import *

# global parameters
global pair_id

pair_lock = threading.Lock()

def initialize_globals():
    with pair_lock:
        global pair_id
        pair_id = 0

class PairFormaterScheduler(threading.Thread):
    def __init__(self, input_queue, output_queue, **pargs):
        super(PairFormaterScheduler, self).__init__()
        self._input_queue = input_queue
        temp_queue = output_queue
        if type(temp_queue) != list:
            temp_queue = [temp_queue]
        self._output_queues = temp_queue
        self._pargs = pargs
        self.start()

    def run(self):
        pair_formater = PairFormater(**self._pargs)
        while True:
            try:
                val = self._input_queue.get(timeout=10)
            except Empty:
                break
            if val == 'DONE':
                break
            
            i, j, distance, circ1, circ2, sep_id = val
            index, formated_pair = pair_formater.format(i, j, distance, circ1, circ2, sep_id)

            for output_queue in self._output_queues:
                output_queue.put((index, formated_pair))
            
class PairFormater():
    def __init__(self, detrending, diag_mask):
        self._detrending = detrending
        self._diag_mask = diag_mask

    def format(self, index1, index2, distance, is_chrom1_circular, is_chrom2_circular, sep_id):
        """Formats a pair of index in the expected matrix extractor format and computes individual index."""
        pair = {
            "Locus1" : index1,
            "Locus2" : index2,
            "Trans": "_trans" in sep_id or sep_id == "trans",
            "Ps": self._detrending == "ps",
            "Diag_mask": self._diag_mask,
            "Distance": distance,
            "Chrom1_circular": is_chrom1_circular,
            "Chrom2_circular": is_chrom2_circular,
            "Sep_id": sep_id
        }

        pair_index = -1
        with pair_lock:
            global pair_id
            pair_index = pair_id
            pair_id += 1

        return pair_index, pair
    