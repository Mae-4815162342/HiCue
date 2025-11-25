from hicue.utils import *

def initialize_globals():
    pass

class ExtracterScheduler(threading.Thread):
    
    def __init__(self, input_queue, output_queues, **extargs):
        super(ExtracterScheduler, self).__init__()

        self._input_queue = input_queue
        self._output_queues = output_queues
        self._extargs = extargs

        self.start()

    def run(self):
        extracter = Extracter(**self._extargs)
        while True:
            try:
                val = self._input_queue.get()
            except Empty:
                break
            if val == 'DONE':
                break

            index, pair = val
            for window, submatrix in extracter.extract_pair(pair):
                for queue in self._output_queues:
                    queue.put((index, window, pair, submatrix))

class Extracter():
    """Class for submatrix extraction from a pair of positions."""
    def __init__(self, cool_file, positions, windows, tracks = None, center = "start", raw = False, random = False):
        self._cool_file = cool_file
        self._positions = positions
        self._tracks = pyBigWig.open(tracks) if tracks is not None and type(tracks) == str else tracks
        self._binning = cool_file.binsize
        self._windows = windows
        self._center = center
        self._raw = raw
        self._random = random
        
    def extract_pair(self, pair):
        """Calls the extract function on the pair for each window size. Yields a square submatrix as a numpy array."""
        for window in self._windows:
            submatrix = extract_window(self._cool_file, 
                                    self._positions.loc[pair['Locus1']], 
                                    self._positions.loc[pair['Locus2']], 
                                    self._binning,
                                    window, 
                                    is_loc1_circ = bool(pair['Chrom1_circular']),
                                    is_loc2_circ = bool(pair['Chrom2_circular']),
                                    center = self._center, 
                                    raw = self._raw)
            if submatrix is None:
                print(f"Pair {pair['Locus1']}:{pair['Locus2']} could not be included in pileup for window {window} due to the following: chromosome ({self._positions.loc[pair['Locus1']]['Chromosome']}) absent from cool.") # TODO write in log
                continue
            if self._tracks:
                subtracks1 = extract_tracks(self._tracks, 
                                    self._positions.loc[pair['Locus1']],
                                    self._binning,
                                    window, 
                                    is_loc_circ = bool(pair['Chrom1_circular']),
                                    center = self._center)
                submatrix = np.concatenate([submatrix, [subtracks1]], axis = 0)
                if pair['Locus1'] != pair['Locus2']:
                    subtracks2 = extract_tracks(self._tracks, 
                                        self._positions.loc[pair['Locus2']],
                                        self._binning,
                                        window, 
                                        is_loc_circ = bool(pair['Chrom2_circular']),
                                        center = self._center)
                    submatrix = np.concatenate([submatrix, [subtracks2]], axis = 0)
            yield window, submatrix
    