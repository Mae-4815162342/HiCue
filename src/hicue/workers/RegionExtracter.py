from hicue.utils import *

def initialize_globals():
    pass

class RegionExtracterScheduler(threading.Thread):
    
    def __init__(self, input_queue, output_queues, **extargs):
        super(RegionExtracterScheduler, self).__init__()

        self._input_queue = input_queue
        self._output_queues = output_queues
        self._extargs = extargs

        self.start()

    def run(self):
        extracter = RegionExtracter(**self._extargs)
        while True:
            try:
                val = self._input_queue.get()
            except Empty:
                break
            if val == 'DONE':
                break

            index, pair = val
            for expected_size, submatrix in extracter.extract_pair(pair):
                for queue in self._output_queues:
                    queue.put((index, expected_size, pair, submatrix))

class RegionExtracter():
    """Class for submatrix extraction from a pair of positions."""
    def __init__(self, cool_file, positions, expected_sizes, tracks = None, raw = False, random = False):
        self._cool_file = cool_file
        self._positions = positions
        self._tracks = pyBigWig.open(tracks) if tracks is not None and type(tracks) == str else tracks
        self._binning = cool_file.binsize
        self._expected_sizes = expected_sizes
        self._raw = raw
        self._random = random
        
    def extract_pair(self, pair):
        """Calls the extract function on the pair for each expected size. Yields a square submatrix as a numpy array."""
        for expected_size in self._expected_sizes:
            submatrix = extract_window_region(self._cool_file, 
                                            self._positions.loc[pair['Locus1']], 
                                            self._positions.loc[pair['Locus2']],
                                            is_loc_circ1 = bool(pair['Chrom1_circular']),
                                            is_loc_circ2 = bool(pair['Chrom2_circular']),
                                            raw = self._raw)
            submatrix = resize_window(submatrix, expected_size = expected_size)
            
            if submatrix is None:
                print(f"Pair {pair['Locus1']}:{pair['Locus2']} could not be included in pileup due to the following: chromosome ({self._positions.loc[pair['Locus1']]['Chromosome']}) absent from cool.") # TODO write in log
                continue
            # if self._tracks: # TODO: to re-write
            #     subtracks1 = extract_tracks(self._tracks, 
            #                         self._positions.loc[pair['Locus1']],
            #                         self._binning,
            #                         resizing, 
            #                         is_loc_circ = bool(pair['Chrom1_circular']),
            #                         center = self._center)
            #     submatrix = np.concatenate([submatrix, [subtracks1]], axis = 0)
            #     if pair['Locus1'] != pair['Locus2']:
            #         subtracks2 = extract_tracks(self._tracks, 
            #                             self._positions.loc[pair['Locus2']],
            #                             self._binning,
            #                             resizing, 
            #                             is_loc_circ = bool(pair['Chrom2_circular']),
            #                             center = self._center)
            #         submatrix = np.concatenate([submatrix, [subtracks2]], axis = 0)
            yield expected_size, submatrix
    