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
            for expected_size, submatrix,tracks in extracter.extract_pair(pair):
                for queue in self._output_queues:
                    queue.put((index, expected_size, pair, submatrix, tracks))

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
            position1 = self._positions.loc[pair['Locus1']]
            position2 = self._positions.loc[pair['Locus2']]
            submatrix = extract_window_region(self._cool_file, 
                                            position1, 
                                            position2,
                                            is_loc_circ1 = bool(pair['Chrom1_circular']),
                                            is_loc_circ2 = bool(pair['Chrom2_circular']),
                                            raw = self._raw)
            
            if submatrix is None:
                print(f"Pair {pair['Locus1']}:{pair['Locus2']} could not be included in pileup due to the following: chromosome ({self._positions.loc[pair['Locus1']]['Chromosome']}) absent from cool.") # TODO write in log
                continue

            subtracks1, subtracks2 = None, None
            if self._tracks:
                subtracks1 = extract_tracks_regions(self._tracks, 
                                    position1,
                                    self._binning,
                                    self._cool_file.chromsizes[position1["Chromosome"]], 
                                    is_loc_circ = bool(pair['Chrom1_circular']))
                
                if pair['Locus1'] != pair['Locus2']:
                    subtracks2 = extract_tracks_regions(self._tracks, 
                                        position2,
                                        self._binning,
                                        self._cool_file.chromsizes[position2["Chromosome"]], 
                                        is_loc_circ = bool(pair['Chrom2_circular']))
                    
            yield expected_size, submatrix, [subtracks1, subtracks2]
    