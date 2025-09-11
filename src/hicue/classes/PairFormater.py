from .utils import *

class PairFormater():
    def __init__(self, positions, separate_by = "", center = "start", overlap = "strict", contact_range = [], separate_regions = "", detrending = "none", diag_mask = 0, min_dist = 0, has_trans = False, circulars = []):
        self._positions = positions
        self._separate_by = separate_by
        self._center = center
        self._overlap = overlap
        self._contact_range = contact_range
        self._separate_regions = separate_regions
        self._detrending = detrending
        self._diag_mask = diag_mask
        self._min_dist = min_dist
        self._has_trans = has_trans
        self._circulars = circulars
    
    def stream_pairs(self, pair_queue, output_queues, threads):
        """Streams the positions and the distance associated with the indexes in the pair_queue."""
        while True:
            try:
                val = pair_queue.get(timeout=10)
            except Empty:
                break
            if val == 'DONE':
                break
            i, j = val
            pos_i = self._positions.loc[i]
            pos_j = self._positions.loc[j]
            is_chrom_i_circular = pos_i["Chromosome"] in self._circulars
            is_chrom_j_circular = pos_j["Chromosome"] in self._circulars
            distance = compute_distance(pos_i, pos_j, self._center) if i != j else 0
            if distance == None:
                for queue in output_queues:
                    queue.put((pos_i, pos_j, i, j, -1, is_chrom_i_circular, is_chrom_j_circular))
            else:
                if distance < 0:
                    i, j = j, i
                    pos_i, pos_j = pos_j, pos_i
                    is_chrom_i_circular, is_chrom_j_circular = is_chrom_j_circular, is_chrom_i_circular
                    distance = abs(distance)
                if distance >= self._min_dist or i == j:
                    for queue in output_queues:
                        queue.put((pos_i, pos_j, i, j, distance, is_chrom_i_circular, is_chrom_j_circular))
        join_queues(output_queues, threads = threads)

    @staticmethod
    def pair_queue_to_df(pair_queue):
        """Creates a pandas DataFrame from a pair queue."""
        pairs = []
        indexes = []
        while True:
            try:
                val = pair_queue.get(timeout=10)
            except Empty:
                break
            if val == "DONE":
                break
            
            index, pair = val
            pairs.append(pair)
            indexes.append(index)
            
        return pd.DataFrame(pairs, index = indexes)
            
    def format_pairs(self, pair_queue, threads = 8):
        """Format pairs into the expected format from the matrix extractor, applying eventual separations."""
        
        to_separate_queue = Queue()
        separated_pairs_queue = Queue()
        formated_pairs_queue = Queue()

        #  initialisation of separators
        separators = schedule_workers(
                    worker_class = "SeparatorScheduler",
                    worker_location = "hicue.workers.Separator",
                    threads = threads,
                    input_queue = to_separate_queue,
                    output_queue = [separated_pairs_queue],
                    separate_by = self._separate_by,
                    center = self._center,
                    overlap = self._overlap, 
                    contact_range = self._contact_range,
                    separate_regions = self._separate_regions,
                    has_trans = self._has_trans
                )        

        # initialisation of pair formaters
        formaters = schedule_workers(
                    worker_class = "PairFormaterScheduler",
                    worker_location = "hicue.workers.PairFormater",
                    threads = threads,
                    input_queue = separated_pairs_queue,
                    output_queue = [formated_pairs_queue],
                    detrending = self._detrending,
                    diag_mask = self._diag_mask
                )   
        
        self.stream_pairs(pair_queue, [to_separate_queue], threads = threads)

        # joining
        join_workers(separators)
        join_queues([separated_pairs_queue], threads =threads)
        join_workers(formaters)
        join_queues([formated_pairs_queue], threads = threads)
        
        formated_pairs_df = self.pair_queue_to_df(formated_pairs_queue)
        
        return formated_pairs_df
