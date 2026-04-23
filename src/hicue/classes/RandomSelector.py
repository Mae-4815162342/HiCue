from hicue.utils import *

class RandomSelector():
    def __init__(self, positions, center = "start", selection_window = 100000, nb_rand_per_pos = 1, random_jitter = 0, is_region = False, padding = None):
        self._positions = positions
        self._center = center
        self._nb_rand_per_pos = nb_rand_per_pos
        self._selection_window = selection_window
        self._jitter = random_jitter
        self._is_region = is_region
        self._padding = padding
        
    @staticmethod
    def stream_pairs(pairs, queue):
        for i, pair in pairs.iterrows():
            queue.put((i, pair))
    
    def select_randoms(self, pairs, threads = 8):
        
        pairs_queue = Queue()
        random_positions_queue = Queue()
        random_pairs_queue = Queue()

        #  initialisation of selectors
        workers = schedule_workers(
                    worker_class = "RandomSelectorScheduler",
                    worker_location = "hicue.workers.RandomSelector",
                    threads = threads,
                    input_queue = pairs_queue,
                    output_queues = [random_positions_queue, random_pairs_queue],
                    nb_pairs = np.max(pairs.index) + 1,
                    positions = self._positions,
                    center = self._center,
                    nb_rand_per_pos = self._nb_rand_per_pos,
                    selection_window = self._selection_window,
                    jitter = self._jitter,
                    is_region = self._is_region,
                    padding = self._padding
                )        

        self.stream_pairs(pairs, pairs_queue)

        # joining
        join_queues([pairs_queue], threads = threads)
        join_workers(workers)
        join_queues([random_positions_queue], threads = threads)
        
        random_positions = position_queue_to_df(random_positions_queue)
        random_pairs = position_queue_to_df(random_pairs_queue)
        
        return random_positions, random_pairs