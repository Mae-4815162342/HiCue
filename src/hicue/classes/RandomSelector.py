from .utils import *

class RandomSelector():
    def __init__(self, center = "start", selection_window = 100000, nb_rand_per_pos = 1):
        self._center = center
        self._nb_rand_per_pos = nb_rand_per_pos
        self._selection_window = selection_window
        
    @staticmethod
    def stream_positions(positions, queue):
        for i, pos in positions.iterrows():
            queue.put((i, pos))
    
    def select_randoms(self, positions, threads = 8):
        
        positions_queue = Queue()
        random_positions_queue = Queue()

        #  initialisation of selectors
        workers = schedule_workers(
                    worker_class = "RandomSelectorScheduler",
                    worker_location = "hicue.workers.RandomSelector",
                    threads = threads,
                    input_queue = positions_queue,
                    output_queues = [random_positions_queue],
                    nb_pos = len(positions),
                    center = self._center,
                    nb_rand_per_pos = self._nb_rand_per_pos,
                    selection_window = self._selection_window,
                )        

        self.stream_positions(positions, positions_queue)

        # joining
        join_queues([positions_queue], threads = threads)
        join_workers(workers)
        join_queues([random_positions_queue], threads = threads)
        
        random_positions = position_queue_to_df(random_positions_queue)
        
        return random_positions