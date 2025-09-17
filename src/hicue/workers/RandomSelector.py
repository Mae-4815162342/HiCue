from .utils import *

def initialize_globals():
    pass

class RandomSelectorScheduler(threading.Thread):
    
    def __init__(self, input_queue, output_queues, **selargs):
        super(RandomSelectorScheduler, self).__init__()

        self._input_queue = input_queue
        self._output_queues = output_queues
        self._selargs = selargs

        self.start()

    def run(self):
        selecter = RandomSelector(self._output_queues, **self._selargs)
        while True:
            try:
                val = self._input_queue.get()
            except Empty:
                break
            if val == 'DONE':
                break

            index, position = val
            selecter.select_positions(index, position)    

class RandomSelector():
    """Class implementing the random selection of positions from a position."""
    def __init__(self, output_queues, nb_pos, center = "start", selection_window = 100000, nb_rand_per_pos = 1):
        self._random_queues = output_queues
        self._nb_pos = nb_pos
        self._center = center
        self._nb_rand_per_pos = nb_rand_per_pos
        self._selection_window = selection_window
        
    def select_positions(self, position_index, position):
        """Randomly creates positions from an original position. 
        Will create nb_rand_per_pos in the selection_window (+- around the start/center/stop).
        If provided, the sign of the original position is kept for the selection. Else, default +.
        The index is computed with the formula k*nb_rand_per_pos + index, permiting the retrival of the original index with % nb_pos_per_rand."""
        
        start = position["Start"]
        stop = position["End"]
        
        center = start
        match self._center:
            case "start":
                center = start
            case "stop":
                center = stop
            case "center":
                center = (start + stop) // 2
                
        interval = [center - self._selection_window, center + self._selection_window]
        
        for k in range(self._nb_rand_per_pos):
            index = k * self._nb_pos + position_index
            random_start = int((random.random() * (interval[1] - interval[0]) + interval[0]) // 1)
            random_position = {
                "Chromosome": position["Chromosome"],
                "Start": random_start,
                "End": random_start,
                "Name": f"rand_pos_{index}",
                "Strand": position["Strand"]
            }
            
            for i in range(len(self._random_queues)):
                self._random_queues[i].put((index, random_position))