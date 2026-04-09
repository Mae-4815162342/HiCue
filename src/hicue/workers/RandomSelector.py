from hicue.utils import *

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

            index, pair = val
            selecter.select_positions(index, pair)

class RandomSelector():
    """Class implementing the random selection of positions from a position."""
    def __init__(self, output_queues, positions, nb_pairs, center = "start", selection_window = 100000, nb_rand_per_pos = 1, jitter = 0):
        self._positions = positions
        self._random_queues = output_queues
        self._nb_pairs = nb_pairs
        self._center = center
        self._nb_rand_per_pos = nb_rand_per_pos
        self._selection_window = selection_window
        self._jitter = jitter

    def _get_center(self, position):
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
        return center
        
    def select_positions(self, pair_index, pair):
        """Randomly creates positions from an original pair of positions. If the contact is diagonal (a position with itself), only one position is created, 2 otherwise.
        Will create nb_rand_per_pos in the selection_window (+- around the start/center/stop).
        If provided, the sign of the original position is kept for the selection. Else, default +.
        The index is computed with the formula k*nb_rand_per_pos + index, permiting the retrieval of the original index with % nb_pos_per_rand."""
        
        locus1 = self._positions.loc[pair['Locus1']]
        locus2 = self._positions.loc[pair['Locus2']]

        center1 = self._get_center(locus1)
        distance = pair["Distance"]

        interval1 = [center1 - self._selection_window, center1 + self._selection_window]

        for k in range(self._nb_rand_per_pos):
            
            # first random position
            index1 =  (k * self._nb_pairs + pair_index) * 2
            random_start = int((random.random() * (interval1[1] - interval1[0]) + interval1[0]) // 1)
            random_position1 = {
                "Chromosome": locus1["Chromosome"],
                "Start": random_start,
                "End": random_start,
                "Name": f"rand_pos_{index1}",
                "Strand": locus1["Strand"]
            }

            self._random_queues[0].put((index1, random_position1))

            random_pair = {
                "Locus1" : index1,
                "Locus2" : index1,
                "Trans": pair["Trans"],
                "Ps": pair["Ps"],
                "Diag_mask": pair["Diag_mask"],
                "Distance": pair["Distance"],
                "Chrom1_circular": pair["Chrom1_circular"],
                "Chrom2_circular": pair["Chrom2_circular"],
                "Sep_id": pair["Sep_id"],
                "Original_pair": f"{pair['Locus1'], pair['Locus2']}"
            }

            # second random position
            if pair['Locus1'] != pair['Locus2']:
                index2 =  (k * self._nb_pairs + pair_index) * 2 + 1

                if locus1["Chromosome"] == locus2["Chromosome"]:
                    random_distance = random.randint(distance - self._jitter, distance + self._jitter)
                    random_start2 = random_start + random_distance

                else:
                    random_distance = None
                    center2 = self._get_center(locus2)
                    interval2 = [center2 - self._selection_window, center2 + self._selection_window]
                    random_start2 = int((random.random() * (interval2[1] - interval2[0]) + interval2[0]) // 1)

                random_position2 = {
                    "Chromosome": locus2["Chromosome"],
                    "Start": random_start2,
                    "End": random_start2,
                    "Name": f"rand_pos_{index2}",
                    "Strand": locus2["Strand"]
                }

                random_pair["Locus2"] = index2
                random_pair["Distance"] = random_distance

                self._random_queues[0].put((index2, random_position2))
                self._random_queues[1].put((index1, random_pair))

            else:
                self._random_queues[1].put((index1, random_pair))