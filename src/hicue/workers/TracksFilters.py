from hicue.imports import *

def initialize_globals():
    pass

class TracksPercentageParser(threading.Thread):
    
    def __init__(self, position_list, tracks, input_queue, **kwargs):
        super(TracksPercentageParser, self).__init__(**kwargs)
        
        self._position_list = position_list
        self._tracks = tracks
        self._input_queue = input_queue

        self.start()

    def run(self):
        while True:
            try:
                val = self._input_queue.get()
            except Empty:
                break
            if val == 'DONE':
                break
            index, chrom, start, stop,  = val
            tracks_value = self._tracks.stats(chrom, start - 1, stop - 1)[0]
            position = {
                "Index": index,
                "Name": f"Locus_{index}",
                "Chromosome": chrom,
                "Start": start,
                "End": stop,
                "Strand": 0,
                "Tracks": tracks_value
            }
            self._position_list.append(position, tracks_value)
            
class TracksThresholdParser(threading.Thread):
    
    def __init__(self, threshold, tracks, input_queue, output_queue, **kwargs):
        super(TracksThresholdParser, self).__init__(**kwargs)
        
        self._min_max, self._threshold = threshold or (None, None)
        self._tracks = tracks
        self._input_queue = input_queue
        self._output_queue = output_queue
        
        self.start()

    def run(self):
        while True:
            try:
                val = self._input_queue.get()
            except Empty:
                break
            if val == 'DONE':
                break
            index, chrom, start, stop = val
            tracks_value = self._tracks.stats(chrom, start - 1, stop - 1)[0]
            if tracks_value is not None:
                if self._threshold is not None:
                    to_keep = (tracks_value >= self._threshold  and self._min_max == "min") or (tracks_value <= self._threshold  and self._min_max == "max")
                else:
                    to_keep = True
                if to_keep:
                    position = {
                        "Name": f"Locus_{index}",
                        "Chromosome": chrom,
                        "Start": start,
                        "End": stop,
                        "Strand": 0,
                        "Tracks": tracks_value
                    }
                    self._output_queue.put((index,position))

class PositionPercentageParser(threading.Thread):
    
    def __init__(self, position_list, tracks, input_queue, **kwargs):
        super(PositionPercentageParser, self).__init__(**kwargs)
        
        self._position_list = position_list
        self._tracks = tracks
        self._input_queue = input_queue

        self.start()

    def run(self):
        while True:
            try:
                val = self._input_queue.get()
            except Empty:
                break
            if val == 'DONE':
                break
            index, position  = val
            tracks_value = float(position["Tracks"])
            self._position_list.append(position, tracks_value)