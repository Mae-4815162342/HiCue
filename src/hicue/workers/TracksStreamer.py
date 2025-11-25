from hicue.imports import *

def initialize_globals():
    pass

class TracksStreamer(threading.Thread):
    def __init__(self, tracks, output_queue, binning = 1000, workers = 8, **kwargs):
        super(TracksStreamer, self).__init__(**kwargs)
        
        self._tracks = pyBigWig.open(tracks) if type(tracks) == str else tracks
        self._output_queue = output_queue
        self._binning = binning
        
        self._threads = workers
        
        self.start()
        
    def run(self):
        index = 0
        for chrom, size in self._tracks.chroms().items():
            for k in range(0, size, self._binning):
                if k + 1 + self._binning < size:
                    self._output_queue.put((index, chrom, k + 1, k + self._binning))
                    index += 1
            if k + self._binning < size:
                self._output_queue.put((index, chrom, k + 1, size))
        for _ in range(self._threads):
            self._output_queue.put("DONE")