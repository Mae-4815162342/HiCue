from .imports import *
from .utils import *

def initialize_globals():
    pass

class FileStreamer(threading.Thread):
    def __init__(self, file, t, output_queues, comments = False):
        super(FileStreamer, self).__init__()
        self._file = file
        self._threads = t
        self._output_queues = output_queues
        self._comments = comments
        self.start()
        
    def run(self):
        """Streams a file into the output queues, line by line, without commented lines if comments is False."""
        try:
            with open(self._file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or (line.startswith('#') and not self._comments):
                        continue
                    for queue in self._output_queues:
                        queue.put(line)
            join_queues(self._output_queues, threads = self._threads)
        except Exception as e:
            print(f"Error in stream_file: {e}")
            join_queues([self._output_queue], threads = self._threads)