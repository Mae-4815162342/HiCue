from .utils import *

class Reader():
    
    def __init__(self, file, file_type, annotation_files={}, overlap = "strict", save_to="", loop = False, record_type = None):
        self._file = file
        self._file_type = file_type
        self._annotation_files = annotation_files
        self._overlap = overlap
        self._save_to = save_to
        self._loop = loop
        self._record_type = record_type

        if len(self._save_to) > 0: # TODO write all outputs with asynchronous functions
            if not os.path.exists(self._save_to):
                os.mkdir(self._save_to)
    
    @staticmethod
    def open_tracks(track_path):
        """Opens a bigwig tracks file."""
        return pyBigWig.open(track_path)
        
    def read_file(self, save_to = None, threads = 8):
        """Reads the file and parses unique positions. Returns the positions Dataframe and the pairing queue."""
        workers = []
        
        # defining queues
        lines_queue = Queue()
        position_queue = Queue()
        pairing_queue = Queue()  
        tmp_queues = []
                
        # launching parsers
        workers += schedule_workers(
                worker_class = "ParserScheduler",
                worker_location = "hicue.workers.Parser",
                threads = threads,
                input_queue = lines_queue,
                output_queues = [position_queue, pairing_queue],
                file_type = self._file_type,
                record_type = self._record_type,
                is_loop = self._loop
            )
        
        # launching anotator if required
        annotators = []
        if len(self._annotation_files) > 0:
            tracks = self.open_tracks(self._annotation_files["tracks"]) if "tracks" in self._annotation_files else None
            gff = split_gff(self._annotation_files["gff"]) if "gff" in self._annotation_files else None
            
            tmp_queue = position_queue
            position_queue = Queue()
            tmp_queues.append(tmp_queue)
            annotators = schedule_workers(
                    worker_class = "AnnotatorScheduler",
                    worker_location = "hicue.workers.Annotator",
                    threads = threads,
                    input_queue = tmp_queue,
                    output_queues = [position_queue],
                    gff = gff,
                    overlap = self._overlap,
                    save_to = self._save_to,
                    gff_type = self._record_type,
                    tracks = tracks,
                )
        
        # reading the file
        workers += schedule_workers(
                    worker_class = "FileStreamer",
                    worker_location = "hicue.workers.FileStreamer",
                    file = self._file,
                    threads = 1,
                    t = threads,
                    output_queues = [lines_queue]
                )

        #joining        
        join_workers(workers)

        join_queues(tmp_queues, threads=threads)

        join_workers(annotators)

        join_queues([position_queue, pairing_queue], threads=threads)
        
        positions = position_queue_to_df(position_queue)
        if save_to:
            positions.to_csv(f"{save_to}/positions_indexed.csv", sep=",")

        return positions, pairing_queue