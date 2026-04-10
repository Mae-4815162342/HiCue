from hicue.utils import *
from hicue.classes.PositionList import PositionList

class TrackReader():
    
    def __init__(self, tracks, threshold = None, percentage = None, min_sep = 1000, positions_file = None, position_type = "bed", overlap = "strict", save_to="", loop = False, record_type = None):
        self._tracks = tracks
        self._threshold = threshold
        self._percentage = percentage
        self._binning = min_sep
        self._positions_file = positions_file
        self._positions_type = position_type
        self._overlap = overlap
        self._save_to = save_to
        self._loop = loop
        self._record_type = record_type

        if len(self._save_to) > 0: # TODO write all outputs with asynchronous functions
            if not os.path.exists(self._save_to):
                os.mkdir(self._save_to)
    
    def read_file(self, save_to = None, threads = 8):
        """Reads the file and parses unique positions. Returns the positions Dataframe and the pairing queue."""

        tracks = pyBigWig.open(self._tracks)

        # building from positions 
        if self._positions_file:
            positions, pairs = self.read_from_positions(tracks, save_to = save_to, threads = threads)
        # building from tracks
        else:
            positions = self.read_from_tracks(tracks, save_to = save_to, threads = threads)
            pairs = None
        return positions, pairs
    
    def read_from_tracks(self, tracks, save_to = None, threads = 8):
        """Reads a track file by bins and applies threshold or percentage filters."""
        tracks_queue = Queue()
        position_queue = Queue()

        workers = []
        position_list = None
        if self._percentage is not None:
            high_low, percentage = self._percentage
            nb_pos = compute_nb_pos_tracks(tracks, percentage, self._binning)
            position_list = PositionList(nb_pos, high_low, binning = self._binning)

            # launching parsers
            workers = schedule_workers(
                worker_class = "TracksPercentageParser",
                worker_location = "hicue.workers.TracksFilters",
                threads = threads,
                position_list = position_list,
                tracks = tracks,
                input_queue = tracks_queue,
            )

        else:
            # launching parsers
            workers = schedule_workers(
                worker_class = "TracksThresholdParser",
                worker_location = "hicue.workers.TracksFilters",
                threads = threads,
                threshold = self._threshold,
                tracks = tracks,
                input_queue = tracks_queue,
                output_queue = position_queue
            )
            
        streamer = schedule_workers(
                worker_class = "TracksStreamer",
                worker_location = "hicue.workers.TracksStreamer",
                threads = 1,
                tracks = tracks,
                output_queue = tracks_queue,
                binning = self._binning,
                workers = threads
            )

        
        #joining        
        join_workers(streamer)

        join_workers(workers)

        join_queues([position_queue], threads=threads)

        positions = position_queue_to_df(position_queue) if self._percentage is None else position_list.get_positions()
        if save_to:
            positions.to_csv(f"{save_to}/positions_indexed.csv", sep=",")

        return positions
            
    def read_from_positions(self, tracks, save_to = None, threads = 8):
        """Annotates a position file with its associated tracks and applies threshold or percentage filters."""

        # defining queues
        lines_queue = Queue()
        position_queue = Queue()
        pairing_queue = Queue()
        annotation_queue = Queue()
                
        # launching parsers
        workers = schedule_workers(
            worker_class = "ParserScheduler",
            worker_location = "hicue.workers.Parser",
            threads = threads,
            input_queue = lines_queue,
            output_queues = [annotation_queue, pairing_queue],
            file_type = self._positions_type,
            record_type = self._record_type,
            is_loop = self._loop
        )
        
        annotators = schedule_workers(
            worker_class = "AnnotatorScheduler",
            worker_location = "hicue.workers.Annotator",
            threads = threads,
            input_queue = annotation_queue,
            output_queues = [position_queue],
            gff = None,
            overlap = self._overlap,
            save_to = self._save_to,
            gff_type = self._record_type,
            tracks = tracks,
            threshold = self._threshold,
        )
        
        percentage_parsers = []
        position_list = None
        if self._percentage is not None:

            # initialazing position_list
            high_low, percentage = self._percentage
            nb_pos = compute_nb_selected_pos(self._positions_file, percentage, self._positions_type, gff_types = [self._record_type]) 
            position_list = PositionList(nb_pos, high_low, binning = self._binning)

            percentage_parsers = schedule_workers(
                worker_class = "PositionPercentageParser",
                worker_location = "hicue.workers.TracksFilters",
                threads = threads,
                input_queue = position_queue,
                tracks = tracks,
                position_list = position_list
            )
    
        # reading the file
        workers += schedule_workers(
                    worker_class = "FileStreamer",
                    worker_location = "hicue.workers.FileStreamer",
                    file = self._positions_file,
                    threads = 1,
                    t = threads,
                    output_queues = [lines_queue]
                )

        #joining       

        join_workers(workers)

        join_queues([annotation_queue], threads=threads)

        join_workers(annotators)

        join_queues([position_queue], threads=threads)

        join_workers(percentage_parsers)

        positions = position_queue_to_df(position_queue) if self._percentage is None else position_list.get_positions()
        
        if save_to:
            positions.to_csv(f"{save_to}/positions_indexed.csv", sep=",")

        return positions, pairing_queue