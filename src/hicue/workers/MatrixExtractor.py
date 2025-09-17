from .utils import *

def initialize_globals():
    pass

class MatrixExtractorScheduler(threading.Thread):
    def __init__(self, input_queue, output_queues, process_threads = 8, **mexargs):
        super(MatrixExtractorScheduler, self).__init__()
        self._cool_queue = input_queue
        self._output_queues = output_queues
        self._threads = process_threads
        self._mexargs = mexargs
        self.start()
        
    def run(self):
        matrixExtractor = MatrixExtractor(**self._mexargs)
        while True:
            try:
                value= self._cool_queue.get()
            except Empty:
                break
            if value == "DONE":
                break
            cool_file, cool_name = value
            pileups = matrixExtractor.extract_from(cool_file, cool_name, threads = self._threads)
            for sep_id in pileups.keys():
                for queue in self._output_queues:
                    queue.put({
                        "pileup":pileups[sep_id],
                        "binning": cool_file.binsize,
                        "sep_id": sep_id
                    })

class MatrixExtractor():
    """Extract the positions submatrices and creates pileups."""
    def __init__(self, formated_pairs, positions, windows, center = "start", raw = False, method = "median", flip = False, randoms = False, nb_rand_per_pos = 1, display_loci = False, display_batch = False, compute_pileups = True, outpath = "", display_args = {}, log = None):
        self._formated_pairs = formated_pairs
        self._positions = positions
        self._compute_pileups = compute_pileups
        self._windows = windows
        self._center = center
        self._raw = raw
        self._method = method
        self._flip = flip
        self._randoms = randoms
        self._nb_rand_per_pos = nb_rand_per_pos
        self._display_loci = display_loci
        self._display_batch = display_batch
        self._outpath = outpath
        self._display_args = display_args
        self._log = log

    @staticmethod
    def stream_pairs(input_queue, formated_pairs, randoms = False, nb_rand_per_pos = 1, threads = 8):
        """Puts each formated pair of sep_id (or random pair keeping the original pair format) in the input queue."""
        nb_pos = len(formated_pairs) // nb_rand_per_pos
        for index, pair in formated_pairs.iterrows():
            if randoms:
                for random_pair in yield_random_pairs(pair, nb_rand_per_pos, nb_pos):
                    input_queue.put((index, random_pair))
            else:
                input_queue.put((index, pair))
        for _ in range(threads):
            input_queue.put("DONE")
        
        
    def extract_from(self, cool_file, cool_name, threads = 8):
        """Extracts the submatrices of each pair of positions and aggregate them in a pileup."""
        unique_sep_ids = np.unique(self._formated_pairs["Sep_id"])

        # for each separation create a pileup
        nb_matrices = len(self._formated_pairs) if not self._randoms else len(self._formated_pairs) * self._nb_rand_per_pos
        pileups = {sep_id : Pileup(nb_matrices = nb_matrices, sep_id = sep_id, mode = self._method, cool_name = cool_name, binning = cool_file.binsize) for sep_id in unique_sep_ids}
            
        # queues initialisation
        input_queue = Queue()
        raw_submatrices_queue = Queue()
        formater_output_queues = [Queue()] if self._compute_pileups else []
        display_submatrices_queue = Queue() if not self._randoms and (self._display_batch or self._display_loci) else None 
        
        # initialisation of the individual extracters
        extracters = schedule_workers(
                worker_class = "ExtracterScheduler",
                worker_location = "hicue.workers.Extracter",
                threads = 1, # to avoid overhead in reading the cool file, only one extracter is alocated
                input_queue = input_queue,
                output_queues = [raw_submatrices_queue],
                cool_file = cool_file,
                positions = self._positions, 
                windows = self._windows, 
                center = self._center, 
                raw = self._raw,
                random = self._randoms
            )
        
        # initialisation of the matrix formaters
        formaters = schedule_workers(
                worker_class = "SubmatrixFormaterScheduler",
                worker_location = "hicue.workers.SubmatrixFormater",
                threads = threads,
                input_queue = raw_submatrices_queue,
                output_queues = formater_output_queues,
                display_queue = display_submatrices_queue,
                outpath = self._outpath,
                cool_file = cool_file,
                cool_name = cool_name,
                positions = self._positions,
                flip = self._flip,
                center = self._center,
                method = self._method,
                raw = self._raw,
                log = self._log
            )
        
        displayers = []
        if not self._randoms:
            displayer_output = [Queue()] if self._display_loci else []
            if self._display_batch:
                #initialisation of the loci displays (individual and batched)
                 # initializing a queue for the the individual display input
                displayers = schedule_workers(
                        worker_class = "DisplayBatch",
                        worker_location = "hicue.classes.AsyncDisplays",
                        threads = 1,
                        input_queue = display_submatrices_queue,
                        output_queues = displayer_output,
                        function = display_batch_submatrices,
                        batch_size = 64,
                        params_to_batch = ["outfolder"],
                        positions = self._positions,
                        chromsizes = cool_file.chromsizes,
                        **self._display_args
                    )

            if self._display_loci:
                display_input = displayer_output[0] if self._display_batch else display_submatrices_queue
                displayers += schedule_workers(
                    worker_class = "Display",
                    worker_location = "hicue.classes.AsyncDisplays",
                    threads = 1,
                    input_queue = display_input,
                    output_queues = [],
                    function = display_submatrix,
                    positions = self._positions,
                    chromsizes = cool_file.chromsizes,
                    **self._display_args
                )

        aggregators = []
        if self._compute_pileups :
        
            # initialisation of the pileup workers
            aggregators = schedule_workers(
                    worker_class = "AggregatorScheduler",
                    worker_location = "hicue.workers.Aggregator",
                    threads = threads,
                    input_queue = formater_output_queues[0],
                    output_queues = [],
                    pileups = pileups
                )
        
        # streaming the pairs in the input_queue
        self.stream_pairs(input_queue, self._formated_pairs, randoms = self._randoms, nb_rand_per_pos = self._nb_rand_per_pos, threads = 1)
        
        # joining         
        join_workers(extracters)
        
        join_queues([raw_submatrices_queue, display_submatrices_queue], threads = threads)

        join_workers(formaters)
        
        join_queues(formater_output_queues, threads = threads)
        
        join_workers(aggregators)

        join_workers(displayers)

        return pileups