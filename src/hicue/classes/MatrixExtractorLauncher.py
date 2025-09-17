from .utils import *

class MatrixExtractorLauncher():
    def __init__(self, cool_files, binnings = [], **extparams):
        self._cool_files = cool_files
        self._binnings = binnings
        self._extparams = extparams
        
    def stream_cools(self, cool_queue, threads = 8):
        for cool_path in self._cool_files:
            cool_name = ".".join(cool_path.split('::')[0].split('/')[-1].split('.')[:-1])
            if cooler.fileops.is_multires_file(cool_path):
                for binning in self._binnings:
                    cool = cooler.Cooler(f"{cool_path}::resolutions/{binning}")
                    cool_queue.put((cool, cool_name))
            else:
                cool = cooler.Cooler(cool_path)
                cool_queue.put((cool, cool_name))
        for _ in range(threads):
            cool_queue.put("DONE")

    def launch_extraction(self, positions, formated_pairs, randoms = False, threads = 8):
        """Launches the multiprocessing of matrix extraction."""
        cool_queue = Queue()
        pileup_queue = Queue()

        #  initialisation of extractors
        extractors = schedule_workers(
                    worker_class = "MatrixExtractorScheduler",
                    worker_location = "hicue.workers.MatrixExtractor",
                    threads = min(threads, len(self._cool_files)), # as one process can only work on one cool at the time, we only allocate at most the number of cool of processes
                    input_queue = cool_queue,
                    output_queues = [pileup_queue],
                    process_threads = threads,
                    formated_pairs = formated_pairs,
                    positions = positions,
                    randoms = randoms,
                    **self._extparams
                )
        
        self.stream_cools(cool_queue, threads = threads)

        join_workers(extractors)

        for _ in range(threads):
            pileup_queue.put("DONE")

        return pileup_queue
        

