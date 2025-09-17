from .utils import *

global chrom_ps
global trans_av

ps_lock = threading.Lock()
trans_lock = threading.Lock()

def initialize_globals():
    with ps_lock:
        global chrom_ps
        chrom_ps = {}
    with trans_lock:
        global trans_av
        trans_av = {}


class SubmatrixFormaterScheduler(threading.Thread):
    
    def __init__(self, input_queue, output_queues, display_queue, outpath, **subfargs):
        super(SubmatrixFormaterScheduler, self).__init__()

        self._input_queue = input_queue
        self._output_queues = output_queues
        self._display_queue = display_queue
        self._outpath = outpath
        self._subfargs = subfargs

        self.start()

    def run(self):
        formater = SubmatrixFormater(**self._subfargs)
        start = time.time()
        while True:
            try:
                val = self._input_queue.get()
            except Empty:
                break
            if val == 'DONE':
                break

            index, window, pair, submatrix = val
            formated_submatrix = formater.format(submatrix, pair)

            if index % 1000 == 0:
                print(index, time.time() - start)

            for queue in self._output_queues:
                queue.put((index, window, pair["Sep_id"], formated_submatrix))

            if self._display_queue is not None:
                self._display_queue.put({
                    "matrix": formated_submatrix,
                    "pair": pair,
                    "window": window,
                    "binning": formater._binning,
                    "outfolder": f"{self._outpath}/{formater._cool_name}/{pair['Sep_id']}/binning_{formater._binning}"
                })

class SubmatrixFormater():
    """Class for submatrix formating (detrending, flipping, masking, ...)"""
    def __init__(self, cool_file, cool_name, positions, flip = False, center = "start", method = "median", raw = False, log = None):
        self._cool_file = cool_file
        self._cool_name = cool_name
        self._binning = cool_file.binsize
        self._positions = positions
        self._flip = flip
        self._center = center
        self._method = method
        self._raw = raw
        self._log = log

    def get_ps(self, chromosome, raw = False):
        """Computes or retrieve ps for the chromosomes"""
        ps = None
        with ps_lock:
            global chrom_ps
            if chromosome in chrom_ps:
                ps = chrom_ps[chromosome]
            else:
                start_time = time.time()
                try:
                    chrom_matrix = self._cool_file.matrix(balance = (not raw)).fetch(chromosome)
                except:
                    chrom_matrix = csr_matrix(self._cool_file.matrix(balance = (not raw), sparse = True).fetch(chromosome))
                finally:
                    ps = distance_law(chrom_matrix, method = self._method)
                if self._log is not None:
                    self._log.write(f"Distance law for chromosome {chromosome} computed in {time.time() - start_time} seconds\n")
                chrom_ps[chromosome] = ps
        return ps

    def get_trans_av(self, chromosome1, chromosome2, method = "median", raw = False):
        """Computes or retrieve trans detrending value for the chromosomes"""
        trans_det = np.nan
        with trans_lock:
            global trans_av
            if (chromosome1, chromosome2) in trans_av:
                trans_det = trans_av[(chromosome1, chromosome2)]
            elif (chromosome2, chromosome1) in trans_av:
                trans_det = trans_av[(chromosome2, chromosome1)]
            else:
                match method:
                    case "mean":
                        trans_det = self._cool_file.matrix(balance = (not raw), sparse = True).fetch(chromosome1, chromosome2).mean()
                    case "median":
                        trans_det = get_sparse_median(self._cool_file.matrix(balance = (not raw), sparse = True).fetch(chromosome1, chromosome2).data, 0)
                trans_av[(chromosome1, chromosome2)] = trans_det
                trans_av[(chromosome2, chromosome1)] = trans_det
        return trans_det

    def detrend(self, submatrix, locus1, locus2, is_trans = False):
        """Detrends cis-contact by the p(s), the median of the inter-chromosomal space otherwise."""
        detrended_submatrix = None
        if is_trans:
            detrending = self.get_trans_av(locus1["Chromosome"], locus2["Chromosome"], method = self._method, raw = self._raw)
            detrended_submatrix = submatrix / detrending
        else:
            ps = self.get_ps(locus1["Chromosome"], locus2["Chromosome"])
            detrended_submatrix = detrend_submatrix(submatrix, locus1, locus2, self._binning, ps, center = self._center)
        return detrended_submatrix

    def format(self, matrix, pair):
        """Applies all formating operations to a submatrix."""
        result_matrix = matrix
        locus1 = self._positions.loc[pair['Locus1']]
        locus2 = self._positions.loc[pair['Locus2']]

        # masking
        result_matrix = mask_diagonal(result_matrix, 
                                      locus1, 
                                      locus2, 
                                      self._binning, 
                                      int(pair["Diag_mask"]), 
                                      center = self._center)
        # detrending
        if bool(pair["Ps"]):
            result_matrix = self.detrend(result_matrix, locus1, locus2, is_trans = pair["Trans"])
        # flipping if required and on diagonal
        if self._flip and pair['Locus1'] == pair['Locus2'] and self._positions.loc[pair['Locus1']]['Strand'] < 0:
            np.flip(result_matrix)
        return result_matrix