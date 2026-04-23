from hicue.utils import *

global chrom_ps
global trans_av

ps_lock = threading.Lock()
diag_mask_lock = threading.Lock()
trans_lock = threading.Lock()

def initialize_globals():
    with ps_lock:
        global chrom_ps
        chrom_ps = {}
    with diag_mask_lock:
        global chrom_diag_mask
        chrom_diag_mask = {}
    with trans_lock:
        global trans_av
        trans_av = {}


class SubmatrixFormaterScheduler(threading.Thread):
    
    def __init__(self, input_queue, output_queues, display_queue, outpath, is_region = False, **subfargs):
        super(SubmatrixFormaterScheduler, self).__init__()

        self._input_queue = input_queue
        self._output_queues = output_queues
        self._display_queue = display_queue
        self._outpath = outpath
        self._is_region = is_region
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

            index, size_metric, pair, submatrix, tracks = val
            formated_submatrix = formater.format(submatrix, tracks, pair, is_region = self._is_region, expected_size = size_metric)

            for queue in self._output_queues:
                queue.put((index, size_metric, pair["Sep_id"], formated_submatrix))

            if self._display_queue is not None:
                to_queue = {
                    "matrix": formated_submatrix,
                    "pair": pair,
                    "binning": formater._binning,
                    "outfolder": f"{self._outpath}/{formater._cool_name}/{pair['Sep_id']}/binning_{formater._binning}",
                    "size_metric": size_metric,
                    "is_region": self._is_region
                }
                self._display_queue.put(to_queue)

class SubmatrixFormater():
    """Class for submatrix formating (detrending, flipping, masking, ...)"""
    def __init__(self, cool_file, cool_name, positions, flip = False, center = "start", method = "median", raw = False, ps_on_all = True, log = None):
        self._cool_file = cool_file
        self._cool_name = cool_name
        self._binning = cool_file.binsize
        self._positions = positions
        self._flip = flip
        self._center = center
        self._method = method
        self._raw = raw
        self._ps_on_all = ps_on_all
        self._log = log

    def get_ps(self, chromosome, is_circular = False, raw = False):
        """Computes or retrieve ps for the chromosomes"""
        ps = None
        chromosome = "All chromsomes" if self._ps_on_all else chromosome
        with ps_lock:
            global chrom_ps
            if chromosome in chrom_ps:
                ps = chrom_ps[chromosome]
            else:
                start_time = time.time()
                if self._ps_on_all:
                    ps = distance_law_genome(self._cool_file, method = self._method)
                    ps = np.concatenate([ps, [np.nan for _ in range(len(ps))]]) if not is_circular else ps
                else:
                    try:
                        chrom_matrix = self._cool_file.matrix(balance = (not raw)).fetch(chromosome)
                    except:
                        chrom_matrix = csr_matrix(self._cool_file.matrix(balance = (not raw), sparse = True).fetch(chromosome))
                    finally:
                        chrom_matrix = np.concatenate([np.concatenate([chrom_matrix, chrom_matrix], axis = 0), np.concatenate([chrom_matrix, chrom_matrix], axis = 0)], axis = 1) if is_circular else chrom_matrix
                        ps = distance_law(chrom_matrix, method = self._method)
                        ps = np.concatenate([ps, [np.nan for _ in range(len(ps))]]) if not is_circular else ps
                if self._log is not None:
                    self._log.write(f"Distance law for chromosome {chromosome} computed in {time.time() - start_time} seconds\n")
                chrom_ps[chromosome] = ps
        return ps
    
    def get_diag_mask(self, diag_mask_size, chromosome, is_circular = False):
        """Computes or retrieve the diagonal mask for a given chromosome and binning size"""
        diag_mask = None
        with diag_mask_lock:
            global chrom_diag_mask
            if self._binning not in chrom_diag_mask:
                chrom_diag_mask[self._binning] = {}
            if chromosome not in chrom_diag_mask[self._binning]:
                chrom_diag_mask[self._binning][chromosome] = compute_chromosome_diag_mask(self._cool_file.chromsizes[chromosome], diag_mask_size, self._binning, is_circular = is_circular)
            diag_mask = chrom_diag_mask[self._binning][chromosome]
        return diag_mask

    def get_trans_av(self, chromosome1, chromosome2):
        """Computes or retrieve trans detrending value for the chromosomes"""
        trans_det = np.nan
        with trans_lock:
            global trans_av
            if (chromosome1, chromosome2) in trans_av:
                trans_det = trans_av[(chromosome1, chromosome2)]
            elif (chromosome2, chromosome1) in trans_av:
                trans_det = trans_av[(chromosome2, chromosome1)]
            else:
                match self._method:
                    case "mean":
                        trans_det = self._cool_file.matrix(balance = (not self._raw), sparse = True).fetch(chromosome1, chromosome2).mean()
                    case "median":
                        trans_det = get_sparse_median(self._cool_file.matrix(balance = (not self._raw), sparse = True).fetch(chromosome1, chromosome2).data, 0)
                if self._log is not None:
                    self._log.write(f"Detrending trans contact between {chromosome1} and {chromosome2} with {self._method} value of contact map: {trans_det}.\n")
                trans_av[(chromosome1, chromosome2)] = trans_det
                trans_av[(chromosome2, chromosome1)] = trans_det
        return trans_det

    def detrend(self, submatrix, locus1, locus2, is_trans = False, is_circular = False):
        """Detrends cis-contact by the p(s), the median of the inter-chromosomal space otherwise."""
        detrended_submatrix = None
        if is_trans:
            detrending = self.get_trans_av(locus1["Chromosome"], locus2["Chromosome"])
            detrended_submatrix = submatrix / detrending
        else:
            ps = self.get_ps(locus1["Chromosome"], is_circular = is_circular)
            detrended_submatrix = detrend_submatrix(submatrix, locus1, locus2, self._binning, ps, center = self._center)
        return detrended_submatrix
    
    def flip(self, submatrix, subtracks, locus1, locus2, is_contact = False):
        """Flips a submatrix in the required axis."""
        # same position, flipped if sense is -1
        result_matrix = submatrix
        result_tracks = subtracks
        if not is_contact:
            if locus1['Strand'] == -1:
                result_matrix = np.flip(result_matrix)
                if result_tracks is not None:
                    result_tracks[0] = np.flip(result_tracks[0])
        
        else:
            # flipping along locus 1 axis
            if locus1['Strand'] == -1:
                result_matrix = np.flip(result_matrix, 0)
                if result_tracks is not None:
                    result_tracks[0] = np.flip(result_tracks[0])

            # flipping along locus 2 axis
            if locus2['Strand'] == -1:
                result_matrix = np.flip(result_matrix, 1)
                if result_tracks is not None:
                    result_tracks[1] = np.flip(result_tracks[1])
        return result_matrix, result_tracks

    def format(self, matrix, tracks, pair, is_region = False, expected_size = 0):
        """Applies all formating operations to a submatrix."""
        result_matrix = matrix

        locus1 = self._positions.loc[pair['Locus1']]
        locus2 = self._positions.loc[pair['Locus2']]

        diagonal_masking = int(pair["Diag_mask"])
        is_trans_contact = bool(pair["Trans"])
        mask = diagonal_masking > 0 and not is_trans_contact
        ps_detrending = bool(pair["Ps"])


        if is_region:
            if mask or ps_detrending:
                # computing the distance to diagonal matrix
                dist_matrix = get_dist_to_diag_matrix(self._cool_file, locus1, locus2)

                # masking diagonal
                if mask:
                    diag_mask = self.get_diag_mask(diagonal_masking, locus1["Chromosome"], is_circular = bool(pair["Chrom1_circular"]))
                    result_matrix = diag_mask[dist_matrix] * result_matrix

                # applying P(s)
                if ps_detrending:
                    if not is_trans_contact:
                        ps = self.get_ps(locus1["Chromosome"], is_circular = bool(pair["Chrom1_circular"]))
                        result_matrix = result_matrix / ps[dist_matrix]
                    else:
                        result_matrix = result_matrix / self.get_trans_av(locus1["Chromosome"], locus2["Chromosome"])
            result_matrix = resize_window(result_matrix, expected_size = expected_size)

        else:
            # masking
            result_matrix = mask_diagonal(result_matrix, 
                                        locus1, 
                                        locus2, 
                                        self._binning, 
                                        diagonal_masking,
                                        center = self._center)
            # detrending
            if ps_detrending:
                result_matrix = self.detrend(result_matrix, locus1, locus2, is_trans = pair["Trans"], is_circular = bool(pair["Chrom1_circular"]))
            
        # flipping if required and on diagonal
        if self._flip:
            result_matrix, tracks = self.flip(result_matrix, tracks, locus1, locus2, is_contact=pair['Locus1']!=pair['Locus2'])

        subtrack1, subtrack2 = tracks[0], tracks[1]

        if subtrack1 is not None:
            result_subtrack1 = resize_tracks(subtrack1, expected_size) if is_region else subtrack1.reshape(1, -1)
            result_matrix = np.concatenate([result_matrix, result_subtrack1], axis = 0)

        if subtrack2 is not None:
            result_subtrack2 = resize_tracks(subtrack2, expected_size) if is_region else subtrack2.reshape(1, -1)
            result_matrix = np.concatenate([result_matrix, result_subtrack2], axis = 0)

        return result_matrix