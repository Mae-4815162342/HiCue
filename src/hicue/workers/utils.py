from .imports import *
import time

def schedule_workers(worker_class, worker_location, threads, **wargs):
    """Creates a worker of class worker_class per thread with the arguments wargs. Returns the instances created in a list."""
    worker_module = importlib.import_module(worker_location)
    worker_module.initialize_globals()
    WorkerClass = getattr(importlib.import_module(worker_location), worker_class)
    workers = []
    for _ in range(threads):
        workers.append(WorkerClass(**wargs))
    return workers

def join_workers(workers):
    """Joins workers in the workers list."""
    for worker in workers:
        worker.join()
        
def join_queues(queues, threads = 1):
    """Signal DONE in each queue of queues"""
    for queue in queues:
        if queue is not None:
            for _ in range(threads):
                queue.put("DONE")

def position_queue_to_df(position_queue):
    """Converts a position_queue output to a dataframe."""
    position_list = []
    index_list = []

    while True:
        try:
            value = position_queue.get(timeout=10)
        except Empty:
            break
        if value == "DONE":
            break
        index, position = value
        index_list.append(index)
        position_list.append(position)
    
    return pd.DataFrame(position_list, index = index_list)

def adjust_locus(locus_position, chromsize, is_circ_chrom = False):
        """Adjust the coordinates of a locus to fit between 0 and the chromosome size. Accounts for circularity."""
        if locus_position < 0:
            locus_position = chromsize + locus_position if is_circ_chrom else 0

        if locus_position >= chromsize:
            locus_position = locus_position - chromsize if is_circ_chrom else chromsize
        return locus_position

def extract_window(cool, locus1, locus2, binning, window, is_loc1_circ = False, is_loc2_circ = False, center="start", raw = False): 
    matrix = cool.matrix(balance=(not raw))

    match center:
        case "start":
            start1 = min(locus1["Start"], locus1["End"]) if locus1["Strand"] == 1 else max(locus1["Start"], locus1["End"])
            start2 = min(locus2["Start"], locus2["End"]) if locus2["Strand"] == 1 else max(locus2["Start"], locus2["End"])
        case "center":
            start1 = (locus1["Start"] + locus1["End"]) // 2
            start2 = (locus2["Start"] + locus2["End"]) // 2
        case "end":
            start1 = max(locus1["Start"], locus1["End"]) if locus1["Strand"] == 1 else min(locus1["Start"], locus1["End"])
            start2 = max(locus2["Start"], locus2["End"]) if locus2["Strand"] == 1 else min(locus2["Start"], locus2["End"])

    chrom_size1 = cool.chromsizes[locus1["Chromosome"]]
    chrom_size2 = cool.chromsizes[locus2["Chromosome"]]

    start1 = adjust_locus(start1, chrom_size1, is_circ_chrom = is_loc1_circ)
    start2 = adjust_locus(start2, chrom_size2, is_circ_chrom = is_loc2_circ)

    start1 = start1 + 1 if start1 % binning == 0 else start1
    start2 = start2 + 1 if start2 % binning == 0 else start2

    # 1. checking overflows
    is_start1_inf = start1 - window < 0
    is_start1_sup = start1 + window > chrom_size1

    is_start2_inf = start2 - window < 0
    is_start2_sup = start2 + window > chrom_size2

    # 2. computing intervales
    pos1 = f"{locus1['Chromosome']}:{start1 - window if not is_start1_inf else 0}-{start1 + window if not is_start1_sup else chrom_size1}"
    pos2 = f"{locus2['Chromosome']}:{start2 - window if not is_start2_inf else 0}-{start2 + window if not is_start2_sup else chrom_size2}"

    # 3. fetching main submatrix
    submatrix = matrix.fetch(pos1, pos2)[:]

    expected_size = (window//binning) * 2 + 1

    start1_overflow = is_start1_inf or is_start1_sup
    start2_overflow = is_start2_inf or is_start2_sup

    # 4. managing overflows
    if (not start1_overflow and not start2_overflow): # no overflow
        return submatrix
    
    bins1_to_fill = expected_size - submatrix.shape[0]
    bins2_to_fill = expected_size - submatrix.shape[1]
    
    fill1, fill2 = [], []
    
    if is_start1_inf: # dim 1 inf
        # computing indexes to fill
        fill1 = np.array([chrom_size1 - (bins1_to_fill * binning), chrom_size1])
        
    if is_start1_sup: # dim 1 sup
        # computing indexes to fill
        fill1 = np.array([0, bins1_to_fill * binning])
    
    fill1_pos = f"{locus1['Chromosome']}:{fill1[0]}-{fill1[1]}" if len(fill1) > 0 else pos1
    len1 = abs(fill1[0] - fill1[1]) // binning if len(fill1) > 0 else 0
        
    if is_start2_inf: # dim 2 inf
        # computing indexes to fill
        fill2 = np.array([chrom_size2 - (bins2_to_fill * binning), chrom_size2])
        
    if is_start2_sup: # dim 2 sup
        # computing indexes to fill
        fill2 = np.array([0, bins2_to_fill * binning])
        
    fill2_pos = f"{locus2['Chromosome']}:{fill2[0]}-{fill2[1]}" if len(fill2) > 0 else pos2
    len2 = abs(fill2[0] - fill2[1]) // binning if len(fill2) > 0 else 0
    
    mat1 = matrix.fetch(pos1, fill2_pos) if is_loc2_circ else np.full((submatrix.shape[0], len2), np.nan)
    if mat1.shape[1] > len2:
        mat1 = mat1[:, 1:] if is_start2_inf else mat1[:, :-1]
    mat2 = matrix.fetch(fill1_pos, pos2) if is_loc1_circ else np.full((len1, submatrix.shape[1]), np.nan)
    if mat2.shape[0] > len1:
        mat2 = mat2[1:] if is_start1_inf else mat2[:-1]
        
    # two dimensions to fill
    if start1_overflow and start2_overflow:

        mat3 = matrix.fetch(fill1_pos, fill2_pos) if is_loc1_circ and is_loc2_circ else np.full((len1,len2), np.nan)
        if mat3.shape[1] > len2:
            mat3 = mat3[:, 1:] if is_start2_inf else mat3[:, :-1]
        if mat3.shape[0] > len1:
            mat3 = mat3[1:] if is_start1_inf else mat3[:-1]

        to_concat1 = [mat1, submatrix] if is_start2_inf else [submatrix, mat1]
        concat1 = np.concatenate(to_concat1, axis = 1)
        
        to_concat2 = [mat3, mat2] if is_start2_inf else [mat2, mat3]
        concat2 = np.concatenate(to_concat2, axis = 1)
        
        to_concat3 = [concat2, concat1] if is_start1_inf else [concat1, concat2]
        submatrix = np.concatenate(to_concat3, axis = 0)
        
    # dim1 to fill
    elif start1_overflow:        
        to_concat = [mat2, submatrix] if is_start1_inf else [submatrix, mat2]
        submatrix = np.concatenate(to_concat, axis = 0)
        
    # dim2 to fill
    elif start2_overflow:
        to_concat = [mat1, submatrix] if is_start2_inf else [submatrix, mat1]
        submatrix = np.concatenate(to_concat, axis = 1)

    return submatrix

def compute_distance(locus1, locus2, center = "start"):
    """Returns the distance in base pairs between two position. None if not in the same chromosome."""
    if locus1["Chromosome"] != locus2["Chromosome"]:
        return None
    match center:
        case "start":
            start1 = min(locus1["Start"], locus1["End"]) if locus1["Strand"] == 1 else max(locus1["Start"], locus1["End"])
            start2 = min(locus2["Start"], locus2["End"]) if locus2["Strand"] == 1 else max(locus2["Start"], locus2["End"])
        case "center":
            start1 = (locus1["Start"] + locus1["End"]) // 2
            start2 = (locus2["Start"] + locus2["End"]) // 2
        case "end":
            start1 = max(locus1["Start"], locus1["End"]) if locus1["Strand"] == 1 else min(locus1["Start"], locus1["End"])
            start2 = max(locus2["Start"], locus2["End"]) if locus2["Strand"] == 1 else min(locus2["Start"], locus2["End"])
    return start2 - start1

def mask_diagonal(submatrix, locus1, locus2, binning, diagonal_mask, center = "start"):
    """Computes the mask to apply to the diagonal from positions 1 and 2"""
    locus_distance = compute_distance(locus1, locus2, center = center)
    if locus_distance is None:
        return submatrix
    
    dist = locus_distance // binning
    if abs(dist) >= len(submatrix):
        return submatrix
    
    if dist == 0: # centered
        for i in range(diagonal_mask//binning):
            np.fill_diagonal(submatrix[i:],  np.nan)
            np.fill_diagonal(submatrix[:,i:],  np.nan)
            
    if dist < 0: # upper diagonal
        for i in range(diagonal_mask//binning):
            dist_i = abs(dist) + i
            np.fill_diagonal(submatrix[:- (abs(dist) + i),  abs(dist) + i:], np.nan)
            np.fill_diagonal(submatrix[:-  (abs(dist) - i),  abs(dist) - i:], np.nan)
            
    if dist > 0: # lower diagonal
        for i in range(diagonal_mask//binning):
            np.fill_diagonal(submatrix[abs(dist) + i:, :-(abs(dist) + i)], np.nan)
            np.fill_diagonal(submatrix[abs(dist) - i:, :- (abs(dist) - i)], np.nan)
            
    return submatrix

def detrend_submatrix(submatrix, locus1, locus2, binning, ps, center="start"):
    """Applies P(s) to a submatrix."""
    dist = compute_distance(locus1, locus2, center = center) // binning
    if dist == 0:
        submatrix_index = [[abs(i - j) for i in range(len(submatrix))] for j in range(len(submatrix))]
    else:
        index_mask =  np.array([
            [abs(i - j) for i in range(len(submatrix) + abs(dist))] 
            for j in range(len(submatrix) + abs(dist))
        ])
        submatrix_index = index_mask[:-dist, dist:] if dist > 0 else index_mask[abs(dist):, :-abs(dist)]
    
    # dealing with overflows: as the submatrix already has NaNs if the chromosomes are not circular, we can return a value for the ps for those bins
    submatrix_index = np.array(submatrix_index) % len(ps)
    submatrix_det = submatrix / ps[submatrix_index]
    return submatrix_det

def yield_random_pairs(pair, nb_rand_per_pos, nb_pos):
    """Yields the formated pair of each computed random pair index"""
    for k in range(nb_rand_per_pos):
        random_pair = pair.copy()
        random_pair["Locus1"] = k * nb_pos + pair["Locus1"] 
        random_pair["Locus2"] = k * nb_pos + pair["Locus2"] 
        yield random_pair

### Distance law adapted from Chromosight (Mathey-Doret et al., 2020)
def distance_law(
    matrix, detectable_bins=None, max_dist=None, smooth=True, method="mean"
):
    """
    Computes genomic distance law by averaging over each diagonal in the upper
    triangle matrix. If a list of detectable bins is provided, pixels in
    missing bins will be excluded from the averages. A maximum distance can be
    specified to define how many diagonals should be computed.

    parameters
    ----------
    matrix: scipy.sparse.csr_matrix
        the input matrix to compute distance law from.
    detectable_bins : numpy.ndarray of ints
        An array of detectable bins indices to consider when computing
        distance law.
    max_dist : int
        Maximum distance from diagonal, in number of bins in which to compute
        distance law
    smooth : bool
        Whether to use isotonic regression to smooth the distance law.
    fun : callable
        A function to apply on each diagonal. Defaults to mean.

    Returns
    -------
    dist: np.ndarray
        the output genomic distance law.

    example
    -------
        >>> m = np.ones((3,3))
        >>> m += np.array([1,2,3])
        >>> m
        array([[2., 3., 4.],
               [2., 3., 4.],
               [2., 3., 4.]])
        >>> distance_law(csr_matrix(m))
        array([3. , 3.5, 4. ])

    """
    mat_n = matrix.shape[0]
    if max_dist is None:
        max_dist = mat_n
    n_diags = min(mat_n, max_dist + 1)
    dist = np.zeros(mat_n)
    if detectable_bins is None:
        detectable_bins = np.array(range(mat_n))
    match method:
        case "mean":
            fun = np.nanmean
        case "median":
            fun = np.nanmedian
        case "sum":
            fun = np.nansum
    
    for diag in range(n_diags):
        # Find detectable which fall in diagonal
        detect_mask = np.zeros(mat_n, dtype=bool)
        detect_mask[detectable_bins] = 1
        # Find bins which are detectable in the diagonal (intersect of
        # hori and verti)
        detect_mask_h = detect_mask[: (mat_n - diag)]
        detect_mask_v = detect_mask[mat_n - (mat_n - diag) :]
        detect_mask_diag = detect_mask_h & detect_mask_v
        detect_diag = matrix.diagonal(diag)[detect_mask_diag]
        
        diag_values = detect_diag[detect_diag > 0]
        if len(diag_values) > 0:
            dist[diag] = fun(diag_values)
        else:
            dist[diag] = np.nan

    # Smooth the curve using isotonic regression: Find closest approximation
    # with the condition that point n+1 cannot be higher than point n.
    # (i.e. contacts can only decrease when increasing distance)
    if smooth and mat_n > 2:
        ir = IsotonicRegression(increasing=False)
        dist[~np.isfinite(dist)] = 0
        dist = ir.fit_transform(range(len(dist)), dist)

    return dist