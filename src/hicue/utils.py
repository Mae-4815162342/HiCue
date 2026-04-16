from hicue.imports import *
path_lock = threading.Lock()

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

def adjust_locus(locus_position, chromsize):
        """Adjust the coordinates of a locus to fit between 0 and the chromosome size. Accounts for circularity."""
        if locus_position < 0:
            locus_position = chromsize + locus_position

        if locus_position >= chromsize:
            locus_position = locus_position - chromsize
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

    if locus1["Chromosome"] not in cool.chromsizes or locus2["Chromosome"] not in cool.chromsizes:
        return None
    chrom_size1 = cool.chromsizes[locus1["Chromosome"]]
    chrom_size2 = cool.chromsizes[locus2["Chromosome"]]

    start1 = adjust_locus(start1, chrom_size1)
    start2 = adjust_locus(start2, chrom_size2)
    
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

def bin_tracks(tracks, chrom, start, stop, binning):
    """Returns the binned extracted region from tracks, delimited by start and stop in chrom."""
    values = []
    k = start
    while k < stop:
        end =  k + binning - 1
        end = end if end < stop else stop
        values.append(tracks.stats(chrom, k, end)[0])
        k += binning
    if k < stop:
        values.append(tracks.values(chrom, k, stop)[0])
    return values
    
def extract_tracks(tracks, locus, binning, window, is_loc_circ = False, center="start"): 
    """Extracts a window from tracks around the locus, binned at binning."""
    # 1. computing central coordinate
    match center:
        case "start":
            coordinate = min(locus["Start"], locus["End"]) if locus["Strand"] == 1 or locus["Strand"] == 0 else max(locus["Start"], locus["End"])
        case "center":
            coordinate = (locus["Start"] + locus["End"]) // 2
        case "end":
            coordinate = max(locus["Start"], locus["End"]) if locus["Strand"] == 1 or locus["Strand"] == 0 else min(locus["Start"], locus["End"])

    chrom_size = tracks.chroms(locus["Chromosome"])

    coordinate = adjust_locus(coordinate, chrom_size) - 1 # re-ajusting at 0 base for bw format

    # 2. computing binned coordinates
    binned_coordinates = (coordinate // binning)*binning
    start = binned_coordinates
    stop = binned_coordinates + binning - 1
    
    # 3. checking overflows
    is_start_inf = start - window < 0
    is_start_sup = stop + window >= chrom_size
    
    # 3. computing intervales
    window_start = start - window if not is_start_inf else 0
    window_stop = stop + window if not is_start_sup else chrom_size - 1

    # 4. fetching main subtracks
    subtracks = bin_tracks(tracks, locus["Chromosome"], window_start, window_stop, binning)

    # 5. managing overflows
    expected_size = (window//binning) * 2 + 1
    start_overflow = is_start_inf or is_start_sup
    
    if not start_overflow: # no overflow
        return np.array(subtracks)
    
    bins_to_fill = expected_size - len(subtracks)
    fill = [np.nan] * bins_to_fill if not is_loc_circ else []
    
    if is_start_inf and is_loc_circ: # dim 1 inf
        start_inf = (chrom_size//binning - bins_to_fill + 1) * binning
        stop_inf = chrom_size
        fill = bin_tracks(tracks, locus["Chromosome"], start_inf, stop_inf, binning)
        
    if is_start_sup and is_loc_circ: # dim 1 sup
        start_sup = 0
        stop_sup =  bins_to_fill * binning - 1
        fill = bin_tracks(tracks, locus["Chromosome"], start_sup, stop_sup, binning)
        
    to_concat = [fill, subtracks] if is_start_inf else [subtracks, fill]
    subtracks = np.concatenate(to_concat)

    return subtracks

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

# def yield_random_pairs(pair, nb_rand_per_pos, nb_pos):
#     """Yields the formated pair of each computed random pair index"""
#     for k in range(nb_rand_per_pos):
#         random_pair = pair.copy()
#         random_pair["Locus1"] = k * nb_pos + pair["Locus1"] 
#         random_pair["Locus2"] = k * nb_pos + pair["Locus2"]
#         yield random_pair

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

def empty_queue_in_dict(queue, keys):
    """Empties a dict queue in a dict, using keys as the dict element key."""
    queue_dict = {}
    while True:
        try:
            value = queue.get(timeout = 10)
        except Empty:
            break
        if value == "DONE":
            break
        key = "_".join([str(value[k]) for k in keys])
        queue_dict[key] = value
    return queue_dict

def create_folder_path(path):
    """If a folder path does not exists, create all the dependencies to this path."""
    global path_lock
    path_list = path.split("/")
    to_add = ""
    if path_list[0] == "":
        path_list = path_list[1:]
        to_add = "/"
    with path_lock:
        for i in range(1, len(path_list) + 1):
            current_path = to_add + "/".join(path_list[:i])
            if not os.path.exists(current_path):
                os.mkdir(current_path)

def compute_nb_pos_tracks(tracks, percentage, binning):
    """Computes the number of position that will be kept from the binned tracks to get the right percentage."""
    chromosomes = tracks.chroms()
    nb_pos_tot = 0
    for chrom in chromosomes:
        nb_pos_tot += chromosomes[chrom] // binning + 1
    return round(nb_pos_tot * percentage / 100)

def compute_nb_pos_gff(gff, percentage, gff_types = ["gene"]):
    """Computes the number of position that will be kept from the gff selected type to get the right percentage."""
    examiner = GFF.GFFExaminer()
    in_handle = open(gff)
    gff_type_counts = examiner.available_limits(in_handle)['gff_type']
    nb_pos_tot = 0
    for gff_type, count in gff_type_counts.items():
        nb_pos_tot += count if gff_type[0] in gff_types else 0
    return round(nb_pos_tot * percentage / 100)

def compute_nb_selected_pos(positions, percentage, positions_type, gff_types = ["gene"]):
    """Computes the number of position that will be kept from the positions file to get the right percentage."""
    in_handle = open(positions)
    if positions_type in ["gff", "gff3", "gtf"]:
        examiner = GFF.GFFExaminer()
        gff_type_counts = examiner.available_limits(in_handle)['gff_type']
        nb_pos_tot = 0
        for gff_type, count in gff_type_counts.items():
            nb_pos_tot += count if gff_type[0] in gff_types else 0
        return round(nb_pos_tot * percentage / 100)

    else:
        return round(len(set(in_handle.readlines())) * percentage / 100)

def split_gff(gff_path, outpath = None):
    """Splits a gff file into an ensemble of gff files, one for each gff id, usually chromosomes.
    Writes in a tmp folder created in the gff directory if no outpath is provided.
    Returns the dictionnary of each id associated to its gff path."""
    outdir = f"{'/'.join(gff_path.split('/')[:-1])}"
    folder_name = f"tmp_{'.'.join(gff_path.split('/')[-1].split('.')[:-1])}"
    outdir = (outdir + ("/" if len(outdir) > 0 else "") + folder_name) if not outpath else outpath + folder_name
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        shutil.rmtree(outdir)
        os.mkdir(outdir)
        
    chrom_gff_path = f"{outdir}/{gff_path.split('/')[-1][:-len('.gff')]}"
    chrom_files = {}
    chrom_files_path = {}
    with open(gff_path, 'r') as file:
        header = ""
        take_header = True
        while True:
            line = file.readline()
            if not line:
                break
            if line[0] == "#": 
                if take_header:
                    header +=line
            else:
                take_header = False
                chrom = line[:line.find("\t")]
                if chrom not in chrom_files:
                    chrom_path = f"{chrom_gff_path}_{chrom}.gff"
                    chrom_files_path[chrom] = chrom_path
                    chrom_files[chrom] = open(chrom_path, 'a')
                    chrom_files[chrom].write(header)
                chrom_files[chrom].write(line)
    for file in chrom_files.values():
        file.close()
    return chrom_files_path


### region mode methods

# extracting individual windows from the list of regions
def extract_window_region(cool, region1, region2, is_loc_circ1 = False, is_loc_circ2 = False, raw = False):
    """Extracts a submatrix of regions interaction zone from cool matrix.
    If required, adds paddings. Takes into account chromosome circularity."""
    matrix = cool.matrix(balance=(not raw))
    binning = cool.binsize

    if region1["Chromosome"] not in cool.chromsizes or region2["Chromosome"] not in cool.chromsizes:
        return None
    chrom_size1 = cool.chromsizes[region1["Chromosome"]]
    chrom_size2 = cool.chromsizes[region2["Chromosome"]]
    
    start1 = adjust_locus(region1["Padded_start"], chrom_size1)
    end1 = adjust_locus(region1["Padded_end"], chrom_size1)
    start2 = adjust_locus(region2["Padded_start"], chrom_size2)
    end2 = adjust_locus(region2["Padded_end"], chrom_size2)
    
    # 1. checking overflows
    region1_overflow = start1 >= end1
    region2_overflow = start2 >= end2

    # 2. computing intervales: assures the original submatrix includes the region center
    lower_overflow1 = chrom_size1 - start1 < end1
    lower_overflow2 = chrom_size2 - start2 < end2
    
    lower_interval1 = f"{region1['Chromosome']}:{start1 if not region1_overflow else 0}-{end1}" 
    higher_interval1 = f"{region1['Chromosome']}:{start1}-{end1 if not region1_overflow else chrom_size1}" 
    pos1 = lower_interval1 if lower_overflow1 else higher_interval1
    
    lower_interval2 = f"{region2['Chromosome']}:{start2 if not region2_overflow else 0}-{end2}" 
    higher_interval2 = f"{region2['Chromosome']}:{start2}-{end2 if not region2_overflow else chrom_size2}" 
    pos2 = lower_interval2 if lower_overflow2 else higher_interval2

    # 3. fetching main submatrix
    submatrix = matrix.fetch(pos1, pos2)[:]
    
    if (not region1_overflow and not region2_overflow): # no overflow
        return submatrix
    
    # 4. managing overflows:
    bins1_to_fill = (chrom_size1 // binning - start1 // binning) + 1 if lower_overflow1 else end1 // binning + (1 if end1 % binning != 0 else 0)
    bins2_to_fill = (chrom_size2 // binning - start2 // binning) + 1 if lower_overflow2 else end2 // binning + (1 if end2 % binning != 0 else 0)

    if bins1_to_fill == 0 and bins2_to_fill == 0:
        return submatrix
        
    if bins1_to_fill > 0:
        fill1 = np.array([chrom_size1 - (bins1_to_fill * binning), chrom_size1]) if lower_overflow1 else np.array([0, end1])
        fill1_pos = f"{region1['Chromosome']}:{fill1[0]}-{fill1[1]}"
    
    if bins2_to_fill > 0:
        fill2 = np.array([chrom_size2 - (bins2_to_fill * binning), chrom_size2]) if lower_overflow2 else np.array([0, end2])
        fill2_pos = f"{region2['Chromosome']}:{fill2[0]}-{fill2[1]}"
        
    mat1, mat2 = None, None
    if bins2_to_fill > 0:
        mat1 = matrix.fetch(pos1, fill2_pos) if is_loc_circ2 else np.full((submatrix.shape[0], bins2_to_fill), np.nan)
            
    if bins1_to_fill > 0:
        mat2 = matrix.fetch(fill1_pos, pos2) if is_loc_circ1 else np.full((bins1_to_fill, submatrix.shape[1]), np.nan)
        
    # two dimensions to fill
    if region1_overflow and region2_overflow:

        mat3 = matrix.fetch(fill1_pos, fill2_pos) if is_loc_circ1 and is_loc_circ2 else np.full((bins1_to_fill, bins2_to_fill), np.nan)

        concat1 = np.concatenate([mat1, submatrix], axis = 1) if lower_overflow2 else np.concatenate([submatrix, mat1], axis = 1)
        concat2 = np.concatenate([mat3, mat2], axis = 1) if lower_overflow2 else np.concatenate([mat2, mat3], axis = 1)
        submatrix = np.concatenate([concat2, concat1], axis = 0) if lower_overflow1 else np.concatenate([concat1, concat2], axis = 0)
        
    # dim1 to fill
    elif region1_overflow:        
        submatrix = np.concatenate([mat2, submatrix], axis = 0) if lower_overflow1 else np.concatenate([submatrix, mat2], axis = 0)
        
    # dim2 to fill
    elif region2_overflow:
        submatrix = np.concatenate([mat1, submatrix], axis = 1) if lower_overflow2 else np.concatenate([submatrix, mat1], axis = 1)
 
    return submatrix

def zoom_array(
    in_array,
    final_shape,
    same_sum=False,
    zoom_function=partial(zoom, order=1),
    **zoom_kwargs
):
    """Rescale an array or image.

    Normally, one can use scipy.ndimage.zoom to do array/image rescaling.
    However, scipy.ndimage.zoom does not coarsegrain images well. It basically
    takes nearest neighbor, rather than averaging all the pixels, when
    coarsegraining arrays. This increases noise. Photoshop doesn't do that, and
    performs some smart interpolation-averaging instead.

    If you were to coarsegrain an array by an integer factor, e.g. 100x100 ->
    25x25, you just need to do block-averaging, that's easy, and it reduces
    noise. But what if you want to coarsegrain 100x100 -> 30x30?

    Then my friend you are in trouble. But this function will help you. This
    function will blow up your 100x100 array to a 120x120 array using
    scipy.ndimage zoom Then it will coarsegrain a 120x120 array by
    block-averaging in 4x4 chunks.

    It will do it independently for each dimension, so if you want a 100x100
    array to become a 60x120 array, it will blow up the first and the second
    dimension to 120, and then block-average only the first dimension.

    (Copied from mirnylib.numutils)

    Parameters
    ----------
    in_array : ndarray
        n-dimensional numpy array (1D also works)
    final_shape : shape tuple
        resulting shape of an array
    same_sum : bool, optional
        Preserve a sum of the array, rather than values. By default, values
        are preserved
    zoom_function : callable
        By default, scipy.ndimage.zoom with order=1. You can plug your own.
    **zoom_kwargs :
        Options to pass to zoomFunction.

    Returns
    -------
    rescaled : ndarray
        Rescaled version of in_array

    """
    in_array = np.asarray(in_array, dtype=np.double)
    in_shape = in_array.shape
    assert len(in_shape) == len(final_shape)
    mults = []  # multipliers for the final coarsegraining
    for i in range(len(in_shape)):
        if final_shape[i] < in_shape[i]:
            mults.append(int(np.ceil(in_shape[i] / final_shape[i])))
        else:
            mults.append(1)
    # shape to which to blow up
    temp_shape = tuple([i * j for i, j in zip(final_shape, mults)])

    # stupid zoom doesn't accept the final shape. Carefully crafting the
    # multipliers to make sure that it will work.
    zoom_multipliers = np.array(temp_shape) / np.array(in_shape) + 0.0000001
    assert zoom_multipliers.min() >= 1

    # applying scipy.ndimage.zoom
    rescaled = zoom_function(in_array, zoom_multipliers, **zoom_kwargs)

    for ind, mult in enumerate(mults):
        if mult != 1:
            sh = list(rescaled.shape)
            assert sh[ind] % mult == 0
            newshape = sh[:ind] + [sh[ind] // mult, mult] + sh[ind + 1 :]
            rescaled.shape = newshape
            rescaled = np.mean(rescaled, axis=ind + 1)
    assert rescaled.shape == final_shape

    if same_sum:
        extra_size = np.prod(final_shape) / np.prod(in_shape)
        rescaled /= extra_size
    return rescaled

def resize_window(submatrix, expected_size = 100):
    """Resizes a submatrix to match the shape expected_size x expected_size."""
    return zoom_array(submatrix, (expected_size, expected_size))

def get_dist_to_diag_matrix(cool, region1, region2):
    """Computes the index matrix to use for p(s) detrending and diagonal masking from regions."""
    binning = cool.binsize
    chrom1_bin = cool.chromsizes[region1["Chromosome"]] // binning
    chrom2_bin = cool.chromsizes[region2["Chromosome"]] // binning

    start1 = adjust_locus(region1["Padded_start"], cool.chromsizes[region1["Chromosome"]])
    end1 = adjust_locus(region1["Padded_end"], cool.chromsizes[region1["Chromosome"]])
    start2 = adjust_locus(region2["Padded_start"], cool.chromsizes[region2["Chromosome"]])
    end2 = adjust_locus(region2["Padded_end"], cool.chromsizes[region2["Chromosome"]])
    
    # Computing the origin bin (up left bin) of the distance matrix value k
    start1 = start1 // binning
    start2 = start2 // binning
    k = abs(start2 - start1)
    
    # Computing matrix dimensions
    end1 = end1 // binning if end1 % binning != 0 else (end1 - 1) // binning
    end2 = end2 // binning if end2 % binning != 0 else (end2 - 1) // binning
    dim0 = end1 - start1 + 1 if start1 <= end1 else (chrom1_bin - start1) + end1 + 2
    dim1 = end2 - start2 + 1 if start2 <= end2 else (chrom2_bin - start2) + end2 + 2
    
    # Filling the dist matrix
    dist_array = np.empty((dim0,  dim1))
    if region1["Chromosome"] != region2["Chromosome"]:
        return dist_array
    dist_array[0, 0] = k
    for i in range(dim0):
        for j in range(dim1):
            dist_array[i, j] = abs(k - i + j)
    return dist_array.astype(int)
    
def compute_chromosome_diag_mask(chromsize, diagonal_masking, binning, is_circular = False):
    """Computes the diagonal mask for a given chromosome, accounting for circularity."""
    diag_mask = np.ones(chromsize // binning + 1)
    masking_bins = diagonal_masking // binning + 1
    diag_mask[:masking_bins] = np.nan
    if is_circular and masking_bins > 1:
        diag_mask[-masking_bins + 1:] = np.nan
    diag_mask = np.concatenate([diag_mask, diag_mask])
    return diag_mask