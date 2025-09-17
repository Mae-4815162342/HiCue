from .cli.imports import *

def get_random_from_locus(cool, locus, nb_pos=2, max_dist=100000):
    """From locus list, computes a list of nb_pos random locus for each loci"""
    random_loci = pd.DataFrame(columns=locus.columns)

    for i in locus.index:
        loci = locus.iloc[i]

        # TODO: include circularity in random picking
        min_pos = max(0, loci["Start"] - max_dist)
        max_pos = min(cool.chromsizes[loci["Chromosome"]], loci["End"] + max_dist)

        for _ in range(nb_pos):

            random_pos = int(np.random.random() * abs(max_pos - min_pos) + min_pos)

            loci["Start"] = random_pos
            loci["End"] = random_pos

            random_loci = random_loci._append(loci, ignore_index = True)

    return random_loci

def get_random_from_locus_2d(cool, locus, nb_pos=2, max_dist=100000):
    """From 2d locus list, computes a list of nb_pos random pair of locus for each loci pair, maintaining the distance between those."""
    random_loci = pd.DataFrame(columns=locus.columns)

    for i in locus.index:
        loci = locus.iloc[i].copy(deep=True)

        min_pos = max(0, loci["Start"] - max_dist)
        max_pos = min(cool.chromsizes[loci["Chromosome"]], loci["End"] + max_dist)

        for _ in range(nb_pos):

            random_pos = int(np.random.random() * abs(max_pos - min_pos) + min_pos)
            loci_dist = abs(random_pos - loci["Start"])

            # Two cases: 
            ## the second location is on the same chromosome, the distance is between the loci physically on the DNA strand;
            ## the second location is on another chromosome, but it is the same operation, as we try to extract the random matrices on the same diagonal as the original positions in the trans contact matrix.
            if random_pos <= loci["Start"]:
                random_pos2 = max(0, loci["Start2"] - loci_dist)
            else: 
                random_pos2 = min(cool.chromsizes[loci["Chromosome2"]], loci["Start2"] + loci_dist)


            loci["Start"] = random_pos
            loci["End"] = random_pos
            loci["Start2"] = random_pos2
            loci["End2"] = random_pos2

            random_loci = random_loci._append(loci, ignore_index = True)

    return random_loci

def detrend(matrix):
    """Applies detrending by p(s) to a chromosome's matrix"""
    y = distance_law(matrix)
    y[np.isnan(y)] = 0.0
    matrix_index = [[abs(i - j) for i in range(len(matrix))] for j in range(len(matrix))]
    matrix = matrix / y[matrix_index]
    return matrix

def extract_window(cool, locus1, locus2, binning, window, circular=[], trans=False, diagonal_mask=0, center="start", detrend_matrix = False, raw = False):
    """Extracts a window from a matix given positions and parameters."""
    if trans and locus1["Chromosome"]!=locus2["Chromosome"]:
        extent1 = cool.extent(locus1["Chromosome"])
        extent2 = cool.extent(locus2["Chromosome"])
        matrix = cool.matrix(balance=(not raw))[extent1[0]: extent1[1], extent2[0]:extent2[1]]
        if detrend_matrix:
            matrix =matrix / np.nanmean(matrix)
    else:
        matrix = detrend(cool.matrix(balance=(not raw)).fetch(locus1["Chromosome"])) if detrend_matrix else cool.matrix(balance=(not raw)).fetch(locus1["Chromosome"])

    if locus1["Chromosome"] == locus2["Chromosome"]:
        for i in range(diagonal_mask//binning):
            np.fill_diagonal(matrix[i:],  np.nan)
            np.fill_diagonal(matrix[:,i:],  np.nan)

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

    pos1 = start1 // binning
    pos2 = start2 // binning
    window_binned = window // binning

    pos1_sub = np.full((window_binned * 2 + 1, len(matrix[0])), np.nan)
    if pos1 - window_binned < 0:
        if locus1["Chromosome"] in circular:
            pos1_sub = np.concatenate([
                matrix[pos1 - window_binned:],
                matrix[:pos1 + window_binned + 1]
            ], axis=0)
        else:
            pos1_sub[-(pos1 + window_binned + 1):] = matrix[:pos1 + window_binned + 1]
    elif pos1 + window_binned + 1 > len(matrix):
        if locus1["Chromosome"] in circular:
            pos1_sub = np.concatenate([
                matrix[pos1 - window_binned:],
                matrix[:pos1 + window_binned + 1 - len(matrix)]
            ], axis=0)
        else:
            pos1_sub[:-(pos1 + window_binned + 1 - len(matrix))] = matrix[pos1 - window_binned:]
    else:
        pos1_sub = matrix[pos1 - window_binned: pos1 + window_binned + 1]


    submatrix = np.full((window_binned * 2 + 1, window_binned * 2 +1), np.nan)
    if pos2 - window_binned < 0:
        if locus2["Chromosome"] in circular:
            submatrix = np.concatenate([
                pos1_sub[:, pos2 - window_binned:],
                pos1_sub[:, :pos2 + window_binned + 1]
            ], axis=1)
        else:
            submatrix[:, -(pos2 + window_binned + 1):] = pos1_sub[:, :pos2 + window_binned + 1]
    elif pos2 + window_binned + 1 > len(matrix[0]):
        if locus2["Chromosome"] in circular:
            submatrix = np.concatenate([
                pos1_sub[:, pos2 - window_binned:],
                pos1_sub[:, :pos2 + window_binned + 1 - len(matrix[0])]
            ], axis=1)
        else:
            submatrix[:, :-(pos2 + window_binned + 1 - len(matrix[0]))] = pos1_sub[:, pos2 - window_binned:]
    else:
        submatrix = pos1_sub[:, pos2 - window_binned: pos2 + window_binned + 1]
    
    return submatrix

def compute_subtracks(bw_tracks, positions, window, center = 'start', circular=[]):
    """Returns a dictionnary with each position's index as key and the corresponding track on the window."""
    chrom_sizes = bw_tracks.chroms()
    subtracks = {}
    for i, locus in positions.iterrows():
        start = locus["Start"]
        match center:
            case "start":
                start = min(locus["Start"], locus["End"]) if locus["Strand"] == 1 else max(locus["Start"], locus["End"])
            case "center":
                start = (locus["Start"] + locus["End"]) // 2
            case "end":
                start = max(locus["Start"], locus["End"]) if locus["Strand"] == 1 else min(locus["Start"], locus["End"])

        chromosome = locus["Chromosome"]
        len_chrom = chrom_sizes[chromosome]

        extracted_track = []
        if start - window < 0:
            bellow_start = abs(start - window)
            extracted_track = np.concatenate([[np.nan] * bellow_start, bw_tracks.values(chromosome, 0, start + window + 1)])
            if chromosome in circular: 
                extracted_track[:bellow_start] = bw_tracks.values(chromosome, len_chrom - bellow_start, len_chrom)

        elif start + window + 1 > len_chrom:
            up_start = start + window + 1 - len_chrom
            extracted_track = np.concatenate([bw_tracks.values(chromosome, start - window, len_chrom), [np.nan] * up_start])
            if chromosome in circular:
                extracted_track[- up_start:] = bw_tracks.values(chromosome, 0, up_start)

        else:
            extracted_track = np.array(bw_tracks.values(chromosome, start - window, start + window + 1))

        subtracks[i] = extracted_track
    return subtracks

def get_dist_positions(positions, index1, index2, center="start"):
     
    locus1 = positions.iloc[index1]
    locus2 = positions.iloc[index2]
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

    if locus1["Chromosome"] == locus2["Chromosome"]:
        return abs(start1 - start2)
    else:
        return np.inf
    
def compute_pairs2d(positions):
    """Separates the 2d lines of a position table into single positions, keeping tab on the pairs indexes. Both the new table and indexation are returned."""
    tmp_positions = pd.DataFrame(columns=["Name", "Chromosome", "Start", "End", "Strand"])
    index_pairs = []
    for i, pos in positions.iterrows():
        tmp_positions = tmp_positions._append({
            "Name":pos["Name"],
            "Chromosome":pos["Chromosome"],
            "Start":pos["Start"],
            "End":pos["End"],
            "Strand":pos["Strand"]
        }, ignore_index = True)
        tmp_positions = tmp_positions._append({
            "Name":pos["Name2"],
            "Chromosome":pos["Chromosome2"],
            "Start":pos["Start2"],
            "End":pos["End2"],
            "Strand":pos["Strand2"]
        }, ignore_index = True)
        index_pairs.append((i * 2, i * 2 + 1))
    return tmp_positions, index_pairs

def compute_submatrices(cool, name, positions, binning, window, circular=[], loops = False, min_dist=0, trans_contact=False, diagonal_mask=0, center="start", sort_contact="None", contact_range="20000:100000:30000", ps_detrend = False, is_2d=False, raw = False):
    locus_pairs = {}
    tmp_locus_pairs = []
    tmp_positions = None
    # if 2d positions, the pairs are pre-computed
    if is_2d:
        tmp_positions, tmp_locus_pairs = compute_pairs2d(positions)
    else:
        tmp_positions = positions.copy()
    if loops:
        tmp_locus_pairs = np.array(list(combinations(tmp_positions.index, r=2))) if not is_2d else tmp_locus_pairs
        tmp_locus_pairs_distances = np.array([get_dist_positions(tmp_positions, i, j, center=center) for i, j in tmp_locus_pairs])
        match sort_contact:
            case "none":
                locus_pairs[name] = tmp_locus_pairs[tmp_locus_pairs_distances >= min_dist]

            case "None":
                locus_pairs[name] = tmp_locus_pairs[tmp_locus_pairs_distances >= min_dist]

            case "":
                locus_pairs[name] = tmp_locus_pairs[tmp_locus_pairs_distances >= min_dist]

            case "distance":
                contact_data = contact_range.split(':')
                contact_min = int(contact_data[0]) if len(contact_data) == 3 else 0
                contact_max = int(contact_data[1]) if len(contact_data) == 3 else 0
                contact_window = int(contact_data[2]) if len(contact_data) == 3 else 0
                tmp_locus_pairs = np.array(list(combinations(tmp_positions.index, r=2)))
                tmp_locus_pairs_distances = [get_dist_positions(tmp_positions, i, j, center=center) for i, j in tmp_locus_pairs]
                for k in range(contact_min, contact_max, contact_window):
                    min_window = k
                    max_window = k + contact_window
                    k_indexes = [i for i in range(len(tmp_locus_pairs_distances)) if min_window <= tmp_locus_pairs_distances[i] < max_window]
                    locus_pairs[f"{name}_{min_window//1000}-{max_window//1000}kb"] = tmp_locus_pairs[k_indexes]

            case "cis_trans":
                cis_list = []
                trans_list = []
                for i in range(len(tmp_positions)):
                    for j in range(i +1, len(tmp_positions)):
                        if get_dist_positions(tmp_positions, i, j, center=center) <= min_dist:
                            continue
                        if tmp_positions.iloc[i]["Chromosome"] == tmp_positions.iloc[j]["Chromosome"]:
                            cis_list.append([i, j])
                        else:
                            trans_list.append([i, j])
                locus_pairs[f"{name}_cis"] = cis_list
                locus_pairs[f"{name}_trans"] = trans_list
    else:
        locus_pairs[name] = list([(i, i) for i in tmp_positions.index]) if not is_2d else tmp_locus_pairs

    all_submatrices = {}
    for locus_name in locus_pairs.keys():

        submatrices = pd.DataFrame(columns=['Loc1', 'Loc2', 'Matrix'])
        for i, j in locus_pairs[locus_name]:
            locus1 = tmp_positions.iloc[i]
            locus2 = tmp_positions.iloc[j]
            if locus1["Chromosome"] != locus2["Chromosome"] and not trans_contact:
                continue
            submatrix= extract_window(cool, locus1, locus2, binning, window, circular = circular, trans = trans_contact, diagonal_mask=diagonal_mask, center=center, detrend_matrix = ps_detrend, raw = raw)
            submatrices = submatrices._append({"Loc1":i, "Loc2": j, "Matrix":submatrix.flatten()}, ignore_index=True)
        all_submatrices[locus_name] = submatrices

    return all_submatrices

def get_windows(matrices, locus, flip, fill=None):
    """Returns windows for pileup flipped if necessary. fill parameter will replace zero and nan values if provided."""
    windows = []
    for _, row in matrices.iterrows():
        i, j, matrix = row["Loc1"], row["Loc2"], row["Matrix"]
        if fill != None:
            matrix[np.isnan(matrix)] = fill
            matrix[matrix <= fill] = fill
        if i == j and flip and locus.iloc[i]['Strand'] == -1:
            windows.append(np.flip(matrix))
        else:
            windows.append(matrix)
    return np.array(windows)

def is_in_region(positions, regions, region, overlap="flex"):
    """Returns the positions from positions list present in the selected region"""
    selected_regions = regions[regions["Id"] == region]
    indexes = np.zeros(len(positions), dtype=bool)
    for i in range(len(positions)):
        position = positions.iloc[i]
        pos1 = min(position["Start"], position["End"])
        pos2 = max(position["Start"], position["End"])
        chromosome = position["Chromosome"]
        for _, reg in selected_regions[selected_regions["Chromosome"] == chromosome].iterrows():
            pos1_reg = min(reg["Start"], reg["End"])
            pos2_reg = max(reg["Start"], reg["End"])
            match overlap:
                case "flex":
                    # included if the position overlaps the regions, even not completely
                    if (pos1 >= pos1_reg and pos1 <= pos2_reg) or (pos2 >= pos1_reg and pos2 <= pos2_reg):
                        indexes[i] = True
                        break
                case "strict":
                    # included only if the position is completely in the region
                    if pos1 >= pos1_reg and pos2 <= pos2_reg:
                        indexes[i] = True
                        break
    return indexes

def sum_up_separate_by(positions_tables, positions, outpath):
    """Recapitulates the separations in a csv file"""
    with open(f"{outpath}/separate_by.csv", "w") as file:
        file.write(f"## Selected positions after separate_by ({sum([len(table) for table in positions_tables.values()])}/{len(positions)} selected)\n")
        for name in positions_tables.keys():
            file.write(f"# {name} {len(positions_tables[name])}\n")
        for name in positions_tables.keys():
            for _, position in positions_tables[name].iterrows():
                file.write(",".join([
                    name,
                    position["Name"],
                    str(position["Start"]),
                    str(position["End"]),
                    str(position["Strand"]),
                    position["Chromosome"]
                ]) + "\n")  

def separate_positions(positions, name, separate_by="", separate_regions="", overlap="flex", outpath=""):
    """Separates a table of positions according to the separate_by parameter. Allows several separations simultaneously."""
    positions_tables = {name:positions}
    sum_up = False
    for separation in separate_by:
        positions_tmp = {}
        match separation:                
            case "direct":
                for pos_name in positions_tables.keys():
                    current_positions = positions_tables[pos_name]
                    forwards = current_positions[current_positions["Strand"] == 1]
                    reverses = current_positions[current_positions["Strand"] == -1]
                    positions_tmp[f"{pos_name}_forward"] = forwards
                    positions_tmp[f"{pos_name}_reverse"] = reverses
                sum_up = True
                    
            case "regions":
                regions = pd.read_csv(separate_regions)
                if len(regions) > 0:
                    for pos_name in positions_tables.keys():
                        current_positions = positions_tables[pos_name]
                        for region in np.unique(regions["Id"]):
                            positions_tmp[f"{pos_name}_{region}"] = current_positions[is_in_region(current_positions, regions, region, overlap=overlap)]
                    sum_up = True
                else:
                    positions_tmp = positions_tables
                        

            case "chroms":
                for pos_name in positions_tables.keys():
                    current_positions = positions_tables[pos_name]
                    for chrom in np.unique(current_positions["Chromosome"]):
                        positions_tmp[f"{pos_name}_{chrom}"] = current_positions[current_positions["Chromosome"] == chrom]
                sum_up = True
            
            # contact functions are implemented in window selection
            case "cis_trans":
                positions_tmp = positions_tables
            case "distance":
                positions_tmp = positions_tables
            case "None":
                positions_tmp = positions_tables
            case "":
                positions_tmp = positions_tables

        for name in positions_tmp:
            positions_tmp[name].index = range(len(positions_tmp[name]))
        positions_tables = positions_tmp

    if len(outpath) > 0 and sum_up:
        sum_up_separate_by(positions_tables, positions, outpath) 
    return positions_tables

def separate_positions_2d(positions, name, separate_by="", separate_regions="", overlap="flex", outpath=""):
    """Separates a table of positions' pairs according to the separate_by parameter. Allows several separations simultaneously."""
    positions_tables = {name:positions}
    sum_up = False
    for separation in separate_by:
        positions_tmp = {}
        match separation:                
            case "direct": # for 2d separation if the strand has been provided, will separate in 4 categories: for/for, for/rev, rev/for, rev/rev
                for pos_name in positions_tables.keys():
                    current_positions = positions_tables[pos_name]
                    forwards = current_positions[current_positions["Strand"] == 1]
                    reverses = current_positions[current_positions["Strand"] == -1]
                    positions_tmp[f"{pos_name}_forward"] = forwards
                    positions_tmp[f"{pos_name}_reverse"] = reverses
                sum_up = True
                    
            case "regions": # will keep pairs in regions if both elements are in the region (with the overlap parameter)
                regions = pd.read_csv(separate_regions)
                if len(regions) > 0:
                    for pos_name in positions_tables.keys():
                        current_positions = positions_tables[pos_name]
                        for region in np.unique(regions["Id"]):
                            positions_tmp[f"{pos_name}_{region}"] = current_positions[is_in_region(current_positions, regions, region, overlap=overlap)]
                    sum_up = True
                else:
                    positions_tmp = positions_tables
                        

            case "chroms":
                for pos_name in positions_tables.keys():
                    current_positions = positions_tables[pos_name]
                    for chrom in np.unique(current_positions["Chromosome"]):
                        positions_tmp[f"{pos_name}_{chrom}"] = current_positions[current_positions["Chromosome"] == chrom]
                sum_up = True
            
            # contact functions are implemented in window selection
            case "cis_trans":
                positions_tmp = positions_tables
            case "distance":
                positions_tmp = positions_tables
            case "None":
                positions_tmp = positions_tables
            case "":
                positions_tmp = positions_tables

        for name in positions_tmp:
            positions_tmp[name].index = range(len(positions_tmp[name]))
        positions_tables = positions_tmp

    if len(outpath) > 0 and sum_up:
        sum_up_separate_by(positions_tables, positions, outpath) 
    return positions_tables

def compile_tracks(positions, tracks, flip = False, method="median"):
    """Compiles tracks according to selected method."""
    selected_tracks = []
    for i, pos in positions.iterrows():
        track = tracks[i]
        if flip and pos["Strand"] == -1:
            track = np.flip(track)
        selected_tracks.append(track)
    selected_tracks = np.array(selected_tracks)
    match method:
        case "median":
            track_pileup = np.nanmedian(selected_tracks, axis=0)
        case "mean":
            track_pileup = np.nanmean(selected_tracks, axis=0)

    return track_pileup