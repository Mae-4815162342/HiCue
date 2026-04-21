from hicue.displays import *
from hicue.classes.Reader import *
from hicue.classes.TrackReader import *
from hicue.classes.PairFormater import *
from hicue.classes.RandomSelector import *
from hicue.classes.MatrixExtractorLauncher import *
from hicue.classes.AsyncDisplays import *

def extract(cool_files, positions, outpath, log = None, **params):

    if not os.path.exists(outpath):
        create_folder_path(outpath)

    format_params = {
        "separate_by" : params["separate_by"],
        "center" : params["center"],
        "contact_range" : params["contact_range"],
        "separate_regions" : params["separation_regions"],
        "min_dist" : params['min_dist'],
        "detrending": params['detrending'],
        "diag_mask": params['diag_mask'],
        "overlap": params['overlap'],
        "has_trans": params["trans"],
        "circulars": params['circulars']
    }

    display_args = {
        "output_format": params["format"],
        "track_unit" : params['track_unit'],
        "display_strand" : params["display_strand"], 
        "display_sense" : params["display_sense"],
        "flipped": params["flip"],
        "cmap": params["indiv_cmap_limits"],
        "color": params["indiv_cmap_color"]
    }

    random_params = {
        "center" : params["center"],
        "selection_window": params['rand_max_dist'],
        "nb_rand_per_pos" : params['nb_pos'],
        "random_jitter": params['random_jitter']
    }

    pileup_display_args = {
        "track_unit" : params['track_unit'],
        "output_format": params["format"],
        "display_strand" : params["display_strand"], 
        "display_sense" : params["display_sense"],
        "flipped": params["flip"],
        "cmap": params["cmap_limits"],
        "cmap_color": params["cmap_color"]
    }

    # checking multiprocessing values
    threads = max(1, params["threads"])

    # reading parameters
    pos_type, pos_file = positions
    data_title = pos_file.split('/')[-1].split('.')[0].replace('/', '_')

    ## Reading file and annotating
    annotation = {}
    if params['gff']:
        annotation['gff'] = params['gff']
    if params['tracks']:
        annotation['tracks'] = params['tracks']
    
    reader = Reader(pos_file, pos_type, annotation_files = annotation, save_to = "", loop = params['loops'], record_type = params['record_type'], overlap = params['overlap'])# TODO: option on verbose annotation
    positions, pairing_queue = reader.read_file(threads = threads)

    ## Formating indexes pairs
    formater = PairFormater(positions, **format_params)
    formated_pairs = formater.format_pairs(pairing_queue, threads = threads)

    if len(formated_pairs) == 0:
        log.write('Empty separation: No position nor pair of positions is matching any possible separation.')   
        return

    if params["save_tmp"]:
        positions.to_csv(f"{outpath}/{data_title}_positions.csv")
        formated_pairs.to_csv(f"{outpath}/{data_title}_formated_pairs.csv")

    ## Random locus selection (for patch only) from positions
    if params['detrending'] == "patch" and params["pileup"]:

        # if random files are provided
        random_pos_path = f"{params['random_path']}_random_positions.csv"
        random_pairs_path = f"{params['random_path']}_random_pairs.csv"
        if os.path.exists(random_pos_path) and os.path.exists(random_pairs_path):
            log.write(f'Using randoms found at {params["random_path"]} for patch detrending.')

            random_positions = pd.read_csv(random_pos_path, header=0, index_col=0)
            random_pairs = pd.read_csv(random_pairs_path, header=0, index_col=0)

        else:
            log.write(f'Generating random positions for patch detrending.')
            selector = RandomSelector(positions, **random_params)
            random_positions, random_pairs = selector.select_randoms(formated_pairs, threads = threads)

            if params["save_tmp"]:
                random_positions.to_csv(f"{outpath}/{data_title}_random_positions.csv")
                random_pairs.to_csv(f"{outpath}/{data_title}_random_pairs.csv")
    
    ## Matrix extraction
    matrix_extractor = MatrixExtractorLauncher(cool_files,
                                               nb_pos = np.max(positions.index),
                                               compute_pileups = params["pileup"],
                                               binnings = params["binnings"], 
                                               windows = params["windows"],
                                               center = params["center"], 
                                               raw = params["raw"], 
                                               method = params["method"], 
                                               flip = params["flip"],
                                               nb_rand_per_pos = params["nb_pos"],
                                               display_loci = params["loci"],
                                               display_batch = params["batch"],
                                               outpath = outpath,
                                               display_args = display_args,
                                               log = log)
    
    pileups = matrix_extractor.launch_extraction(positions, formated_pairs, threads=threads)

    if params["pileup"]:
        pileups_random = {}
        if params['detrending'] == "patch":
            pileups_random_queue = matrix_extractor.launch_extraction(random_positions, random_pairs, randoms = True, threads=threads)
            pileups_random = empty_queue_in_dict(pileups_random_queue, keys = ["sep_id", "binning", "cool_name"]) # exporting the patch detrending as an dict for access

        ## Pileup detrending and display
        pileup_display_args["display_strand"] = pileup_display_args["display_strand"] and (0 not in positions["Strand"])
        pileup_display = Display(
            input_queue = pileups,
            output_queues = [],
            function = display_pileup,
            patch_detrending = pileups_random,
            outpath = outpath,
            title = data_title,
            size_metrics = params["windows"],
            is_contact = (pos_type == "bed2d" or params["loops"]),
            **pileup_display_args
        )
        pileup_display.join()
        
def tracks(cool_files, tracks, outpath, log = None, **params):

    if not os.path.exists(outpath):
        create_folder_path(outpath)

    # tracks_params
    tracks_params = {
        "threshold": params['threshold'],
        "percentage": params['percentage'],
        "min_sep": params['min_sep']
    }

    format_params = {
        "separate_by" : params["separate_by"],
        "center" : params["center"],
        "contact_range" : params["contact_range"],
        "separate_regions" : params["separation_regions"],
        "min_dist" : params['min_dist'],
        "detrending": params['detrending'],
        "diag_mask": params['diag_mask'],
        "overlap": params['overlap'],
        "has_trans": params["trans"],
        "circulars": params['circulars'],
        "loop": params['loops']
    }

    display_args = {
        "track_unit" : params['track_unit'],
        "output_format": params["format"],
        "display_strand" : params["display_strand"], 
        "display_sense" : params["display_sense"],
        "flipped": params["flip"],
        "cmap": params["indiv_cmap_limits"],
        "color": params["indiv_cmap_color"],
        "display_tracks": True
    }

    random_params = {
        "center" : params["center"],
        "selection_window": params['rand_max_dist'],
        "nb_rand_per_pos" : params['nb_pos'],
        "random_jitter": params['random_jitter']
    }

    pileup_display_args = {
        "track_unit" : params['track_unit'],
        "output_format": params["format"],
        "display_strand" : params["display_strand"], 
        "display_sense" : params["display_sense"],
        "flipped": params["flip"],
        "cmap": params["cmap_limits"],
        "cmap_color": params["cmap_color"]
    }

    pos_type, pos_file = params['positions'] or (None, None)

    # checking multiprocessing values
    threads = max(1, params["threads"])

    # reading parameters
    data_title = tracks.split('/')[-1].split('.')[0].replace('/', '_')

    ## Reading tracks file    
    reader = TrackReader(tracks, positions_file = pos_file, position_type = pos_type, save_to = "", loop = params['loops'], record_type = params['record_type'], overlap = params['overlap'], **tracks_params) # TODO: option on verbose annotation
    positions, pair_queue = reader.read_file(threads = threads)

    ## Formating indexes pairs
    formater = PairFormater(positions, **format_params)
    formated_pairs = formater.format_pairs(pair_queue=pair_queue, threads = threads)

    if len(formated_pairs) == 0:
        log.write('Empty separation: No position nor pair of positions is matching any possible separation.')   
        return

    if params["save_tmp"]:
        positions.to_csv(f"{outpath}/{data_title}_positions.csv")
        formated_pairs.to_csv(f"{outpath}/{data_title}_formated_pairs.csv")

    ## Random locus selection (for patch only) from positions
    if params['detrending'] == "patch" and params["pileup"]:

        # if random files are provided
        random_pos_path = f"{params['random_path']}_random_positions.csv"
        random_pairs_path = f"{params['random_path']}_random_pairs.csv"
        if os.path.exists(random_pos_path) and os.path.exists(random_pairs_path):
            log.write(f'Using randoms found at {params["random_path"]} for patch detrending.\n')

            random_positions = pd.read_csv(random_pos_path, header=0, index_col=0)
            random_pairs = pd.read_csv(random_pairs_path, header=0, index_col=0)

        else:
            log.write(f'Generating random positions for patch detrending.\n')
            selector = RandomSelector(positions, **random_params)
            random_positions, random_pairs = selector.select_randoms(formated_pairs, threads = threads)

            if params["save_tmp"]:
                random_positions.to_csv(f"{outpath}/{data_title}_random_positions.csv")
                random_pairs.to_csv(f"{outpath}/{data_title}_random_pairs.csv")
    
    ## Matrix extraction
    matrix_extractor = MatrixExtractorLauncher(cool_files,
                                               tracks = tracks,
                                               nb_pos = np.max(positions.index) + 1,
                                               compute_pileups = params["pileup"],
                                               binnings = params["binnings"], 
                                               windows = params["windows"],
                                               center = params["center"], 
                                               raw = params["raw"], 
                                               method = params["method"], 
                                               flip = params["flip"],
                                               nb_rand_per_pos = params["nb_pos"],
                                               display_loci = params["loci"],
                                               display_batch = params["batch"],
                                               outpath = outpath,
                                               display_args = display_args,
                                               log = log)
    
    pileups = matrix_extractor.launch_extraction(positions, formated_pairs, threads=threads)

    if params["pileup"]:
        pileups_random = {}
        if params['detrending'] == "patch":
            pileups_random_queue = matrix_extractor.launch_extraction(random_positions, random_pairs, randoms = True, threads=threads)
            pileups_random = empty_queue_in_dict(pileups_random_queue, keys = ["sep_id", "binning", "cool_name"]) # exporting the patch detrending as an dict for access

        ## Pileup detrending and display
        pileup_display_args["display_strand"] = pileup_display_args["display_strand"] and (0 not in positions["Strand"])
        pileup_display = Display(
            input_queue = pileups,
            output_queues = [],
            function = display_pileup,
            patch_detrending = pileups_random,
            outpath = outpath,
            title = data_title,
            size_metrics = params["windows"],
            is_contact = params["loops"] or pos_type == "bed2d",
            **pileup_display_args
        )
        pileup_display.join()

def regions(cool_files, positions, outpath, log = None, **params):

    if not os.path.exists(outpath):
        create_folder_path(outpath)

    format_params = {
        "separate_by" : params["separate_by"],
        "center" : params["center"],
        "contact_range" : params["contact_range"],
        "separate_regions" : params["separation_regions"],
        "min_dist" : params['min_dist'],
        "detrending": params['detrending'],
        "diag_mask": params['diag_mask'],
        "overlap": params['overlap'],
        "has_trans": params["trans"],
        "circulars": params['circulars']
    }

    display_args = {
        "output_format": params["format"],
        "track_unit" : params['track_unit'],
        "display_strand" : params["display_strand"], 
        "display_sense" : params["display_sense"],
        "display_log": params['display_log'],
        "flipped": params["flip"],
        "cmap": params["indiv_cmap_limits"],
        "color": params["indiv_cmap_color"]
    }

    random_params = {
        "center" : params["center"],
        "selection_window": params['rand_max_dist'],
        "nb_rand_per_pos" : params['nb_pos'],
        "random_jitter": params['random_jitter'],
    }

    pileup_display_args = {
        "track_unit" : params['track_unit'],
        "output_format": params["format"],
        "display_strand" : params["display_strand"], 
        "display_sense" : params["display_sense"],
        "display_log": params['display_log'],
        "flipped": params["flip"],
        "cmap": params["cmap_limits"],
        "cmap_color": params["cmap_color"]
    }

    # checking multiprocessing values
    threads = max(1, params["threads"])

    # reading parameters
    pos_type, pos_file = positions
    data_title = pos_file.split('/')[-1].split('.')[0].replace('/', '_')

    ## Reading file and annotating
    annotation = {}
    if params['gff']:
        annotation['gff'] = params['gff']
    if params['tracks']:
        annotation['tracks'] = params['tracks']

    reader = Reader(pos_file, pos_type, annotation_files = annotation, save_to = "", padding = params['padding'], loop = params['loops'], record_type = params['record_type'], min_region_size = params['min_region_size'], overlap = params['overlap'])# TODO: option on verbose annotation
    positions, pairing_queue = reader.read_file(threads = threads)

    ## Formating indexes pairs
    formater = PairFormater(positions, **format_params)
    formated_pairs = formater.format_pairs(pairing_queue, threads = threads)

    if len(formated_pairs) == 0:
        log.write('Empty separation: No position nor pair of positions is matching any possible separation.')   
        return

    if params["save_tmp"]:
        positions.to_csv(f"{outpath}/{data_title}_positions.csv")
        formated_pairs.to_csv(f"{outpath}/{data_title}_formated_pairs.csv")

    # Random locus selection (for patch only) from positions
    if params['detrending'] == "patch" and params["pileup"]:

        # if random files are provided
        random_pos_path = f"{params['random_path']}_random_positions.csv"
        random_pairs_path = f"{params['random_path']}_random_pairs.csv"
        if os.path.exists(random_pos_path) and os.path.exists(random_pairs_path):
            log.write(f'Using randoms found at {params["random_path"]} for patch detrending.')

            random_positions = pd.read_csv(random_pos_path, header=0, index_col=0)
            random_pairs = pd.read_csv(random_pairs_path, header=0, index_col=0)

        else:
            log.write(f'Generating random positions for patch detrending.')
            selector = RandomSelector(positions, is_region = True, padding = params["padding"], **random_params)
            random_positions, random_pairs = selector.select_randoms(formated_pairs, threads = threads)

            if params["save_tmp"]:
                random_positions.to_csv(f"{outpath}/{data_title}_random_positions.csv")
                random_pairs.to_csv(f"{outpath}/{data_title}_random_pairs.csv")
    
    ## Matrix extraction
    matrix_extractor = MatrixExtractorLauncher(cool_files,
                                               nb_pos = np.max(positions.index),
                                               compute_pileups = params["pileup"],
                                               binnings = params["binnings"], 
                                               windows = params["windows"],
                                               center = params["center"], 
                                               raw = params["raw"], 
                                               method = params["method"], 
                                               flip = params["flip"],
                                               nb_rand_per_pos = params["nb_pos"],
                                               display_loci = params["loci"],
                                               display_batch = params["batch"],
                                               outpath = outpath,
                                               display_args = display_args,
                                               resizing = params["expected_sizes"],
                                               padding = params['padding'],
                                               extract_regions = True,
                                               log = log)
    
    pileups = matrix_extractor.launch_extraction(positions, formated_pairs, threads=threads) #TODO: here add the multi-sized extraction

    if params["pileup"]:
        pileups_random = {}
        if params['detrending'] == "patch":
            pileups_random_queue = matrix_extractor.launch_extraction(random_positions, random_pairs, randoms = True, threads=threads)
            pileups_random = empty_queue_in_dict(pileups_random_queue, keys = ["sep_id", "binning", "cool_name"]) # exporting the patch detrending as an dict for access

        ## Pileup detrending and display
        pileup_display_args["display_strand"] = pileup_display_args["display_strand"] and (0 not in positions["Strand"])
        pileup_display = Display(
            input_queue = pileups,
            output_queues = [],
            function = display_pileup,
            patch_detrending = pileups_random,
            outpath = outpath,
            title = data_title,
            size_metrics = params["expected_sizes"],
            is_contact = (pos_type == "bed2d" or params["loops"]),
            is_region = True,
            padding = params['padding'],
            **pileup_display_args
        )
        pileup_display.join()

# def compare(cool_pair, positions, outpath, params, log = None):
#     if not os.path.exists(outpath):
#         os.mkdir(outpath)

#     # parsing parameters
#     cool_pair_list = np.array(cool_pair)
#     cool_files = np.unique(cool_pair_list.reshape(-1, 1))
#     gff = params['gff']
#     windows = params['windows']
#     detrending = params['detrending']
#     nb_pos = params['nb_pos']
#     max_dist = params['random_max_dist']
#     loops = params['loops']
#     raw = params['raw']
#     min_dist = params['min_dist']
#     diagonal_mask = params['diagonal_mask']
#     trans_contact = params['trans_contact']
#     circular_chromosomes = params['circular_chromosomes']
#     display_strand = params['display_strand']
#     output_format = params['output_formats']
#     compute_pileup = params['pileup']
#     method = params['method']
#     flip = params['flip']
#     cmap = params['cmap_pileup']
#     cmap_color = params['cmap_color']
#     display_sense = params['display_sense']
#     center = params["center"]
#     separate_by = params["separate_by"]
#     separate_regions = params["separation_regions"]
#     overlap = params["overlap"]
#     contact_separation = params["contact_separation"]
#     contact_range = params["contact_range"]
#     ps_detrending = params["detrending"] == "ps"

#     # parsing gff file
#     pos_type, pos_file = positions
#     match pos_type:
#         case 'gff':
#             positions_parsed = parse_gff(pos_file) 
#         case 'bed':
#             # if gff file is provided, annotating each line of the bed file with its genes. The submatrices will be computed on such positions.
#             if gff != None:
#                 positions_parsed = parse_bed_annotated(pos_file, gff, overlap=overlap)
#                 # TODO write position file
#                 # write_positions(positions_parsed, outpath)
#             else:
#                 positions_parsed = parse_bed(pos_file, default_strand=1)
#         case 'global':
#             positions_parsed = None

#     data_title = pos_file.split('/')[-1].split('.')[0].replace('/', '_')
#     cools_meta = {}
#     global_matrices = {}

#     selected_positions = separate_positions(positions_parsed, data_title, separate_by=separate_by, separate_regions=separate_regions, overlap=overlap, outpath=outpath) if pos_type != 'global' else {'global':None}
#     random_locus = []
#     for position_name, positions in selected_positions.items():

#         # retrieving submatrices for each unique cool file in the cool pairs
#         for cool_path in cool_files:
#             cool = cooler.Cooler(cool_path)
#             bins = cool.binsize

#             matrix_outfolder = f"{outpath}/{cool.filename.split('/')[-1].split('.')[0]}"

#             cools_meta[cool_path] = (bins, cool.chromsizes, matrix_outfolder)

#             if pos_type == "global":
#                 global_matrices[cool_path] = cool.matrix(balance=(not raw))[:]

#             else:
                
#                 if detrending == "patch":
#                     random_locus = get_random_from_locus(cool, positions_parsed, nb_pos=nb_pos, max_dist=max_dist) if len(random_locus) == 0 else random_locus
            
#                 if not os.path.exists(matrix_outfolder):
#                     os.mkdir(matrix_outfolder)

#                 matrix_tmp = matrix_outfolder + "/tmp"
#                 if not os.path.exists(matrix_tmp):
#                     os.mkdir(matrix_tmp)


#                 for window in windows:
#                     # retrieving submatrices
#                     submatrices = compute_submatrices(cool, 
#                                                     position_name, 
#                                                     positions, 
#                                                     bins, 
#                                                     window, 
#                                                     loops=loops, 
#                                                     min_dist=min_dist,
#                                                     circular=circular_chromosomes, 
#                                                     trans_contact=trans_contact, 
#                                                     diagonal_mask=diagonal_mask, 
#                                                     center=center, 
#                                                     sort_contact=contact_separation, 
#                                                     contact_range = contact_range,
#                                                     ps_detrend = ps_detrending,
#                                                     raw=raw)

#                     if detrending == "patch":
#                         random_submatrices = compute_submatrices(cool, 
#                                                                 position_name, 
#                                                                 random_locus, 
#                                                                 bins, 
#                                                                 window, 
#                                                                 loops=loops, 
#                                                                 min_dist=min_dist,
#                                                                 circular=circular_chromosomes, 
#                                                                 trans_contact=False, 
#                                                                 center = center, 
#                                                                 sort_contact=contact_separation, 
#                                                                 contact_range = contact_range,
#                                                                 raw=raw)
                    

#                     # writting all submatrices
#                     for name in submatrices.keys():
#                         if len(submatrices[name]) == 0:
#                             continue
#                         for _, submatrix in submatrices[name].iterrows():
#                             matrix = submatrix.Matrix
#                             size = int(np.sqrt(len(matrix)))
#                             outname = f"{name}_{submatrix.Loc1}_{submatrix.Loc2}_{window}_{bins}" if loops else f"{name}_{submatrix.Loc1}_{window}_{bins}"
#                             pd.DataFrame(matrix.reshape((size, size))).to_csv(f"{matrix_tmp}/{outname}.csv")
                        
#                         # computing pileup
#                         if compute_pileup:
#                             ## aggregating the matrices
#                             pileup_matrices = get_windows(submatrices[name], positions, flip)
#                             match method:
#                                 case "median":
#                                     pileup = np.apply_along_axis(np.nanmedian, 0, pileup_matrices)
#                                 case "mean":
#                                     pileup = np.apply_along_axis(np.nanmean, 0, pileup_matrices)

#                             ## detrending
#                             if name[-len("trans"):] != "trans":
#                                 match detrending:
#                                     case "patch":
#                                         if nb_pos >= 1:
#                                             random_pileup_matrices = get_windows(random_submatrices[name], random_locus, flip)
#                                             match method:
#                                                 case "median":
#                                                     pileup_null = np.apply_along_axis(np.nanmedian, 0, random_pileup_matrices)
#                                                 case "mean":
#                                                     pileup_null = np.apply_along_axis(np.nanmean, 0, random_pileup_matrices)
#                                             pileup = pileup / pileup_null
                            
#                             size = int(np.sqrt(len(pileup)))
#                             pd.DataFrame(pileup.reshape((size, size))).to_csv(f"{matrix_tmp}/{name}_pileup_{window}_{bins}.csv")

#         # building pairs of submatrices and displaying
#         for cools in cool_pair_list:

#             # retrieving the outfolders
#             bins, chromsizes, matrix_outfolder1 = cools_meta[cools[0]]
#             _, _, matrix_outfolder2 = cools_meta[cools[1]]

#             if pos_type == "global":
#                 global_matrices[cool_path] = cool.matrix(balance=(not raw))[:]

#                 mat1 = global_matrices[cools[0]]
#                 mat2 = global_matrices[cools[1]]
#                 mat_name1 = matrix_outfolder1.split('/')[-1]
#                 mat_name2 = matrix_outfolder2.split('/')[-1]

#                 # display
#                 save_to = matrix_outfolder1 + "_vs_" + matrix_outfolder2.split('/')[-1]
#                 if not os.path.exists(save_to):
#                     os.mkdir(save_to)
#                 display_compare(mat1, mat2, mat_name1, mat_name2, bins, None, None, None, None, None, chromsizes=chromsizes, display_sense=display_sense, output_format=output_format, cmap=cmap, cmap_color=cmap_color, outfolder=save_to, is_global=True)

#             else:

#                 # listing the submatrices in each folder and keeping the commons
#                 submatrices1 = os.listdir(matrix_outfolder1 + "/tmp")
#                 submatrices2 = os.listdir(matrix_outfolder2 + "/tmp")

#                 submatrices_commons =  np.intersect1d(submatrices1, submatrices2)

#                 for sub in submatrices_commons:

#                     # retrieving matrices
#                     mat1 = np.array(pd.read_csv(matrix_outfolder1 + "/tmp/" + sub, index_col=0))
#                     mat2 = np.array(pd.read_csv(matrix_outfolder2 + "/tmp/" + sub, index_col=0))
#                     mat_name1 = matrix_outfolder1.split('/')[-1]
#                     mat_name2 = matrix_outfolder2.split('/')[-1]

#                     # retrieving matrices positions
#                     is_pileup = False
#                     if "pileup" in sub:
#                         if not compute_pileup:
#                             continue
#                         is_pileup = True
#                         window = int(sub[len(position_name) + 1:].split('_')[1])
#                         bins = int(sub[len(position_name) + 1:].split('_')[2].split('.')[0])
#                     else:
#                         if loops:
#                             pos1 = int(sub[len(position_name) + 1:].split('_')[0])
#                             pos2 = int(sub[len(position_name) + 1:].split('_')[1])
#                             window = int(sub[len(position_name) + 1:].split('_')[2])
#                             bins = int(sub[len(position_name) + 1:].split('_')[3].split('.')[0])
#                         else:
#                             pos1 = int(sub[len(position_name) + 1:].split('_')[0])
#                             pos2 = None
#                             window = int(sub[len(position_name) + 1:].split('_')[1])
#                             bins = int(sub[len(position_name) + 1:].split('_')[2].split('.')[0])

#                     # display
#                     save_to = matrix_outfolder1 + "_vs_" + matrix_outfolder2.split('/')[-1]
#                     if not os.path.exists(save_to):
#                         os.mkdir(save_to)
#                     display_compare(mat1, mat2, mat_name1, mat_name2, bins, window, pos1, pos2, positions, position_name, chromsizes=chromsizes, display_sense=display_sense, display_strand=display_strand, circular=circular_chromosomes, output_format=output_format, cmap=cmap, cmap_color=cmap_color, is_pileup=is_pileup, is_contact=loops, outfolder=save_to)
#                     if log:
#                         log.write(f"Displaying {mat_name1} vs {mat_name2} ({bins//1000}kb binning, {window//1000}kb window)\n")

#     # removing tmp files
#     for cool_path in cool_files:
#         _, _, matrix_outfolder = cools_meta[cool_path]
#         if os.path.exists(matrix_outfolder):
#             rmtree(matrix_outfolder + "/tmp")
#             os.rmdir(matrix_outfolder)