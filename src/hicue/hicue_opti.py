from .utils_opti import *
from .classes.Reader import *
from .classes.PairFormater import *

def extract_opti(cool_files, positions, outpath, log = None, **params):

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    # parsing parameters

    ## used
    # loops = params['loops']
    # min_dist = params['min_dist']
    # separate_by = params["separate_by"]
    # separate_regions = params["separation_regions"]
    # overlap = params["overlap"]
    # contact_separation = params["contact_separation"]
    # contact_range = params["contact_range"]
    # has_trans = params['trans']
    # center = params["center"]

    ## passed
    # ps_detrending = params["detrending"] == "ps"

    # circular_chromosomes = params['circular_chromosomes']


    # windows = params['windows']
    # detrending = params['detrending']
    # nb_pos = params['nb_pos']
    # max_dist = params['random_max_dist']
    # raw = params['raw']
    # diagonal_mask = params['diagonal_mask']
    # display_strand = params['display_strand']
    # output_format = params['output_formats']
    # compute_pileup = params['pileup']
    # plot_loci = params['loci']
    # method = params['method']
    # flip = params['flip']
    # cmap = params['cmap_pileup']
    # cmap_color = params['cmap_color']
    # display_sense = params['display_sense']

    # checking multiprocessing values
    threads = max(1, params["threads"])
    cpus = max(1, min(params["cpu"], cpu_count() - 1)) # allocates the required CPUs without overflowing the machine cores

    # reading parameters
    pos_type, pos_file = positions
    data_title = pos_file.split('/')[-1].split('.')[0].replace('/', '_')

    ## Reading file and annotating
    annotation = {}
    if params['gff']:
        annotation['gff'] = params['gff']
    if params['tracks']:
        annotation['tracks'] = params['tracks']
    
    reader = reader = Reader(pos_file, pos_type, annotation_files = annotation, save_to = "", loop = params['loops'])# TODO: option on verbose annotation
    positions, pairing_queue = reader.read_file(threads = threads)

    ## Formating indexes pairs
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
    formater = PairFormater(positions, **format_params)
    formated_pairs = formater.format_pairs(pairing_queue, threads = threads)

    print(formated_pairs.head())
    


    # random_locus = []
    # for position_name, positions in selected_positions.items():
    #     for cool_path in cool_files:
    #         cool = cooler.Cooler(cool_path)
    #         bins = cool.binsize
    #         if detrending == "patch":
    #             random_locus = get_random_from_locus(cool, positions_parsed, nb_pos=nb_pos, max_dist=max_dist) if len(random_locus) == 0 else random_locus

    #         matrix_outfolder = f"{outpath}/{cool.filename.split('/')[-1].split('.')[0]}"
    #         if not os.path.exists(matrix_outfolder):
    #             os.mkdir(matrix_outfolder)


    #         for window in windows:
    #             outfolder = f"{matrix_outfolder}/window_{window//1000}kb"
    #             if not os.path.exists(outfolder):
    #                 os.mkdir(outfolder)

    #             # retrieving submatrices
    #             submatrices = compute_submatrices_opti(cool, 
    #                                               position_name, 
    #                                               positions, 
    #                                               bins, 
    #                                               window, 
    #                                               loops=loops, 
    #                                               min_dist=min_dist,
    #                                               circular=circular_chromosomes, 
    #                                               trans_contact=trans_contact, 
    #                                               diagonal_mask=diagonal_mask, 
    #                                               center=center, 
    #                                               sort_contact=contact_separation, 
    #                                               contact_range = contact_range,
    #                                               ps_detrend = ps_detrending,
    #                                               raw=raw)

    #             if detrending == "patch":
    #                 random_submatrices = compute_submatrices_opti(cool, 
    #                                                         position_name, 
    #                                                         random_locus, 
    #                                                         bins, 
    #                                                         window, 
    #                                                         loops=loops, 
    #                                                         min_dist=min_dist,
    #                                                         circular=circular_chromosomes, 
    #                                                         trans_contact=False, 
    #                                                         center = center, 
    #                                                         sort_contact=contact_separation, 
    #                                                         contact_range = contact_range,
    #                                                         raw=raw)
                 
    #             for name in submatrices.keys():
    #                 if len(submatrices[name]) == 0:
    #                     continue
    #                 # displaying submatrices
    #                 if plot_loci:
    #                     display_submatrices(submatrices[name], positions, window, outfolder=outfolder + f"/{name}", circular=circular_chromosomes, chromsizes = cool.chromsizes, output_format=output_format, display_strand=display_strand, display_sense=display_sense, binning = bins)
    #                 display_all_submatrices(submatrices[name], positions, window, outfolder=outfolder + f"/{name}", circular=circular_chromosomes, chromsizes = cool.chromsizes, output_format=output_format, display_strand=display_strand, display_sense=display_sense, binning = bins)
                    
    #                 # computing pileup
    #                 if compute_pileup:
    #                     ## aggregating the matrices
    #                     pileup_matrices = get_windows(submatrices[name], positions, flip)
    #                     match method:
    #                         case "median":
    #                             pileup = np.apply_along_axis(np.nanmedian, 0, pileup_matrices)
    #                         case "mean":
    #                             pileup = np.apply_along_axis(np.nanmean, 0, pileup_matrices)

    #                     ## detrending
    #                     if name[-len("trans"):] != "trans":
    #                         match detrending:
    #                             case "patch":
    #                                 if nb_pos >= 1:
    #                                     random_pileup_matrices = get_windows(random_submatrices[name], random_locus, flip)
    #                                     match method:
    #                                         case "median":
    #                                             pileup_null = np.apply_along_axis(np.nanmedian, 0, random_pileup_matrices)
    #                                         case "mean":
    #                                             pileup_null = np.apply_along_axis(np.nanmean, 0, random_pileup_matrices)
    #                                     pileup = pileup / pileup_null
                        
    #                     title = f"{name.replace('_', ' ')} pileup ({len(pileup_matrices)} matrices)"
    #                     pileup_outpath = f"{outfolder}/{name}_pileup" if len(outfolder) > 0 else ""
    #                     if len(pileup_outpath) > 0:
    #                         if not os.path.exists(f"{outfolder}/matrices_tables"):
    #                             os.mkdir(f"{outfolder}/matrices_tables")
    #                         size = int(np.sqrt(len(pileup)))
    #                         pd.DataFrame(pileup.reshape((size, size))).to_csv(f"{outfolder}/matrices_tables/{name}_pileup.csv")
    #                     display_pileup(pileup, window, cmap=cmap, cmap_color=cmap_color, title=title, outpath=pileup_outpath, output_format=output_format, display_strand=flip, display_sense=display_sense, binning = bins)
