from .utils_opti import *
from .classes.Reader import *
from .classes.PairFormater import *
from .classes.RandomSelector import *
from .classes.MatrixExtractorLauncher import *
from .classes.AsyncDisplays import *

def extract(cool_files, positions, outpath, log = None, **params):

    if not os.path.exists(outpath):
        os.mkdir(outpath)

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
        "display_strand" : params["display_strand"], 
        "display_sense" : params["display_sense"],
        "cmap": params["indiv_cmap_limits"],
        "color": params["indiv_cmap_color"]
    }

    random_params = {
        "center" : params["center"],
        "selection_window": params['rand_max_dist'],
        "nb_rand_per_pos" : params['nb_pos']
    }

    pileup_display_args = {
        "output_format": params["format"],
        "display_strand" : params["display_strand"], 
        "display_sense" : params["display_sense"],
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
    
    reader = reader = Reader(pos_file, pos_type, annotation_files = annotation, save_to = "", loop = params['loops'])# TODO: option on verbose annotation
    positions, pairing_queue = reader.read_file(threads = threads)

    ## Formating indexes pairs
    formater = PairFormater(positions, **format_params)
    formated_pairs = formater.format_pairs(pairing_queue, threads = threads)

    if params["save_tmp"]:
        positions.to_csv(f"{outpath}/{data_title}_positions.csv")
        formated_pairs.to_csv(f"{outpath}/{data_title}_formated_pairs.csv")

    ## Random locus selection (for patch only) from positions
    random_selection = None
    if params['detrending'] == "patch" and params["pileup"]:
        selector = RandomSelector(**random_params)
        random_selection = selector.select_randoms(positions, threads = threads)

        if params["save_tmp"]:
            random_selection.to_csv(f"{outpath}/{data_title}_random_patch.csv") # TODO: add method to re-use randoms for reproducibility
    
    ## Matrix extraction
    matrix_extractor = MatrixExtractorLauncher(cool_files,
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

    # works until then
    if params["pileup"]:
        pileups_random = {}
        if params['detrending'] == "patch":
            pileups_random_queue = matrix_extractor.launch_extraction(random_selection, formated_pairs, randoms = True, threads=threads)
            pileups_random = empty_queue_in_dict(pileups_random_queue, keys = ["sep_id", "binning"]) # exporting the patch detrending as an dict for access

        ## Pileup detrending and display
        # seems to crash after here
        pileup_display = Display(
            input_queue = pileups,
            output_queues = [],
            function = display_pileup,
            patch_detrending = pileups_random,
            outpath = outpath,
            title = data_title,
            windows = params["windows"],
            is_contact = (pos_type == "bed2d" or params["loops"]),
            **pileup_display_args
        )
        pileup_display.join()
        