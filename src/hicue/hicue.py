from .utils import *
from .displays import *
from .parser import *

def extract(cool_files, positions, outpath, params):

    if not os.path.exists(outpath):
        os.mkdir(outpath)

    # parsing parameters
    gff = params['gff']
    windows = params['windows']
    detrending = params['detrending']
    nb_pos = params['nb_pos']
    max_dist = params['random_max_dist']
    loops = params['loops']
    min_dist = params['min_dist']
    diagonal_mask = params['diagonal_mask']
    trans_contact = params['trans_contact']
    circular_chromosomes = params['circular_chromosomes']
    display_strand = params['display_strand']
    output_format = params['output_formats']
    compute_pileup = params['pileup']
    plot_loci = params['loci']
    method = params['method']
    flip = params['flip']
    cmap = params['cmap_pileup']
    display_sense = params['display_sense']
    center = params["center"]
    separate_by = params["separate_by"]
    separate_regions = params["separation_regions"]
    overlap = params["overlap"]
    contact_separation = params["contact_separation"]
    contact_range = params["contact_range"]
    ps_detrending = params["detrending"] == "ps"

    # parsing gff file
    pos_type, pos_file = positions
    match pos_type:
        case 'gff':
            positions_parsed = parse_gff(pos_file) 
        case 'bed':
            # if gff file is provided, annotating each line of the bed file with its genes. The submatrices will be computed on such positions.
            if gff != None:
                positions_parsed = parse_bed_annotated(pos_file, gff, overlap=overlap)
                # TODO write position file
                # write_positions(positions_parsed, outpath)
            else:
                positions_parsed = parse_bed(pos_file, default_strand=1)

    data_title = pos_file.split('/')[-1].split('.')[0].replace('/', '_')

    selected_positions = separate_positions(positions_parsed, data_title, separate_by=separate_by, separate_regions=separate_regions, overlap=overlap, outpath=outpath)
    random_locus = []
    for position_name, positions in selected_positions.items():
        for cool_path in cool_files:
            cool = cooler.Cooler(cool_path)
            bins = cool.binsize
            if detrending == "patch":
                random_locus = get_random_from_locus(cool, positions_parsed, nb_pos=nb_pos, max_dist=max_dist) if len(random_locus) == 0 else random_locus

            matrix_outfolder = f"{outpath}/{cool.filename.split('/')[-1].split('.')[0]}"
            if not os.path.exists(matrix_outfolder):
                os.mkdir(matrix_outfolder)


            for window in windows:
                outfolder = f"{matrix_outfolder}/window_{window//1000}kb"
                if not os.path.exists(outfolder):
                    os.mkdir(outfolder)

                # retrieving submatrices
                submatrices = compute_submatrices(cool, 
                                                  position_name, 
                                                  positions, 
                                                  bins, 
                                                  window, 
                                                  loops=loops, 
                                                  min_dist=min_dist,
                                                  circular=circular_chromosomes, 
                                                  trans_contact=trans_contact, 
                                                  diagonal_mask=diagonal_mask, 
                                                  center=center, 
                                                  sort_contact=contact_separation, 
                                                  contact_range = contact_range,
                                                  ps_detrend = ps_detrending)

                if detrending == "patch":
                    random_submatrices = compute_submatrices(cool, 
                                                            position_name, 
                                                            random_locus, 
                                                            bins, 
                                                            window, 
                                                            loops=loops, 
                                                            min_dist=min_dist,
                                                            circular=circular_chromosomes, 
                                                            trans_contact=False, 
                                                            center = center, 
                                                            sort_contact=contact_separation, 
                                                            contact_range = contact_range)
                 
                for name in submatrices.keys():
                    if len(submatrices[name]) == 0:
                        continue
                    # displaying submatrices
                    if plot_loci:
                        display_submatrices(submatrices[name], positions, window, outfolder=outfolder + f"/{name}", circular=circular_chromosomes, chromsizes = cool.chromsizes, output_format=output_format, display_strand=display_strand, display_sense=display_sense)
                    display_all_submatrices(submatrices[name], positions, window, outfolder=outfolder + f"/{name}", circular=circular_chromosomes, chromsizes = cool.chromsizes, output_format=output_format, display_strand=display_strand, display_sense=display_sense)
                    
                    # computing pileup
                    if compute_pileup:
                        ## aggregating the matrices
                        pileup_matrices = get_windows(submatrices[name], positions, flip)
                        match method:
                            case "median":
                                pileup = np.apply_along_axis(np.nanmedian, 0, pileup_matrices)
                            case "mean":
                                pileup = np.apply_along_axis(np.nanmean, 0, pileup_matrices)

                        ## detrending
                        if name[-len("trans"):] != "trans":
                            match detrending:
                                case "patch":
                                    if nb_pos >= 1:
                                        random_pileup_matrices = get_windows(random_submatrices[name], random_locus, flip)
                                        match method:
                                            case "median":
                                                pileup_null = np.apply_along_axis(np.nanmedian, 0, random_pileup_matrices)
                                            case "mean":
                                                pileup_null = np.apply_along_axis(np.nanmean, 0, random_pileup_matrices)
                                        pileup = pileup / pileup_null
                        
                        title = f"{name.replace('_', ' ')} pileup ({len(pileup_matrices)} matrices)"
                        pileup_outpath = f"{outfolder}/{name}_pileup" if len(outfolder) > 0 else ""
                        # TODO add display strand True for separated sense of transcription
                        display_strand = False
                        if len(pileup_outpath) > 0:
                            if not os.path.exists(f"{outfolder}/matrices_tables"):
                                os.mkdir(f"{outfolder}/matrices_tables")
                            pd.DataFrame(pileup).to_csv(f"{outfolder}/matrices_tables/{name}_pileup.csv")
                        display_pileup(pileup, window, cmap=cmap, title=title, outpath=pileup_outpath, output_format=output_format, display_strand=display_strand, display_sense=display_sense)
                

def tracks(cool_files, tracks, outpath, params):
    if not os.path.exists(outpath):
        os.mkdir(outpath)
    # parsing parameters
    threshold = params['threshold']
    percentage = params['percentage']
    gff = params['gff']
    overlap = params['overlap']
    binning = params['binning']
    windows = params['windows']
    detrending = params['detrending']
    nb_pos = params['nb_pos']
    max_dist = params['random_max_dist']
    loops = params['loops']
    min_dist = params['min_dist']
    diagonal_mask = params['diagonal_mask']
    trans_contact = params['trans_contact']
    circular_chromosomes = params['circular_chromosomes']
    display_strand = params['display_strand']
    output_format = params['output_formats']
    compute_pileup = params['pileup']
    plot_loci = params['loci']
    method = params['method']
    flip = params['flip']
    cmap = params['cmap_pileup']
    display_sense = params['display_sense']
    center = params["center"]
    contact_separation = params["contact_separation"]
    contact_range = params["contact_range"]
    ps_detrending = params["detrending"] == "ps"

    # parsing tracks
    tracks_parsed, bw_tracks = parse_tracks(tracks, threshold = threshold, percentage = percentage, gff=gff, binning=binning)
    data_title = tracks.split('/')[-1].split('.')[0].replace('/', '_')

    selected_positions = separate_positions(tracks_parsed, data_title, outpath=outpath)
    random_locus = []
    for position_name, positions in selected_positions.items():
        for cool_path in cool_files:
            cool = cooler.Cooler(cool_path)
            bins = cool.binsize
            if detrending == "patch":
                random_locus = get_random_from_locus(cool, tracks_parsed, nb_pos=nb_pos, max_dist=max_dist) if len(random_locus) == 0 else random_locus

            matrix_outfolder = f"{outpath}/{cool.filename.split('/')[-1].split('.')[0]}"
            if not os.path.exists(matrix_outfolder):
                os.mkdir(matrix_outfolder)


            for window in windows:
                outfolder = f"{matrix_outfolder}/window_{window//1000}kb"
                if not os.path.exists(outfolder):
                    os.mkdir(outfolder)

                # retrieving submatrices
                submatrices = compute_submatrices(cool, 
                                                  position_name, 
                                                  positions, 
                                                  bins, 
                                                  window, 
                                                  loops=loops, 
                                                  min_dist=min_dist,
                                                  circular=circular_chromosomes, 
                                                  trans_contact=trans_contact, 
                                                  diagonal_mask=diagonal_mask, 
                                                  center=center, 
                                                  sort_contact=contact_separation, 
                                                  contact_range = contact_range,
                                                  ps_detrend = ps_detrending)

                if detrending == "patch":
                    random_submatrices = compute_submatrices(cool, 
                                                            position_name, 
                                                            random_locus, 
                                                            bins, 
                                                            window, 
                                                            loops=loops, 
                                                            min_dist=min_dist,
                                                            circular=circular_chromosomes, 
                                                            trans_contact=False, 
                                                            center = center, 
                                                            sort_contact=contact_separation, 
                                                            contact_range = contact_range)
                 
                for name in submatrices.keys():
                    if len(submatrices[name]) == 0:
                        continue
                    # displaying submatrices
                    if plot_loci:
                        display_submatrices(submatrices[name], positions, window, outfolder=outfolder + f"/{name}", circular=circular_chromosomes, chromsizes = cool.chromsizes, output_format=output_format, display_strand=display_strand, display_sense=display_sense)
                    display_all_submatrices(submatrices[name], positions, window, outfolder=outfolder + f"/{name}", circular=circular_chromosomes, chromsizes = cool.chromsizes, output_format=output_format, display_strand=display_strand, display_sense=display_sense)
                    
                    # computing pileup
                    if compute_pileup:
                        ## aggregating the matrices
                        pileup_matrices = get_windows(submatrices[name], positions, flip)
                        match method:
                            case "median":
                                pileup = np.apply_along_axis(np.nanmedian, 0, pileup_matrices)
                            case "mean":
                                pileup = np.apply_along_axis(np.nanmean, 0, pileup_matrices)

                        ## detrending
                        if name[-len("trans"):] != "trans":
                            match detrending:
                                case "patch":
                                    if nb_pos >= 1:
                                        random_pileup_matrices = get_windows(random_submatrices[name], random_locus, flip)
                                        match method:
                                            case "median":
                                                pileup_null = np.apply_along_axis(np.nanmedian, 0, random_pileup_matrices)
                                            case "mean":
                                                pileup_null = np.apply_along_axis(np.nanmean, 0, random_pileup_matrices)
                                        pileup = pileup / pileup_null
                        
                        title = f"{name.replace('_', ' ')} pileup ({len(pileup_matrices)} matrices)"
                        pileup_outpath = f"{outfolder}/{name}_pileup" if len(outfolder) > 0 else ""
                        # TODO add display strand True for separated sense of transcription
                        display_strand = False
                        if len(pileup_outpath) > 0:
                            if not os.path.exists(f"{outfolder}/matrices_tables"):
                                os.mkdir(f"{outfolder}/matrices_tables")
                            pd.DataFrame(pileup).to_csv(f"{outfolder}/matrices_tables/{name}_pileup.csv")
                        display_pileup(pileup, window, cmap=cmap, title=title, outpath=pileup_outpath, output_format=output_format, display_strand=display_strand, display_sense=display_sense)

def annotate(gff, bed):
    pass