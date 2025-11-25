from hicue.utils import *

def adjust_extents(ax, chrom1, chrom2, is_chrom1_circ, is_chrom2_circ, chromsizes={}):
    min, max = np.min(ax.get_xlim()), np.max(ax.get_xlim())
    extent_x = [item.get_text() for item in ax.get_xticklabels() if item._x >= min and item._x < max]
    x_ticks = [tick for tick in ax.get_xticks() if tick >= min and tick < max]
    for i in range(len(extent_x)):
        if extent_x[i][0] == '−':
            extent_x[i] = str(chromsizes[chrom2]//1000 - int(extent_x[i][1:])) if is_chrom2_circ else " "
        elif int(extent_x[i]) > chromsizes[chrom2]//1000:
            extent_x[i] = str(int(extent_x[i]) - chromsizes[chrom2]//1000) if is_chrom2_circ else " "
    ax.set(xticks=x_ticks, xticklabels=extent_x)

    min, max = np.min(ax.get_ylim()), np.max(ax.get_ylim())
    extent_y = [item.get_text() for item in ax.get_yticklabels() if item._y > min and item._y <= max]
    y_ticks = [tick for tick in ax.get_yticks() if tick > min and tick <= max]
    for i in range(len(extent_y)):
        if extent_y[i][0] == '−':
            extent_y[i] = str(chromsizes[chrom1]//1000 - int(extent_y[i][1:])) if is_chrom1_circ else " "
        elif int(extent_y[i]) > chromsizes[chrom1]//1000:
            extent_y[i] = str(int(extent_y[i]) - chromsizes[chrom1]//1000) if is_chrom1_circ else " "
    ax.set(yticks=y_ticks, yticklabels=extent_y)

def plot_strands(ax, locus1, locus2, window, is_contact = False, display_sense = "forward", flip = False, strand_level = 1.2, adjustment = 0.9):
    """Plots the direction of the positions strand around the matrix axis."""
    pos1 = min(locus1["Start"], locus1["End"]) if locus1["Strand"] == 1 else max(locus1["Start"], locus1["End"])
    pos2 = min(locus2["Start"], locus2["End"]) if locus2["Strand"] == 1 else max(locus2["Start"], locus2["End"])

    flip1, flip2 = locus1["Strand"] == -1 and not flip, locus2["Strand"] == -1 and not flip
    x_right = [ARROW_RIGHT, "left"]
    x_left = [ARROW_LEFT, "right"]
    y_up = [ARROW_UP, "bottom"]
    y_down = [ARROW_DOWN, "top"]

    if not is_contact:
        if locus1["Strand"] != 0:
            arrow = x_right
            match display_sense:
                case "forward":
                    arrow = x_left if locus1["Strand"] == -1 and not flip else x_right
                    to_add = - window if flip and locus1["Strand"] == -1 else window
                case "reverse":
                    arrow = x_right if locus1["Strand"] == -1 and not flip else x_left
                    to_add = window if flip and locus1["Strand"] == -1 else - window
            ax.text(pos2//1000, (pos1 + to_add * strand_level)//1000, arrow[0], horizontalalignment=arrow[1], fontsize=20)

    else:
        if locus2["Strand"] != 0:
            arrow = x_right
            match display_sense:
                case "forward":
                    arrow = x_left if locus2["Strand"] == -1 and not flip else x_right
                    to_add = - window if flip and locus1["Strand"] == -1 else window
                case "reverse":
                    arrow = x_right if locus2["Strand"] == -1 and not flip else x_left
                    to_add = window if flip and locus1["Strand"] == -1 else -window
            ax.text(pos2//1000, (pos1 + to_add * strand_level)//1000, arrow[0], horizontalalignment=arrow[1], fontsize=20)

        if locus1["Strand"] != 0:
            arrow = y_down
            match display_sense:
                case "forward":
                    arrow = y_up if locus1["Strand"] == -1 and not flip else y_down
                    to_add = - window if flip and locus2["Strand"] == -1 else window
                case "reverse":
                    arrow = y_down if locus1["Strand"] == -1 and not flip else y_up
                    to_add = window if flip and locus2["Strand"] == -1 else - window
            ax.text((pos2 + to_add * (strand_level * adjustment) )//1000, pos1//1000 , arrow[0], verticalalignment=arrow[1], fontsize=20)

def plot_map(ax, matrix, loc1, loc2, window, locus1, locus2, is_chrom1_circ, is_chrom2_circ, title="", display_sense="forward", chromsizes={}, display_strand=False, flipped = False, cmap=None, color="afmhot_r", adjust=True, show_title=True, log=True, strand_level = 1.2, adjustment = 0.9):
    """Plots a single matrix on the provided axis"""
    is_contact = loc1 != loc2
    pos1 = min(locus1["Start"], locus1["End"]) if locus1["Strand"] == 1 else max(locus1["Start"], locus1["End"])
    pos2 = min(locus2["Start"], locus2["End"]) if locus2["Strand"] == 1 else max(locus2["Start"], locus2["End"])

    name = f"{locus1['Name'].replace('/', '_')}" if not is_contact else f"{locus1['Name'].replace('/', '_')}-{locus2['Name'].replace('/', '_')}"
    title = name if len(title) == 0 else title
    if show_title:
        ax.set_title(title)

    with np.errstate(divide='ignore', invalid='ignore'): # cancelling the divide by 0 warning, as those value are replaced by np.nan and won't affect the final result
        display_matrix = np.log10(matrix) if log else matrix
    vmin = cmap[0] if cmap != None else None
    vmax = cmap[1] if cmap != None else None

    match display_sense:
        case "forward":
            y_extent = [(pos1 + window)//1000, (pos1 - window)//1000] if locus1["Strand"] != -1 or not flipped else [(pos1 - window)//1000, (pos1 + window)//1000]
            x_extent = [(pos2 - window)//1000, (pos2 + window)//1000] if locus2["Strand"] != -1 or not flipped else [(pos2 + window)//1000, (pos2 - window)//1000]
            mat = ax.imshow(display_matrix, extent=x_extent + y_extent, cmap=color, vmin=vmin, vmax=vmax)
        case "reverse":
            y_extent = [(pos1 - window)//1000, (pos1 + window)//1000] if locus1["Strand"] != -1 or not flipped else [(pos1 + window)//1000, (pos1 - window)//1000]
            x_extent = [(pos2 + window)//1000, (pos2 - window)//1000] if locus2["Strand"] != -1 or not flipped else [(pos2 - window)//1000, (pos2 + window)//1000]
            mat = ax.imshow(np.flip(display_matrix), extent=x_extent + y_extent, cmap=color, vmin=vmin, vmax=vmax)
    
    chrom1 = locus1["Chromosome"]
    chrom2 = locus2["Chromosome"]
    if adjust:
        adjust_extents(ax, chrom1, chrom2, is_chrom1_circ, is_chrom2_circ, chromsizes)

    if display_strand:
        plot_strands(ax, locus1, locus2, window, is_contact=pos1!=pos2, display_sense=display_sense, flip = flipped, strand_level = strand_level, adjustment = adjustment)

    return mat

def plot_tracks(tracks, ax_tracks, start, stop, axis = "horizontal", xlabel = "", ylabel = "", flip = False):
    """Plots tracks on the track axis aligned between start and stop, on the horizontal or vertical axis."""
    index = np.array([i/len(tracks) * (stop - start) + start for i in range(len(tracks))])

    match axis:
        case "horizontal":
            index = np.flip(index) if flip else index
            ax_tracks.plot(index, tracks)
            if len(ylabel) > 0:
                ax_tracks.set_ylabel(ylabel)
            if len(xlabel) > 0:
                ax_tracks.set_xlabel(xlabel)
        case "vertical":
            index = np.flip(index) if not flip else index
            ax_tracks.plot(tracks, index)
            ax_tracks.yaxis.tick_right()
            if len(ylabel) > 0:
                ax_tracks.set_xlabel(ylabel)
            if len(xlabel) > 0:
                ax_tracks.set_ylabel(xlabel)

async def display_submatrix(matrix, pair, window, binning = 1000, outfolder="", positions = None, output_format=['pdf'], chromsizes = {}, display_strand=False, flipped = False, display_sense="forward", cmap = None, color = "afmhot_r", display_tracks = False, track_unit = ""):
        """Displays a single submatrix for a pair of positions."""

        # checking outpath
        individual_outfolder = f"{outfolder}/individual_{window//1000}kb_window"
        create_folder_path(individual_outfolder)

        i, j = pair["Locus1"], pair["Locus2"]
        pos1, pos2 = positions.loc[i], positions.loc[j]
        is_contact = i!=j
        title = f"Window centered on\n{pos1['Name']} vs {pos2['Name']}" if is_contact else f"Window centered on\n{positions.loc[i]['Name']}"
        y_label = f"{pos1['Name']}" if is_contact else ""
        x_label = f"{pos2['Name']}" if is_contact else ""
        outpath = "" if len(individual_outfolder) == 0 else f"{individual_outfolder}/{pos1['Name'].replace('/', '_')}" if not is_contact else f"{individual_outfolder}/{pos1['Name'].replace('/', '_')}-{pos2['Name'].replace('/', '_')}"

        # checking tracks
        submatrix = matrix[:len(matrix[0])]
        subtracks = matrix[len(matrix[0]):] if len(matrix) > len(matrix[0]) else []
        
        if len(subtracks) == 0 or not display_tracks:

            plt.figure(figsize=(6,6))

            mat = plot_map(plt.gca(), submatrix, i, j, window, pos1, pos2, pair["Chrom1_circular"], pair["Chrom2_circular"], title=title, chromsizes = chromsizes, display_sense=display_sense, display_strand=display_strand, flipped = flipped, strand_level=1.2, cmap=cmap, color=color)

            plt.colorbar(mat, fraction=0.01)
            plot_xlabel =  "\n" + f"{pos2['Chromosome']} Genomic coordinates (in kb)"
            plot_xlabel += "\n" + x_label if len(x_label) > 0 else ""
            plt.xlabel(plot_xlabel)
            plot_ylabel =  f"{pos1['Chromosome']} Genomic coordinates (in kb)"
            plot_ylabel = y_label + "\n" + plot_ylabel if len(y_label) > 0 else plot_ylabel
            plt.ylabel(plot_ylabel)
        
        else:

            width = 3 if is_contact else 2
            wspace = 0.4 if is_contact else 0.1
            hspace = 0.7 if is_contact else 0.5
            ratios = [6, 1, 0.1] if is_contact else [6, 0.1]
            figwidth = 8 if is_contact else 6
            figheight = 7 if is_contact else 8
            plt.figure(figsize=(figwidth,figheight))
            gs = grid.GridSpec(5, width, height_ratios = [1,1,1,1,1], width_ratios = ratios, wspace=wspace, hspace=hspace) 

            # matrix ax
            ax = plt.subplot(gs[:4, 0])
            mat = plot_map(ax, submatrix, i, j, window, pos1, pos2, pair["Chrom1_circular"], pair["Chrom2_circular"], title=title, show_title=(not is_contact), chromsizes = chromsizes, display_sense=display_sense, display_strand=display_strand, flipped=flipped, strand_level=1.2, adjust=False, cmap=cmap, color=color)

            # colorbar ax
            ax_cb = plt.subplot(gs[1, width - 1])
            plt.colorbar(mat, fraction=0.01, cax=ax_cb)
            
            flip1, flip2 = pos1["Strand"] == -1, pos2["Strand"] == -1

            track_labelling = track_unit
            ax_tracks = plt.subplot(gs[4, 0], sharex=ax)
            plot_xlabel =  "\n" + f"{pos2['Chromosome']} Genomic coordinates (in kb)"
            plot_xlabel += "\n" + x_label if len(x_label) > 0 else ""

            if not is_contact:
                xstart, xstop = ax.get_xlim()
                plot_tracks(subtracks[0], ax_tracks, xstart, xstop, axis = "horizontal", xlabel = plot_xlabel, ylabel = track_labelling, flip = (display_sense == "reverse"))

            else:
                # first track
                ax_tracks1 = plt.subplot(gs[:4, 1],sharey = ax)
                ystart, ystop = ax.get_ylim()
                plot_tracks(subtracks[0], ax_tracks1, ystart, ystop, axis = "vertical", ylabel = track_labelling, flip = (display_sense == "reverse"))

                # second track
                if len(subtracks) > 1:
                    xstart, xstop = ax.get_xlim()
                    plot_tracks(subtracks[1], ax_tracks, xstart, xstop, axis = "horizontal", xlabel = plot_xlabel, ylabel = track_labelling, flip = (display_sense == "reverse"))

            adjust_extents(ax, pos1["Chromosome"], pos2["Chromosome"], pair["Chrom1_circular"], pair["Chrom2_circular"], chromsizes)

            plot_ylabel =  f"{pos1['Chromosome']} Genomic coordinates (in kb)"
            plot_ylabel = y_label + "\n" + plot_ylabel if len(y_label) > 0 else plot_ylabel
            ax.set_ylabel(plot_ylabel)

            if is_contact:
                plt.suptitle(title)

        if len(outpath) > 0 :
            for format in output_format:
                plt.savefig(outpath + f".{format}", bbox_inches="tight")
        plt.close()

async def display_batch_submatrices(submatrices, positions, window, title = "", batch_size = 64, outfolder="", output_format=['pdf'], chromsizes = {}, display_strand=False, flipped = False, display_sense="forward", display_tracks = False, track_unit="", cmap=None, color="afmhot_r"):
    """Displays a batch of submatrices in a single figure. References each submatrix in a csv file."""
    # checking outpath
    batched_outfolder = f"{outfolder}/batched_{window//1000}kb_window"
    create_folder_path(batched_outfolder)

    cols = math.ceil(math.sqrt(batch_size))
    rows = math.ceil(batch_size / cols)
    nb_matrices = len(submatrices)
    
    plt.figure(figsize=(20,20) if not display_tracks else (23, 23))
    gs = grid.GridSpec(rows, cols, height_ratios=[1]*rows, width_ratios=[1]*cols, wspace=0.4 if display_tracks else 0.2, hspace=0.5)

    index = []
    n = 1
    for k in range(nb_matrices):
        matrix = submatrices[k]["matrix"]

        # checking tracks
        submatrix = matrix[:len(matrix[0])]
        subtracks = matrix[len(matrix[0]):] if len(matrix) > len(matrix[0]) else []

        i = k // rows
        j = k % rows

        loc1, loc2 = submatrices[k]["pair"]["Locus1"], submatrices[k]["pair"]["Locus2"]
        is_chrom1_circ, is_chrom2_circ = submatrices[k]["pair"]["Chrom1_circular"], submatrices[k]["pair"]["Chrom2_circular"]
        pos1, pos2 = positions.loc[loc1], positions.loc[loc2]
        is_contact = loc1 != loc2

        name = f"{pos1['Name'].replace('/', '_')}" if not is_contact else f"{pos1['Name'].replace('/', '_')}-{pos2['Name'].replace('/', '_')}"
        map_title = f"Position {n}" if not is_contact else f"Contact {n}"
        n += 1
        index.append([map_title, name])
        
        if not display_tracks:
            ax = plt.subplot(gs[i, j])
            plot_map(ax, submatrix, loc1, loc2, window, pos1, pos2, is_chrom1_circ, is_chrom2_circ, title=map_title, display_sense=display_sense, chromsizes = chromsizes, display_strand=display_strand, flipped = flipped, strand_level=1.6, adjustment= 0.6, cmap=cmap, color=color)

            if k >= nb_matrices - cols :
                ax.set_xlabel('\nGenomic\ncoordinates in kb')
            if j == 0:
                ax.set_ylabel('Genomic\ncoordinates in kb\n')
                
        else:
            subgs = gs[i, j].subgridspec(2, 3 if is_contact and len(subtracks) > 1 else 2,
                                         width_ratios = [9, 2, 1] if is_contact else [4, 2],
                                         height_ratios = [4, 1],
                                         wspace = 0.5 if is_contact else 0.1,
                                         hspace = 0.6
                                        )
            ax_mat = plt.subplot(subgs[0, 0])
            plot_map(ax_mat, submatrix, loc1, loc2, window, pos1, pos2, is_chrom1_circ, is_chrom2_circ, title=map_title, display_sense=display_sense, chromsizes = chromsizes, display_strand=display_strand, flipped = flipped, strand_level=1.8 if is_contact else 1.75, adjustment= 0.55, cmap=cmap, color=color)
            
            ax_tracks1 = plt.subplot(subgs[1, 0], sharex = ax_mat)
            start, stop  = ax_mat.get_xlim()
            tracks = subtracks[0] if not is_contact else subtracks[1]
            plot_tracks(tracks, ax_tracks1, start, stop, axis = "horizontal",  ylabel = track_unit, flip = (display_sense == "reverse"))
            
            if k >= nb_matrices - cols :
                ax_tracks1.set_xlabel('\nGenomic\ncoordinates in kb')
            if j == 0:
                ax_mat.set_ylabel('Genomic\ncoordinates in kb\n')
            
            if is_contact and len(subtracks) > 1:
                ax_tracks2 = plt.subplot(subgs[0, 1], sharey = ax_mat)
                start, stop  = ax_mat.get_ylim()
                plot_tracks(subtracks[0], ax_tracks2, start, stop, axis = "vertical",  ylabel = track_unit, flip = (display_sense == "reverse"))
                
    outpath = batched_outfolder + f"/{title}"
    if len(outpath) > 0 :
        for format in output_format:
            plt.savefig(outpath + f".{format}", bbox_inches="tight")
        pd.DataFrame(index, columns=['Reference', 'Name']).to_csv(batched_outfolder + f"/{title}_references.csv")
    
    plt.close()

async def display_pileup(pileup, sep_id, cool_name = "", patch_detrending = {}, windows = [], binning = 1000, cmap=None, cmap_color="seismic", title="", outpath="", output_format=['.pdf'], display_strand=True, flipped = False, display_sense="forward", is_contact = False, track_label="Average Track", track_unit=""):
    """Displays a pileup with or without tracks."""
    vmin = None if cmap == None else cmap[0]
    vmax = None if cmap == None else cmap[1]
    xlabel = "\nGenomic coordinates (in kb)"
    ylabel = "Genomic coordinates (in kb)"

    outfolder = f"{outpath}/{pileup.get_cool_name()}/{sep_id}/binning_{pileup.get_binning()}"
    create_folder_path(outfolder)

    for window in windows:
        pileup_matrix = pileup.get_matrix(window)
        if pileup_matrix is None:
            print(f"Empty pileup for {sep_id} with window size {window}bp.") # TODO add to log
            continue
        track_pileup = []
        if pileup.has_tracks(window):
            square_pileup_matrix = pileup_matrix[:pileup.get_size(window)]
            track_pileup = pileup_matrix[pileup.get_size(window):] if len(pileup_matrix) > pileup.get_size(window) else []
            pileup_matrix = square_pileup_matrix

        # applying patch detrending
        identificator = f"{sep_id}_{binning}_{cool_name}"
        if identificator in patch_detrending:
            detrending = patch_detrending[identificator]["pileup"].get_matrix(window)
            detrending_size = patch_detrending[identificator]["pileup"].get_size(window)
            pileup_detrending = detrending[:detrending_size]
            tracks_detrending = detrending[detrending_size:] if len(detrending) > detrending_size else []

            pileup_matrix = pileup_matrix / pileup_detrending
            track_pileup = track_pileup / tracks_detrending

        pileup_title = f"{title} pileup in {pileup.get_cool_name()} \n ({pileup.get_nb_matrices(window)} matrices)"

        pileup_sense = np.flip(pileup_matrix) if display_sense == "reverse" else pileup_matrix
        
        if len(track_pileup) == 0:
            plt.figure(figsize=(6,6))
            plt.title(pileup_title)
            match display_sense:
                case "forward":
                    mat = plt.imshow(np.log10(pileup_sense), extent=[-window//1000, window//1000, window//1000, -window//1000], cmap=cmap_color, vmin=vmin, vmax=vmax)
                case "reverse":
                    mat = plt.imshow(np.log10(pileup_sense), extent=[window//1000, -window//1000, -window//1000, window//1000], cmap=cmap_color, vmin=vmin, vmax=vmax)
            plt.colorbar(mat, fraction=0.01)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)

            if display_strand and pileup.is_directed(axis=1):
                transcription_sens = ARROW_LEFT if display_sense == "reverse" else ARROW_RIGHT
                arrow_alignment = "right" if display_sense == "reverse" else "left"
                to = -window//1000 * 1.2 if display_sense == "reverse" else window//1000 * 1.2
                plt.text(0, to, transcription_sens, horizontalalignment=arrow_alignment, fontsize=20)

        else:
            has_second_exp = is_contact and len(track_pileup) > 1

            width = 3 if has_second_exp else 2
            wspace = 0.4 if has_second_exp else 0.1
            hspace = 0.7 if has_second_exp else 0.5
            ratios = [6, 1, 0.1] if has_second_exp else [6, 0.1]
            figwidth = 8 if has_second_exp else 6
            figheight = 7 if has_second_exp else 8
            plt.figure(figsize=(figwidth,figheight))
            gs = grid.GridSpec(5, width, height_ratios = [1,1,1,1,1], width_ratios = ratios, wspace=wspace, hspace=hspace) 

            # matrix ax
            ax = plt.subplot(gs[:4, 0])
            match display_sense:
                case "forward":
                    mat = plt.imshow(np.log10(pileup_sense), extent=[-window//1000, window//1000, window//1000, -window//1000], cmap=cmap_color, vmin=vmin, vmax=vmax)
                case "reverse":
                    mat = plt.imshow(np.log10(pileup_sense), extent=[window//1000, -window//1000, -window//1000, window//1000], cmap=cmap_color, vmin=vmin, vmax=vmax)

            if display_strand:
                pos = {"Start": 0, "End": 0, "Strand":1}
                plot_strands(ax, pos, pos, window, is_contact=is_contact, display_sense=display_sense, flip = flipped, strand_level = 1.2, adjustment = 0.9)

            # colorbar ax
            ax_cb = plt.subplot(gs[1, width - 1])
            plt.colorbar(mat, fraction=0.01, cax=ax_cb)

            # tracks axes
            track_labelling = track_label + f"\n (in {track_unit})" if track_unit != "" else track_label
            ax_tracks = plt.subplot(gs[4, 0], sharex=ax)

            if not is_contact:
                xstart, xstop = ax.get_xlim()
                plot_tracks(track_pileup[0], ax_tracks, xstart, xstop, axis = "horizontal", xlabel = xlabel, ylabel = track_labelling, flip = (display_sense == "reverse"))

            else:
                # first track
                ax_tracks1 = plt.subplot(gs[:4, 1],sharey = ax)
                ystart, ystop = ax.get_ylim()
                plot_tracks(track_pileup[0], ax_tracks1, ystart, ystop, axis = "vertical", ylabel = track_labelling, flip = (display_sense == "reverse"))

                # second track
                if len(track_pileup) > 1:
                    xstart, xstop = ax.get_xlim()
                    plot_tracks(track_pileup[1], ax_tracks, xstart, xstop, axis = "horizontal", xlabel = xlabel, ylabel = track_labelling, flip = (display_sense == "reverse"))
            
            ax.set_ylabel(ylabel)
            if is_contact:
                plt.suptitle(title)

            if has_second_exp:
                plt.suptitle(pileup_title)
            else:
                ax.set_title(pileup_title)

        if len(outpath) > 0 :
            for format in output_format:
                plt.savefig(outfolder + f"/pileup_{window // 1000}kb_window.{format}", bbox_inches="tight")
        
        plt.close()

# def display_compare(matrix1, matrix2, mat_name1, mat_name2, binning, window, pos1, pos2, positions, position_name, chromsizes={}, output_format=['pdf'], is_pileup=False, outfolder="", is_contact=False, display_sense="forward", display_strand=False, circular=[], cmap=None, cmap_color="afmhot_r", is_global=False):
#     plt.figure(figsize=(16, 5))
#     gs = grid.GridSpec(5, 4, width_ratios=[1, 1, 1, 0.01])

#     ax_mat1 = plt.subplot(gs[:, 0])
#     ax_mat2 = plt.subplot(gs[:, 1])
#     ax_ratio = plt.subplot(gs[:, 2])

#     ax_mat_colorbar = plt.subplot(gs[1, 3])
#     ax_ratio_colorbar = plt.subplot(gs[3, 3])

#     secondary_pos = pos2 if is_contact and pos2 != None else pos1
#     xlabel = "\nGenomic coordinates (in kb)"
#     ylabel = "Genomic coordinates (in kb)"

#     # plotting matrices
#     opti_vmin, opti_vmax = opti_limits([np.log10(matrix1), np.log10(matrix2)])
#     cmap_submat = [opti_vmin, opti_vmax] if cmap == None else cmap

#     # centering log ratio on 0
#     log_ratio = np.log2(matrix1/ matrix2)
#     log_ratio[np.isinf(log_ratio)] = np.nan
#     vmin, vmax = np.nanmin(log_ratio), np.nanmax(log_ratio)
#     cmap_ratio = np.max([abs(vmin), abs(vmax)])

#     if not is_pileup:
#         if not is_global:
#             plot_map(ax_mat1, matrix1, pos1, secondary_pos, window, positions, chromsizes=chromsizes, show_title=False, display_sense=display_sense, display_strand=display_strand, circular=circular, cmap=cmap_submat)
#             im_matrix = plot_map(ax_mat2, matrix2, pos1, secondary_pos, window, positions, chromsizes=chromsizes, show_title=False, display_sense=display_sense, display_strand=display_strand, circular=circular, cmap=cmap_submat)
#             im_ratio = plot_map(ax_ratio, log_ratio, pos1, secondary_pos, window, positions, chromsizes=chromsizes, show_title=False, display_sense=display_sense, display_strand=display_strand, circular=circular, log=False, color = "bwr", cmap=[-cmap_ratio, cmap_ratio])
#         else:
#             plot_global_map(ax_mat1, matrix1, chromsizes, display_sense=display_sense)
#             im_matrix = plot_global_map(ax_mat2, matrix2, chromsizes, display_sense=display_sense)
#             im_ratio = plot_global_map(ax_ratio, log_ratio, chromsizes, display_sense=display_sense, log=False, color = "bwr")

#     else:
#         match display_sense:
#             case "forward":
#                 ax_mat1.imshow(np.log10(matrix1), extent=[-window//1000, window//1000, window//1000, -window//1000], cmap=cmap_color, vmin=cmap_submat[0], vmax=cmap_submat[1])
#                 im_matrix = ax_mat2.imshow(np.log10(matrix2), extent=[-window//1000, window//1000, window//1000, -window//1000], cmap=cmap_color, vmin=cmap_submat[0], vmax=cmap_submat[1])
#                 im_ratio = ax_ratio.imshow(log_ratio, extent=[-window//1000, window//1000, window//1000, -window//1000], cmap="bwr", vmin=-cmap_ratio, vmax=cmap_ratio)
#             case "reverse":
#                 ax_mat1.imshow(np.log10(matrix1), extent=[window//1000, -window//1000, -window//1000, window//1000], cmap=cmap_color, vmin=cmap_submat[0], vmax=cmap_submat[1])
#                 im_matrix = ax_mat2.imshow(np.log10(matrix2), extent=[window//1000, -window//1000, -window//1000, window//1000], cmap=cmap_color, vmin=cmap_submat[0], vmax=cmap_submat[1])
#                 im_ratio = ax_ratio.imshow(log_ratio, extent=[window//1000, -window//1000, -window//1000, window//1000], cmap="bwr", vmin=-cmap_ratio, vmax=cmap_ratio)

#         if display_strand:
#             transcription_sens = ARROW_LEFT if display_sense == "reverse" else ARROW_RIGHT
#             arrow_alignment = "right" if display_sense == "reverse" else "left"
#             to = -window//1000 * 1.2 if display_sense == "reverse" else window//1000 * 1.2
#             ax_mat1.text(0, to, transcription_sens, horizontalalignment=arrow_alignment, fontsize=20)
#             ax_mat2.text(0, to, transcription_sens, horizontalalignment=arrow_alignment, fontsize=20)
#             ax_ratio.text(0, to, transcription_sens, horizontalalignment=arrow_alignment, fontsize=20)

#     # plotting axis and colorbars
#     ax_mat1.set_ylabel(ylabel)
#     ax_mat1.set_xlabel(xlabel)
#     ax_mat2.set_xlabel(xlabel)
#     ax_ratio.set_xlabel(xlabel)

#     plt.colorbar(im_matrix, cax=ax_mat_colorbar)
#     plt.colorbar(im_ratio, cax=ax_ratio_colorbar)

#     # titles
#     if not is_global:
#         name = f"{positions.iloc[pos1]['Name'].replace('/', '_')} submatrices" if not is_contact else f"{positions.iloc[pos1]['Name'].replace('/', '_')}-{positions.iloc[pos2]['Name'].replace('/', '_')} submatrices"
#         if is_pileup:
#             name = position_name + " pileups"
#     else:
#         name = "Complete genome contact matrices"
#     binning_title = f" ({binning//1000}kb binning)"
#     plt.suptitle(name + binning_title)
#     ax_mat1.set_title(mat_name1)
#     ax_mat2.set_title(mat_name2)
#     ax_ratio.set_title(f"{mat_name1}/{mat_name2}")
#     ax_mat_colorbar.set_title("Normalized\ncontact\n(in log10)", fontsize=9)
#     ax_ratio_colorbar.set_title("Log2 ratio", fontsize=9)

#     outpath = outfolder + f"/{name.replace(' ','_')}"
#     if len(outfolder) > 0 :
#         for format in output_format:
#             if not is_global:
#                 plt.savefig(outpath + f".{binning // 1000}kb.{window // 1000}kb_window.{format}", bbox_inches="tight")
#             else:
#                 plt.savefig(outpath + f".{binning // 1000}kb.{format}", bbox_inches="tight")
#     else:
#         plt.show()