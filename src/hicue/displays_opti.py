from .imports_opti import *

path_lock = threading.Lock()

def create_folder_path(path):
    """If a folder path does not exists, create all the dependencies to this path."""
    global path_lock
    path_list = path.split("/")
    if path_list[0] == "":
        path_list = path_list[1:]
    with path_lock:
        for i in range(1, len(path_list) + 1):
            current_path = "/".join(path_list[:i])
            if not os.path.exists(current_path):
                os.mkdir(current_path)

def adjust_extents(ax, chrom1, chrom2, is_chrom1_circ, is_chrom2_circ, chromsizes={}):
    min, max = ax.get_xlim()
    extent_x = [item.get_text() for item in ax.get_xticklabels() if item._x >= min and item._x < max]
    x_ticks = [tick for tick in ax.get_xticks() if tick >= min and tick < max]
    for i in range(len(extent_x)):
        if extent_x[i][0] == '−':
            extent_x[i] = str(chromsizes[chrom2]//1000 - int(extent_x[i][1:])) if is_chrom2_circ else " "
        elif int(extent_x[i]) > chromsizes[chrom2]//1000:
            extent_x[i] = str(int(extent_x[i]) - chromsizes[chrom2]//1000) if is_chrom2_circ else " "
    ax.set(xticks=x_ticks, xticklabels=extent_x)

    max, min = ax.get_ylim()
    extent_y = [item.get_text() for item in ax.get_yticklabels() if item._y > min and item._y <= max]
    y_ticks = [tick for tick in ax.get_yticks() if tick > min and tick <= max]
    for i in range(len(extent_y)):
        if extent_y[i][0] == '−':
            extent_y[i] = str(chromsizes[chrom1]//1000 - int(extent_y[i][1:])) if is_chrom1_circ else " "
        elif int(extent_y[i]) > chromsizes[chrom1]//1000:
            extent_y[i] = str(int(extent_y[i]) - chromsizes[chrom1]//1000) if is_chrom1_circ else " "
    ax.set(yticks=y_ticks, yticklabels=extent_y)

def plot_map(ax, matrix, loc1, loc2, window, locus1, locus2, is_chrom1_circ, is_chrom2_circ, title="", display_sense="forward", chromsizes={}, display_strand=False, strand_level=1.2, cmap=None, color="afmhot_r", adjust=True, show_title=True, log=True):
    """Plots a single matrix on the provided axis"""
    is_contact = loc1 != loc2
    strand = not is_contact and locus1["Strand"] == -1
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
            mat = ax.imshow(display_matrix, extent=[(pos2 - window)//1000, (pos2 + window)//1000, (pos1 + window)//1000, (pos1 - window)//1000], cmap=color, vmin=vmin, vmax=vmax)
        case "reverse":
            mat = ax.imshow(np.flip(display_matrix), extent=[(pos2 + window)//1000, (pos2 - window)//1000, (pos1 - window)//1000, (pos1 + window)//1000], cmap=color, vmin=vmin, vmax=vmax)
    
    chrom1 = locus1["Chromosome"]
    chrom2 = locus2["Chromosome"]
    if adjust:
        adjust_extents(ax, chrom1, chrom2, is_chrom1_circ, is_chrom2_circ, chromsizes)

    if display_strand:
        match display_sense:
            case "forward":
                transcription_sens = ARROW_LEFT if strand else ARROW_RIGHT
                arrow_alignment = "right" if strand else "left"
                pos_up = (pos1 + window * strand_level)//1000
            case "reverse":
                transcription_sens = ARROW_RIGHT if strand else ARROW_LEFT
                arrow_alignment = "left" if strand else "right"
                pos_up = (pos1 - window * strand_level)//1000
        ax.text(pos2//1000, pos_up, transcription_sens, horizontalalignment=arrow_alignment, fontsize=20)
    return mat

async def display_submatrix(matrix, pair, window, binning = 1000, outfolder="", positions = None, output_format=['pdf'], chromsizes = {}, display_strand=False, display_sense="forward", cmap = None, color = "afmhot_r", tracks = None, track_label = ""):
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

        if tracks == None:

            plt.figure(figsize=(6,6))

            mat = plot_map(plt.gca(), matrix, i, j, window, pos1, pos2, pair["Chrom1_circular"], pair["Chrom2_circular"], title=title, chromsizes = chromsizes, display_sense=display_sense, display_strand=display_strand, strand_level=1.2, cmap=cmap, color=color)

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
            mat = plot_map(ax, matrix, i, j, window, pos1, pos2, pair["Chrom1_circular"], pair["Chrom2_circular"], title=title, show_title=(not is_contact), chromsizes = chromsizes, display_sense=display_sense, display_strand=display_strand, strand_level=1.2, adjust=False, cmap=cmap, color=color)

            # colorbar ax
            ax_cb = plt.subplot(gs[1, width - 1])
            plt.colorbar(mat, fraction=0.01, cax=ax_cb)
            
            flip_tracks = display_sense == "reverse"

            # first track
            ax_track1 = plt.subplot(gs[4, 0], sharex=ax)
            track = tracks[j]
            xstart, xstop = ax.get_xlim()
            index1 = [i/len(track) * (xstop - xstart) + xstart for i in range(len(track))]
            index1 = np.flip(index1) if flip_tracks else index1
            ax_track1.plot(index1, track)
            ax_track1.set_ylabel(track_label)

            # second track
            if is_contact:
                ax_track2 = plt.subplot(gs[:4, 1],sharey = ax)
                track = tracks[i]
                xstart, xstop = ax.get_ylim()
                index2 = [i/len(track) * (xstop - xstart) + xstart for i in range(len(track))]
                index2 = index2 if flip_tracks else np.flip(index2)
                ax_track2.plot(track, index2)
                ax_track2.yaxis.tick_right()
                ax_track2.set_xlabel(track_label)

            adjust_extents(ax, pos1["Chromosome"], pos2["Chromosome"], pair["Chrom1_circular"], pair["Chrom2_circular"], chromsizes)

            plot_xlabel =  "\n" + f"{pos2['Chromosome']} Genomic coordinates (in kb)"
            plot_xlabel += "\n" + x_label if len(x_label) > 0 else ""
            ax_track1.set_xlabel(plot_xlabel)
            plot_ylabel =  f"{pos1['Chromosome']} Genomic coordinates (in kb)"
            plot_ylabel = y_label + "\n" + plot_ylabel if len(y_label) > 0 else plot_ylabel
            ax.set_ylabel(plot_ylabel)

            if is_contact:
                plt.suptitle(title)

        if len(outpath) > 0 :
            for format in output_format:
                plt.savefig(outpath + f".{format}", bbox_inches="tight")
        plt.close()

async def display_batch_submatrices(submatrices, positions, window, title = "", batch_size = 64, outfolder="", output_format=['pdf'], circular=[], chromsizes = [], display_strand=False, display_sense="forward", cmap=None, color="afmhot_r"):
    """Displays a batch of submatrices in a single figure. References each submatrix in a csv file."""
    # checking outpath
    batched_outfolder = f"{outfolder}/batched_{window//1000}kb_window"
    create_folder_path(batched_outfolder)

    cols = math.ceil(math.sqrt(batch_size))
    rows = math.ceil(batch_size / cols)
    nb_matrices = len(submatrices)
    
    plt.figure(figsize=(20,20))
    gs = grid.GridSpec(rows, cols, height_ratios=[1]*rows, width_ratios=[1]*cols, wspace=0.2, hspace=0.5)

    index = []
    n = 1
    for k in range(nb_matrices):
        matrix = submatrices[k]["matrix"]

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
        
        ax = plt.subplot(gs[i, j])
        plot_map(ax, matrix, loc1, loc2, window, pos1, pos2, is_chrom1_circ, is_chrom2_circ, title=map_title, display_sense=display_sense, chromsizes = chromsizes, display_strand=display_strand, strand_level=1.6, cmap=cmap, color=color)

        if i == rows - 1 or i >= nb_matrices // rows:
            ax.set_xlabel('\nGenomic coordinates in kb')
        if j == 0:
            ax.set_ylabel('Genomic coordinates in kb\n')


    outpath = batched_outfolder + f"/{title}"
    if len(outpath) > 0 :
        for format in output_format:
            plt.savefig(outpath + f".{format}", bbox_inches="tight")
        pd.DataFrame(index, columns=['Reference', 'Name']).to_csv(batched_outfolder + f"/{title}_references.csv")

    plt.close()

async def display_pileup(pileup, sep_id, patch_detrending = {}, windows = [], track_pileup=[], cmap=None, cmap_color="seismic", title="", track_title="", outpath="", output_format=['.pdf'], display_strand=True, display_sense="forward", is_contact = False, track_label="Average Track"):
    """Displays a pileup with or without tracks."""
    vmin = None if cmap == None else cmap[0]
    vmax = None if cmap == None else cmap[1]
    xlabel = "\nGenomic coordinates (in kb)"
    ylabel = "Genomic coordinates (in kb)"

    outfolder = f"{outpath}/{pileup.get_cool_name()}/{sep_id}/binning_{pileup.get_binning()}"
    create_folder_path(outfolder)

    for window in windows:
        pileup_matrix = pileup.get_matrix(window)
        # applying patch detrending
        
        if sep_id in patch_detrending:
            pileup_matrix = pileup_matrix / patch_detrending[sep_id]["pileup"].get_matrix(window)
        
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

            if display_strand:
                transcription_sens = ARROW_LEFT if display_sense == "reverse" else ARROW_RIGHT
                arrow_alignment = "right" if display_sense == "reverse" else "left"
                to = -window//1000 * 1.2 if display_sense == "reverse" else window//1000 * 1.2
                plt.text(0, to, transcription_sens, horizontalalignment=arrow_alignment, fontsize=20)

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
            match display_sense:
                case "forward":
                    mat = plt.imshow(np.log10(pileup_sense), extent=[-window//1000, window//1000, window//1000, -window//1000], cmap=cmap_color, vmin=vmin, vmax=vmax)
                case "reverse":
                    mat = plt.imshow(np.log10(pileup_sense), extent=[window//1000, -window//1000, -window//1000, window//1000], cmap=cmap_color, vmin=vmin, vmax=vmax)

            if display_strand:
                transcription_sens = ARROW_LEFT if display_sense == "reverse" else ARROW_RIGHT
                arrow_alignment = "right" if display_sense == "reverse" else "left"
                to = -window//1000 * 1.2 if display_sense == "reverse" else window//1000 * 1.2
                ax.text(0, to, transcription_sens, horizontalalignment=arrow_alignment, fontsize=20)

            # colorbar ax
            ax_cb = plt.subplot(gs[1, width - 1])
            plt.colorbar(mat, fraction=0.01, cax=ax_cb)
            
            flip_tracks = display_sense == "reverse"

            # first track
            ax_track1 = plt.subplot(gs[4, 0], sharex=ax)
            xstart, xstop = ax.get_xlim()
            index1 = [i/len(track_pileup) * (xstop - xstart) + xstart for i in range(len(track_pileup))]
            index1 = np.flip(index1) if flip_tracks else index1
            ax_track1.plot(index1, track_pileup)
            ax_track1.set_ylabel(track_label)

            # second track
            if is_contact:
                ax_track2 = plt.subplot(gs[:4, 1],sharey = ax)
                xstart, xstop = ax.get_ylim()
                index2 = [i/len(track_pileup) * (xstop - xstart) + xstart for i in range(len(track_pileup))]
                index2 = index2 if flip_tracks else np.flip(index2)
                ax_track2.plot(track_pileup, index2)
                ax_track2.yaxis.tick_right()
                ax_track2.set_xlabel(track_label)

            if len(track_title) > 0:
                ax_track1.set_title(track_title, fontsize=11)
                if is_contact:
                    ax_track2.set_ylabel(track_title, fontsize=11)
            ax_track1.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)

            if is_contact:
                plt.suptitle(pileup_title)
            else:
                ax.set_title(pileup_title)

        if len(outpath) > 0 :
            for format in output_format:
                plt.savefig(outfolder + f"/pileup_{window // 1000}kb_window.{format}", bbox_inches="tight")
        
        plt.close()