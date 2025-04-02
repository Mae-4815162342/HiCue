from .cli.imports import *

ARROW_LEFT = "←"
ARROW_RIGHT = "→"

def adjust_extents(ax, chrom1, chrom2, circular=[], chromsizes={}):
    extent_x = [item.get_text() for item in ax.get_xticklabels()]
    for i in range(len(extent_x)):
        if extent_x[i][0] == '−':
            extent_x[i] = str(chromsizes[chrom2]//1000 - int(extent_x[i][1:])) if chrom2 in circular else " "
        elif int(extent_x[i]) > chromsizes[chrom2]//1000:
            extent_x[i] = str(int(extent_x[i]) - chromsizes[chrom2]//1000) if chrom2 in circular else " "
    ax.set_xticklabels(extent_x)

    extent_y = [item.get_text() for item in ax.get_yticklabels()]
    for i in range(len(extent_y)):
        if extent_y[i][0] == '−':
            extent_y[i] = str(chromsizes[chrom1]//1000 - int(extent_y[i][1:])) if chrom1 in circular else " "
        elif int(extent_y[i]) > chromsizes[chrom1]//1000:
            extent_y[i] = str(int(extent_y[i]) - chromsizes[chrom1]//1000) if chrom1 in circular else " "
    ax.set_yticklabels(extent_y)

def plot_map(ax, matrix, loc1, loc2, window, locus, title="", display_sense="forward", circular = [], chromsizes={}, display_strand=False, strand_level=1.2, cmap=None, adjust=True, show_title=True):
    """Plots a single matrix on the provided axis"""
    is_contact = loc1 != loc2
    strand = not is_contact and locus.iloc[loc1]["Strand"] == -1
    pos1 = min(locus.iloc[loc1]["Start"], locus.iloc[loc1]["End"]) if locus.iloc[loc1]["Strand"] == 1 else max(locus.iloc[loc1]["Start"], locus.iloc[loc1]["End"])
    pos2 = min(locus.iloc[loc2]["Start"], locus.iloc[loc2]["End"]) if locus.iloc[loc2]["Strand"] == 1 else max(locus.iloc[loc2]["Start"], locus.iloc[loc2]["End"])

    name = f"{locus.iloc[loc1]['Name'].replace('/', '_')}" if not is_contact else f"{locus.iloc[loc1]['Name'].replace('/', '_')}-{locus.iloc[loc2]['Name'].replace('/', '_')}"
    title = name if len(title) == 0 else title
    if show_title:
        ax.set_title(title)

    log_matrix = np.log10(matrix)
    vmin = cmap[0] if cmap != None else None
    vmax = cmap[1] if cmap != None else None

    match display_sense:
        case "forward":
            mat = ax.imshow(log_matrix, extent=[(pos2 - window)//1000, (pos2 + window)//1000, (pos1 + window)//1000, (pos1 - window)//1000], cmap="afmhot_r", vmin=vmin, vmax=vmax)
        case "reverse":
            mat = ax.imshow(np.flip(log_matrix), extent=[(pos2 + window)//1000, (pos2 - window)//1000, (pos1 - window)//1000, (pos1 + window)//1000], cmap="afmhot_r", vmin=vmin, vmax=vmax)
    
    chrom1 = locus.iloc[loc1]["Chromosome"]
    chrom2 = locus.iloc[loc2]["Chromosome"]
    if adjust:
        adjust_extents(ax, chrom1, chrom2, circular, chromsizes)

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

def display_submatrices(submatrices, locus, window, outfolder="", output_format=['pdf'], circular=[], chromsizes = {}, display_strand=False, display_sense="forward", tracks = None, track_label = ""):
    if len(outfolder) > 0 and not os.path.exists(outfolder):
        os.mkdir(outfolder)
    if len(outfolder) > 0 and not os.path.exists(outfolder + "/indiviudal_displays"):
        os.mkdir(outfolder  + "/indiviudal_displays")
    
    for _, row in submatrices.iterrows():
        i, j, matrix = row["Loc1"], row["Loc2"], row["Matrix"]
        is_contact = i!=j
        reshape_size = int(np.sqrt(len(matrix)))
        matrix = matrix.reshape((reshape_size, reshape_size))
        title = f"Window centered on\n{locus.iloc[i]['Name']} vs {locus.iloc[j]['Name']}" if is_contact else f"Window centered on\n{locus.iloc[i]['Name']}"
        y_label = f"{locus.iloc[i]['Name']}" if is_contact else ""
        x_label = f"{locus.iloc[j]['Name']}" if is_contact else ""
        outpath = "" if len(outfolder) == 0 else f"{outfolder}/indiviudal_displays/{locus.iloc[i]['Name'].replace('/', '_')}" if not is_contact else f"{outfolder}/indiviudal_displays/{locus.iloc[i]['Name'].replace('/', '_')}-{locus.iloc[j]['Name'].replace('/', '_')}"

        if tracks == None:

            plt.figure(figsize=(6,6))

            mat = plot_map(plt.gca(), matrix, i, j, window, locus, title=title, circular=circular, chromsizes = chromsizes, display_sense=display_sense, display_strand=display_strand, strand_level=1.2)

            plt.colorbar(mat, fraction=0.01)
            plot_xlabel =  "\n" + f"{locus.iloc[j]['Chromosome']} Genomic coordinates (in kb)"
            plot_xlabel += "\n" + x_label if len(x_label) > 0 else ""
            plt.xlabel(plot_xlabel)
            plot_ylabel =  f"{locus.iloc[i]['Chromosome']} Genomic coordinates (in kb)"
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
            mat = plot_map(ax, matrix, i, j, window, locus, title=title, show_title=(not is_contact), circular=circular, chromsizes = chromsizes, display_sense=display_sense, display_strand=display_strand, strand_level=1.2, adjust=False)

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

            adjust_extents(ax, locus.iloc[i]["Chromosome"], locus.iloc[j]["Chromosome"], circular, chromsizes)

            plot_xlabel =  "\n" + f"{locus.iloc[j]['Chromosome']} Genomic coordinates (in kb)"
            plot_xlabel += "\n" + x_label if len(x_label) > 0 else ""
            ax_track1.set_xlabel(plot_xlabel)
            plot_ylabel =  f"{locus.iloc[i]['Chromosome']} Genomic coordinates (in kb)"
            plot_ylabel = y_label + "\n" + plot_ylabel if len(y_label) > 0 else plot_ylabel
            ax.set_ylabel(plot_ylabel)

            if is_contact:
                plt.suptitle(title)

        if len(outpath) > 0 :
            for format in output_format:
                plt.savefig(outpath + f".{format}", bbox_inches="tight")
        else:
            plt.show()

def display_pileup(pileup, window, track_pileup=[], cmap=None, title="", outpath="", output_format=['.pdf'], display_strand=True, display_sense="forward", is_contact = False, track_label="Average Track"):
    vmin = None if cmap == None else cmap[0]
    vmax = None if cmap == None else cmap[1]
    xlabel = "\nGenomic coordinates (in kb)"
    ylabel = "Genomic coordinates (in kb)"
    pileup_sense = np.flip(pileup) if display_sense == "reverse" else pileup
    if len(pileup.shape) < 2:
        reshape_size = int(np.sqrt(len(pileup)))
        pileup_sense = pileup_sense.reshape((reshape_size, reshape_size))
    if len(track_pileup) == 0:
        plt.figure(figsize=(6,6))
        plt.title(title)
        match display_sense:
            case "forward":
                mat = plt.imshow(np.log10(pileup_sense), extent=[-window//1000, window//1000, window//1000, -window//1000], cmap="seismic", vmin=vmin, vmax=vmax)
            case "reverse":
                mat = plt.imshow(np.log10(pileup_sense), extent=[window//1000, -window//1000, -window//1000, window//1000], cmap="seismic", vmin=vmin, vmax=vmax)
        plt.colorbar(mat, fraction=0.01)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

        if display_strand:
            transcription_sens = ARROW_LEFT if display_sense == "reverse" else ARROW_RIGHT
            arrow_alignment = "right" if display_sense == "reverse" else "left"
            to = -window//1000 * 1.2 if display_sense == "reverse" else window//1000 * 1.2
            plt.text(0, to, horizontalalignment=arrow_alignment, fontsize=20)

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
                mat = plt.imshow(np.log10(pileup_sense), extent=[-window//1000, window//1000, window//1000, -window//1000], cmap="seismic", vmin=vmin, vmax=vmax)
            case "reverse":
                mat = plt.imshow(np.log10(pileup_sense), extent=[window//1000, -window//1000, -window//1000, window//1000], cmap="seismic", vmin=vmin, vmax=vmax)

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

        ax_track1.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)

        if is_contact:
            plt.suptitle(title)
        else:
            ax.set_title(title)

    if len(outpath) > 0 :
        for format in output_format:
            plt.savefig(outpath + f".{format}", bbox_inches="tight")
    else:
        plt.show()

def display_all_submatrices(submatrices, locus, window, outfolder="", output_format=['pdf'], circular=[], chromsizes = [], display_strand=False, display_sense="forward"):
    if len(outfolder) > 0 and not os.path.exists(outfolder):
        os.mkdir(outfolder)
    
    indexes = {}
    for index, row in submatrices.iterrows():
        i, j = row["Loc1"], row["Loc2"]
        indexes[(i,j)] = index
    keys = np.array(list(indexes.keys()))
    # displaying in batches of 64 matrices
    max_value = len(submatrices)//64 if len(submatrices) % 64 == 0 else len(submatrices)//64 + 1
    index = []
    n = 1
    for k in range(0, max_value):
        current_keys = keys[k * 64: (k+1)*64] if len(keys) > (k+1)*64 else keys[k * 64:]
        plt.figure(figsize=(20,20))
        gs = grid.GridSpec(8, 8, height_ratios=[1]*8, width_ratios=[1]*8, wspace=0.2, hspace=0.5)

        i, j = 0, 0
        for key in current_keys:
            matrix = submatrices.iloc[indexes[tuple(key)]]["Matrix"]
            reshape_size = int(np.sqrt(len(matrix)))
            matrix = matrix.reshape((reshape_size, reshape_size))
            loc1, loc2 = key
            is_contact = loc1 != loc2
            name = f"{locus.iloc[loc1]['Name'].replace('/', '_')}" if not is_contact else f"{locus.iloc[loc1]['Name'].replace('/', '_')}-{locus.iloc[loc2]['Name'].replace('/', '_')}"
            title = f"Position {n}" if not is_contact else f"Contact {n}"
            n += 1
            index.append([title, name])
            
            ax = plt.subplot(gs[j, i])
            plot_map(ax, matrix, loc1, loc2, window, locus, title=title, display_sense=display_sense,  circular=circular, chromsizes = chromsizes, display_strand=display_strand, strand_level=1.6)

            if j == 8 - 1:
                ax.set_xlabel('\nGenomic coordinates in kb')
            if i % 8 == 0:
                ax.set_ylabel('Genomic coordinates in kb\n')

            i += 1
            if i % 8 == 0:
                j += 1
                i = 0

        outpath = outfolder + f"/all_submatrix_batch_{(k + 1)}_out_of_{max_value}"
        if len(outpath) > 0 :
            for format in output_format:
                plt.savefig(outpath + f".{format}", bbox_inches="tight")
        else:
            plt.show()
    
    if len(outpath) > 0 :
        pd.DataFrame(index, columns=['Reference', 'Name']).to_csv(outfolder + "/all_submatrix_references.csv")
            