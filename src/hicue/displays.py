"""
displays.py
===========
Matplotlib rendering functions for HiCue submatrices, pileups and tracks.

All public functions in this module are ``async def`` coroutines intended to be
called from an :class:`~hicue.classes.AsyncDisplays.Display` or
:class:`~hicue.classes.AsyncDisplays.DisplayBatch` worker thread via
``loop.run_until_complete()``.

Because HiCue uses the non-interactive ``Agg`` matplotlib backend (set in
``imports.py``), all figures are saved to disk and never shown interactively.

Functions
---------
adjust_extents
    Correct tick labels on circular chromosomes.
plot_strands
    Overlay transcription-strand arrows on a matrix axis.
plot_map
    Render a single Hi-C submatrix onto a given ``Axes``.
plot_tracks
    Plot a 1-D genomic track aligned to a matrix axis.
display_submatrix
    (async) Full single-submatrix figure with optional tracks.
display_batch_submatrices
    (async) Multi-panel batch figure for many submatrices.
display_pileup
    (async) Pileup figure with optional tracks and detrending.
"""

from hicue.utils import *


# ---------------------------------------------------------------------------
# Axis helpers
# ---------------------------------------------------------------------------

def adjust_extents(ax, chrom1, chrom2, is_chrom1_circ, is_chrom2_circ, chromsizes={}):
    """Correct tick labels that overflow chromosome boundaries on circular chromosomes.

    When a window extends past the end of a circular chromosome the genomic
    coordinates are negative or exceed the chromosome length.  This function
    replaces those spurious labels with the correctly wrapped coordinates.
    Non-circular overflows are blanked (replaced with a single space).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The matrix axes whose tick labels should be corrected.
    chrom1 : str
        Name of the chromosome on the Y axis (locus 1).
    chrom2 : str
        Name of the chromosome on the X axis (locus 2).
    is_chrom1_circ : bool
        Whether *chrom1* is circular.
    is_chrom2_circ : bool
        Whether *chrom2* is circular.
    chromsizes : dict[str, int], optional
        Mapping of chromosome name → size in base-pairs.
    """
    # --- X axis (chrom2) ---
    min_x, max_x = np.min(ax.get_xlim()), np.max(ax.get_xlim())
    extent_x = [item.get_text() for item in ax.get_xticklabels()
                if item._x >= min_x and item._x < max_x]
    x_ticks = [tick for tick in ax.get_xticks() if tick >= min_x and tick < max_x]

    for i in range(len(extent_x)):
        if extent_x[i][0] == '−':
            # Negative coordinate: wrap around chromosome end if circular.
            extent_x[i] = (
                str(chromsizes[chrom2] // 1000 - int(extent_x[i][1:]))
                if is_chrom2_circ else " "
            )
        elif int(extent_x[i]) > chromsizes[chrom2] // 1000:
            # Coordinate past chromosome end: subtract chromosome length if circular.
            extent_x[i] = (
                str(int(extent_x[i]) - chromsizes[chrom2] // 1000)
                if is_chrom2_circ else " "
            )
    ax.set(xticks=x_ticks, xticklabels=extent_x)

    # --- Y axis (chrom1) ---
    min_y, max_y = np.min(ax.get_ylim()), np.max(ax.get_ylim())
    extent_y = [item.get_text() for item in ax.get_yticklabels()
                if item._y > min_y and item._y <= max_y]
    y_ticks = [tick for tick in ax.get_yticks() if tick > min_y and tick <= max_y]

    for i in range(len(extent_y)):
        if extent_y[i][0] == '−':
            extent_y[i] = (
                str(chromsizes[chrom1] // 1000 - int(extent_y[i][1:]))
                if is_chrom1_circ else " "
            )
        elif int(extent_y[i]) > chromsizes[chrom1] // 1000:
            extent_y[i] = (
                str(int(extent_y[i]) - chromsizes[chrom1] // 1000)
                if is_chrom1_circ else " "
            )
    ax.set(yticks=y_ticks, yticklabels=extent_y)


def plot_strands(ax, locus1, locus2, window, is_contact=False,
                 display_sense="forward", flip=False,
                 strand_level=1.2, adjustment=0.9):
    """Overlay transcription-direction arrows around a matrix axis.

    Arrows are drawn just outside the matrix extent so they do not overlap the
    data.  The arrow direction and position depend on the strand, display sense
    and whether the view is flipped.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    locus1 : pandas.Series
        Row from the positions DataFrame for locus 1 (Y axis).
        Must have ``"Start"``, ``"End"`` and ``"Strand"`` keys.
    locus2 : pandas.Series
        Row from the positions DataFrame for locus 2 (X axis).
    window : int
        Half-window size in base-pairs used to determine arrow placement.
    is_contact : bool, optional
        ``True`` when the plot shows an inter-locus contact (locus1 ≠ locus2).
    display_sense : {"forward", "reverse"}, optional
        Orientation of the genomic axis.
    flip : bool, optional
        Whether the submatrix was flipped (strand-normalised orientation).
    strand_level : float, optional
        Multiplier controlling how far outside the window the arrow is placed.
        Default 1.2.
    adjustment : float, optional
        Secondary multiplier for the vertical arrow (contact mode only).
        Default 0.9.
    """
    pos1 = (min(locus1["Start"], locus1["End"])
            if locus1["Strand"] == 1
            else max(locus1["Start"], locus1["End"]))
    pos2 = (min(locus2["Start"], locus2["End"])
            if locus2["Strand"] == 1
            else max(locus2["Start"], locus2["End"]))

    x_right = [ARROW_RIGHT, "left"]
    x_left  = [ARROW_LEFT,  "right"]
    y_up    = [ARROW_UP,    "bottom"]
    y_down  = [ARROW_DOWN,  "top"]

    if not is_contact:
        # Single-locus view: only the X arrow is drawn.
        if locus1["Strand"] != 0:
            arrow = x_right
            match display_sense:
                case "forward":
                    arrow = x_left if locus1["Strand"] == -1 and not flip else x_right
                    to_add = -window if flip and locus1["Strand"] == -1 else window
                case "reverse":
                    arrow = x_right if locus1["Strand"] == -1 and not flip else x_left
                    to_add = window if flip and locus1["Strand"] == -1 else -window
            ax.text(pos2 // 1000,
                    (pos1 + to_add * strand_level) // 1000,
                    arrow[0],
                    horizontalalignment=arrow[1],
                    fontsize=20)
    else:
        # Contact view: X arrow for locus2, Y arrow for locus1.
        if locus2["Strand"] != 0:
            arrow = x_right
            match display_sense:
                case "forward":
                    arrow = x_left if locus2["Strand"] == -1 and not flip else x_right
                    to_add = -window if flip and locus1["Strand"] == -1 else window
                case "reverse":
                    arrow = x_right if locus2["Strand"] == -1 and not flip else x_left
                    to_add = window if flip and locus1["Strand"] == -1 else -window
            ax.text(pos2 // 1000,
                    (pos1 + to_add * strand_level) // 1000,
                    arrow[0],
                    horizontalalignment=arrow[1],
                    fontsize=20)

        if locus1["Strand"] != 0:
            arrow = y_down
            match display_sense:
                case "forward":
                    arrow = y_up if locus1["Strand"] == -1 and not flip else y_down
                    to_add = -window if flip and locus2["Strand"] == -1 else window
                case "reverse":
                    arrow = y_down if locus1["Strand"] == -1 and not flip else y_up
                    to_add = window if flip and locus2["Strand"] == -1 else -window
            ax.text((pos2 + to_add * (strand_level * adjustment)) // 1000,
                    pos1 // 1000,
                    arrow[0],
                    verticalalignment=arrow[1],
                    fontsize=20)


def plot_map(ax, matrix, loc1, loc2, window, locus1, locus2,
             is_chrom1_circ, is_chrom2_circ,
             title="", display_sense="forward", chromsizes={},
             display_strand=False, flipped=False,
             cmap=None, color="afmhot_r",
             adjust=True, show_title=True, log=True,
             strand_level=1.2, adjustment=0.9):
    """Render a single Hi-C submatrix onto *ax*.

    The matrix is displayed using ``imshow`` with genomic coordinates (in kb)
    as axis extents.  Optionally applies a log10 transform, adjusts circular-
    chromosome tick labels, and overlays strand arrows.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Target axes.
    matrix : numpy.ndarray
        2-D contact matrix (already square and formatted).
    loc1 : int
        Index of locus 1 in the positions DataFrame (used to detect
        intra-locus vs. inter-locus view).
    loc2 : int
        Index of locus 2 in the positions DataFrame.
    window : int
        Half-window size in base-pairs.
    locus1 : pandas.Series
        Metadata row for locus 1 (``Chromosome``, ``Start``, ``End``,
        ``Strand``, ``Name``).
    locus2 : pandas.Series
        Metadata row for locus 2.
    is_chrom1_circ : bool
        Whether chromosome 1 is circular.
    is_chrom2_circ : bool
        Whether chromosome 2 is circular.
    title : str, optional
        Axes title.  Defaults to the locus name(s).
    display_sense : {"forward", "reverse"}, optional
        Orientation of the genomic axis.
    chromsizes : dict[str, int], optional
        Chromosome name → size mapping for tick correction.
    display_strand : bool, optional
        Whether to overlay strand arrows.
    flipped : bool, optional
        Whether the matrix was strand-normalised (flipped).
    cmap : list[float, float] or None, optional
        ``[vmin, vmax]`` colour scale limits.  ``None`` uses auto-scaling.
    color : str, optional
        Matplotlib colormap name.  Default ``"afmhot_r"``.
    adjust : bool, optional
        Whether to call :func:`adjust_extents` for circular chromosomes.
    show_title : bool, optional
        Whether to set the axes title.
    log : bool, optional
        Apply ``log10`` transform before rendering.  Default ``True``.
    strand_level : float, optional
        Passed to :func:`plot_strands`.
    adjustment : float, optional
        Passed to :func:`plot_strands`.

    Returns
    -------
    matplotlib.image.AxesImage
        The image object (useful for adding a colour bar).
    """
    is_contact = loc1 != loc2
    pos1 = (min(locus1["Start"], locus1["End"])
            if locus1["Strand"] == 1
            else max(locus1["Start"], locus1["End"]))
    pos2 = (min(locus2["Start"], locus2["End"])
            if locus2["Strand"] == 1
            else max(locus2["Start"], locus2["End"]))

    name = (f"{locus1['Name'].replace('/', '_')}"
            if not is_contact
            else f"{locus1['Name'].replace('/', '_')}-{locus2['Name'].replace('/', '_')}")
    title = name if len(title) == 0 else title
    if show_title:
        ax.set_title(title)

    # Log-transform (suppressing divide-by-zero warnings on zero contacts).
    with np.errstate(divide='ignore', invalid='ignore'):
        display_matrix = np.log10(matrix) if log else matrix

    vmin = cmap[0] if cmap is not None else None
    vmax = cmap[1] if cmap is not None else None

    match display_sense:
        case "forward":
            y_extent = (
                [(pos1 + window) // 1000, (pos1 - window) // 1000]
                if locus1["Strand"] != -1 or not flipped
                else [(pos1 - window) // 1000, (pos1 + window) // 1000]
            )
            x_extent = (
                [(pos2 - window) // 1000, (pos2 + window) // 1000]
                if locus2["Strand"] != -1 or not flipped
                else [(pos2 + window) // 1000, (pos2 - window) // 1000]
            )
            mat = ax.imshow(display_matrix,
                            extent=x_extent + y_extent,
                            cmap=color, vmin=vmin, vmax=vmax)
        case "reverse":
            y_extent = (
                [(pos1 - window) // 1000, (pos1 + window) // 1000]
                if locus1["Strand"] != -1 or not flipped
                else [(pos1 + window) // 1000, (pos1 - window) // 1000]
            )
            x_extent = (
                [(pos2 + window) // 1000, (pos2 - window) // 1000]
                if locus2["Strand"] != -1 or not flipped
                else [(pos2 - window) // 1000, (pos2 + window) // 1000]
            )
            mat = ax.imshow(np.flip(display_matrix),
                            extent=x_extent + y_extent,
                            cmap=color, vmin=vmin, vmax=vmax)

    chrom1 = locus1["Chromosome"]
    chrom2 = locus2["Chromosome"]
    if adjust:
        adjust_extents(ax, chrom1, chrom2, is_chrom1_circ, is_chrom2_circ, chromsizes)

    if display_strand:
        plot_strands(ax, locus1, locus2, window,
                     is_contact=(pos1 != pos2),
                     display_sense=display_sense,
                     flip=flipped,
                     strand_level=strand_level,
                     adjustment=adjustment)

    return mat

def plot_map_region(ax, matrix, region1, region2, padding, is_contact = False,
             title="", display_sense="forward",
             display_strand=False, flipped=False,
             cmap=None, color="afmhot_r",
             show_title=True, log=True,strand_level=1.2, adjustment=0.9):

    name = (f"{region1['Name'].replace('/', '_')}"
            if not is_contact
            else f"{region1['Name'].replace('/', '_')}-{region2['Name'].replace('/', '_')}")
    title = name if len(title) == 0 else title
    if show_title:
        ax.set_title(title)

    # Log-transform (suppressing divide-by-zero warnings on zero contacts).
    with np.errstate(divide='ignore', invalid='ignore'):
        display_matrix = np.log10(matrix) if log else matrix

    vmin = cmap[0] if cmap is not None else None
    vmax = cmap[1] if cmap is not None else None

    forward_labels = ["Start", "End"]
    reverse_labels = ["End", "Start"]
    padding_size = padding * matrix.shape[0] /(2 * padding + 1)
    
    match display_sense:
        case "forward":
            axis_labels = forward_labels if region2["Strand"] != -1 or not flipped else reverse_labels
            mat = ax.imshow(display_matrix, cmap=color, vmin=vmin, vmax=vmax, interpolation = None)
            
        case "reverse":
            axis_labels = reverse_labels if region2["Strand"] != -1 or not flipped else forward_labels
            mat = ax.imshow(np.flip(display_matrix), cmap=color, vmin=vmin, vmax=vmax, interpolation = None)
            
    ax_start, ax_end = ax.get_xlim()
    labels_indexes = [ax_start + padding_size, ax_end - padding_size]
    
    ax.set_xticks(labels_indexes, axis_labels)
    ax.set_yticks(labels_indexes, axis_labels)

    # TODO: to re-write for region mode
    # if display_strand:
    #     plot_strands(ax, region1, region2, window,
    #                  is_contact=(region1['Name'] != region2['Name']),
    #                  display_sense=display_sense,
    #                  flip=flipped,
    #                  strand_level=strand_level,
    #                  adjustment=adjustment)

    return mat


def plot_tracks(tracks, ax_tracks, start, stop,
                axis="horizontal", xlabel="", ylabel="", flip=False):
    """Plot a 1-D genomic track on *ax_tracks* aligned to the matrix axis.

    Parameters
    ----------
    tracks : array-like
        1-D array of binned track values (e.g. ChIP-seq signal).
    ax_tracks : matplotlib.axes.Axes
        Target axes for the track.
    start : float
        Start coordinate (in kb) of the matrix axis to align to.
    stop : float
        Stop coordinate (in kb) of the matrix axis to align to.
    axis : {"horizontal", "vertical"}, optional
        Whether the track runs horizontally (below the matrix) or vertically
        (to the right of the matrix, for contact plots).  Default horizontal.
    xlabel : str, optional
        Label for the X axis.
    ylabel : str, optional
        Label for the Y axis (or the track value axis for vertical mode).
    flip : bool, optional
        Whether to flip the index array so the track aligns with a
        reverse-strand display.  Default ``False``.
    """
    index = np.array([i / len(tracks) * (stop - start) + start
                      for i in range(len(tracks))])

    match axis:
        case "horizontal":
            index = np.flip(index) if flip else index
            ax_tracks.plot(index, tracks)
            if len(ylabel) > 0:
                ax_tracks.set_ylabel(ylabel)
            if len(xlabel) > 0:
                ax_tracks.set_xlabel(xlabel)
        case "vertical":
            # For vertical tracks the index maps to the Y axis.
            index = np.flip(index) if not flip else index
            ax_tracks.plot(tracks, index)
            ax_tracks.yaxis.tick_right()
            # Swap label assignment: ylabel goes on X, xlabel on Y.
            if len(ylabel) > 0:
                ax_tracks.set_xlabel(ylabel)
            if len(xlabel) > 0:
                ax_tracks.set_ylabel(xlabel)


# ---------------------------------------------------------------------------
# Async display functions
# ---------------------------------------------------------------------------

async def display_submatrix(
    matrix, pair, size_metric,
    binning=1000, outfolder="", positions=None,
    output_format=['pdf'], chromsizes={},
    display_strand=False, flipped=False, display_sense="forward",
    cmap=None, color="afmhot_r",
    display_tracks=False, track_unit="",
    is_region = False, padding = None, display_log = True
):
    """Render and save a single submatrix figure.

    Produces either a simple matrix plot or a combined matrix+tracks layout
    depending on *display_tracks* and whether track data is present in
    *matrix*.

    The figure is saved under *outfolder*/individual_{window//1000}kb_window/
    in each of the requested *output_format*\\ s.

    Parameters
    ----------
    matrix : numpy.ndarray
        Combined matrix/track array.  Rows ``0 … N-1`` form the square contact
        matrix (where ``N = len(matrix[0])``); any additional rows are binned
        track values.
    pair : pandas.Series
        Formatted pair row.  Must contain ``"Locus1"``, ``"Locus2"``,
        ``"Chrom1_circular"``, ``"Chrom2_circular"``.
    window : int
        Half-window size in base-pairs.
    binning : int, optional
        Bin size in base-pairs.  Default 1000.
    outfolder : str, optional
        Root output directory.  Skips saving if empty.
    positions : pandas.DataFrame, optional
        Full positions DataFrame indexed by locus index.
    output_format : list[str], optional
        List of file-format suffixes (e.g. ``["pdf", "png"]``).
    chromsizes : dict[str, int], optional
        Chromosome name → size mapping.
    display_strand : bool, optional
        Whether to overlay strand arrows.
    flipped : bool, optional
        Whether the matrix was strand-normalised.
    display_sense : {"forward", "reverse"}, optional
        Axis orientation.
    cmap : list[float, float] or None, optional
        ``[vmin, vmax]`` colour scale limits.
    color : str, optional
        Matplotlib colormap name.
    display_tracks : bool, optional
        Whether to render the track rows.
    track_unit : str, optional
        Unit label for the track Y axis.
    """
    # Resolve output sub-folder.
    individual_outfolder = f"{outfolder}/individual_{size_metric // 1000}kb_window" if not is_region else f"{outfolder}/individual_{size_metric}px"
    create_folder_path(individual_outfolder)

    i, j = pair["Locus1"], pair["Locus2"]
    pos1, pos2 = positions.loc[i], positions.loc[j]
    is_contact = i != j

    title = (f"Window centered on\n{pos1['Name']} vs {pos2['Name']}"
             if is_contact
             else f"Window centered on\n{positions.loc[i]['Name']}")
    y_label = f"{pos1['Name']}" if is_contact else ""
    x_label = f"{pos2['Name']}" if is_contact else ""
    outpath = (
        ""
        if len(individual_outfolder) == 0
        else (
            f"{individual_outfolder}/{pos1['Name'].replace('/', '_')}"
            if not is_contact
            else f"{individual_outfolder}/{pos1['Name'].replace('/', '_')}-{pos2['Name'].replace('/', '_')}"
        )
    )

    # Separate the square matrix from any appended track rows.
    submatrix = matrix[: len(matrix[0])]
    subtracks = matrix[len(matrix[0]) :] if len(matrix) > len(matrix[0]) else []

    if len(subtracks) == 0 or not display_tracks:
        # --- Simple matrix-only layout ---
        plt.figure(figsize=(6, 6))
        if is_region:
            mat = plot_map_region(plt.gca(), submatrix, pos1, pos2, padding, 
                is_contact=is_contact,
                title=title, display_sense=display_sense,
                display_strand=display_strand, flipped=flipped,
                cmap=cmap, color=color, log=display_log, 
                strand_level=1.2)
        else:
            mat = plot_map(
                plt.gca(), submatrix, i, j, size_metric, pos1, pos2,
                pair["Chrom1_circular"], pair["Chrom2_circular"],
                title=title, chromsizes=chromsizes,
                display_sense=display_sense, display_strand=display_strand,
                flipped=flipped, strand_level=1.2, cmap=cmap, color=color,
                log = display_log
            )
        plt.colorbar(mat, fraction=0.01)
        plot_xlabel = "\n" + f"{pos2['Chromosome']} Genomic coordinates (in kb)" if not is_region else "\n" + f"{pos2['Chromosome']} Resized genomic coordinates"
        plot_xlabel += "\n" + x_label if len(x_label) > 0 else ""
        plt.xlabel(plot_xlabel)
        plot_ylabel = f"{pos1['Chromosome']} Genomic coordinates (in kb)" if not is_region else "\n" + f"{pos1['Chromosome']} Resized genomic coordinates"
        plot_ylabel = y_label + "\n" + plot_ylabel if len(y_label) > 0 else plot_ylabel
        plt.ylabel(plot_ylabel)

    else:
        # --- Matrix + tracks layout ---
        width   = 3 if is_contact else 2
        wspace  = 0.4 if is_contact else 0.1
        hspace  = 0.7 if is_contact else 0.5
        ratios  = [6, 1, 0.1] if is_contact else [6, 0.1]
        figwidth  = 8 if is_contact else 6
        figheight = 7 if is_contact else 8

        plt.figure(figsize=(figwidth, figheight))
        gs = grid.GridSpec(
            5, width,
            height_ratios=[1, 1, 1, 1, 1],
            width_ratios=ratios,
            wspace=wspace, hspace=hspace,
        )

        # Matrix sub-plot (top 4 rows, left column).
        ax = plt.subplot(gs[:4, 0])

        if is_region:
            mat = plot_map_region(ax, submatrix, pos1, pos2, padding, 
                is_contact=is_contact,
                title=title, display_sense=display_sense,
                display_strand=display_strand, flipped=flipped,
                cmap=cmap, color=color, log=display_log, 
                strand_level=1.8 if is_contact else 1.75,
                adjust=False, show_title=(not is_contact))
        else:
        
            mat = plot_map(
                ax, submatrix, i, j, size_metric, pos1, pos2,
                pair["Chrom1_circular"], pair["Chrom2_circular"],
                title=title, show_title=(not is_contact),
                chromsizes=chromsizes,
                display_sense=display_sense, display_strand=display_strand,
                flipped=flipped, strand_level=1.2, adjust=False,
                cmap=cmap, color=color, log = display_log
            )

        # Colour bar (small sub-plot, right column).
        ax_cb = plt.subplot(gs[1, width - 1])
        plt.colorbar(mat, fraction=0.01, cax=ax_cb)

        track_labelling = track_unit
        ax_tracks = plt.subplot(gs[4, 0], sharex=ax)
        plot_xlabel = "\n" + f"{pos2['Chromosome']} Genomic coordinates (in kb)" if not is_region else "\n" + f"{pos2['Chromosome']} Resized genomic coordinates"
        plot_xlabel += "\n" + x_label if len(x_label) > 0 else ""

        if not is_contact:
            # Single-locus: one horizontal track below the matrix.
            xstart, xstop = ax.get_xlim()
            plot_tracks(
                subtracks[0], ax_tracks, xstart, xstop,
                axis="horizontal",
                xlabel=plot_xlabel,
                ylabel=track_labelling,
                flip=(display_sense == "reverse"),
            )
        else:
            # Contact: vertical track on the right + optional horizontal track.
            ax_tracks1 = plt.subplot(gs[:4, 1], sharey=ax)
            ystart, ystop = ax.get_ylim()
            plot_tracks(
                subtracks[0], ax_tracks1, ystart, ystop,
                axis="vertical",
                ylabel=track_labelling,
                flip=(display_sense == "reverse"),
            )
            if len(subtracks) > 1:
                xstart, xstop = ax.get_xlim()
                plot_tracks(
                    subtracks[1], ax_tracks, xstart, xstop,
                    axis="horizontal",
                    xlabel=plot_xlabel,
                    ylabel=track_labelling,
                    flip=(display_sense == "reverse"),
                )

        adjust_extents(
            ax,
            pos1["Chromosome"], pos2["Chromosome"],
            pair["Chrom1_circular"], pair["Chrom2_circular"],
            chromsizes,
        )

        plot_ylabel = f"{pos1['Chromosome']} Genomic coordinates (in kb)" if not is_region else "\n" + f"{pos2['Chromosome']} Resized genomic coordinates"
        plot_ylabel = y_label + "\n" + plot_ylabel if len(y_label) > 0 else plot_ylabel
        ax.set_ylabel(plot_ylabel)

        if is_contact:
            plt.suptitle(title)

    if len(outpath) > 0:
        for fmt in output_format:
            plt.savefig(outpath + f".{fmt}", bbox_inches="tight")
    plt.close()


async def display_batch_submatrices(
    submatrices, positions, size_metric,
    title="", batch_size=64, outfolder="",
    output_format=['pdf'], chromsizes={},
    display_strand=False, flipped=False, display_sense="forward",
    display_tracks=False, track_unit="",
    cmap=None, color="afmhot_r", is_region = False, padding = None, display_log = True
):
    """Render and save a multi-panel batch figure.

    Lays out up to *batch_size* submatrices in a grid (≈ square) and writes an
    accompanying CSV index mapping panel labels to locus names.

    The figure is saved under *outfolder*/batched_{window//1000}kb_window/.

    Parameters
    ----------
    submatrices : list[dict]
        Each dict must contain ``"matrix"`` (numpy array) and ``"pair"``
        (formatted pair Series with ``"Locus1"``, ``"Locus2"``,
        ``"Chrom1_circular"``, ``"Chrom2_circular"``).
    positions : pandas.DataFrame
        Full positions DataFrame.
    window : int
        Half-window size in base-pairs.
    title : str, optional
        Figure/file title prefix (used as file stem).
    batch_size : int, optional
        Maximum number of panels.  Controls grid dimensions.  Default 64.
    outfolder : str, optional
        Root output directory.
    output_format : list[str], optional
        File format suffixes.
    chromsizes : dict[str, int], optional
        Chromosome name → size mapping.
    display_strand : bool, optional
        Whether to overlay strand arrows on each panel.
    flipped : bool, optional
        Whether matrices were strand-normalised.
    display_sense : {"forward", "reverse"}, optional
        Axis orientation.
    display_tracks : bool, optional
        Whether to render appended track rows.
    track_unit : str, optional
        Unit label for track axes.
    cmap : list[float, float] or None, optional
        ``[vmin, vmax]`` colour limits.
    color : str, optional
        Matplotlib colormap name.
    """
    batched_outfolder = f"{outfolder}/batched_{size_metric // 1000}kb_window" if not is_region else f"{outfolder}/batched_{size_metric}px"
    create_folder_path(batched_outfolder)

    cols = math.ceil(math.sqrt(batch_size))
    rows = math.ceil(batch_size / cols)
    nb_matrices = len(submatrices)

    plt.figure(figsize=(20, 20) if not display_tracks else (23, 23))
    gs = grid.GridSpec(
        rows, cols,
        height_ratios=[1] * rows,
        width_ratios=[1] * cols,
        wspace=0.4 if display_tracks else 0.2,
        hspace=0.5,
    )

    index = []  # CSV reference: [[label, name], …]
    n = 1       # Panel counter for human-readable labels.

    for k in range(nb_matrices):
        matrix = submatrices[k]["matrix"]

        submatrix = matrix[: len(matrix[0])]
        subtracks = matrix[len(matrix[0]) :] if len(matrix) > len(matrix[0]) else []

        # GridSpec is row-major; map flat index to (row, col).
        i = k // rows
        j = k % rows

        loc1 = submatrices[k]["pair"]["Locus1"]
        loc2 = submatrices[k]["pair"]["Locus2"]
        is_chrom1_circ = submatrices[k]["pair"]["Chrom1_circular"]
        is_chrom2_circ = submatrices[k]["pair"]["Chrom2_circular"]
        pos1, pos2 = positions.loc[loc1], positions.loc[loc2]
        is_contact = loc1 != loc2

        name = (
            f"{pos1['Name'].replace('/', '_')}"
            if not is_contact
            else f"{pos1['Name'].replace('/', '_')}-{pos2['Name'].replace('/', '_')}"
        )
        map_title = f"Position {n}" if not is_contact else f"Contact {n}"
        n += 1
        index.append([map_title, name])

        if not display_tracks:
            # --- Simple grid cell (matrix only) ---
            ax = plt.subplot(gs[i, j])
            if is_region:
                plot_map_region(ax, submatrix, pos1, pos2, padding, 
                            is_contact=is_contact,
                            title=map_title, display_sense=display_sense,
                            display_strand=display_strand, flipped=flipped,
                            cmap=cmap, color=color, log=display_log, 
                            strand_level=1.6, adjustment=0.6,)
                if k >= nb_matrices - cols:
                    ax.set_xlabel('\nResized\ngenomic\ncoordinates')
                if j == 0:
                    ax.set_ylabel('Resized\ngenomic\ncoordinates\n')
            else:
                plot_map(
                    ax, submatrix, loc1, loc2, size_metric, pos1, pos2,
                    is_chrom1_circ, is_chrom2_circ,
                    title=map_title, display_sense=display_sense,
                    chromsizes=chromsizes, display_strand=display_strand,
                    flipped=flipped, strand_level=1.6, adjustment=0.6,
                    cmap=cmap, color=color, log = display_log
                )
                if k >= nb_matrices - cols:
                    ax.set_xlabel('\nGenomic\ncoordinates in kb')
                if j == 0:
                    ax.set_ylabel('Genomic\ncoordinates in kb\n')

        else:
            # --- Grid cell with embedded track sub-plots ---
            subgs = gs[i, j].subgridspec(
                2,
                3 if is_contact and len(subtracks) > 1 else 2,
                width_ratios=[9, 2, 1] if is_contact else [4, 2],
                height_ratios=[4, 1],
                wspace=0.5 if is_contact else 0.1,
                hspace=0.6,
            )
            ax_mat = plt.subplot(subgs[0, 0])
            if is_region:
                plot_map_region(ax, submatrix, pos1, pos2, padding, 
                            is_contact=is_contact,
                            title=map_title, display_sense=display_sense,
                            display_strand=display_strand, flipped=flipped,
                            cmap=cmap, color=color, log=display_log, 
                            strand_level=1.8 if is_contact else 1.75,
                            adjustment=0.55,)
            else:
                plot_map(
                    ax_mat, submatrix, loc1, loc2, size_metric, pos1, pos2,
                    is_chrom1_circ, is_chrom2_circ,
                    title=map_title, display_sense=display_sense,
                    chromsizes=chromsizes, display_strand=display_strand,
                    flipped=flipped,
                    strand_level=1.8 if is_contact else 1.75,
                    adjustment=0.55,
                    cmap=cmap, color=color, log = display_log
                )

            ax_tracks1 = plt.subplot(subgs[1, 0], sharex=ax_mat)
            start, stop = ax_mat.get_xlim()
            # For contact plots the horizontal track is index 1 (X axis locus).
            tracks = subtracks[0] if not is_contact else subtracks[1]
            plot_tracks(
                tracks, ax_tracks1, start, stop,
                axis="horizontal",
                ylabel=track_unit,
                flip=(display_sense == "reverse"),
            )
    
            if k >= nb_matrices - cols:
                if is_region:
                    ax_tracks1.set_xlabel('\nResized\ngenomic\ncoordinates')
                else:
                    ax_tracks1.set_xlabel('\nGenomic\ncoordinates in kb')
            if j == 0:
                if is_region:
                    ax_mat.set_ylabel('Resized\ngenomic\ncoordinates\n')
                else:
                    ax_mat.set_ylabel('Genomic\ncoordinates in kb\n')

            if is_contact and len(subtracks) > 1:
                # Vertical track to the right (Y axis locus).
                ax_tracks2 = plt.subplot(subgs[0, 1], sharey=ax_mat)
                start, stop = ax_mat.get_ylim()
                plot_tracks(
                    subtracks[0], ax_tracks2, start, stop,
                    axis="vertical",
                    ylabel=track_unit,
                    flip=(display_sense == "reverse"),
                )

    outpath = batched_outfolder + f"/{title}"
    if len(outpath) > 0:
        for fmt in output_format:
            plt.savefig(outpath + f".{fmt}", bbox_inches="tight")
        pd.DataFrame(index, columns=["Reference", "Name"]).to_csv(
            batched_outfolder + f"/{title}_references.csv"
        )
    plt.close()


async def display_pileup(
    pileup, sep_id,
    cool_name="", patch_detrending={}, size_metrics=[],
    binning=1000, cmap=None, cmap_color="seismic",
    title="", outpath="", output_format=['.pdf'],
    display_strand=True, flipped=False, display_sense="forward",
    is_contact=False, track_label="Average Track", track_unit="", is_region = False, padding = None,
):
    """Render and save pileup figures (one per window size).

    Optionally applies patch detrending (divides the pileup by a null-model
    pileup) and overlays averaged track data.  One figure is produced per
    entry in *size_metrics*.

    The figures are saved under
    *outpath*/{cool_name}/{sep_id}/binning_{binning}/.

    Parameters
    ----------
    pileup : hicue.classes.Pileup.Pileup
        Aggregated pileup object.
    sep_id : str
        Separation group identifier used in the output path and figure title.
    cool_name : str, optional
        Name of the source cool file (used in the output path).
    patch_detrending : dict, optional
        Mapping of ``"{sep_id}_{binning}_{cool_name}"`` → ``{"pileup": Pileup}``
        for patch detrending.  Empty dict disables detrending.
    size_metrics : list[int], optional
        List of expected or window size in pixels or base-pairs respectively.  One figure per size_metrics.
    binning : int, optional
        Bin size in base-pairs.  Default 1000.
    cmap : list[float, float] or None, optional
        ``[vmin, vmax]`` colour scale limits.
    cmap_color : str, optional
        Matplotlib colormap name.  Default ``"seismic"``.
    title : str, optional
        Figure title prefix.
    outpath : str, optional
        Root output directory.  Skips saving if empty.
    output_format : list[str], optional
        File format suffixes.
    display_strand : bool, optional
        Whether to overlay strand arrows on the pileup.
    flipped : bool, optional
        Whether matrices were strand-normalised.
    display_sense : {"forward", "reverse"}, optional
        Axis orientation.
    is_contact : bool, optional
        Whether the pileup represents inter-locus contacts.
    track_label : str, optional
        Label for the averaged track.
    track_unit : str, optional
        Unit string appended to *track_label*.
    """
    vmin = None if cmap is None else cmap[0]
    vmax = None if cmap is None else cmap[1]
    xlabel = "\nGenomic coordinates (in kb)" if not is_region else "\nResized genomic coordinates"
    ylabel = "Genomic coordinates (in kb)" if not is_region else "Resized genomic coordinates"

    outfolder = f"{outpath}/{pileup.get_cool_name()}/{sep_id}/binning_{pileup.get_binning()}"
    create_folder_path(outfolder)

    for size_metric in size_metrics:
        pileup_matrix = pileup.get_matrix(size_metric)
        has_tracks = pileup.has_tracks(size_metric)

        if pileup_matrix is None:
            sentence = f"Empty pileup for {sep_id} with window size {size_metric}bp." if not is_region else f"Empty pileup for {sep_id} with expected size of {size_metric}px."
            print(sentence)
            continue

        track_pileup = []
        if has_tracks:
            square_pileup_matrix = pileup_matrix[: pileup.get_size(size_metric)]
            track_pileup = (
                pileup_matrix[pileup.get_size(size_metric):]
                if len(pileup_matrix) > pileup.get_size(size_metric)
                else []
            )
            pileup_matrix = square_pileup_matrix

        # Apply patch detrending if a null-model pileup is available.
        identificator = f"{sep_id}_{binning}_{cool_name}"
        if identificator in patch_detrending:
            detrending = patch_detrending[identificator]["pileup"].get_matrix(size_metric)
            detrending_size = patch_detrending[identificator]["pileup"].get_size(size_metric)
            pileup_detrending = detrending[:detrending_size]
            tracks_detrending = (
                detrending[detrending_size:]
                if len(detrending) > detrending_size
                else []
            )
            pileup_matrix = pileup_matrix / pileup_detrending
            if has_tracks:
                track_pileup = track_pileup / tracks_detrending

        pileup_title = (
            f"{title} pileup in {pileup.get_cool_name()} \n"
            f"({pileup.get_nb_matrices(size_metric)} matrices)"
        )

        pileup_sense = np.flip(pileup_matrix) if display_sense == "reverse" else pileup_matrix

        if len(track_pileup) == 0:
            # --- Simple pileup (matrix only) ---
            plt.figure(figsize=(6, 6))
            plt.title(pileup_title)
            if is_region:
                mat = plt.imshow(
                    np.log10(pileup_sense),
                    cmap=cmap_color, vmin=vmin, vmax=vmax,
                    interpolation = None
                )
                padding_size = padding * pileup_sense.shape[0] /(2 * padding + 1)

                axis_labels = ["Start", "End"] if display_sense == "forward" else ["End", "Start"]
                ax_start, ax_end = plt.gca().get_xlim()
                labels_indexes = [ax_start + padding_size, ax_end - padding_size]
                plt.xticks(labels_indexes, axis_labels)
                plt.yticks(labels_indexes, axis_labels)

            else:
                match display_sense:
                    case "forward":
                        mat = plt.imshow(
                            np.log10(pileup_sense),
                            extent=[-size_metric // 1000, size_metric // 1000,
                                    size_metric // 1000, -size_metric // 1000],
                            cmap=cmap_color, vmin=vmin, vmax=vmax,
                        )
                    case "reverse":
                        mat = plt.imshow(
                            np.log10(pileup_sense),
                            extent=[size_metric // 1000, -size_metric // 1000,
                                    -size_metric // 1000,  size_metric // 1000],
                            cmap=cmap_color, vmin=vmin, vmax=vmax,
                        )
            plt.colorbar(mat, fraction=0.01)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)

            if display_strand:
                transcription_sens = ARROW_LEFT if display_sense == "reverse" else ARROW_RIGHT
                arrow_alignment = "right" if display_sense == "reverse" else "left"
                to = -size_metric // 1000 * 1.2 if display_sense == "reverse" else size_metric // 1000 * 1.2
                plt.text(0, to, transcription_sens,
                         horizontalalignment=arrow_alignment, fontsize=20)

        else:
            # --- Pileup + tracks layout ---
            has_second_exp = is_contact and len(track_pileup) > 1
            width     = 3 if has_second_exp else 2
            wspace    = 0.4 if has_second_exp else 0.1
            hspace    = 0.7 if has_second_exp else 0.5
            ratios    = [6, 1, 0.1] if has_second_exp else [6, 0.1]
            figwidth  = 8 if has_second_exp else 6
            figheight = 7 if has_second_exp else 8

            plt.figure(figsize=(figwidth, figheight))
            gs = grid.GridSpec(
                5, width,
                height_ratios=[1, 1, 1, 1, 1],
                width_ratios=ratios,
                wspace=wspace, hspace=hspace,
            )

            ax = plt.subplot(gs[:4, 0])
            if is_region:
                mat = ax.imshow(
                    np.log10(pileup_sense),
                    cmap=cmap_color, vmin=vmin, vmax=vmax,
                    interpolation = None
                )
                padding_size = padding * pileup_sense.shape[0] /(2 * padding + 1)

                axis_labels = ["Start", "End"] if display_sense == "forward" else ["End", "Start"]
                ax_start, ax_end = ax.get_xlim()
                labels_indexes = [ax_start + padding_size, ax_end - padding_size]
                ax.set_xticks(labels_indexes, axis_labels)
                ax.set_yticks(labels_indexes, axis_labels)

            else:
                match display_sense:
                    case "forward":
                        mat = ax.imshow(
                            np.log10(pileup_sense),
                            extent=[-size_metric // 1000, size_metric // 1000,
                                    size_metric // 1000, -size_metric // 1000],
                            cmap=cmap_color, vmin=vmin, vmax=vmax,
                        )
                    case "reverse":
                        mat = ax.imshow(
                            np.log10(pileup_sense),
                            extent=[size_metric // 1000, -size_metric // 1000,
                                    -size_metric // 1000,  size_metric // 1000],
                            cmap=cmap_color, vmin=vmin, vmax=vmax,
                        )

                if display_strand:
                    pos = {"Start": 0, "End": 0, "Strand": 1}
                    plot_strands(ax, pos, pos, size_metric,
                                is_contact=is_contact,
                                display_sense=display_sense,
                                flip=flipped,
                                strand_level=1.2, adjustment=0.9)

            ax_cb = plt.subplot(gs[1, width - 1])
            plt.colorbar(mat, fraction=0.01, cax=ax_cb)

            track_labelling = (
                track_label + f"\n (in {track_unit})"
                if track_unit != ""
                else track_label
            )
            ax_tracks = plt.subplot(gs[4, 0], sharex=ax)

            if not is_contact:
                xstart, xstop = ax.get_xlim()
                plot_tracks(
                    track_pileup[0], ax_tracks, xstart, xstop,
                    axis="horizontal",
                    xlabel=xlabel,
                    ylabel=track_labelling,
                    flip=(display_sense == "reverse"),
                )
            else:
                ax_tracks1 = plt.subplot(gs[:4, 1], sharey=ax)
                ystart, ystop = ax.get_ylim()
                plot_tracks(
                    track_pileup[0], ax_tracks1, ystart, ystop,
                    axis="vertical",
                    ylabel=track_labelling,
                    flip=(display_sense == "reverse"),
                )
                if len(track_pileup) > 1:
                    xstart, xstop = ax.get_xlim()
                    plot_tracks(
                        track_pileup[1], ax_tracks, xstart, xstop,
                        axis="horizontal",
                        xlabel=xlabel,
                        ylabel=track_labelling,
                        flip=(display_sense == "reverse"),
                    )

            ax.set_ylabel(ylabel)
            if is_contact:
                plt.suptitle(title)

            if has_second_exp:
                plt.suptitle(pileup_title)
            else:
                ax.set_title(pileup_title)

        if len(outpath) > 0:
            for fmt in output_format:
                plt.savefig(
                    outfolder + f"/pileup_{size_metric // 1000}kb_window.{fmt}" if not is_region else outfolder + f"/pileup_{size_metric}px.{fmt}",
                    bbox_inches="tight",
                )
        plt.close()