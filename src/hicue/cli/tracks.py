import click

from .imports import *

from .custom_types import COOL, INT_LIST, STR_LIST, TRACK_FILE, GFF_FILE

import hicue.hicue as h

@click.command("tracks")
@click.argument("outpath", type=click.Path(file_okay=False))
@click.argument('tracks', type=TRACK_FILE)
@click.argument("cool_files", type=COOL)
@click.option('-t', '--threshold', type=(click.Choice(['min', 'max']), float), help="Threshold applied to the tracks for selection in the tracks' unit. If 'min', the threshold is considered as the minimum value for tracks. If 'max' is selected, the threshold is the maximum value.")
@click.option('-p', '--percentage', type=(click.Choice(['high', 'low']), click.IntRange(0, 100)), help="Threshold applied to the tracks for selection in percent. The first parameter indicates whether to take the percentage of high or low values. The second parameter is the percentage between 0 and 100. Will be applied after the threshold if --threshold is provided.")
@click.option('--gff', type=GFF_FILE, help="Gff file provided for the position file automatic annotation if the file is a bed2d. The positions considered for pileup are all the genes contained in the bed2d files. For more options, use the hicue annotate command.")
@click.option('--track_unit', type=str, default="Track Unit", help="Label of the tracks unit axis in display. Default value: Track Unit")
@click.option('--overlap', type=click.Choice(['strict', 'flex'], case_sensitive=False), default='strict', help="When evaluating the belonging of a position to an interval, sets the severity of the discrimination: 'strict' will only allow position having start and end positions within the interval, whereas 'flex' will consider the position even if the overlap is not complete. Default value: strict.")
@click.option('--binning', type=int, default=1000, help='Size of binning in case gff is not provided in bp. All tracks values will therefore be computed for bins of size --binning. Default value: 1000.')
# @click.option('--on_genes', is_flag=True, default=False, help="Project the tracks on genes for selection. The threshold will therefore be applied on the average track value found on each gene.")
@click.option('-w', '--windows', type=INT_LIST, default="30000", help="Window size for sub-matrices extraction in bp. Several window sizes can be provided as a comma-separated list. Default value: 30000.")
@click.option('-d', '--detrending',type=click.Choice(['patch', 'ps', 'none'], case_sensitive=False), default='none', help='Detrending option. Default value: none.')
@click.option('-m', '--pileup_method', type=click.Choice(['median', 'mean'], case_sensitive=False), default='median', help="Aggregation method. If the selected detrending is patch, the method will also be used to aggregate the random sub-matrices. Default value: median.")
@click.option('-f', '--flip', is_flag=True, help="Enables sub-matrices flipping depending on their sense of transcription in the pileups. Requires the strand annotation of provided positions. If not provided, will consider all position in forward.")
@click.option('-r', '--raw', is_flag=True, default=False, help="Use the raw matrices in the cool files (sets balance to False). Default value: False")
@click.option('--nb_pos', type=int, default=2, help="Number of random positions selected for patch detrending. Default value: 2.")
@click.option('--rand_max_dist', type=int, default="100000", help="Maximum distance in bp between provided positions and the random positions selected for patch detrending. Default value: 100000.")
@click.option('--format', type=STR_LIST, default="pdf", help="Figures saving formats. Default value: pdf")
@click.option('--circulars', type=STR_LIST, default="none", help="Coma-separated list of the chromosomes to treat as circular. By default, chromosomes are not considered circular.")
@click.option('--loops', is_flag=True, help="Centers the sub-matrices on pairs of positions instead of single position.")
@click.option('--trans', is_flag=True, help="Enables trans-chromosomal contacts in extracted sub-matrices when in --loops option.")
@click.option('--min_dist', type=int, default="30000", help="Minimal distance in bp between two positions before the pair is used in the --loops option. Default value: 30000.")
@click.option('--diag_mask', type=int, default="0", help="Distance from the diagonal in bp to which the matrix are set to NaN. Is applied only if superior to the bin size. Default value: 0")
@click.option('--pileup/--no-pileup', default=True, help="Compute and display pileups.")
@click.option('--loci/--no-loci', default=False, help="Display single loci as individual figures.")
@click.option('--display_strand', is_flag=True, help="Display strands on the single matrices and pileup. Requires the strand annotation of provided positions.")
@click.option('--cmap_limits', type=(float, float), help="Min and Max value for matrix display. Usage: --cmap_limits MIN MAX.")
@click.option('--cmap_color', type=click.Choice(list(colormaps)), default="seismic", help="Colormap used for pileup. Must be a valid matplotlib colormap. Default: seismic")
@click.option('--display_sense', type=click.Choice(['forward', 'reverse'], case_sensitive=False), default='forward', help="Sense of display. In 'forward' mode, the matrices are represented with the forward sense going from left to right, and from right to left in 'reverse' mode.")
@click.option('--center', type=click.Choice(['start', 'center', 'end'], case_sensitive=False), default='start', help="Defines the positional parameter of each position chosen as the window center. 'start' for the start site, 'end' for the end site, and 'center' for the average between those last two.")
@click.option('--separate_by', type=STR_LIST, default="None", help="Comma-separated list of the separation operations. Allowed operations: ['direction', 'regions', 'chroms']. As the option is an inclusion of the separate command, enter: hicue separate --help for more information. ")
@click.option('--separation_regions', type=click.Path(exists=True, dir_okay=False, readable=True, path_type='csv'), help="Path to csv file providing the regions when --separate_by regions mode is selected. As the option is an inclusion of the separate command, enter: hicue separate --help for more information.")# Allows discountinuous interval if provided with the same ID. Csv format: Id,Chromosome,Start,End.")
@click.option('--contact_separation', type=click.Choice(['cis_trans','distance','none'], case_sensitive=False), default='none', help="Separation option for the --loops.\ncis_trans: cis contacts and trans contact will be plotted separatly. Recomended parameters: --detrending ps --method_pileup mean.\ndistance: using the --contact_range parameter, will compute a selection for each distance in the range. Only applicable on cis contact (if the --trans option is enabled, will not compute the trans contacts). Default value: None.")
@click.option('--contact_range', type=(int, int, int), default=["20000", "200000", "30000"], help="Provides MIN MAX STEP in bp as a range for the distance separation on contacts. Overrides the --min_dist option. Default value: (20000,200000,30000).")
@click.pass_context
def tracks(ctx, 
            outpath,
            tracks,
            cool_files, 
            threshold,
            percentage,
            gff,
            track_unit,
            overlap,
            binning,
            windows, 
            detrending, 
            pileup_method,
            flip,
            raw,
            nb_pos,
            rand_max_dist,
            format,
            circulars,
            loops,
            trans,
            min_dist,
            diag_mask,
            pileup,
            loci,
            display_strand,
            cmap_limits,
            cmap_color,
            display_sense,
            center, 
            separate_by,
            separation_regions,
            contact_separation,
            contact_range
            ):
        
        # oppening log
        if not os.path.exists(outpath):
                os.mkdir(outpath)
        log = open(f"{outpath}/{datetime.datetime.now()}_log.txt", 'w')
        log.write(f"Tracks mode.\nExecuting command: hicue {' '.join(sys.argv[1:])}\n")
        log.write(f"""Extracting from {cool_files}
                tracks file: {tracks}
                outpath: {outpath}

                Options:
                --threshold = {threshold}
                --percentage = {percentage}
                --gff = {gff}
                --track_unit = {track_unit}
                --overlap = {overlap}
                --binning = {binning}
                --windows = {windows}
                --detrending = {detrending}
                --pileup_method = {pileup_method}
                --flip = {flip}
                --raw = {raw}
                --nb_pos = {nb_pos}
                --rand_max_dist = {rand_max_dist}
                --format = {format}
                --circulars = {circulars}
                --loops = {loops}
                --trans = {trans}
                --min-dist = {min_dist}
                --diag_mask = {diag_mask}
                --pileup = {pileup}
                --loci = {loci}
                --display_strand = {display_strand}
                --cmap_limits = {cmap_limits}
                --cmap_color = {cmap_color}
                --display_sense = {display_sense}
                --center = {center}
                --separate_by = {separate_by}
                --separation_regions = {separation_regions}
                --contact_separation = {contact_separation}
                --contact_range = {contact_range}
        """)

        params = {
                "threshold":threshold,
                "percentage":percentage,
                "gff":gff,
                "track_unit":track_unit,
                "overlap":overlap,
                "binning": binning,
                "windows":windows,
                "detrending":detrending,
                "nb_pos":nb_pos,
                "random_max_dist":rand_max_dist,
                "loops":loops,
                "raw":raw,
                "min_dist":min_dist,
                "diagonal_mask":diag_mask,
                "trans_contact":trans,
                "circular_chromosomes":circulars,
                "display_strand":display_strand,
                "output_formats":format,
                "pileup":pileup,
                "loci":loci,
                "method":pileup_method,
                "flip":flip,
                "cmap_pileup": cmap_limits,
                "cmap_color": cmap_color,
                "display_sense":display_sense,
                "center":center,
                "separate_by":separate_by,
                "separation_regions":separation_regions,
                "contact_separation":contact_separation,
                "contact_range":contact_range
        }

        h.tracks(cool_files, tracks, outpath, params, log = log)