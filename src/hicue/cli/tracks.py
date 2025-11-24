import click

from .imports import *
from hicue.utils_opti import *

from .custom_types import COOL, INT_LIST, STR_LIST, GFF_FILE, SEPARATOR_LIST, TRACK_FILE

import hicue.hicue_opti as h

@click.command("tracks")
@click.argument("outpath", type=str)
@click.argument('tracks', type=TRACK_FILE)
@click.argument("cool_files", type=COOL)
@click.option('-t', '--threshold', type=(click.Choice(['min', 'max']), float), help="Threshold applied to the tracks for selection in the tracks' unit. If 'min', the threshold is considered as the minimum value for tracks. If 'max' is selected, the threshold is the maximum value.")
@click.option('-p', '--percentage', type=(click.Choice(['high', 'low']), click.IntRange(0, 100)), help="Threshold applied to the tracks for selection in percent. The first parameter indicates whether to take the percentage of high or low values. The second parameter is the percentage between 0 and 100.")# Will be applied after the threshold if --threshold is provided.")
@click.option('--gff', type=GFF_FILE, help="Gff file provided for the position file automatic annotation if the file is a bed2d. The positions considered for pileup are all the genes contained in the bed2d files. For more options, use the hicue annotate command.")
@click.option('--track_unit', type=str, default="", help="Label of the tracks unit axis in display. Default value: Tracks")
@click.option('-b', '--binnings', type=INT_LIST, default="1000", help="Bin size in bp. Used only if the provided cool files are in mcool format. Several bin sizes can be provided as a comma-separated list. Default value: 1000.")
@click.option('-w', '--windows', type=INT_LIST, default="30000", help="Window size for sub-matrices extraction in bp. Several window sizes can be provided as a comma-separated list. Default value: 30000.")
@click.option('-d', '--detrending',type=click.Choice(['patch', 'ps', 'none'], case_sensitive=False), default='none', help='Detrending option. Default value: none.')
@click.option('-m', '--method', type=click.Choice(['median', 'mean', 'sum'], case_sensitive=False), default='median', help="Aggregation method. If the selected detrending is patch, the method will also be used to aggregate the random sub-matrices. Default value: median.")
@click.option('-f', '--flip', is_flag=True, default=False, help="Enables sub-matrices flipping depending on their sense of transcription in the pileups. Requires the strand annotation of provided positions. If not provided, will consider all position in forward.")
@click.option('-r', '--raw', is_flag=True, default=False, help="Use the raw matrices in the cool files (sets balance to False). Default value: False")
@click.option('-e', '--threads', type=int, default=8, help="Number of threads used by each multithreaded worker type. Default: 8.")
@click.option('--nb_pos', type=int, default=2, help="Number of random positions selected for patch detrending. Default value: 2.")
@click.option('--rand_max_dist', type=int, default="100000", help="Maximum distance in bp between provided positions and the random positions selected for patch detrending. Default value: 100000.")
@click.option('--format', type=STR_LIST, default="pdf", help="Figures saving formats. Default value: pdf")
@click.option('--circulars', type=STR_LIST, default="none", help="Coma-separated list of the chromosomes to treat as circular. By default, chromosomes are not considered circular.")
@click.option('--loops', is_flag=True, help="Centers the sub-matrices on pairs of positions instead of single position.")
@click.option('--trans', is_flag=True, help="Enables trans-chromosomal contacts in extracted sub-matrices when in --loops option.")
@click.option('--min_sep', type=int, default="1000", help="Minimal distance separating two selected positions when calling the threshold on tracks. Default value: 1000.")
@click.option('--min_dist', type=int, default="30000", help="Minimal distance in bp between two positions before the pair is used in the --loops option. Default value: 30000.")
@click.option('--diag_mask', type=int, default=0, help="Distance from the diagonal in bp to which the matrix are set to NaN. Is applied only if superior to the bin size. Default value: 0")
@click.option('--pileup/--no-pileup', default=True, help="Compute and display pileups.")
@click.option('--loci/--no-loci', default=False, help="Display single loci as individual figures.")
@click.option('--batch/--no-batch', default=False, help="Display batched loci figures.")
@click.option('--display_strand', is_flag=True, help="Display strands on the single matrices and pileup. Requires the strand annotation of provided positions.")
@click.option('--cmap_limits', type=(float, float), help="Min and Max value for matrix display. Usage: --cmap_limits MIN MAX.")
@click.option('--indiv_cmap_limits', type=(float, float), help="Min and Max value for matrix display. Usage: --cmap_limits MIN MAX.")
@click.option('--cmap_color', type=click.Choice(list(colormaps)), default="seismic", help="Colormap used for pileup. Must be a valid matplotlib colormap. Default: seismic")
@click.option('--indiv_cmap_color', type=click.Choice(list(colormaps)), default="afmhot_r", help="Colormap used for individual displays. Must be a valid matplotlib colormap. Default: afmhot_r")
@click.option('--display_sense', type=click.Choice(['forward', 'reverse'], case_sensitive=False), default='forward', help="Sense of display. In 'forward' mode, the matrices are represented with the forward sense going from left to right, and from right to left in 'reverse' mode.")
@click.option('--center', type=click.Choice(['start', 'center', 'end'], case_sensitive=False), default='start', help="Defines the positional parameter of each position chosen as the window center. 'start' for the start site, 'end' for the end site, and 'center' for the average between those last two.")
@click.option('--overlap', type=click.Choice(['strict', 'flex'], case_sensitive=False), default='strict', help="When evaluating the belonging of a position to an interval, sets the severity of the discrimination: 'strict' will only allow position having start and end positions within the interval, whereas 'flex' will consider the position even if the overlap is not complete. Default value: strict.")
@click.option('--separate_by', type=SEPARATOR_LIST, help=f"Comma-separated list of the separation operations. Allowed operations: {AUTHORIZED_SEPARATORS}. As the option is an inclusion of the separate command, enter: hicue separate --help for more information. ")
@click.option('--separation_regions', type=click.Path(exists=True, dir_okay=False, readable=True), help="Path to csv file providing the regions when --separate_by regions mode is selected. As the option is an inclusion of the separate command, enter: hicue separate --help for more information.")# Allows discountinuous interval if provided with the same ID. Csv format: Id,Chromosome,Start,End.")
@click.option('--contact_range', type=(int, int, int), help="Provides MIN MAX STEP in bp as a range for the distance separation on contacts. Overrides the --min_dist option. Default value: (20000,200000,30000).")
@click.option('--record_type', default="gene", type=click.Choice(['gene']), help="GFF selected records type. If not provided, all records type are selected in the annotation.")
@click.option('--save_tmp', is_flag=True, help="Save the temporary files to outpath")
@click.pass_context
def tracks(ctx, outpath, tracks, cool_files, **params):
        # oppening log
        if not os.path.exists(outpath):
               create_folder_path(outpath)
        log = open(f"{outpath}/{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}_log.txt", 'w')
        log.write(f"Tracks mode.\nExecuting command: hicue {' '.join(sys.argv[1:])}\n")
        log.write(f"""Extracting from {cool_files}
                tracks file: {tracks}
                outpath: {outpath}
        """)
        log.write(f"\nOptions:")
        for param in params:
                log.write(f"""
        {param} = {params[param]}""")
        log.write("\n")

        start_time = time.time()

        h.tracks(cool_files, tracks, outpath, log = log, **params)
        
        end_time = time.time()

        log.write(f"Total time: {end_time - start_time} seconds ({(end_time - start_time)/60} min)")
        if log != None:
                log.close()