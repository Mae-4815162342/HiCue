# import logging
import click

from .imports import *

from .custom_types import COOL, INT_LIST, STR_LIST, POSITION_FILE, GFF_FILE

import hicue.hicue_opti as h

@click.command("extract")
@click.argument("outpath", type=click.Path(file_okay=False))
@click.argument('positions', type=POSITION_FILE)
@click.argument("cool_files", type=COOL)
@click.option('--gff', type=GFF_FILE, help="Gff file provided for the position file automatic annotation if the file is a bed2d. The positions considered for pileup are all the genes contained in the bed2d files. For more options, use the hicue annotate command.")
@click.option('--tracks', type=GFF_FILE, help="Track file provided for the position file automatic annotation if the file is a bed2d. The positions considered for pileup are all the genes contained in the bed2d files. For more options, use the hicue annotate command.") # TODO re-write
@click.option('-b', '--binnings', type=INT_LIST, default="1000", help="Bin size in bp. Used only if the provided cool files are in mcool format. Several bin sizes can be provided as a comma-separated list. Default value: 1000.")
@click.option('-w', '--windows', type=INT_LIST, default="30000", help="Window size for sub-matrices extraction in bp. Several window sizes can be provided as a comma-separated list. Default value: 30000.")
@click.option('-d', '--detrending',type=click.Choice(['patch', 'ps', 'none'], case_sensitive=False), default='none', help='Detrending option. Default value: none.')
@click.option('-m', '--method', type=click.Choice(['median', 'mean', 'sum'], case_sensitive=False), default='median', help="Aggregation method. If the selected detrending is patch, the method will also be used to aggregate the random sub-matrices. Default value: median.")
@click.option('-f', '--flip', is_flag=True, help="Enables sub-matrices flipping depending on their sense of transcription in the pileups. Requires the strand annotation of provided positions. If not provided, will consider all position in forward.")
@click.option('-r', '--raw', is_flag=True, default=False, help="Use the raw matrices in the cool files (sets balance to False). Default value: False")
@click.option('-t', '--threads', type=int, default=8, help="Number of threads used by each multithreaded worker type. Default: 8.")
@click.option('--nb_pos', type=int, default=2, help="Number of random positions selected for patch detrending. Default value: 2.")
@click.option('--rand_max_dist', type=int, default="100000", help="Maximum distance in bp between provided positions and the random positions selected for patch detrending. Default value: 100000.")
@click.option('--format', type=STR_LIST, default="pdf", help="Figures saving formats. Default value: pdf")
@click.option('--circulars', type=STR_LIST, default="none", help="Coma-separated list of the chromosomes to treat as circular. By default, chromosomes are not considered circular.")
@click.option('--loops', is_flag=True, help="Centers the sub-matrices on pairs of positions instead of single position.")
@click.option('--trans', is_flag=True, help="Enables trans-chromosomal contacts in extracted sub-matrices when in --loops option.")
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
@click.option('--separate_by', type=STR_LIST, help="Comma-separated list of the separation operations. Allowed operations: ['direction', 'regions', 'chroms']. As the option is an inclusion of the separate command, enter: hicue separate --help for more information. ")
@click.option('--separation_regions', type=click.Path(exists=True, dir_okay=False, readable=True, path_type='csv'), help="Path to csv file providing the regions when --separate_by regions mode is selected. As the option is an inclusion of the separate command, enter: hicue separate --help for more information.")# Allows discountinuous interval if provided with the same ID. Csv format: Id,Chromosome,Start,End.")
@click.option('--contact_range', type=(int, int, int), help="Provides MIN MAX STEP in bp as a range for the distance separation on contacts. Overrides the --min_dist option. Default value: (20000,200000,30000).")
@click.option('--save_tmp', is_flag=True, help="Save the temporary files to outpath")
@click.pass_context
def extract(ctx, outpath, positions, cool_files, **params):
        
        # oppening log
        if not os.path.exists(outpath):
                os.mkdir(outpath)
        log = open(f"{outpath}/{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}_log.txt", 'w')
        log.write(f"Extract mode.\nExecuting command: hicue {' '.join(sys.argv[1:])}\n")
        log.write(f"""Extracting from {cool_files}
                positions file: {positions}
                outpath: {outpath}
                
                Options:
                {params}
        """)

        start_time = time.time()

        h.extract(cool_files, positions, outpath, log = log, **params)
        
        end_time = time.time()

        log.write(f"Total time: {end_time - start_time}")
        if log != None:
                log.close()