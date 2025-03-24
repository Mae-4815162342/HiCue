import click

# from hicue import hicue as h

# default_params = {
#     "windows":[30000],
#     "detrending":"patch",
#     "nb_pos":2,
#     "random_max_dist":100000,
#     "loops":False,
#     "min_dist":60000,
#     "diagonal_mask":0,
#     "trans_contact":False,
#     "circular_chromosomes":[],#["chrI"],
#     "display_strand":True,
#     "output_formats":['pdf'],
#     "pileup":True,
#     "loci":False,
#     "method":"median",
#     "flip":True,
#     "cmap_pileup": [-0.25, 0.25],
#     "display_sense":"reverse",
#     "center":"center",
#     "overlap":"strict",
#     "separate_by":"",
#     "separation_regions":"./test_data/TYs.csv",
#     "contact_separation":"cis_trans",
#     "contact_range":"20000:200000:30000"
# }

from . import (
    extract,
    tracks,
    separate
)

@click.group()
@click.pass_context
def cli(ctx):
    pass

cli.add_command(extract.extract)
cli.add_command(tracks.tracks)
cli.add_command(separate.separate)

if __name__ == "__main__":
    cli()