from hicue.imports import *

from . import (
    extract,
    tracks,
    regions,
    # compare
)

@click.group()
@click.pass_context
def cli(ctx):
    pass

cli.add_command(extract.extract)
cli.add_command(tracks.tracks)
cli.add_command(regions.regions)
# cli.add_command(compare.compare)

if __name__ == "__main__":
    cli()