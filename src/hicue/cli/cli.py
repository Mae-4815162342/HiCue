from hicue.imports import *

from . import (
    extract,
    tracks,
    # compare
)

@click.group()
@click.pass_context
def cli(ctx):
    pass

cli.add_command(extract.extract)
cli.add_command(tracks.tracks)
# cli.add_command(compare.compare)

if __name__ == "__main__":
    cli()