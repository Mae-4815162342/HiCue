import click

from . import (
    extract,
    tracks,
    extract2d, 
    # separate,
    # annotate
)

@click.group()
@click.pass_context
def cli(ctx):
    pass

cli.add_command(extract.extract)
cli.add_command(tracks.tracks)
cli.add_command(extract2d.extract2d)
# cli.add_command(separate.separate)
# cli.add_command(annotate.annotate)

if __name__ == "__main__":
    cli()