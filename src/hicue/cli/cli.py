import click

from hicue import hicue as h

@click.command()
@click.argument('filename')
@click.option("--name", default="World", help="Say hello to NAME.")
def cli(name, filename):
    print(filename)
    """Simple program that greets NAME."""
    h.say_hello(f"Goodbye, {name}!")
    h.extract()

# from . import (
#     extract,
#     tracks
# )

# __all__ = [
#     "extract",
#     "tracks"
# ]

if __name__ == "__main__":
    cli()