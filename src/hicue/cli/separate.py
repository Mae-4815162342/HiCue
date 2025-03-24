# import logging
import click

from .imports import *

@click.command("separate")
@click.pass_context
def separate(ctx):
    click.echo('Separating')