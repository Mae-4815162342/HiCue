import click

import hicue.hicue as h

@click.command("annotate")
@click.argument('gff_file', type=click.Path(exists=True, dir_okay=False))
@click.argument('bed_file', type=click.Path(exists=True, dir_okay=False))
@click.option('--overlap', type=click.Choice(['strict', 'flex'], case_sensitive=False), default='flex', help="When evaluating the belonging of a position to an interval, sets the severity of the discrimination: 'strict' will only allow position having start and end positions within the interval, whereas 'flex' will consider the position even if the overlap is not complete. Default value: strict.")
def annotate(gff_file, bed_file):
    h.annotate(gff_file, bed_file)
    
