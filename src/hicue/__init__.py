"""
HiCue: A visualization tool for Hi-C datasets.

HiCue provides tools for extracting, analyzing, and visualizing Hi-C chromatin
interaction data with support for various genomic file formats.
"""

from .version import __version__
from . import hicue
from . import parser
from . import displays
from . import utils

__all__ = ["__version__", "hicue", "parser", "displays", "utils"]
