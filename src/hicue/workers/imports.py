# basic imports
import numpy as np
import pandas as pd
import importlib
import random
import math

# system
import sys
import os

# threading & multiprocessing
from multiprocessing import Process
from multiprocessing import Queue as MultiQueue
from queue import Queue, Empty
import threading

# scikit-learn
from sklearn.utils.sparsefuncs import _get_median as get_sparse_median
from sklearn.isotonic import IsotonicRegression
from scipy.sparse import csr_matrix

# parsing
from BCBio import GFF
import pyBigWig

# Hi-C librairies
import cooler

# local
from ..classes.Pileup import Pileup
from ..displays_opti import *

# display
default_display = {
    "output_format": ["pdf"],
    "display_strand" : False, 
    "display_sense"  :"forward",
    "indiv_cmap_limits": None,
    "indiv_cmap_color": "afmhot_r"
}
