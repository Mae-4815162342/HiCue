# data management
import pandas as pd
import numpy as np

# UI
import click

# errors managment
pd.options.mode.chained_assignment = None 
np.seterr(all="ignore")

# computation
from collections import defaultdict
from itertools import combinations
from functools import partial
import random
import math

# system
from shutil import rmtree
import importlib
import datetime
import shutil
import time
import sys
import os

# multithreading
from multiprocessing import SimpleQueue
from multiprocessing import Queue as MultiQueue
from multiprocessing import Process
from queue import Queue, Empty
import threading
import asyncio

# matplotlib setting
import matplotlib
matplotlib.use('Agg')

# plotting
from matplotlib import gridspec as grid
from matplotlib import colormaps
import matplotlib.pyplot as plt

# biological data 
from BCBio import GFF
import pyBigWig
import cooler

# scikit-learn & scipy
from sklearn.utils.sparsefuncs import _get_median as get_sparse_median
from sklearn.isotonic import IsotonicRegression
from scipy.sparse import csr_matrix
from scipy.ndimage import zoom

# constants
ARROW_LEFT = "←"
ARROW_RIGHT = "→"
ARROW_UP = "↑"
ARROW_DOWN = "↓"
AUTHORIZED_SEPARATORS = ["direction", "regions", "chroms"]
GTF_FIELDS = {
    "Name":"gene_name",
    "ID": "gene_id"
}
default_display = {
    "output_format": ["pdf"],
    "display_strand" : False, 
    "display_sense"  :"forward",
    "indiv_cmap_limits": None,
    "indiv_cmap_color": "afmhot_r"
}