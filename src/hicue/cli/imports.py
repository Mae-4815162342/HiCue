# basic librairies
import pandas as pd
import numpy as np

# tools
from itertools import combinations

# plotting
import matplotlib
matplotlib.use('Agg')
from matplotlib import gridspec as grid
from matplotlib import colormaps
import matplotlib.pyplot as plt

# system
from shutil import rmtree
import datetime
import time
import sys
import os

# multiprocessing
from queue import Empty

# biological data 
from chromosight.utils.preprocessing import distance_law # TODO copy the method with reference to avoid costly importation
from BCBio import GFF
import pyBigWig
import cooler

# errors managment
pd.options.mode.chained_assignment = None 
np.seterr(all="ignore")