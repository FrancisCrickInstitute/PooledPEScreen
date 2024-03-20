"""
[19/03/2024, Michael Herger] data_processing.py
Description: Processing of pegRNA-ST and endogenous target data, including scoring and correlation analyses
"""

###----- IMPORTS ----------------------------------------------------------------------------------------------------###
import os
import sys
import subprocess
import time
import argparse
import pandas as pd
import math
import numpy as np
from functools import reduce
import matplotlib.pyplot as plt
import seaborn as sns
# import auxiliary functions
import surrogate_target_analysis as pooledpe_sta

###----- MAIN -------------------------------------------------------------------------------------------------------###
if __name__ == '__main__':
    # ----- Define input parameters -----#
