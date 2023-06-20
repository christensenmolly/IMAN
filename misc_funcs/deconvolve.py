#! /usr/bin/env python

import sys
from os.path import exists
from os import remove
from scipy.ndimage.filters import convolve
from scipy import signal
import numpy as np
from numpy import max, sum
import pyfits
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
import subprocess
import os
from scipy.optimize import curve_fit
FNULL = open(os.devnull, 'w')

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
sys.path.append(PATH_TO_PACKAGE+'/Decomposition')





def deconvolution_1d(array, psf_array):
    '''
    Function that performs deconvolution of the given array (signal) with the divisor
    '''
    recovered, remainder = signal.deconvolve(array, psf_array)
    return recovered
