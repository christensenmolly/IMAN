from astropy.io import fits as pyfits
import numpy as np
import math
import itertools
import matplotlib.pyplot as plt
from scipy import ndimage
import sys
from matplotlib.path import Path
from itertools import product
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import os
import subprocess
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats
import glob
import warnings
from astropy.stats import sigma_clipped_stats


def find_sky(input_image, mask_image=None):
    try:
        hdulist = pyfits.open(input_image)
        data = hdulist[0].data
        ny,nx = np.shape(data)
        if mask_image is not None:
            hdulist_mask = pyfits.open(mask_image)
            mask_data = hdulist_mask[0].data    
        I = []
        if mask_image is not None:
            for k in range(ny):
                for i in range(nx):
                    if mask_data[k,i]==0:
                        I.append(data[k,i])
        else:
            for k in range(ny):
                for i in range(nx):
                        I.append(data[k,i])
        
        mean, median, std = sigma_clipped_stats(I, sigma=3.0, iters=5)
        return mean, median, std
    except:
        return float('nan'),float('nan'),float('nan')