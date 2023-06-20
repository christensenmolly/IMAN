import sys
import math
import numpy as np
#from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import astropy.io.fits as pyfits
from astropy.stats import sigma_clipped_stats
import argparse
import glob
import warnings
warnings.filterwarnings("ignore")




def crea_one_plot(fig_plot=0, xsize=5, ysize=5):
        f = figure(fig_plot, figsize=(xsize,ysize))
        ax = f.add_subplot(111)
        
        return f,ax

def crea_two_row_grids(fig_plot=0, xsize=4, ysize=6.7, left=0.25, right=3.0, hspace=0.0, ratios=[4,1]):
        # Create plots:
        f = figure(fig_plot, figsize=(xsize,ysize))
        gs = gridspec.GridSpec(2, 1, height_ratios=ratios)
        gs.update(left=left, right=right, hspace=hspace)
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        
        return f,ax1,ax2


def crea_three_row_grids(fig_plot=0, xsize=4, ysize=10, left=0.25, right=3.0, hspace=0.0):
        # Create plots:
        f = figure(fig_plot, figsize=(xsize,ysize))
        gs = gridspec.GridSpec(3, 1)
        gs.update(left=left, right=right, hspace=hspace)
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        ax3 = plt.subplot(gs[2])
        
        return f,ax1,ax2,ax3
    

def save_plot(f, output_file, DPI=300, show=False):
        if show:
            plt.show()
        f.savefig(output_file, transparent = False, dpi=DPI, bbox_inches='tight', pad_inches=0.05)
        f.clf()
        plt.close(f)   
        
        
def merge_two_row_grids(ax1, ax2):
        xticklabels = ax1.get_xticklabels()
        setp(xticklabels, visible=False)

def merge_three_row_grids(ax1, ax2, ax3):
        xticklabels = ax1.get_xticklabels()
        setp(xticklabels, visible=False) 
        xticklabels = ax2.get_xticklabels()
        setp(xticklabels, visible=False)         