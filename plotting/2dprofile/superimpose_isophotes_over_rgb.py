#! /usr/bin/env python


import sys
import math
import numpy as np
from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import NullFormatter
from mpl_toolkits.axes_grid1 import AxesGrid
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
import re
import glob
from matplotlib.colors import LogNorm
import astropy.io.fits as pyfits
import argparse
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.stats import sigma_clipped_stats
import scipy.ndimage as ndimage
import warnings
from PIL import Image



matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
#matplotlib.use('agg',warn=False, force=True)   ####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

warnings.filterwarnings("ignore")




def convert_polygonline_to_shapleyregion(line):
    coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
    pol = []
    for kk in range(0,len(coords)-1,2):
        pol.append((int(float(coords[kk])),int(float(coords[kk+1]))))
    polygon = Polygon(pol)
    return polygon

def main(input_rgb_image, region_file, output_image, xlim=None, ylim=None, text=None):
    image = Image.open(input_rgb_image)
    
    fig = plt.figure(0, (8, 8))
    ax = fig.add_subplot(111)
    
    imgplot = ax.imshow(np.flipud(image), origin="lower", aspect='equal')
    
    f = open(region_file, 'r')
    lines = f.readlines()
    for line in lines:
        if 'polygon(' in line:
            poly = convert_polygonline_to_shapleyregion(line)
            x,y = poly.exterior.xy
            ax.plot(x, y, color='red', alpha=0.7, linewidth=0.5, solid_capstyle='round', zorder=2)
            
            #exit()
    #plt.show()

    if xlim is not None:
        ax.set_xlim(xlim[0],xlim[1])

    if ylim is not None:
        ax.set_ylim(ylim[0],ylim[1])

    if text is not None:
        ax.text(0.03, 0.9, text, fontsize=17, color='black',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline')

    ax.set_frame_on(False)
    ax.set_xticks([]); ax.set_yticks([])
    plt.axis('off')
    

    plt.savefig(output_image, bbox_inches='tight', pad_inches=0)
    plt.clf()
    plt.close()


#main('fitscut.jpg', 'isophotes.reg',xlim=[27,779],ylim=[297,502], text='NGC7332')
#main('fitscut.jpg', 'isophotes.reg',xlim=[166,636],ylim=[316,486], text='NGC3869')

#main('fitscut_inv.jpg', 'isophotes.reg', 'NGC7332_isophotes.png',xlim=[27,779],ylim=[297,502], text='NGC7332')
main('fitscut_inv.jpg', 'isophotes.reg', 'NGC3869_isophotes.png', xlim=[166,636],ylim=[316,486], text='NGC3869')