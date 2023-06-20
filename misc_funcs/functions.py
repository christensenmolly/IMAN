#! /usr/bin/env python
import pylab
import random
from random import gauss
import sys
import os
import math
import numpy as np
from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import NullFormatter
from numpy import *
from pylab import *
tmp_out = sys.stdout
from scipy import stats
import re
from scipy.optimize import fsolve
import pyfits
import subprocess
from scipy.odr.odrpack import *
import argparse
import heapq
from scipy.interpolate import splprep, splev
from scipy.interpolate import interp1d

from photutils import datasets
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import EllipticalAperture
import shutil


def replace_nans(array, max_iter, tol, kernel_size=1, method='localmean'):
    # Initialize arrays
    filled = np.empty( [array.shape[0], array.shape[1]], dtype=np.float64)
    kernel = np.empty( (2*kernel_size+1, 2*kernel_size+1), dtype=np.float64)

    # Indices where array is NaN
    inans, jnans = np.nonzero( np.isnan(array) )
    
    # Number of NaN elements
    n_nans = len(inans)
    
    # Arrays which contain replaced values to check for convergence
    replaced_new = np.zeros( n_nans, dtype=np.float64)
    replaced_old = np.zeros( n_nans, dtype=np.float64)

    # Depending on kernel type, fill kernel array
    if method == 'localmean':
        for i in range(2*kernel_size+1):
            for j in range(2*kernel_size+1):
                kernel[i,j] = 1.

    elif method == 'idw':
        kernel = np.array([[0, 0.5, 0.5, 0.5,0],
                  [0.5,0.75,0.75,0.75,0.5],
                  [0.5,0.75,1,0.75,0.5],
                  [0.5,0.75,0.75,0.5,1],
                  [0, 0.5, 0.5 ,0.5 ,0]])
        #print kernel, 'kernel'

    else:
        raise ValueError( 'method not valid. Should be one of `localmean`.')
    
    # Fill new array with input elements
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            filled[i,j] = array[i,j]

    # Make several passes
    # until we reach convergence
    for it in range(max_iter):
        
        # for each NaN element
        for k in range(n_nans):
            i = inans[k]
            j = jnans[k]
            
            # Initialize to zero
            filled[i,j] = 0.0
            n = 0
            
            # Loop over the kernel
            for I in range(2*kernel_size+1):
                for J in range(2*kernel_size+1):
                   
                    # if we are not out of the boundaries
                    if i+I-kernel_size < array.shape[0] and i+I-kernel_size >= 0:
                        if j+J-kernel_size < array.shape[1] and j+J-kernel_size >= 0:
                                                
                            # if the neighbour element is not NaN itself.
                            if filled[i+I-kernel_size, j+J-kernel_size] == filled[i+I-kernel_size, j+J-kernel_size] :
                                
                                # do not sum itself
                                if I-kernel_size != 0 and J-kernel_size != 0:
                                    
                                    # convolve kernel with original array
                                    filled[i,j] = filled[i,j] + filled[i+I-kernel_size, j+J-kernel_size]*kernel[I, J]
                                    n = n + 1

            # Divide value by effective number of added elements
            if n != 0:
                filled[i,j] = filled[i,j] / n
                replaced_new[k] = filled[i,j]
            else:
                filled[i,j] = 0.0
                
        # Check if mean square difference between values of replaced
        # elements is below a certain tolerance
        if np.mean( (replaced_new-replaced_old)**2 ) < tol:
            break
        else:
            for l in range(n_nans):
                replaced_old[l] = replaced_new[l]
    
    return filled

def split_seq(seq, p):
    newseq = []
    n = len(seq) / p    # min items per subsequence
    r = len(seq) % p    # remaindered items
    b,e = 0, n + min(1, r)  # first split
    for i in range(p):
        newseq.append(seq[b:e])
        r = max(0, r-1)  # use up remainders
        b,e = e, e + n + min(1, r)  # min(1,r) is always 0 or 1

    return newseq  
  
def crea_av_vert_cuts(image_data,x0,y0,radius_x_in,radius_x_out,radius_y,n_bins):
    # Function to create vertical cuts of the image in the given frame
    ny,nx = image_data.shape
    
    x0 = int(x0+0.5)
    radius_x_in = int(radius_x_in)
    radius_x_out = int(radius_x_out)
    n_bins = min([radius_x_out-radius_x_in,n_bins])
    
    y0 = int(y0+0.5)
    radius_y = int(radius_y)
    
    left_grid = range(x0-radius_x_out,x0-radius_x_in,1)
    right_grid = range(x0+radius_x_in,x0+radius_x_out,1)
    
    left_grid = split_seq(left_grid, n_bins)
    right_grid = split_seq(right_grid, n_bins)
    
    RR = []
    for k in range(n_bins):
      RR.append(np.mean(right_grid[k])-x0)
    
    Final_ar = []
    for k in range(n_bins):
	# Left part:
	Array = []
	for i in left_grid[k]:
	  Array.append(image_data[y0-radius_y:y0+radius_y,i])
	Array = np.array(Array)
	left_ar = Array.sum(axis=0)/len(Array)
	
	# Right part:
	Array = []
	for i in right_grid[k]:
	  Array.append(image_data[y0-radius_y:y0+radius_y,i])
	Array = np.array(Array)
	right_ar = Array.sum(axis=0)/len(Array)
	
	final_ar = (left_ar + right_ar)/2.
	Final_ar.append(final_ar)
    
    Final_ar = np.array(Final_ar)
    y = range(-radius_y,radius_y,1)
    '''
    for k in range(len(Final_ar)):
      plt.plot(y,Final_ar[k])
    plt.show()
    '''
    
    return Final_ar,y,RR



'''
imageHDU = pyfits.open('cropped_i.fits')[0]
image_data = imageHDU.data
x0 = 300
y0 = 100
radius_x_in = 50
radius_x_out = 100
radius_y = 25.
n_bins = 10

Final_ar,y = crea_av_vert_cuts(image_data,x0,y0,radius_x_in,radius_x_out,radius_y,n_bins)
'''
