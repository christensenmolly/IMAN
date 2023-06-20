#!/usr/bin/python
# -*- coding:  cp1251 -*-
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
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
import re
import glob
import pyfits
import collections
from skimage import measure
import pyfits
from skimage.morphology import skeletonize
from photutils import find_peaks
from skimage.transform import (hough_line, hough_line_peaks,
                               probabilistic_hough_line)
from skimage.feature import canny
from skimage import data
from matplotlib import cm
from astropy.stats import sigma_clipped_stats
import argparse

import imp_masking
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

def find_best_lines(angles):
    indices = []; total_angles = []

    for k in range(len(angles)-1):
      for i in range(k+1,len(angles)):
	angle1=angles[k]
	angle2=angles[i]
	if angle1*angle2<0.:  
	  #print 'here'
	  total_angle = abs(angle1)+abs(angle2)
	  total_angles.append(total_angle)
	  indices.append([k,i])
    idx = find_nearest(np.array(total_angles, float),90.)
    [ind1,ind2] = indices[idx]
    #print angles[ind1], angles[ind2],total_angles[idx]
    if total_angles[idx]<=88. or total_angles[idx]>=92.:
      return float('nan')
    else:
      print angles[ind1],angles[ind2]
      #return (np.min([abs(angles[ind1]),abs(angles[ind2])])+90.-np.max([abs(angles[ind1]),abs(angles[ind2])]))/2.
      return np.max([angles[ind1],angles[ind2]])

def find_lines(image):
    # Constructing test image
    '''
    image = np.zeros((100, 100))
    idx = np.arange(25, 75)
    image[idx[::-1], idx] = 255
    image[idx, idx] = 255
    '''
    # Classic straight-line Hough transform
    h, theta, d = hough_line(image)

    # Generating figure 1
    '''
    fig, axes = plt.subplots(1, 3, figsize=(15, 6),
			    subplot_kw={'adjustable': 'box-forced'})
    ax = axes.ravel()

    ax[0].imshow(image, cmap=cm.gray)
    ax[0].set_title('Input image')
    ax[0].set_axis_off()

    ax[1].imshow(np.log(1 + h),
		extent=[np.rad2deg(theta[-1]), np.rad2deg(theta[0]), d[-1], d[0]],
		cmap=cm.gray, aspect=1/1.5)
    ax[1].set_title('Hough transform')
    ax[1].set_xlabel('Angles (degrees)')
    ax[1].set_ylabel('Distance (pixels)')
    ax[1].axis('image')

    ax[2].imshow(image, cmap=cm.gray)
    '''
    angles = []
    for _, angle, dist in zip(*hough_line_peaks(h, theta, d)):
	#print np.degrees(angle), dist
	y0 = (dist - 0 * np.cos(angle)) / np.sin(angle)
	y1 = (dist - image.shape[1] * np.cos(angle)) / np.sin(angle)
	#ax[2].plot((0, image.shape[1]), (y0, y1), '-r')
	angles.append(np.degrees(angle))
    #ax[2].set_xlim((0, image.shape[1]))
    #ax[2].set_ylim((image.shape[0], 0))
    #ax[2].set_axis_off()
    #ax[2].set_title('Detected lines')

    #plt.tight_layout()
    #plt.show()
    '''
    # Line finding using the Probabilistic Hough Transform
    image = data.camera()
    edges = canny(image, 2, 1, 25)
    lines = probabilistic_hough_line(edges, threshold=10, line_length=5,
				    line_gap=3)
    '''
    try:
      angle = find_best_lines(angles)
    except:
      angle = float('nan')
    #print 'PSF Angle is', angle
    return angle




def detect_spike(segm_image, XC, YC):
  hdulist = pyfits.open(segm_image) # open FITS file
  Data = hdulist[0].data
  ny,nx = np.shape(Data)

  X_min = []; Y_min = []; X_max = []; Y_max = []; Angle = []
  for k in range(len(XC)):
    xc = int(math.ceil(XC[k]))
    yc = int(math.ceil(YC[k]))
    (yy,xx) = np.where(Data==Data[yc, xc])
    X_min.append(min(xx))
    Y_min.append(min(yy))

    X_max.append(max(xx))
    Y_max.append(max(yy))

    data_star = Data[min(yy):max(yy),min(xx):max(xx)]
    ny_star,nx_star = np.shape(data_star)
    for kk in range(ny_star):
      for ii in range(nx_star):
	if data_star[kk,ii]==Data[yc, xc]:
	  data_star[kk,ii]=1.
	else:
	  data_star[kk,ii]=0.
    '''
    contours = measure.find_contours(Data, 1.)
    # Display the image and plot all contours found
    fig, ax = plt.subplots()
    ax.imshow(Data, interpolation='nearest', cmap=plt.cm.gray)

    for n, contour in enumerate(contours):
	ax.plot(contour[:, 1], contour[:, 0], linewidth=2)

    ax.axis('image')
    ax.set_xticks([])
    ax.set_yticks([])
    plt.show()
    '''
    #plt.imshow(Data, interpolation='nearest', cmap=plt.cm.gray)
    #plt.show()
    #exit()
    # perform skeletonization
    
    skeleton = skeletonize(data_star)
    angle = find_lines(skeleton)
    print angle,xc,yc
    if np.isnan(angle)==False:
      if angle<0.:
	angle = angle + 90.
      Angle.append(angle)
    '''
    exit()
    # display results
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(8, 4),
			    sharex=True, sharey=True,
			    subplot_kw={'adjustable': 'box-forced'})

    ax = axes.ravel()

    ax[0].imshow(Data, cmap=plt.cm.gray)
    ax[0].axis('off')
    ax[0].set_title('original', fontsize=20)

    ax[1].imshow(skeleton, cmap=plt.cm.gray)
    ax[1].axis('off')
    ax[1].set_title('skeleton', fontsize=20)

    fig.tight_layout()
    plt.show()
    '''

  try:
    mean, median, std = sigma_clipped_stats(Angle, sigma=3.0, iters=5)
    return median, std
  except:
    try:
      return np.mean(Angle), np.std(Angle)
    except:
      try:
	return Angle[0], float('nan')
      except:
	return float('nan'), float('nan')

def main(input_image, npeaks=10, sextr_setup='wise.sex'):
  imp_masking.auto_mask(input_image,'mask_tmp.reg',sextr_setup=sextr_setup, models=False)
  os.remove('mask_tmp.reg')
  segm_image = 'segm.fits'

  # Detect the brightest peaks:
  hdulist = pyfits.open(input_image) # open FITS file  
  Data = hdulist[0].data
  mean, median, std = sigma_clipped_stats(Data, sigma=3.0, iters=5)
  threshold = median + 3.*std
  
  tbl = find_peaks(Data, threshold, box_size=30, npeaks=npeaks)
  x = []; y = []
  for k in range(len(tbl)):
    x.append(float(tbl['x_peak'][k]))
    y.append(float(tbl['y_peak'][k]))

  Angle,Angle_std = detect_spike(segm_image, x, y)
  print 'PSF Angle is', Angle,' +/- ', Angle_std
  os.remove('segm.fits')
  os.remove('field.cat')
  return Angle,Angle_std

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Detect spikes")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("npeaks", nargs='?', const=1, help="Optional: Number of peaks to study",type=int, default=10) 
    args = parser.parse_args()

    input_image = args.input_image
    npeaks = int(args.npeaks)

    Angle,Angle_std = main(input_image, npeaks=npeaks)
