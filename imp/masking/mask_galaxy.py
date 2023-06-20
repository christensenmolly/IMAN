#!/usr/bin/python

# Import the necessary modules
import pyfits
import numpy as np
import math
import itertools
import matplotlib
matplotlib.use('TkAgg')
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
from photutils import data_properties, properties_table
import numpy.ma as ma


from photutils import detect_threshold
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources



from skimage.draw import ellipse_perimeter
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats

from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources

import shapely
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon

sys.path.append('/Users/mosenkov/MEGA/MyPrograms/IMAN/IMP_NEW')
import convert_segm_to_region_new
import imp_masking
import mask_indiv
import random




def inner_mask(input_image, x0, y0, max_radius, q, PA, sigma=2., enlarge_masks=1.2, enlarge_masks_type='multiply', min_radius=5., step = 1., FWHM=3., npixels=5, output_mask_file='inner_mask.fits', output_reg_file='inner_mask.reg'):
    '''
    DESCRIPTION:
    To create an inner mask we draw ellipses with increasing radius from min_radius
    to max_radius, with some step. We collect all pixels which intersect these ellipse pixels
    and use sigma clipping to mask pixels which have intensities larger than sigma*mean.
    After that we detect sources using photutils to create a segmentaion map from the created mask,
    by removing masked pixels which are not in the group (npixels).
    Then we smooth data with a kernel and again detect sources which are grouped.
    Finally, we make the masked area larger (proportionally, by some factor), with a final region file
    which can be revisited.
    
    enlarge_masks_type can be 'multiply' (the poplygons will made larger by the factor of enlarge_masks) or 'add' (value of enlarge_masks will be added or subtracted).  
    '''
    
    print 'Inner masking ...'

    if (enlarge_masks!=1. and enlarge_masks_type=='multiply') or (enlarge_masks!=0. and enlarge_masks_type=='add'):
        if enlarge_masks_type=='multiply':
            offset = [enlarge_masks, 1.]
        else:
            offset = [enlarge_masks]
    else:
        offset = None

    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    header = hdulist[0].header
    
    ny,nx = np.shape(data)
    
    PA = np.radians(180. - PA) # It should be orientation of the Major axis in clockwise direction (in radians)
    
    radii = np.arange(min_radius, max_radius, step) # Create a set of radii, use a smaller step (at least 1 pix) 
    # The minimual radius should not be small, otherwise the procedure may mask the centre.
    
    mask = np.zeros(np.shape(data))
    
    for radius in radii:
        #print '\tRadius:\t',radius
        # Draw ellipse in pixel coordinates (using Bresenham's approach):
        rr, cc = ellipse_perimeter(int(x0), int(y0), int(radius), int(radius*q), PA, shape=None)
        
        # Remove those pixels which are outside the image
        RR = []; CC = []
        for kk in range(len(rr)):
            if rr[kk]>=0 and cc[kk]>=0 and rr[kk]<nx and cc[kk]<ny:
                RR.append(rr[kk])
                CC.append(cc[kk])
        try:
            mean, median, std = sigma_clipped_stats(data[CC,RR], sigma=3., iters=5)        
            
            for k in range(len(RR)):
                if data[CC[k],RR[k]]>sigma*mean:
                    mask[CC[k],RR[k]]=1.
        except:
            z=1 #TODO: This happend if we are outside of the image borders, but I have solved this above (check)
    
    # Create segmneation which is above a threshold, grouped by npixels - this is important to remove single masked pixels (actually, many such pixels can be found after this masking)!
    segm = detect_sources(mask, 0.1, npixels=npixels, filter_kernel=None)
            
    # Create a kernel for averaging - this important to mask some missing areas which apperantly should be masked
    Sigma = FWHM * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(Sigma, x_size=int(math.ceil(FWHM))*2, y_size=int(math.ceil(FWHM))*2)
    kernel.normalize()
    
    # Create new segmentation map for the averaged data - to get more roundy masks
    segm_new = detect_sources(segm, 0.1, npixels=npixels, filter_kernel=kernel)
    
    outHDU = pyfits.PrimaryHDU(segm_new, header=header)
    outHDU.writeto(output_mask_file, clobber=True)       

    # Convert segmentation mask to region file with some extension of the masked areas, defined by factor_enlarge_masks (this is important to mask the outskirts of the masked objects as well). This mask then can be changed by hand.
    
    #TODO as below:
    if offset is not None:  
        convert_segm_to_region_new.main(output_mask_file, output_reg_file, output_mask_image=None, offset=offset,xc=x0,yc=y0)
    print '\t Done!'



def outer_mask(input_image, x0, y0, sma, q, PA, snr=2., enlarge_masks=1.2, enlarge_masks_type='multiply', output_mask_file='outer_mask.fits', output_reg_file='outer_mask.reg'):
    '''
    DESCRIPTION:
    Function to create an outer mask around the galaxy which is surrounded by the given ellipse.
    '''
    
    if (enlarge_masks!=1. and enlarge_masks_type=='multiply') or (enlarge_masks!=0. and enlarge_masks_type=='add'):
        if enlarge_masks_type=='multiply':
            offset = [enlarge_masks, 1.]
        else:
            offset = [enlarge_masks]
    else:
        offset = None
    
    # Create outer mask
    print 'Outer masking ...'
    circle = Point(x0,y0).buffer(1)
    ellipse = shapely.affinity.scale(circle,sma,sma*q)
    rot_ellipse = shapely.affinity.rotate(ellipse, PA, origin='center', use_radians=False)    
    imp_masking.auto_mask(input_image, output_reg_file, sextr_setup=None, models=False, galaxy_polygon=rot_ellipse, snr=snr)
    
    
    # Do offsetting, if needed
    if offset is not None:        
        convert_segm_to_region_new.do_offseting(output_reg_file, 'tmp.reg', offset=offset, xc=None, yc=None)
        shutil.move('tmp.reg', output_reg_file)
    
    
    # Create final fits mask
    if output_mask_file is not None:
        mask_indiv.mask(input_image, output_reg_file, factor=1, output_image=None, output_mask=output_mask_file, show_running=False)

    #convert_segm_to_region_new.main('mask_tmp.fits', 'outer_mask.reg', output_mask_image='outer_mask.fits', offset=[factor_enlarge_masks,1.],xc=x0,yc=y0)
    print '\t Done!'

def convert_segm_to_boolean(mask):
    return np.ma.make_mask(mask, shrink=False)



def fill_masked_pixels(input_image, mask_image, x0, y0, max_radius, q, PA, min_radius=5., step = 1., output_filled_file='galaxy_filled.fits', outer_mask_filled_value=0., outer_mask_filled_value_std=0.):
    '''
    DESCRIPTION:
    Function to fill in the masked pixels with the average value for each ellipse which is generated (the same as in inner_mask). The values which are used here should be the same as were used for inner_mask.
    '''
    print 'Filling masked pixels ...' 
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    header = hdulist[0].header
    
    ny,nx = np.shape(data)

    hdulist_mask = pyfits.open(mask_image)
    mask = hdulist_mask[0].data


    mask = convert_segm_to_boolean(mask)

    PA = np.radians(180. - PA)
    if outer_mask_filled_value is None or outer_mask_filled_value_std is None:
        # In this case we will run ellipse to the image corners
        max_radius = np.max([math.sqrt(x0**2+(ny-y0)**2),math.sqrt((nx-x0)**2+(ny-y0)**2),math.sqrt((nx-x0)**2+y0**2),math.sqrt(x0**2+y0**2)])
    
    radii = np.arange(min_radius, max_radius, step)

    ellipse_data = np.empty((ny,nx,))
    ellipse_data[:] = np.nan    

    filled_data = np.copy(data)
    
    for radius in radii:
        #print '\tRadius:\t',radius
        # Draw ellipse in pixel coordinates:
        rr, cc = ellipse_perimeter(int(x0), int(y0), int(radius), int(radius*q), PA, shape=None)

        # Remove those pixels which are outside the image
        RR = []; CC = []
        for kk in range(len(rr)):
            if rr[kk]>=0 and cc[kk]>=0 and rr[kk]<nx and cc[kk]<ny:
                RR.append(rr[kk])
                CC.append(cc[kk])
        try:
            mean, median, std = sigma_clipped_stats(data[CC,RR], mask=mask[CC,RR], sigma=3., iters=5)        
        except:
            #print 'here'
            continue
        
        # Fill in masked pixels with a value determined by the normal distribution.
        for k in range(len(RR)):
            ellipse_data[CC[k],RR[k]]=random.normalvariate(mean, std)
    
    kernel = Gaussian2DKernel(1)
    ellipse_data = convolve(ellipse_data, kernel)

    if outer_mask_filled_value is None or outer_mask_filled_value_std is None:
        for i in range(ny):
            for k in range(nx):
                if mask[i,k] == True:
                    filled_data[i,k] = ellipse_data[i,k]

    else:
        for i in range(ny):
            for k in range(nx):
                if mask[i,k] == True:
                    filled_data[i,k] = ellipse_data[i,k]
                if np.isnan(filled_data[i,k]):
                    filled_data[i,k] = random.normalvariate(outer_mask_filled_value, outer_mask_filled_value_std)

    outHDU = pyfits.PrimaryHDU(filled_data, header=header)
    outHDU.writeto(output_filled_file, clobber=True)
    print '\t Done!'




def combine_masks(input_region_files, output_region_file):
	    with open(output_region_file, 'w') as outfile:
		outfile.write('image\n')
		for file in input_region_files:
		    with open(file) as infile:
			Lines = infile.readlines()
			for Line in Lines:
			  if 'image' not in Line:
			    outfile.write(Line)
 
    
           
            


def main(input_image, x0, y0, sma, q, PA, snr_outer=2., sigma_inner=2., enlarge_masks=1.2, enlarge_masks_type='multiply'):
    # Create outer mask
    outer_mask(input_image, x0, y0, sma, q, PA, snr=snr_outer, enlarge_masks=enlarge_masks, enlarge_masks_type=enlarge_masks_type, output_mask_file=None, output_reg_file='outer_mask.reg')

    # Create inner mask
    inner_mask(input_image, x0, y0, sma, q, PA, sigma=sigma_inner, enlarge_masks=enlarge_masks, enlarge_masks_type=enlarge_masks_type, output_mask_file='inner_mask.fits', output_reg_file='inner_mask.reg')
    
    # Combine two regions with masks 
    combine_masks(['inner_mask.reg','outer_mask.reg'], 'final_mask.reg')    

    # Create final mask 
    print 'Final masking ...'
    mask_indiv.mask(input_image, 'final_mask.reg', factor=1, output_image='galaxy_clean.fits', output_mask='final_mask.fits', show_running=False)        
    print '\t Done!'

    
    # Create image with filled pixels
    fill_masked_pixels(input_image, 'final_mask.fits', x0, y0, sma, q, PA, output_filled_file='galaxy_filled.fits', outer_mask_filled_value=None, outer_mask_filled_value_std=None) 
    

    


'''
input_image = 'NGC5529_SDSS_r.fits'
x0 = 1997.
y0 = 2010.
sma = 368.024235545 / 0.45
q = 1./4.97486278969
PA = 23.4680577545
sigma=1.2
factor_enlarge_masks=1.2
'''
'''
#PGC039799,185.08845,0.36679,82.3406007134,1.02407394974,130.150563872,,
input_image = 'PGC039799_SDSS_r.fits'
x0 = 366.
y0 = 366.
sma = 82.3406007134 / 0.45
q = 0.34
PA = 36.9
sigma = 1.2
factor_enlarge_masks = 1.2
'''

'''
#NGC3381,162.10335,34.7114,125.194613662,1.15009367979,129.087538274,,0.00454522060347,0.000207957998613,
input_image = 'frame-r-004504-2-0047.fits'
x0 = 947.
y0 = 698.
sma = 125.194613662 / 0.45
q = 1./1.15
PA = 129.087538274+90.
sigma = 2.0
factor_enlarge_masks = 1.0
'''

#NGC0628,24.174,15.78332,586.892119573,1.02094094781,-3.13354644609,,



'''
input_image = 'NGC0628_SDSS_r_crop.fits'
x0 = 1400
y0 = 1400
sma = 586.892119573 / 0.45
q = 1./1.02094094781
PA = -3.13354644609
sigma = 2.0
factor_enlarge_masks = 1.4

main(input_image, x0, y0, sma, q, PA, snr_outer=2., sigma_inner=sigma, enlarge_masks=factor_enlarge_masks, enlarge_masks_type='multiply')
'''
