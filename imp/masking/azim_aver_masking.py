#!/usr/bin/python
'''
    DESCRIPTION:
    To create an inner mask we draw ellipses with increasing radius from min_radius
    to max_radius, with some step. We collect all pixels which intersect these ellipse pixels
    and use sigma clipping to mask pixels which have intensities larger than mean+sigmas*std.
    After that we detect sources using photutils to create a segmentaion map from the created mask,
    by removing masked pixels which are not in the group (npixels).
    Then we smooth data with a kernel and again detect sources which are grouped.
    Finally, we make the masked area larger (proportionally, by some factor), with a final region file
    which can be revisited.
    MINIMAL USAGE: python azim_aver_masking.py [input_image]
'''
# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import math
import sys

from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import os
from astropy.stats import sigma_clipped_stats
from skimage.draw import ellipse_perimeter
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from astropy.stats import gaussian_fwhm_to_sigma
from photutils import detect_sources
import shapely
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


import convert_segm_to_region


def read_region(reg_file):
    f = open(reg_file, "r")    
    for line in f:
        if "ellipse" in line:
            params = line.split(",")
            cen = [float(params[0].split('(')[1]),float(params[1])]
            ellA = float(params[2])
            ellB = float(params[3])
            ellPA = float(params[4].split(')')[0])
            if ellA < ellB:
                ellA, ellB = ellB, ellA
                ellPA += 90    
            break
    f.close()
    return cen[0],cen[1],ellA,ellB,ellPA

def main(input_image, input_mask=None, x0=None, y0=None, ell=0., PA=0., max_radius=None, sigmas=3., offset_size=1., offset_pix=0., min_radius=5., step = 1., FWHM=3., npixels=5, output_mask_file='inner_mask.fits', output_reg_file='inner_mask.reg', galaxy_region=None):
    print('Azimuthally averaged masking ...')

    if galaxy_region is not None:
        x0,y0,sma,smb,PA = read_region(galaxy_region)
        ell = 1. - smb/sma

    q = 1. - ell
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    header = hdulist[0].header

    if input_mask is not None:
            hdulist_mask = pyfits.open(input_mask)
            input_mask = hdulist_mask[0].data
    else:
        input_mask = np.zeros(np.shape(data))
        
    ny,nx = np.shape(data)
    if x0 is None:
        x0 = nx/2.
    if y0 is None:
        y0 = ny/2.
    PA = PA - 90.
    PA = np.radians(180. - PA) # It should be orientation of the Major axis in clockwise direction (in radians)
    
    if max_radius is None:
        max_radius = 10000
    
    radii = np.arange(min_radius, max_radius, step) # Create a set of radii, use a smaller step (at least 1 pix) 
    # The minimual radius should not be small, otherwise the procedure may mask the centre.
    
    mask = np.zeros(np.shape(data))
    
    for radius in radii:
        # Draw ellipse in pixel coordinates (using Bresenham's approach):
        rr, cc = ellipse_perimeter(int(x0), int(y0), int(radius), int(radius*q), PA, shape=None)
        
        # Remove those pixels which are outside the image
        RR = []; CC = []; intensity = []
        for kk in range(len(rr)):
            if rr[kk]>=0 and cc[kk]>=0 and rr[kk]<nx and cc[kk]<ny:
                if input_mask[cc[kk],rr[kk]]==0.:
                    intensity.append(data[cc[kk],rr[kk]])
                    RR.append(rr[kk])
                    CC.append(cc[kk])
        if RR == [] and radius>1:
            break

        try:
            mean, median, std = sigma_clipped_stats(intensity, sigma=3., iters=5)  #sigma_clipped_stats(data[CC,RR], sigma=3., iters=5)        
            
            for k in range(len(RR)):
                if data[CC[k],RR[k]]>mean + sigmas*std:   #sigma*mean:
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
    outHDU.writeto(output_mask_file, overwrite=True)       

    # Convert segmentation mask to region file with some extension of the masked areas, defined by factor_enlarge_masks (this is important to mask the outskirts of the masked objects as well). This mask then can be changed by hand.
    convert_segm_to_region.main(output_mask_file, output_reg_file, output_mask_image=None, offset_size=offset_size, offset_pix=offset_pix,xc=x0,yc=y0)
    print('Done!')
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky background estimation")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("--galaxy_region", help="Galaxy region", type=str, default=None)     
    
    parser.add_argument("--outputRegion", help="Output region file with mask",type=str,default='inner_mask.reg') 
    parser.add_argument("--outputMask", help="Output mask file",type=str,default='inner_mask.fits') 
    parser.add_argument("--inputMask", help="Input mask image (zero pixels are good).",
                        type=str, default=None)
    parser.add_argument("--xcen", help="x-coordinate of the object center (image center by default).",
                        type=float, default=None)
    parser.add_argument("--ycen", help="y-coordinate of the object center (image center by default).",
                        type=float, default=None)
    parser.add_argument("--ell", help="Ellipticity. Default is ell=0.0", default=0.0, type=float)
    parser.add_argument("--posang", help="Position angle of the ellipse in degrees. Up=90, Right=0, counterclockwise. Default is posang=0.0",
                        default=0.0, type=float) # Now in degrees!
    parser.add_argument("--sma_max", help="Maximal major axis of ellipse. Default: fitted to the image size.",
                        default=None, type=float)
    parser.add_argument("--sma_min", help="Minimal major axis of ellipse. Default: 1.",
                        default=5., type=float)
    parser.add_argument("--step", help="Linear step in semi-major axis length between successive ellipses.",
                        default=1., type=float)
    parser.add_argument("--fwhm", help="FWHM for smoothing the mask.",
                        default=1., type=float)
    parser.add_argument("--sigmas", nargs='?', const=1, help="Number of sigmas above the mean which should be masked",type=float,default=3.)
    parser.add_argument("--min_pix", nargs='?', const=1, help="Number of joint pixels",type=int,default=5)
    parser.add_argument("--offset_pix", nargs='?', const=1, help="Offset to make the regions larger or smaller, in pixels. Default 0..",type=float, default=0.) 
    parser.add_argument("--offset_size", nargs='?', const=1, help="Offset to make the regions larger or smaller, in units of the region size (multiplication). Default 1.5",type=float, default=1.5) 
    
    args = parser.parse_args()

    input_image = args.inputImage
    galaxy_region = args.galaxy_region
    
    output_region = args.outputRegion
    outputMask = args.outputMask
    input_mask = args.inputMask
    xcen = args.xcen
    ycen = args.ycen
    ell = args.ell
    posang = args.posang
    sma_min = args.sma_min
    sma_max = args.sma_max
    step = args.step
    
    sigmas = args.sigmas
    min_pix = args.min_pix
    offset_pix = args.offset_pix
    offset_size = args.offset_size
    fwhm = args.fwhm

    main(input_image, input_mask, xcen, ycen, ell, posang, max_radius=sma_max, sigmas=sigmas, offset_size=offset_size, offset_pix=0., min_radius=sma_min, step = step, FWHM=fwhm, npixels=min_pix, output_mask_file=outputMask, output_reg_file=output_region, galaxy_region=galaxy_region)



