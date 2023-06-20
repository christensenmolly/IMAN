#! /usr/bin/env python
import numpy as np
import gzip
import shutil
from joblib import Parallel, delayed
import astropy.io.fits as pyfits
import sys
import os
import shutil
import time
import subprocess
import glob
import math
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from scipy import special
import scipy.ndimage as ndimage
import argparse
import collections
from skimage import data, filters, measure, morphology
from skimage.measure import label, regionprops, regionprops_table
from PIL import Image


def clean_from_small_objects(img, min_size=5, area_threshold=10, connectivity=1):
    from skimage import data, filters, measure, morphology
    from skimage.measure import label, regionprops, regionprops_table
    mask = np.copy(img)
    mask[img<0.] = 0.
    mask[img>0.] = 1
    mask = img.astype(bool)
    
    mask = morphology.remove_small_objects(mask, min_size=min_size, connectivity=connectivity)
    mask = morphology.remove_small_holes(mask, area_threshold=area_threshold, connectivity=connectivity)
    mask = mask.astype(int)
    
    img[mask==0] = 0.
    
    #hdu = pyfits.PrimaryHDU(img)
    #hdu.writeto('test.fits', overwrite=True)
    return img,mask


def main(input_image, I_min, filter_size_img=None, filter_size_mask=10, output_mask = None, output_image = None):
    # Background should be subtracted!!!
    # Remove pixels with I<I_min (make them 0)
    
    hdulist = pyfits.open(input_image)
    img = hdulist[0].data
    ySize, xSize = np.shape(img)
    
    #cleaned_img = np.zeros(shape=(ySize,xSize))
    cleaned_img = np.copy(img)
    
    cleaned_img[cleaned_img<I_min] = 0.
    cleaned_img,mask = clean_from_small_objects(cleaned_img)
    
    if filter_size_img is not None:
        cleaned_img = ndimage.uniform_filter(cleaned_img, size=filter_size_img)

    if filter_size_mask is not None:
        mask = ndimage.uniform_filter(mask, size=filter_size_mask)
    
    cleaned_mask = measure.label(mask)
    
    if output_mask is not None:
        hdu = pyfits.PrimaryHDU(cleaned_mask)
        hdu.writeto(output_mask, overwrite=True)

    if output_image is not None:
        hdu = pyfits.PrimaryHDU(cleaned_img)
        hdu.writeto(output_image, overwrite=True)

    return cleaned_img, cleaned_mask
