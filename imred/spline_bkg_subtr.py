
from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from astropy.io import fits
from astropy.stats import sigma_clipped_stats, SigmaClip
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint
from photutils.segmentation import (detect_sources,
                                    make_2dgaussian_kernel)
import numpy as np
from astropy.convolution import convolve
import sys
import os
LOCAL_DIR = "/imred"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/sky_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/masking'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))

import photometric_calibration
import remove_cosmic_rays
import rebin_image
import determine_sky
import auto_masking
import convert_reg_to_mask
import backg_subtr
import crop_image
import subprocess

import platform
opsyst = platform.system()

def convert_segm_to_boolean(mask):
    return np.ma.make_mask(mask, shrink=False)

def mask_sources_segm(file):

    hdulist = fits.open(file, do_not_scale_image_data=True, mode='update')
    header = hdulist[0].header
    mask = hdulist[0].data
   
    
    mask_astropy = convert_segm_to_boolean(mask)
    
    return mask_astropy



def main(input_image, input_mask, output_image=None, sigma=3.0, box_size=50, filter_size=3):
    hdulist = fits.open(input_image, do_not_scale_image_data=True, mode='update')
    header = hdulist[0].header
    data = hdulist[0].data
    

    mask = mask_sources_segm(input_mask)


    sigma_clip = SigmaClip(sigma=sigma)#, maxiters=10)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (box_size, box_size), filter_size=(filter_size, filter_size),
                    sigma_clip=sigma_clip, bkg_estimator=bkg_estimator, mask=mask)
    

    outHDU = fits.PrimaryHDU(bkg.background, header=header)
    outHDU.writeto('bkg.fits', overwrite=True)  

    outHDU = fits.PrimaryHDU(data-bkg.background, header=header)
    
    if output_image is None:
        output_image = input_image.split('/')[-1].split('.fits')[0]+'_skysub.fits'
        
    outHDU.writeto(output_image, overwrite=True)
    
    if opsyst=='Linux':
        ds9Proc = subprocess.Popen(["ds9", output_image,
                                            "-scale", "histequ"])
        ds9Proc.wait() 
    elif opsyst=='Darwin':
        #subprocess.call(["open","-W","-n","-a","/Applications/SAOImageDS9.app",output_image,"-scale", "histequ"])
        ds9Proc = subprocess.Popen(["/Applications/SAOImageDS9.app/Contents/MacOS/ds9", output_image,
                                            "-scale", "histequ"])
        ds9Proc.wait()
