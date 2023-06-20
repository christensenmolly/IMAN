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

sys.path.append(os.path.join(IMAN_DIR, 'misc'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/sky_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/masking'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/ds9_regions'))

import photometric_calibration
import remove_cosmic_rays
import rebin_image
import determine_sky
import auto_masking
import convert_reg_to_mask
import backg_subtr
import crop_image
import subprocess
import spline_bkg_subtr
import arithm_operations
import convert_fk5_to_image

def convert_segm_to_boolean(mask):
    return np.ma.make_mask(mask, shrink=False)

def mask_sources(data):
    sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
    threshold = detect_threshold(data, nsigma=2.0, sigma_clip=sigma_clip)
    
    kernel = make_2dgaussian_kernel(3.0, size=3)  # FWHM = 3.
    convolved_data = convolve(data, kernel)
    
    segment_img = detect_sources(convolved_data, threshold, npixels=3)

    mask = segment_img.make_source_mask()
    return mask


def mask_sources_segm(file):

    hdulist = fits.open(file, do_not_scale_image_data=True, mode='update')
    header = hdulist[0].header
    segm = hdulist[0].data

    
    
    

    hdulist1 = fits.open('mask.fits', do_not_scale_image_data=True, mode='update')
    mask = hdulist1[0].data    
    
    mask_astropy = convert_segm_to_boolean(segm+mask)
    
    return mask_astropy



def main(input_image, output_image, sigma=3.0, box_size=50, filter_size=3, polynomial_degree=1):

    auto_masking.main(input_image, output_region_file='general_mask.reg', snr=2., min_pix=5, region_type='polygon', sextr_setup=None, galaxy_ellipse=None, offset_size=1.5, offset_pix=0., verbosity=True)

        

    ds9Proc = subprocess.Popen(["ds9", input_image,
                                            "-regions", 'general_mask.reg',
                                            "-scale", "histequ"])
    ds9Proc.wait() 
    
    convert_reg_to_mask.mask(input_image, 'general_mask.reg', output_image=None, output_mask='mask.fits', mask_value=1, mask_DN=None, verbosity=True)
    
    if os.path.exists('unified_mask.reg'):
        use_mask = str(input('\n Would you like to use the unified mask? (yes):') or 'yes')
        if use_mask=='yes':
            print("Using the unified_mask.reg ...")
            convert_fk5_to_image.main('unified_mask.reg', input_image, 'mask_tmp.reg')
            
            convert_reg_to_mask.mask(input_image, 'mask_tmp.reg', output_image=None, output_mask='mask_tmp.fits', mask_value=1, mask_DN=None, verbosity=True)        
            
            arithm_operations.main('mask_tmp.fits', 'mask.fits', 'add', 'mask.fits')
            if os.path.exists('mask_tmp.fits'):
                os.remove('mask_tmp.fits')
            if os.path.exists('mask_tmp.reg'):    
                os.remove('mask_tmp.reg')        
        else:
            print("Using the general_mask.reg ...")
        
    arithm_operations.main('mask.fits', 'segm.fits', 'add', 'final_mask.fits')
    

    determine_sky.sky_subtraction(input_image, 'final_mask.fits', polynomial_degree=polynomial_degree, output_image='sky_subtr.fits', output_sky=None, hdu_inp=0, sampling_factor=1., sigma_smooth=3, verbosity=True, sky_value='mode')



    ds9Proc = subprocess.Popen(["ds9", 'sky_subtr.fits',
                                        "-scale", "histequ"])
    ds9Proc.wait() 


    spline_bkg_subtr.main('sky_subtr.fits', 'final_mask.fits', output_image, sigma=sigma, box_size=box_size, filter_size=filter_size)
    
#input_image = '/media/mosav/MY_DATA_DRIVE/APO_observations_red/UT211009/reduced/UGC10043/UGC10043.g.0003_reduced_crop_rayremoved.fits'
#main(input_image)
