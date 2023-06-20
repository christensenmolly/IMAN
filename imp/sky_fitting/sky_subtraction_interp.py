# img_interp.py
import os
import sys
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from PIL import Image
from astropy.io import fits as pyfits

from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
import random
import scipy.ndimage as ndimage
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
import argparse

LOCAL_DIR = "/imp/sky_fitting"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))
sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))
import rebin_image
import arithm_operations

def remove_files(files):
    for file in files:
        if os.path.exists(file):
            os.remove(file)


def main(input_image, mask_image, output_image='sky_sub_interp.fits', sky_image=None, sigma_smooth=30., sampling_factor=5.):
    hdulist0 = pyfits.open(input_image)
    data0 = hdulist0[0].data
    header0 = hdulist0[0].header 
    ny0,nx0 = np.shape(data0)

    hdulist_mask = pyfits.open(mask_image)
    mask0 = hdulist_mask[0].data
    input_image0 = input_image

    
    if sampling_factor is not None:
        rebin_image.downsample(input_image, 5., output_image='tmp_rebin.fits', set_wcs=True, print_mes=True, norm=False)
        rebin_image.downsample(mask_image, 5., output_image='tmp_mask_rebin.fits', set_wcs=True, print_mes=True, norm=False, no_interp=True)
        input_image = 'tmp_rebin.fits'
        mask_image = 'tmp_mask_rebin.fits'
    
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    header = hdulist[0].header
    ny,nx = np.shape(data)

    hdulist1 = pyfits.open(mask_image)
    mask = hdulist1[0].data    
    
    data[mask > 0.] = np.nan
    
    kernel = Gaussian2DKernel(x_stddev=sigma_smooth)
    
    astropy_conv = convolve(data, kernel)
    
    outHDU = pyfits.PrimaryHDU(astropy_conv, header)
    outHDU.writeto('interp_tmp_rebin.fits', overwrite=True)    
    
    if nx!=nx0:
        #rebin_image.downsample('interp_tmp_rebin.fits', float(nx)/float(nx0), output_image='interp_tmp.fits', set_wcs=True, print_mes=True, norm=False)
        rebin_image.rebin(input_image0, 'interp_tmp_rebin.fits', output_image='interp_tmp.fits', hdu_ref=0, hdu_inp=0, preserve_bad_pixels=True)
    else:
        os.rename('interp_tmp_rebin.fits', 'interp_tmp.fits')
    
    arithm_operations.main(input_image0, 'interp_tmp.fits', 'sub', output_image)
    if sky_image is not None:
        os.rename('interp_tmp.fits', sky_image)
    
    
    remove_files(['tmp_rebin.fits','tmp_mask_rebin.fits','interp_tmp_rebin.fits','interp_tmp.fits'])
    print('Done!')
    


#main('35840_r_trim.fits', 'general_mask_r.fits', sky_image='sky_interp.fits', sigma_smooth=30., sampling_factor=5.)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky background interpolation. The input image should have WCS!")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("inputMask", help="Input mask image",type=str) 
    parser.add_argument("--output_image", help="Optional: Output image", type=str, default='sky_sub_interp.fits') 
    parser.add_argument("--output_sky", help="Optional: Output sky", type=str, default=None)
    parser.add_argument("--sigma", help="Optional: Sigma for Gaussian smoothing", type=float, default=None) 
    parser.add_argument("--factor", help="Optional: Sampling factor (>=1)", type=float, default=None) 
    
    args = parser.parse_args()

    input_image = args.inputImage
    mask_image = args.inputMask
    output_image = args.output_image
    output_sky = args.output_sky
    sigma_smooth = args.sigma
    factor = args.factor

    main(input_image, mask_image, output_image=output_image, sky_image=output_sky, sigma_smooth=sigma_smooth, sampling_factor=factor)