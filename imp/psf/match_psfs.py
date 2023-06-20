import matplotlib.pyplot as plt
from astropy.visualization import LogStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.datasets import load_irac_psf
from photutils import CosineBellWindow, create_matching_kernel
import sys
from os.path import exists
from os import remove
from scipy.ndimage.filters import convolve
from scipy import signal
from numpy import max, sum
import pyfits
import numpy as np

from astropy import wcs
import os
import shutil

sys.path.append('/home/amosenko/MyPrograms/IMAN/IMP_NEW')
import imp_rebin

def add_wcs(input_image, output_image):
    hdulist_inp = pyfits.open(input_image)
    data_inp = hdulist_inp[0].data
    header = hdulist_inp[0].header
    
    ny,nx = np.shape(data_inp)
    xc = nx/2.
    yc = ny/2.

    header.append(('WCSAXES', 2))
    
    header.append(('CRPIX1', xc))
    header.append(('CRPIX2', yc))

    header.append(('CDELT1', -0.066667))
    header.append(('CDELT1',  0.066667))

    header.append(('CUNIT1', 'deg'))
    header.append(('CUNIT2', 'deg'))
    
    header.append(('CTYPE1', 'RA---TAN'))
    header.append(('CTYPE2', 'DEC--TAN'))

    header.append(('CRVAL1', 0.0))
    header.append(('CRVAL2', -90.0))   

    header.append(('LONPOLE', 180.0))
    header.append(('LATPOLE', 67.9122158365))    

    outHDU = pyfits.PrimaryHDU(data_inp, header=header)
    outHDU.writeto(output_image, clobber=True)
    return nx,ny



def main(input_images, input_psf, match_psf, output_images):
    # ATTENTION: input_psf and match_psf images should have the same pixel scale but may have different shape!
    # input_images is a list. If just one image should be processed, then use [input_image].
    
    
    # 1. Add wcs to PSF images
    nx_inp_psf,ny_inp_psf = add_wcs(input_psf, 'inp_psf_tmp.fits')
    nx_match_psf,ny_match_psf = add_wcs(match_psf, 'match_psf_tmp.fits')    
    #exit()
    # 2. Rebin tmp PSF images to the same shape (with minimal size)
    if nx_inp_psf>nx_match_psf:
        imp_rebin.rebin('match_psf_tmp.fits', 'inp_psf_tmp.fits', 'inp_psf_tmp_tmp.fits', preserve_bad_pixels=True)
        shutil.move('inp_psf_tmp_tmp.fits', 'inp_psf_tmp.fits')
    else:
        imp_rebin.rebin('inp_psf_tmp.fits', 'match_psf_tmp.fits', 'match_psf_tmp_tmp.fits', preserve_bad_pixels=True)
        shutil.move('match_psf_tmp_tmp.fits', 'match_psf_tmp.fits')        
    

    # 3. Read in all input files
    hdulist_inp_psf = pyfits.open('inp_psf_tmp.fits')
    data_inp_psf = hdulist_inp_psf[0].data
    ny_inp_psf,nx_inp_psf = np.shape(data_inp_psf)

    hdulist_match_psf = pyfits.open('match_psf_tmp.fits')
    data_match_psf = hdulist_match_psf[0].data
    ny_match_psf,nx_match_psf = np.shape(data_match_psf)
    
    # 4. Create kernel for matching
    window = CosineBellWindow(alpha=0.35)
    kernel = create_matching_kernel(data_inp_psf, data_match_psf, window=window)
    kernel = kernel / sum(kernel)

    for k in range(len(input_images)):
        input_image = input_images[k]
        output_image = output_images[k]
        
        # 5. Convolve input image to match match_psf
        hdulist_inp = pyfits.open(input_image)
        data_inp = hdulist_inp[0].data
        header = hdulist_inp[0].header
        
        convolvedImage = signal.fftconvolve(data_inp, kernel, mode='same')

        outHDU = pyfits.PrimaryHDU(convolvedImage, header=header)
        outHDU.writeto(output_image, clobber=True)  
    print 'Done!'
    

#main(['galaxy_100_galf.fits','sigma_100.fits'], 'Combined_PACS100_V20_100_100.fits', 'SPIRE_250_100.fits', ['galaxy_100_galf_match_250.fits','sigma_100_match.fits'])