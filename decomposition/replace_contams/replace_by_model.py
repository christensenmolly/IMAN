import numpy as np
import matplotlib.pyplot as plt

import astropy.io.fits as pyfits
from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve
import os
import sys


LOCAL_DIR = "/decomposition/replace_contams"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/masking/mto-lib'))

import mto_func



def main(input_image, xc, yc, mask_image=None, backg_mean=0.0, backg_std = -1, verbosity=True, sigma=5):

    galaxyHDU = pyfits.open(input_image)
    galaxyData = galaxyHDU[0].data
    galaxyHeader = galaxyHDU[0].header
    ny,nx = np.shape(galaxyData)    
    if mask_image is None:
        # mto
        mto_func.main(input_image, bg_mean=backg_mean, bg_variance=backg_std, verbosity=verbosity)
        mask_image = 'segm_mto.fits'
    
    segmHDU = pyfits.open(mask_image)
    segmData = segmHDU[0].data
    
    I_galaxy = segmData[int(yc), int(xc)]
    
    nandata = np.copy(galaxyData)
    
    for k in range(ny):
        for i in range(nx):
            if segmData[k,i]!=I_galaxy and segmData[k,i]>0:
                nandata[k,i] = float('nan')
    
    kernel = Gaussian2DKernel(x_stddev=sigma)
    nandata_conv = convolve(nandata, kernel)
    

    for k in range(ny):
        for i in range(nx):
            if segmData[k,i]!=I_galaxy and segmData[k,i]>0:
                galaxyData[k,i] = nandata_conv[k,i] 

    hdu = pyfits.PrimaryHDU(galaxyData, galaxyHeader)
    hdu.writeto('replaced.fits', clobber=True)      
    
#main('galaxy_rot_crop.fits', 638, 432, backg_mean=0.0, backg_std = 89.9, verbosity=True) # NGC2683
main('galaxy_rot_crop.fits', 1653, 1252, backg_mean=0.0, backg_std = 146.54, verbosity=True, sigma=25., mask_image='mask_rot_crop.fits')  # NGC4594

    
#main('NGC2683.phot.1_nonan.fits', 1074, 665, backg_mean=0.0, backg_std = float(4.64599e-05), verbosity=True)
#main('NGC2683.phot.2_nonan.fits', 895, 1175, backg_mean=0.0, backg_std = float(0.00013), verbosity=True)
#main('combined.fits', 1045, 1182, backg_mean=0.0, backg_std = float(0.00013), verbosity=True)