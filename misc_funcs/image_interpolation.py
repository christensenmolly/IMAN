import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans








def main(input_image, mask_image, output_image):
    # Load the data from data.astropy.org
    #filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')

    hdu = fits.open(input_image)[0]
    img = hdu.data

    hdu1 = fits.open(mask_image)[0]
    mask = hdu1.data    
    
    
    # This example is intended to demonstrate how astropy.convolve and
    # scipy.convolve handle missing data, so we start by setting the brightest
    # pixels to NaN to simulate a "saturated" data set
    img[mask !=0] = np.nan

    # We smooth with a Gaussian kernel with x_stddev=1 (and y_stddev=1)
    # It is a 9x9 array
    kernel = Gaussian2DKernel(x_stddev=1, y_stddev=1)

    # create a "fixed" image with NaNs replaced by interpolated values
    fixed_image = interpolate_replace_nans(img, kernel)
    
    maskedHDU = fits.PrimaryHDU(data=fixed_image)
    maskedHDU.writeto(output_image, clobber=True)    
    

main('galaxy.fits', 'mask.fits', 'galaxy_filled.fits')