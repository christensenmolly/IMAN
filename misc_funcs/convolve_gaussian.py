import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astropy.convolution import Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve
import argparse


def main(input_image, FWHM, output_image):
        hdulist = fits.open(input_image)
        img = hdulist[0].data
        header = hdulist[0].header
        
        sigma_gaus = FWHM/2.35482

        # We smooth with a Gaussian kernel
        kernel = Gaussian2DKernel(x_stddev=sigma_gaus)
        astropy_conv = convolve(img, kernel)
        
        hdu = fits.PrimaryHDU(astropy_conv, header)
        if output_image is None:
            output_image = input_image.split('.fits')[0] + '_conv.fits'
        
        hdu.writeto(output_image, overwrite=True)        



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arithmetical operations")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("FWHM", help="FWHM (in pixels)", type=float, default=1.)
    parser.add_argument("--output_image", help="Output fits image", type=str, default=None)  
    args = parser.parse_args()

    input_image = args.input_image
    FWHM = args.FWHM
    output_image = args.output_image

    main(input_image, FWHM, output_image)
