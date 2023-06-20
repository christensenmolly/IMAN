#!/usr/bin/python
# DESCRIPTION:
# Script to rebin image to a reference image.
# factor = pixelscale_new/pixelscale_old

# Import standard modules
import os
import os.path
import math
import numpy as np
import shutil
import astropy.io.fits as pyfits
from astropy import wcs
import glob
import argparse
from astropy import log
import warnings
from scipy import ndimage
warnings.filterwarnings("ignore")
import sys

LOCAL_DIR = "/imp/rebinning"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))

from hcongrid import hcongrid
import rebin_image

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Sampling")
    parser.add_argument("input_image", help="Input image")    
    parser.add_argument("factor", help="Input factor (<1 - upsample, >1 downsample)",type=float,default=1.)
    parser.add_argument("--o", help="Optional: Output image",type=str,default=None)
    parser.add_argument("--norm", help="Norm images",
                        default=False, action="store_true")   
    args = parser.parse_args()


    input_image = args.input_image
    output_image = args.o
    factor = args.factor
    norm = args.norm
    
    rebin_image.downsample(input_image, factor, output_image=output_image, set_wcs=True, print_mes=True, norm=norm)
