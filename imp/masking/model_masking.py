#!/usr/bin/python
# WARNING: SHOULD BE DONE! Began on Jan. 4, 2019
# Import standard modules
import sys
import math
import numpy as np
from numpy import *
import os
import shutil
import subprocess
import warnings
from astropy.io import fits as pyfits
import argparse
import scipy.ndimage as ndimage
warnings.filterwarnings("ignore")


LOCAL_DIR = "/imp/masking"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'sextractor'))

import run_SExtractor


def main(input_image, model_image, output_mask, mask_image=None, resid_lower_limit=0.25, resid_upper_limit=1.0, sigma_aver=3.0, hdu_model=0):      
        hdulist_image = pyfits.open(input_image)
        frame_image = hdulist_image[0].data

        hdulist_model = pyfits.open(model_image)
        frame_model = hdulist_model[hdu_model].data
        
        frame_resid = (frame_image-frame_model)/frame_image
        
        
        header_resid = hdulist_image[0].header


        (dimy,dimx) = frame_resid.shape
        mask_data = np.zeros(shape=(dimy, dimx))

        for k in range(dimy):
             for i in range(dimx):
                if frame_resid[k,i]>resid_lower_limit and frame_resid[k,i]<resid_upper_limit:
                    mask_data[k,i] = 1.
                else:
                    mask_data[k,i] = 0.
        mask_data = ndimage.gaussian_filter(np.array(mask_data,float), sigma=sigma_aver, order=0)
        mask_data[mask_data<=0.5]=0
        mask_data[mask_data>0.5]=1

        hdu = pyfits.PrimaryHDU(mask_data, header_resid)
        hdu.writeto(output_mask,clobber=True)

        run_SExtractor.call_SE(output_mask, snr=2., min_pix=10, sextr_dir=None, sextr_setup='cold.sex', sextr_param='default.param', output_cat='field.cat', checkimage_type='SEGMENTATION',checkimage_name='segm.fits', sextr_add_string=None)
        hdulist_mask = pyfits.open('segm.fits')

        frame_mask = hdulist_mask[0].data
        mask_data = ndimage.gaussian_filter(np.array(frame_mask,float), sigma=sigma_aver, order=0)
        mask_data[mask_data<=0.5]=0
        mask_data[mask_data>0.5]=1
        hdu = pyfits.PrimaryHDU(mask_data, header_resid)
        hdu.writeto('%s_conv.fits' % (output_mask.split('.fits')[0]),clobber=True)
        os.remove('field.cat')
        os.remove('segm.fits')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create mask of spiral arms")
    parser.add_argument("inputImage", help="Input image of a galaxy")
    parser.add_argument("modelImage", help="Model image of a galaxy")
    parser.add_argument("--outputMask", help="Output mask", type=str, default='mask_spirals.fits') 
    parser.add_argument("--inputMask", nargs='?', const=1, help="Optional: Input mask image",type=str,default=None)
    parser.add_argument("--hdu_model", nargs='?', const=1, help="Optional: Layer of model image",type=int,default=0)

    args = parser.parse_args()
    input_image = args.inputImage
    model_image = args.modelImage
    output_mask = args.outputMask
    input_mask = args.inputMask
    hdu_model = args.hdu_model
    main(input_image, model_image, output_mask, mask_image=input_mask, hdu_model = hdu_model)
  
