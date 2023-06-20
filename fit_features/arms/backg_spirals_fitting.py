#!/usr/bin/python

# Import standard modules
import sys
import math
import numpy as np
from numpy import *
import os
import shutil
import subprocess
import warnings
import pyfits
import argparse
import scipy.ndimage as ndimage
warnings.filterwarnings("ignore")

path_to_azimProfile = '/home/amosenko/MEGA/MyPrograms/Ellipse_photometry' #### CHANGE HERE!
path_to_Ellipse_photometry = '/home/amosenko/MEGA/MyPrograms/Ellipse_photometry/ellipse_phorometry' #### CHANGE HERE!

def subtract(first_image, second_image, output_image):
    hdulist_first = pyfits.open(first_image)
    frame_data1 = hdulist_first[0].data
    header_data1 = hdulist_first[0].header

    hdulist_second = pyfits.open(second_image)
    frame_data2 = hdulist_second[0].data
    header_data2 = hdulist_second[0].header

    final_image = frame_data1 - frame_data2
    hdu = pyfits.PrimaryHDU(np.array(final_image,float), header_data1)
    hdu.writeto(output_image, clobber=True)
    return final_image

def create_background(input_image, spiral_image, spiral_mask_image):
    hdulist_first = pyfits.open(input_image)
    frame_data1 = hdulist_first[0].data
    header_data1 = hdulist_first[0].header

    hdulist_second = pyfits.open(spiral_image)
    frame_data2 = hdulist_second[0].data
    header_data2 = hdulist_second[0].header    

    hdulist_mask = pyfits.open(spiral_mask_image)
    frame_mask = hdulist_mask[0].data
    
    spirals_data = np.copy(frame_data1)
    (dimy,dimx) = frame_data1.shape
    for k in range(dimy):
        for i in range(dimx):
            if frame_mask[k,i]!=0.:
                spirals_data[k,i] = frame_data1[k,i] - frame_data2[k,i]
    hdu = pyfits.PrimaryHDU(spirals_data, header_data1)
    hdu.writeto('galaxy_background.fits', clobber=True)   


def main(input_image, spiral_mask_image, output_image, xcen=None, ycen=None, sma_max=None, sma_min=1., ell=0., PA=0.,method='ellipse'):      
    hdulist = pyfits.open(input_image)
    frame_data = hdulist[0].data
    (dimy,dimx) = frame_data.shape
    
    if xcen is None:
        xcen = dimx / 2.0
    
    if ycen is None:
        ycen = dimy / 2.0        

    if sma_max is None:
        sma_max = max(xcen, ycen, dimx-xcen, dimy-ycen)
        
    if method=='ellipse':
        # Create image of background (with masked spirals)
        subprocess.call('python3 %s/ellipse_photometry.py %s --x-cen %.3f --y-cen %.3f --min-aper %.3f --max-aper %.3f --aper-step 1.0 --mag-zpt 20.0 --make-model --scale 0.396' % (path_to_Ellipse_photometry,input_image, xcen, ycen, sma_min, sma_max), shell=True)
        
        # Subtract model of background from galaxy image -> spirals:
        subtract(input_image, './result/model.fits', output_image)
        #shutil.copy('./result/filled.fits', 'galaxy_background.fits')
        shutil.rmtree('./result')
        os.remove('field.cat')
        os.remove('seg.fits')
    else:
        subprocess.call('python3 %s/azimProfile.py %s --mask %s --model iso_model.fits --ell %.3f --posang %.2f --xcen %.2f --ycen %.2f --sma-max %.3f --sma-min %.3f' % (path_to_azimProfile,input_image,spiral_mask_image,ell,PA,xcen,ycen,sma_max,sma_min), shell=True)
        
        # Subtract model of background from galaxy image -> spirals:
        subtract(input_image, 'iso_model.fits', output_image)
        os.remove('azim_model.txt')
        os.remove('iso_model.fits')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract spiral arms")
    parser.add_argument("inputImage", help="Input fits image")
    parser.add_argument("inputMask", help="Input mask image with masked spirals")
    
    parser.add_argument("--outputImage", help="Input fits image with spirals left", type=str, default='spirals.fits') 
    

    parser.add_argument("--xcen", nargs='?', const=1, help="Optional: Galaxy center x",type=float,default=None)  
    parser.add_argument("--ycen", nargs='?', const=1, help="Optional: Galaxy center y",type=float,default=None)      
    parser.add_argument("--ell", nargs='?', const=1, help="Optional: Galaxy ellipticity",type=float,default=0.0)    
    parser.add_argument("--posang", nargs='?', const=1, help="Optional: Galaxy position angle in degrees. Up=0, Left=90. Default is posang=0.0",type=float,default=0.0)      
    parser.add_argument("--sma_min", nargs='?', const=1, help="Optional: Minimum radius where to mask",type=float,default=1.)  
    parser.add_argument("--sma_max", nargs='?', const=1, help="Optional: Maximum radius where to mask",type=float,default=None)  
    parser.add_argument("--method", help="Method to fit ellipses: ellipse or circle", type=str, default='ellipse') 

    args = parser.parse_args()

    input_image = args.inputImage
    output_image = args.outputImage
    input_mask = args.inputMask
    xcen = args.xcen
    ycen = args.ycen
    ell = args.ell
    posang = args.posang
    sma_min = args.sma_min
    sma_max = args.sma_max
    method = args.method

    main(input_image, input_mask, output_image, xcen=xcen, ycen=ycen, sma_max=sma_max, sma_min=sma_min, ell=ell, PA=posang,method=method)



#input_image = '/home/amosenko/CurrentWork/MOL_A/Arm_masking/PGC52017/test3/image.fits'
#spiral_mask_image = '/home/amosenko/CurrentWork/MOL_A/Arm_masking/PGC52017/test3/mask_spirals.fits'
#main(input_image, spiral_mask_image)
#main(input_image, input_mask, xcen=143., ycen=137., sma_max=None, sma_min=1., ell=0., PA=0.,method='circle')


























