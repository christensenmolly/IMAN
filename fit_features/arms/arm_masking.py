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

# Additional module
import run_SExtractor

path_to_sextractor_dir = '/home/amosenko/MEGA/MyPrograms/IMAN/IMP_NEW/Sextractor' ##########CHANGE HERE!!!
path_to_azimProfile = '/home/amosenko/MEGA/MyPrograms/Ellipse_photometry'  ##########CHANGE HERE!!!

def combine_masks(mask_images):
    final_mask = np.copy(mask_images[0])
    
    for k in range(len(mask_images)):
        mask = mask_images[k]
        hdulist_data = pyfits.open(mask)
        if k==0:
            final_mask = hdulist_data[0].data
        else:
            final_mask = final_mask + hdulist_data[0].data
    return final_mask
        


def main(input_image, input_mask, output_mask, xcen, ycen, ell, posang, sma_min, sma_max, resid_lower_limit=0.25, resid_upper_limit=1.0, sigma_aver=3.0):

    
    for kk in range(2):
        hdulist_data = pyfits.open(input_image)
        frame_data = hdulist_data[0].data
        header_data = hdulist_data[0].header
        (dimy,dimx) = frame_data.shape
        
        if xcen is None:
            xcen = dimx/2.
        if ycen is None:
            ycen = dimy/2.
        if sma_min is None:
            sma_min = 1.

        if sma_max is None:
            sma_max = max(xcen, ycen, dimx-xcen, dimy-ycen)
 
        if input_mask is not None:
            subprocess.call('python3 %s/azimProfile.py %s --mask %s --model iso_model.fits --ell %.3f --posang %.2f --xcen %.2f --ycen %.2f --sma-max %.3f --sma-min %.3f' % (path_to_azimProfile,input_image,input_mask,ell,posang,xcen,ycen,sma_max,sma_min), shell=True)
        else:
            subprocess.call('python3 %s/azimProfile.py %s --model iso_model.fits --ell %.3f --posang %.2f --xcen %.2f --ycen %.2f --sma-max %s --sma-min %.3f' % (path_to_azimProfile,input_image,ell,posang,xcen,ycen,sma_max,sma_min), shell=True)            

        hdulist_model = pyfits.open('iso_model.fits')
        frame_model = hdulist_model[0].data
        

        resid_data = np.copy(frame_model)
        mask_data = np.zeros(shape=(dimy, dimx))


        for k in range(dimy):
             for i in range(dimx):
                resid_data[k,i] = (frame_data[k,i] - frame_model[k,i]) / frame_data[k,i]
                if resid_data[k,i]>resid_lower_limit and resid_data[k,i]<resid_upper_limit:
                    mask_data[k,i] = 1.
                else:
                    mask_data[k,i] = 0.

        hdu = pyfits.PrimaryHDU(resid_data, header_data)
        hdu.writeto('resid.fits',clobber=True)      
        #mask_data = ndimage.gaussian_filter(mask_data, sigma=sigma_aver, order=0)
        #mask_data[mask_data<=0.5]=0
        #mask_data[mask_data>0.5]=1
        hdu = pyfits.PrimaryHDU(np.array(mask_data,float), header_data)
        hdu.writeto(output_mask, clobber=True)

        run_SExtractor.call_SE(output_mask, snr=2., min_pix=10, sextr_dir=path_to_sextractor_dir, sextr_setup='cold.sex', sextr_param='default.param', output_cat='field.cat', checkimage_type='SEGMENTATION',checkimage_name='segm.fits', sextr_add_string=None)
        hdulist_mask = pyfits.open('segm.fits')
        frame_mask = hdulist_mask[0].data
        mask_data = frame_mask#ndimage.gaussian_filter(frame_mask, sigma=sigma_aver, order=0)
        #mask_data[mask_data<=0.5]=0
        #mask_data[mask_data>0.5]=1
        hdu = pyfits.PrimaryHDU(mask_data, header_data)
        hdu.writeto(output_mask,clobber=True)

	shutil.copy(output_mask, '%s_%i.fits' % (output_mask.split('.fits')[0],kk+1))
	shutil.copy('resid.fits', 'resid_%i.fits' % (kk+1))
	if k==0:
            if input_mask is None:
                input_mask = output_mask
            else:
                combined_mask = combine_masks(['%s_1.fits' % (output_mask.split('.fits')[0]), input_mask])
                hdu = pyfits.PrimaryHDU(combined_mask, header_data)
                hdu.writeto(output_mask, clobber=True)            

    combined_mask = combine_masks(['%s_1.fits' % (output_mask.split('.fits')[0]), '%s_2.fits' % (output_mask.split('.fits')[0])])
    hdu = pyfits.PrimaryHDU(combined_mask, header_data)
    hdu.writeto(output_mask, clobber=True)      
    
    
    mask_data = ndimage.gaussian_filter(np.array(combined_mask,float), sigma=sigma_aver, order=0)
    mask_data[mask_data<0.5]=0
    mask_data[mask_data>=0.5]=1
    hdu = pyfits.PrimaryHDU(mask_data, header_data)
    hdu.writeto('%s_conv.fits' % (output_mask.split('.fits')[0]), clobber=True)    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create mask of spiral arms")
    parser.add_argument("inputImage", help="Input fits image")
    parser.add_argument("--outputMask", help="Output mask", type=str, default='mask_spirals.fits') 
    
    parser.add_argument("--inputMask", nargs='?', const=1, help="Optional: Input mask image",type=str,default=None)

    parser.add_argument("--xcen", nargs='?', const=1, help="Optional: Galaxy center x",type=float,default=None)  
    parser.add_argument("--ycen", nargs='?', const=1, help="Optional: Galaxy center y",type=float,default=None)      
    parser.add_argument("--ell", nargs='?', const=1, help="Optional: Galaxy ellipticity",type=float,default=0.0)    
    parser.add_argument("--posang", nargs='?', const=1, help="Optional: Galaxy position angle in degrees. Up=0, Left=90. Default is posang=0.0",type=float,default=0.0)      
    parser.add_argument("--sma_min", nargs='?', const=1, help="Optional: Minimum radius where to mask",type=float,default=None)  
    parser.add_argument("--sma_max", nargs='?', const=1, help="Optional: Maximum radius where to mask",type=float,default=None)  


    args = parser.parse_args()

    input_image = args.inputImage
    output_mask = args.outputMask
    input_mask = args.inputMask
    xcen = args.xcen
    ycen = args.ycen
    ell = args.ell
    posang = args.posang
    sma_min = args.sma_min
    sma_max = args.sma_max


    main(input_image, input_mask, output_mask, xcen, ycen, ell, posang, sma_min, sma_max)


    
'''    
input_image = 'image.fits'
#ellipse_pars = [143.,137.,150.,0.,0.] #52017
#ellipse_pars = [389.,279.,284.,0.,0.] #33371
ellipse_pars = [169.,162.,180.,0.,0.] #33720
main(input_image, ellipse_pars, mask_image=None)

hdulist_mask = pyfits.open('mask_spirals.fits')
frame_mask = hdulist_mask[0].data
header_mask = hdulist_mask[0].header
mask_data = ndimage.gaussian_filter(np.array(frame_mask,float), sigma=3., order=0)
mask_data[mask_data<0.5]=0
mask_data[mask_data>=0.5]=1
hdu = pyfits.PrimaryHDU(mask_data, header_mask)
hdu.writeto('mask_spirals_conv.fits',clobber=True)
'''