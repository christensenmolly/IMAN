#! /usr/bin/env python
import numpy as np
import gzip
import shutil
from joblib import Parallel, delayed
import astropy.io.fits as pyfits
import sys
import os
import shutil
import time
import subprocess
import glob
import math
import matplotlib.pyplot as plt
import glob
import pickle
import collections
from joblib import Parallel, delayed
import tarfile
import argparse

LOCAL_DIR = "/stacking"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/masking'))
sys.path.append(os.path.join(IMAN_DIR, 'plotting/2dprofile'))
#sys.path.append('/home/amosenko/MyGit/IMAN/misc_funcs')
#sys.path.append('/home/amosenko/MyGit/IMAN/imp/masking')
#sys.path.append('/home/amosenko/MyGit/IMAN/plotting/2dprofile')


FNULL = open(os.devnull, 'w')


#import sdss_downloader
import imcombine
import convolve
import plot_smoothed_image

def swarp_name(verbosity):
        # Check what name has SWarp package on this system
        rCode = subprocess.call("which swarp >/dev/null", shell=True)
        if rCode == 0:
            swarpName = "swarp"
        else:
            rCode = subprocess.call("which SWarp >/dev/null", shell=True)
            if rCode == 0:
                swarpName = "SWarp"
            else:
                if verbosity: print("\033[31m Error: SWarp was not found on your system.\033[0m")
                if verbosity: print("\033[31m The command has to be either 'swarp' or 'SWarp'\033[0m")
                if verbosity: print("\033[31m Intall SWarp package or try to run this script without -s option.\033[0m")
                exit(1)
        return swarpName

def combine_images(input_images, output_combined_image):
    # WARNING: Works worse than imcombine!
    hdulist1 = pyfits.open(input_images[0])
    frame1 = hdulist1[0].data
    header = hdulist1[0].header
    
    final_frame = np.copy(frame1)
    
    for file in input_images[1:]:
        hdulist = pyfits.open(file)
        frame = hdulist1[0].data
        final_frame = final_frame + frame
    

    hdu = pyfits.PrimaryHDU(final_frame/float(len(input_images)), header)
    hdu.writeto(output_combined_image, overwrite=True)    
        


def combine_images_swarp(input_images, output_image, verbosity=True):
                        # SWarping the image to reduce all image to a similar wcs:
                        if verbosity: print('SWarping to fix the image...')
                        swarpName = swarp_name(verbosity)
                        callSt = "%s -verbose_type quiet -BACK_TYPE MANUAL -COMBINE_TYPE AVERAGE -VERBOSE_TYPE QUIET " % (swarpName)
                        callSt += " ".join(["%s[0]" % (s) for s in input_images])
                        subprocess.call(callSt, shell="True")
                        if verbosity: print('Done!')
                        
                        # Rename and remove SWarp tmp files 
                        shutil.move('coadd.fits', output_image)
                        os.remove('coadd.weight.fits')
                        os.remove('swarp.xml')    
    
    

def main(input_images, psf_images=None, mask_image=None, scale=0.396, m0=None, psf_number_to_match=1, output_combined_image='combined.fits', deep_picture='combined.png', SB_bright=24., SB_faint=27., sigma_smooth=2., add_axes=False, galaxy_name='', method='imcombine'):
        output_combined_image = os.path.abspath(output_combined_image)
        
        if psf_images is not None:
            for psf_image in psf_images:
                # add pixelscale to psf fits
                subprocess.call("addpixscl %s %f" % (psf_image, scale), shell=True, stdout=FNULL, stderr=subprocess.STDOUT) 

            
            
            for psf_file in ['match_psf_ui.fits','match_psf_gi.fits','match_psf_ri.fits','match_psf_zi.fits']:
                if os.path.exists(psf_file):
                    os.remove(psf_file)
            
            
            
            # match psf to i band
            matched_images = []
            for k in range(len(psf_images)):
                if k+1!=psf_number_to_match:
                    subprocess.call("pypher %s %s match_psf_%i%i.fits -r 1.e-5" % (psf_images[k], psf_images[psf_number_to_match-1], k+1, psf_number_to_match), shell=True, stdout=FNULL, stderr=subprocess.STDOUT) 
            
            
                    #subprocess.call("pypher psf_g.fits psf_i.fits match_psf_gi.fits -r 1.e-5", shell=True, stdout=FNULL, stderr=subprocess.STDOUT) 
                    #subprocess.call("pypher psf_r.fits psf_i.fits match_psf_ri.fits -r 1.e-5", shell=True, stdout=FNULL, stderr=subprocess.STDOUT) 
                    #subprocess.call("pypher psf_z.fits psf_i.fits match_psf_zi.fits -r 1.e-5", shell=True, stdout=FNULL, stderr=subprocess.STDOUT) 
                
                
                    convolve.convolution(input_images[k], 'match_psf_%i%i.fits' % (k+1, psf_number_to_match), 'galaxy_%i%i.fits' % (k+1, psf_number_to_match))
                    #convolve.convolution('galaxy_g.fits', 'match_psf_gi.fits', 'galaxy_gi.fits')
                    #convolve.convolution('galaxy_r.fits', 'match_psf_ri.fits', 'galaxy_ri.fits')
                    #convolve.convolution('galaxy_z.fits', 'match_psf_zi.fits', 'galaxy_zi.fits')
                    matched_images.append(os.path.abspath('galaxy_%i%i.fits' % (k+1, psf_number_to_match)))
            
            matched_images.append(input_images[psf_number_to_match-1])
        else:
            matched_images = input_images
        
        # Combine bands
        if os.path.exists(output_combined_image):
            os.remove(output_combined_image)
        
        if method=='imcombine':
            imcombine.imcombine(matched_images, output_combined_image, 0.1, 0.1)
            hdulist0 = pyfits.open(input_images[0])
            frame0 = hdulist0[0].data
            header0 = hdulist0[0].header

            hdulist1 = pyfits.open(output_combined_image)
            frame1 = hdulist1[0].data
            ny,nx = np.shape(frame1)
            
            header0['NAXIS1'] = nx 
            header0['NAXIS2'] = ny
            
            hdu = pyfits.PrimaryHDU(frame1, header0)
            hdu.writeto('tmp.fits', overwrite=True) 
            os.rename('tmp.fits', output_combined_image)
        else:
            #combine_images(input_images, output_combined_image) # WARNING: Works worse than imcombine!                 
            combine_images_swarp(matched_images, output_combined_image, verbosity=True)
        


        #'''
        
        if m0 is None:
            m0 = 28.
        else:
            m0 = np.mean(np.array(m0,float))
        
        if sigma_smooth is not None:
            # plot smoothed image
            png_image = plot_smoothed_image.main(output_combined_image, mask_image, galaxy_name, m0, scale, SB_bright=SB_bright, SB_faint=SB_faint, sigma_smooth=sigma_smooth, add_axes=add_axes)
            if png_image!=deep_picture:
                shutil.copy(png_image, deep_picture)
        #'''
        return m0



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Stacking images")
    parser.add_argument("input_images", help="Input fits images separated by comma")
    parser.add_argument("--psf_images", nargs='?', const=1, help="Optional: Input psf images separated by comma. In this case psf matching will be done.",type=str, default=None)     
    parser.add_argument("--mask_image", nargs='?', const=1, help="Optional: Input mask image.",type=str, default=None)     
   

    parser.add_argument("--m0", nargs='?', const=1, help="Optional: Zero-points separated by comma", type=str, default=None)  
    parser.add_argument("--scale", nargs='?', const=1, help="Optional: Pixel scale", type=float, default=1.0) 

    parser.add_argument("--psf_number", nargs='?', const=1, help="Optional: Number of image (srating from 1) in your list to which all images will be comvolved",type=int,default=1) 
    
    parser.add_argument("--o", nargs='?', const=1, help="Optional: Output combined image.",type=str, default='combined.fits')  
    parser.add_argument("--p", nargs='?', const=1, help="Optional: Output combined picture.",type=str, default='combined.png') 
    
    parser.add_argument("--SBbright", nargs='?', const=1, help="Optional: Brightest SB level", type=float, default=24.0) 
    parser.add_argument("--SBfaint", nargs='?', const=1, help="Optional: Fintest SB level", type=float, default=27.0) 
    parser.add_argument("--sigma", nargs='?', const=1, help="Optional: Smoothing sigma for Gaussian filter", type=float, default=5.0) 
    
    parser.add_argument("--name", nargs='?', const=1, help="Optional: Galaxy name.",type=str, default='') 
    parser.add_argument("--axes", action="store_true", default=False, help="Add axes")
    
    args = parser.parse_args()
    
    input_images = args.input_images.split(',')

    mask_image = args.mask_image
    
    if args.psf_images is not None:
        psf_images = args.psf_images.split(',')
    else:
        psf_images = args.psf_images
    scale = args.scale
    
    if args.m0 is not None:
        m0 = args.m0.split(',')
    else:
        m0=28.
    
    psf_number = int(args.psf_number)
    
    output_combined_image = args.o
    deep_picture = args.p
    SB_bright = args.SBbright
    SB_faint = args.SBfaint
    sigma_smooth = args.sigma
    add_axes = args.axes
    galaxy_name = args.name    

    main(input_images, psf_images=psf_images, mask_image=mask_image, scale=scale, m0=m0, psf_number_to_match=psf_number, output_combined_image=output_combined_image, deep_picture=deep_picture, SB_bright=SB_bright, SB_faint=SB_faint, sigma_smooth=sigma_smooth, add_axes=add_axes, galaxy_name=galaxy_name)
