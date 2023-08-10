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

from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from scipy import special
import scipy.ndimage as ndimage
import argparse
import collections

from PIL import Image

#LOCAL_DIR = "/plotting/2dprofile"
#IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

#sys.path.append(os.path.join(IMAN_DIR, 'plotting/2dprofile'))

import merge_regb_images
import plot_image



def convert_segm_to_boolean(mask):
    return np.ma.make_mask(mask, shrink=False)


def roundup(x, round_number=10):
    return int(math.ceil(x / float(round_number))) * round_number




def add_scale_bar(ax, nx, ny, pix2sec, kpc_per_arc=None, L_bar=60.):
    if L_bar is not None:
        #if Distance!=None:
        #kpc_per_arc = 4.84*float(Distance)/1000.
        #L_bar = roundup((nx/5.)*pix2sec*kpc_per_arc, round_number=1)

        #plot_image.crea_scale_bar(ax,nx,ny,pix2secs[0],l_bar=L_bar,kpc_per_arc=kpc_per_arc)
        #else:
        #L_bar = 60.#roundup((nx/5.)*pix2sec)
        if L_bar>=60.:
            L_bar = roundup(L_bar/60., round_number=1)
            angle_units='arcmin'
        else:
            angle_units='arcsec'
        if kpc_per_arc is None:
            kpc_per_arc = float('nan')
        
        plot_image.crea_scale_bar(ax,nx,ny,pix2sec,l_bar=L_bar,kpc_per_arc=kpc_per_arc,angle_units=angle_units, plot_scale_size=False)        


def resize_deep_images(deep_image, output_image):
    basewidth = 600
    img = Image.open(deep_image)
    
    wpercent = (basewidth/float(img.size[0]))
    hsize = int((float(img.size[1])*float(wpercent)))
    img = img.resize((basewidth,hsize), Image.ANTIALIAS)
    if not '.png' in output_image:
        img = img.convert("RGB") # Remove this if you want png!!!!
    img.save(output_image) 



def main(input_image, mask_image=None, galaxy_name='', m0=28., pix2sec=1., SB_bright=24., SB_faint=27., sigma_smooth=1., add_axes=False, output_image=None, show_mask=False, verbosity=True, L_bar=60., kpc_per_arc=None):
    if verbosity: print('Plotting smoothed image...')

    galaxyHDU = pyfits.open(input_image)
    galaxyData = galaxyHDU[0].data
    ny,nx = np.shape(galaxyData)
    
    if mask_image is not None:
        # Create astropy mask (boolean)
        maskHDU = pyfits.open(mask_image)
        mask = maskHDU[0].data
    else:
        mask = np.zeros([ny,nx])
    
    mask_astropy = convert_segm_to_boolean(mask)    
    
    # Smooth data
    galaxyData = ndimage.gaussian_filter(galaxyData, sigma=sigma_smooth, order=0)
    galaxyData_clean = galaxyData * ~mask_astropy
    if show_mask:
        galaxyData = galaxyData_clean
    

    fig =  plt.figure(0)
    ax = plt.subplot()
    vmax = np.nanmax(galaxyData_clean)
    vmin = np.nanmin(galaxyData_clean)
    vmax = 10**(0.4*(m0-SB_faint))*(pix2sec**2)
    
    plt.imshow(galaxyData,cmap='gist_heat_r', vmin=0., vmax=vmax, interpolation='none', origin="lower", aspect='equal')
    ax.set_frame_on(False)
    ax.set_xticks([]); ax.set_yticks([])
    plt.axis('off')    
    ax.text(0.03, 0.9, galaxy_name, fontsize=17, color='black',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline',backgroundcolor='whitesmoke')
    add_scale_bar(ax, nx, ny, pix2sec, L_bar=L_bar, kpc_per_arc=kpc_per_arc)
    plt.draw() 

    plt.savefig('1.png', bbox_inches='tight', pad_inches=0.01, dpi = 300)
    plt.clf()
    plt.close()       


    fig =  plt.figure(1)
    ax = plt.subplot()
    vmax = np.nanmax(galaxyData_clean) #np.mean([np.max(data_ski),np.max(data_ref)])
    vmin = np.nanmin(galaxyData_clean) #np.mean([np.min(data_ski),np.min(data_ref)])
    vmax = 10**(0.4*(m0-SB_bright))*(pix2sec**2)
    
    plt.imshow(galaxyData,cmap='gist_heat_r', vmin=0., vmax=vmax, interpolation='none', origin="lower", aspect='equal')
    ax.set_frame_on(False)
    ax.set_xticks([]); ax.set_yticks([])
    plt.axis('off')    
    ax.text(0.03, 0.9, galaxy_name, fontsize=17, color='black',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline',backgroundcolor='whitesmoke')
    add_scale_bar(ax, nx, ny, pix2sec, L_bar=L_bar, kpc_per_arc=kpc_per_arc)
    plt.draw() 

    plt.savefig('2.png', bbox_inches='tight', pad_inches=0.01, dpi = 300)
    plt.clf()
    plt.close()   


    fig =  plt.figure(2)
    ax = plt.subplot()
    vmax = np.nanmax(galaxyData_clean) #np.mean([np.max(data_ski),np.max(data_ref)])
    vmin = np.nanmin(galaxyData_clean) #np.mean([np.min(data_ski),np.min(data_ref)])
    #vmax = 10**(0.4*(m0-SB_bright))*(pix2sec**2)
    
    plt.imshow(galaxyData,cmap='gist_stern', vmin=0., vmax=vmax, interpolation='none', origin="lower", aspect='equal')
    ax.set_frame_on(False)
    ax.set_xticks([]); ax.set_yticks([])
    plt.axis('off')    
    ax.text(0.03, 0.9, galaxy_name, fontsize=17, color='black',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline',backgroundcolor='whitesmoke')
    add_scale_bar(ax, nx, ny, pix2sec, L_bar=L_bar, kpc_per_arc=kpc_per_arc)
    plt.draw() 

    plt.savefig('3.png', bbox_inches='tight', pad_inches=0.01, dpi = 300)
    plt.clf()
    plt.close() 
    
    merge_regb_images.main('1.png', '2.png', output_image='merged.png',verbosity=verbosity)
    merge_regb_images.main('merged.png', '3.png', output_image='%s.png' % ('combined1'),verbosity=verbosity)
    
    for file in ['1.png','2.png','3.png','merged.png']:
        os.remove(file)
    
    if add_axes:
        plot_image.axes_direct(input_image, '%s.png' % ('combined1'), 'tmp.png')
        os.rename('tmp.png', '%s.png' % ('combined1'))
    
    if output_image is not None:
        #shutil.move('%s.png' % ('combined1'), output_image)
        resize_deep_images('%s.png' % ('combined1'), output_image)
        if os.path.exists('%s.png' % ('combined1')):
            os.remove('%s.png' % ('combined1'))
        if verbosity: print('Done!')
        return output_image
    else:
        #os.rename('combined1.png', '%s.png' % ('combined'))
        resize_deep_images('%s.png' % ('combined1'), '%s.png' % ('combined'))

        if os.path.exists('%s.png' % ('combined1')):
            os.remove('%s.png' % ('combined1'))
        if verbosity: print('Done!')
        return '%s.png' % ('combined')
    

    
    









if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Plot smoothed image to highlight faint details")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("m0", help="Zero-point")
    parser.add_argument("pix2sec", help="Scale")
    parser.add_argument("--mask_image", help="Input fits mask", default=None)
    parser.add_argument("--name", default='',
                        help="Galaxy name")


    parser.add_argument("--SBbright", default=24.,
                        help="Maximum SB level")
    parser.add_argument("--SBfaint", default=27.,
                        help="Minimum SB level")
    parser.add_argument("--sigma", default=1.,
                        help="Sigma for gaussian smooth")
    parser.add_argument("--axes", action="store_true", default=False,
                        help="Add axes to the output image")

    args = parser.parse_args()
    
    
    input_image = args.input_image
    mask_image = args.mask_image
    
    
    
    m0 = float(args.m0)
    pix2sec = float(args.pix2sec)
    axes= args.axes
    sigma = float(args.sigma)
    galaxy_name = args.name
    
    SBbright = float(args.SBbright)
    SBfaint = float(args.SBfaint)
    
    
    
    
    
    main(input_image, mask_image, galaxy_name, m0, pix2sec, SB_bright=SBbright, SB_faint=SBfaint, sigma_smooth=sigma, add_axes=axes)


