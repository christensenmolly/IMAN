#!/usr/bin/python

# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from photutils import detect_sources, segmentation
import argparse
import pyregion
import os
import astropy.units as u
import scipy.ndimage as ndi
import sys
from skimage.draw import ellipse
from skimage.measure import label, regionprops, regionprops_table
from skimage.transform import rotate
from photutils import make_source_mask
from astropy.stats import sigma_clipped_stats
import shutil
'''
image = np.zeros((600, 600))

rr, cc = ellipse(300, 350, 100, 220)
image[rr, cc] = 1

image = rotate(image, angle=15, order=0)

rr, cc = ellipse(100, 100, 60, 50)
image[rr, cc] = 1

print(image)
label_img = label(image)
regions = regionprops(label_img)

exit()
'''


LOCAL_DIR = "/imp/cropping"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'detect_objects'))
sys.path.append(os.path.join(IMAN_DIR, 'decomposition/simple_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/masking'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/sky_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))
sys.path.append(os.path.join(IMAN_DIR, 'Ellipse_photometry'))
sys.path.append(os.path.join(IMAN_DIR, 'plotting/1dprofile'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/psf'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/add_astrometry'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/phot_calibration'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/rotate'))
sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))
sys.path.append(os.path.join(IMAN_DIR, 'extinction_maps'))
sys.path.append(os.path.join(IMAN_DIR, 'iraf_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/masking/mto-lib'))
sys.path.append(os.path.join(IMAN_DIR, 'sextractor'))
sys.path.append(os.path.join(IMAN_DIR, 'plotting/2dprofile'))


import image_interpolation_astropy
import convert_segm_to_region
import remove_nonobject_isophotes
import convert_reg_to_mask

def create_ellipse_region(input_image, xc, yc, sma, smb, PA, file, coord_format='image'):
    if coord_format=='image':
        f_galaxy = open(file, 'w')
        f_galaxy.write('%s\n' % ('image') )
        f_galaxy.write('ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # color=red\n' % (xc, yc, sma, smb, PA))
        f_galaxy.close()
    else:
        region = 'image;ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f)' % (xc, yc, sma, smb, PA)
        r = pyregion.parse(region)
        
        f = pyfits.open(input_image)
        r2 = pyregion.parse(r).as_imagecoord(f[0].header)
        
        f_galaxy = open(file, 'w')
        f_galaxy.write('%s\n' % ('fk5') )
        f_galaxy.write('ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # color=red\n' % (xc, yc, sma, smb, PA))
        f_galaxy.close()        

def convert_segm_to_boolean(mask):
    return np.ma.make_mask(mask)

def main(input_image, mask_image, output_region='galaxy_polygon.reg', xc=None, yc=None, I_min=None, min_radius=10., verbosity=True):
    '''
    Function to determine the ellipse which fits the outermost galaxy isophote at signal-to-noise ratio equal snr.
    The galaxy has a center of xc, yc. If they are None, then the center of the image is taken as the galaxy center.
    min_radius says that the galaxy semi-major axis should be larger than this when searching for the outermost ellipse.
    '''
    input_image_old = input_image
    
    if mask_image is not None:
        image_interpolation_astropy.astropy_smoothing(input_image, mask_image, output_image='interpolated_tmp.fits', sigma_smooth=10., sampling_factor=3, sigma_back=None)
        input_image = 'interpolated_tmp.fits'
    else:
        shutil.copy(input_image, 'interpolated_tmp.fits')

    
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data 

    #hdulist_segm = pyfits.open(mask_image)
    #mask = hdulist_segm[0].data    
    
    
    ny,nx = np.shape(data)
    if xc is None or yc is None:
        xc = nx/2.
        yc = ny/2.
    
    if I_min is None:
        mask_I = make_source_mask(data, nsigma=2, npixels=5, dilate_size=11)
        mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask_I)
        I_min = 3.*std
        if verbosity: print('We estimate the outer isophote at the 3-sigma level: %f' % (I_min))
    
    galaxy_segm = np.zeros((ny,nx))
    galaxy_segm[data>I_min]=1 # Segm for objects with I > I_min
    
    # Smooth the segmentation
    segmap_float = ndi.uniform_filter(np.float64(galaxy_segm), size=3)
  #WARNING: NEW
    galaxy_segm = segmap_float > 0.5 #WARNING: NEW
    
    outHDU = pyfits.PrimaryHDU(galaxy_segm.astype(int))
    outHDU.writeto('segm_tmp.fits', overwrite=True)    
    
    
    # convert to region
    convert_segm_to_region.main('segm_tmp.fits', 'segm_tmp.reg', output_mask_image=None, fits_slice = 0, offset_size=1., offset_pix=0., xc=None, yc=None, system='image', ignore_value=None, verbosity=True)
    
    remove_nonobject_isophotes.main('segm_tmp.reg', output_region_file=output_region, xc=xc, yc=yc)
    
    os.remove('interpolated_tmp.fits')
    os.remove('segm_tmp.fits')
    os.remove('segm_tmp.reg')
    
    convert_reg_to_mask.mask(input_image_old, output_region, output_image=None, output_mask=output_region.split('.reg')[0]+'.fits', mask_value=1, show_running=True, mask_DN=None, verbosity=True)

    hdulist = pyfits.open(output_region.split('.reg')[0]+'.fits')
    label_img = hdulist[0].data     

    regions = regionprops(label_img.astype(int))

    os.remove(output_region.split('.reg')[0]+'.fits')    
    
    for props in regions:
        y0, x0 = props.centroid
        PA = np.degrees(props.orientation)
        smb = 0.5*props.minor_axis_length
        sma = 0.5*props.major_axis_length



    if verbosity: print('Galaxy ellipse:') 
    if verbosity: print('\txc,yc: %.1f, %.1f' % (x0, y0)) 
    if verbosity: print('\tsma [pix]: %.1f' % (sma))
    if verbosity: print('\tsmb [pix]: %.1f' % (smb))
    if verbosity: print('\tell: %.2f' % (1.-smb/sma))    
    if verbosity: print('\tPA [Degrees: Up=90, Right=0, counterclockwise]: %.1f' % (PA)) 

    create_ellipse_region(input_image_old, x0, y0, sma, smb, PA, output_region.split('.reg')[0]+'_ell.reg') 


  
    mask_astropy = convert_segm_to_boolean(label_img.astype(int))
    Flux = np.nansum(data * mask_astropy)
    
    if verbosity: print('\tTotal flux within %f DN: %.1f' % (I_min, Flux)) 

    return x0,y0,sma,smb,PA,Flux,mask_astropy




    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Determine galaxy polygon and ellipse")
    parser.add_argument("inputImage", help="Input image")
    parser.add_argument("mask_image", help="Mask image")
    parser.add_argument("--reg", help="Output region file", type=str, default='galaxy_polygon.reg') 
    parser.add_argument("--xcen", nargs='?', const=1, help="Optional: Galaxy center x",type=float, default=None)  
    parser.add_argument("--ycen", nargs='?', const=1, help="Optional: Galaxy center y",type=float, default=None)  
    parser.add_argument("--I_min", nargs='?', const=1, help="Optional: SB of the outer isophote (DN)", type=float, default=None)  
    parser.add_argument("--min_radius", nargs='?', const=1, help="Optional: Minumum radius of the galaxy to be determined",type=float,default=10.)  
    args = parser.parse_args()

    input_image = args.inputImage
    mask_image = args.mask_image
    output_region = args.reg
    xcen = args.xcen
    ycen = args.ycen
    I_min = args.I_min
    min_radius = args.min_radius    
    

    main(input_image, mask_image, output_region=output_region, xc=xcen, yc=ycen, I_min=I_min, min_radius=min_radius)
