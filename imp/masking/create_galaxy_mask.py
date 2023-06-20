from astropy.io import fits as pyfits
import numpy as np
import math
import itertools
import matplotlib.pyplot as plt
from scipy import ndimage
import sys
from matplotlib.path import Path
from itertools import product
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import os
import subprocess
from astropy.modeling import models, fitting
from astropy.stats import sigma_clipped_stats
#from photutils import data_properties, properties_table
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from astropy import wcs
from astroquery import ned
from astropy import units as u
from astropy.coordinates import SkyCoord

LOCAL_DIR = "/imp/masking"
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

import run_SExtractor
import convert_segm_to_region
import convert_reg_to_mask
import auto_masking
import determine_sky
import get_galaxy_ellipse
import azimProfile
import plot_azim_averaged_profile
import crop_galaxy_image
import azim_aver_masking
import get_galaxy_center_ned
import model_masking
import merge_masks
import inner_masking
import crop_psf_star
import add_keyw_to_header
import rebin_image
import sersic_fitting
import add_wcs
#import photometric_calibration
import sky_around_galaxy
import find_good_psf_stars
import test_for_deepness
import rotate_galaxy
import find_iraf_ellipse_magnitude
import iraf_ellipse
import schlafly_extinction
import mark_objects_to_unmask
import fix_mask
import mto_func
import show_segm_over_image
import arithm_operations
import check_background
import plot_smoothed_image
import mark_mto_regions_of_galaxy

def remove_file(file):
    if os.path.exists(file):
        os.remove(file)


def main(input_image, mask_image, region_file, sky_std, masking_tool='combined', user_interact=True, verbosity=True):
        if masking_tool=='sextractor':
            xc,yc,sma,smb,PA = crop_galaxy_image.read_region(region_file)  
            output_region_file = auto_masking.main(input_image, output_region_file='mask_outside.reg', snr=2., min_pix=5, region_type='polygon', sextr_setup=None, galaxy_ellipse=[[xc,yc], sma/4., smb/4., PA], offset_size=1.2, offset_pix=0.)
        
            inner_masking.detect_stars(input_image, psf_image, output_region='mask_inside.reg', median_sky=0., std_sky=sky_std, sigma_numb=5., min_distance=15, offset_size=1.5, min_pix=5, galaxy_ellipse=[[xc,yc], sma, smb, PA])
            
            merge_masks.main(['mask_outside.reg','mask_inside.reg'], 'galaxy_mask.reg')
            if user_interact:
                subprocess.call("ds9 %s -scale histequ -regions galaxy_mask.reg" % (input_image), shell=True)
                answer='no'
                while answer=='no':
                    convert_reg_to_mask.mask(input_image, 'galaxy_mask.reg', output_image=output_image, output_mask=mask_image, mask_value=1, show_running=verbosity, mask_DN=None)
                    subprocess.call("ds9 %s -scale histequ -regions galaxy_mask.reg" % (input_image.split('.fits')[0] + '_clean.fits'), shell=True)
                    answer = 'yes'
                    answer = str(input('Are you happy with this (YES/no)?') or 'yes')
                    if answer!='no':
                        convert_reg_to_mask.mask(input_image, 'galaxy_mask.reg', output_image=output_image, output_mask=mask_image, mask_value=1, show_running=verbosity, mask_DN=None)

            else:
                convert_reg_to_mask.mask(input_image, 'galaxy_mask.reg', output_image=output_image, output_mask=mask_image, mask_value=1, show_running=verbosity, mask_DN=None)
        
        elif masking_tool=='mto':
            #mto_func.main(input_image)
            if verbosity:
                    subprocess.call("python3 %s/mto.py %s -out %s -bg_mean 0. -verbosity 1" % (os.path.join(IMAN_DIR, 'imp/masking/mto-lib'),input_image,'segm_mto.fits'), shell=True)
            else:
                    subprocess.call("python3 %s/mto.py %s -out %s -bg_mean 0. -verbosity 0" % (os.path.join(IMAN_DIR, 'imp/masking/mto-lib'),input_image,'segm_mto.fits'), shell=True)   


            xc,yc,sma,smb,PA = crop_galaxy_image.read_region(region_file)  
            
            coordinates_to_unmask = mark_mto_regions_of_galaxy.main('segm_mto.fits',xc, yc, 'parameters.csv')
            mark_objects_to_unmask.main(input_image, coordinates=coordinates_to_unmask, output_region='objects_to_unmask.reg', coord_system='pix')
            fix_mask.main(input_image, 'segm_mto.fits', 'objects_to_unmask.reg', new_segm_image=mask_image, clean_image=input_image.split('.fits')[0]+'_clean.fits',verbosity=verbosity)
            
            if user_interact:
                show_segm_over_image.main(input_image, mask_image, out_region_file = 'objects_to_unmask.reg')
                #subprocess.call("ds9 %s -scale histequ -regions objects_to_unmask.reg" % (input_image.split('.fits')[0] + '_clean.fits'), shell=True)
                answer = 'yes'
                answer = str(input('Are you happy with this (YES/no)?') or 'yes')

                while answer=='no':
                    fix_mask.main(input_image, mask_image, 'objects_to_unmask.reg', new_segm_image=mask_image, clean_image=input_image.split('.fits')[0]+'_clean.fits',verbosity=verbosity)
                    #show_segm_over_image.main(input_image, mask_image, out_region_file = 'objects_to_unmask.reg')
                    subprocess.call("ds9 %s -scale histequ -regions objects_to_unmask.reg" % (input_image.split('.fits')[0] + '_clean.fits'), shell=True)
                    answer = 'yes'
                    answer = str(input('Are you happy with this (YES/no)?') or 'yes')
       
            os.remove('objects_to_unmask.reg')
       
        elif masking_tool=='combined':
            answer = 'yes'
            if os.path.exists('galaxy_mask.reg') and os.path.exists(mask_image) and user_interact:
                answer = str(input('File with the galaxy mask already exists! Do you want to remove it and start masking again? (yes/NO)?') or 'no')
            
            if answer=='no':
                if verbosity: print('Editing the existing mask...')
                show_segm_over_image.main(input_image, mask_image, out_region_file = 'galaxy_mask.reg')
                #fix_mask.main(input_image, 'segm_mto.fits','galaxy_mask.reg', region_file, new_segm_image=mask_image, clean_image=input_image.split('.fits')[0]+'_clean.fits')
                
                #convert_reg_to_mask.mask(input_image, 'galaxy_mask.reg', output_image=None, output_mask='tmp_mask.fits', mask_value=1, show_running=verbosity, mask_DN=None)
                #merge_masks.main(['tmp_mask.fits','mask.fits'], 'tmp_tmp_mask.fits')
                #shutil.move('tmp_tmp_mask.fits', 'mask.fits')
                #remove_file('tmp_mask.fits')

                fix_mask.main_repeat(input_image, 'segm_mto.fits', 'galaxy_mask.reg', new_segm_image=mask_image, clean_image=input_image.split('.fits')[0]+'_clean.fits',verbosity=verbosity)
            else:
                # Outer Sextractor mask
                output_region_file = auto_masking.main(input_image, output_region_file='mask_outside.reg', snr=2., min_pix=5, region_type='polygon', sextr_setup=None, galaxy_ellipse=region_file, offset_size=1.2, offset_pix=0., verbosity=verbosity)

                
                # MTO mask for the inner part
                #gain = mto_func.main(input_image, out='segm_mto.fits', bg_mean= 0., bg_variance=sky_std**2, verbosity=verbosity) #### WARNING: Often crashes!
                if verbosity:
                    subprocess.call("python3 %s/mto.py %s -out %s -bg_mean 0. -bg_variance %f -verbosity 1" % (os.path.join(IMAN_DIR, 'imp/masking/mto-lib'),input_image,'segm_mto.fits', sky_std**2), shell=True)
                else:
                    subprocess.call("python3 %s/mto.py %s -out %s -bg_mean 0. -bg_variance %f -verbosity 0" % (os.path.join(IMAN_DIR, 'imp/masking/mto-lib'),input_image,'segm_mto.fits', sky_std**2), shell=True)                
                
                xc,yc,sma,smb,PA = crop_galaxy_image.read_region(region_file)
                coordinates_to_unmask = mark_mto_regions_of_galaxy.main('segm_mto.fits',xc, yc, 'parameters.csv')
                
                mark_objects_to_unmask.main(input_image, coordinates=coordinates_to_unmask, output_region='objects_to_unmask.reg', coord_system='pix')
                merge_masks.main(['mask_outside.reg','objects_to_unmask.reg'], 'galaxy_mask.reg', verbosity=verbosity)
   
                fix_mask.main(input_image, 'segm_mto.fits', 'galaxy_mask.reg', region_file, new_segm_image=mask_image, clean_image=input_image.split('.fits')[0]+'_clean.fits',verbosity=verbosity)            

            if user_interact:
                show_segm_over_image.main(input_image, mask_image, out_region_file = 'galaxy_mask.reg')
                answer = 'yes'
                answer = str(input('Are you happy with this (YES/no)?') or 'yes')

                while answer=='no':
                    #fix_mask.main(input_image, 'segm_mto.fits','galaxy_mask.reg', region_file, new_segm_image=mask_image, clean_image=input_image.split('.fits')[0]+'_clean.fits')
                    fix_mask.main_repeat(input_image, 'segm_mto.fits', 'galaxy_mask.reg', new_segm_image=mask_image, clean_image=input_image.split('.fits')[0]+'_clean.fits',verbosity=verbosity)
                    show_segm_over_image.main(input_image, mask_image, out_region_file = 'galaxy_mask.reg')
                    #subprocess.call("ds9 %s -scale histequ -regions galaxy_mask.reg" % (input_image.split('.fits')[0] + '_clean.fits'), shell=True)
                    answer = 'yes'
                    answer = str(input('Are you happy with this (YES/no)?') or 'yes')            

            
            remove_file('segm.fits')
            remove_file('mask_tmp.reg')
            remove_file('parameters.csv')
            remove_file('inner_mask_mto.reg')
            remove_file('mask_outside.reg')
            remove_file('mask_segm.reg')
            remove_file('mask_segm.fits')
            remove_file('mask_segm_mto.fits')
            remove_file('mask_segm_mto.reg') 
            remove_file('objects_to_unmask.reg')
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Test for deepness of the image")
    parser.add_argument("input_image", help="Input image")
    parser.add_argument("mask", help="Input mask image")
    parser.add_argument("region_file", help="Region file with an ellipse which outlines the galaxy")
    parser.add_argument("skystd", help="RMS of the sky", type=float)


    args = parser.parse_args()

    input_image = args.input_image
    mask_image = args.mask
    region_file = args.region_file
    sky_std = args.skystd

    
    main(input_image, mask_image, region_file, sky_std, masking_tool='combined', user_interact=True, verbosity=True)
