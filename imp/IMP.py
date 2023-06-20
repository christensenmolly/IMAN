#!/usr/bin/python
# DESCRIPTION:
# Script to prepare galaxy images
# MINIMAL USAGE: python IMP.py [input_image]

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

LOCAL_DIR = "/imp"
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
import photometric_calibration
import sky_around_galaxy
import find_good_psf_stars
import test_for_deepness
import rotate_galaxy
import find_iraf_ellipse_magnitude
import iraf_ellipse
import schlafly_extinction
import extinction
import Arenou_model
import mark_objects_to_unmask
import fix_mask
import mto_func
import show_segm_over_image
import arithm_operations
import check_background
import plot_smoothed_image
import mark_mto_regions_of_galaxy


def retrieve_l_b(name, RA=None, DEC=None):
    if RA is None or DEC is None:
        result_table = ned.Ned.query_object(name)
        #print(result_table)
        try:
            RA,DEC = float(result_table['RA(deg)']),float(result_table['DEC(deg)']) # an astropy.table.Table
        except:
            RA,DEC = float(result_table['RA']),float(result_table['DEC']) # an astropy.table.Table
    c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='fk5').galactic
    return c.l.degree,c.b.degree


def ext_func(Av, Rv, wavelength):
    # Wavelength in um
    return extinction.ccm89(np.array([wavelength*10000.]), Av, Rv)[0]

def remove_file(file):
    if os.path.exists(file):
        os.remove(file)
        


def IMP_main(input_image,
             action,
             output_image=None,
             psf_image='psf.fits',
             xc=None, yc=None,
             mask_image=None,
             region_file=None,
             user_interact = True,
             sma=None, smb=None, PA=None,
             azim_tab=None,
             galaxy_ellipse=None,
             number=1,
             m0=20.,scale=1.0,
             band='R',catalogue='NOMAD',
             out_directory=None,
             annulus_width=None,
             sky_degree=5, sky_mean=0., sky_std=None, sky_med_unc=None,
             offset_size=1.,
             crop_box=False,
             step=1.,
             wavelength=None,
             RA=None, DEC=None, name=None,
             ellipse_file=None,
             Distance=None,
             sigma_cal=0.,
             SB_lim=28.,
             linear='yes',
             rmin=1., rmax=None,
             hcenter='no',
             masking_tool='combined',
             exptime = 1.0,
             Ncombine = 1,
             gain = None,
             ron = None,
             verbosity=True,
             sigma_smooth=1.,
             galaxy_name=None,
             ext_map='Green', 
             factor=None):


    
    if action=='add_astrometry':
        add_wcs.main(input_image, scale, band, region_file, output_image=None, define_ZP=False)
    
    if action=='find_good_stars':
        fwhm_med, ell_med, theta_med = find_good_psf_stars.main(input_image, None, band, catalogue, mag_min=15.6, mag_max=18.3, sex_cat='field.cat', verbosity=verbosity)
        os.remove('mask.fits')
        return fwhm_med, ell_med, theta_med
    
    if action=='phot_calibration':
        #band = ['r','R','r','r']
        #catalogue = ['I/345/gaia2','NOMAD','V/139','II/349/ps1']
        results = photometric_calibration.main(input_image, region_file, mask_image, bands=band, catalogues=catalogue, star_number=None, obs_wavelength=wavelength, ext_map=ext_map, verbosity=verbosity, N_phot_stars=5)
        return results
    
    if action=='ima_stats':
        scale, note = rebin_image.resolution(input_image)
        scale = float(scale)
        
        limit, median_MEAN, median_MEDIAN, std_MEDIAN, median_STD = test_for_deepness.sky_in_boxes(input_image, m0, scale, mask_image=mask_image, box_size_arcsec=10., Nboxes=1000, n_sigma=3, units='mag', upper=False, verbosity=verbosity)
        return limit, median_MEAN, median_MEDIAN, std_MEDIAN, median_STD, float(scale)

    if action == 'image_masking':
        output_region_file = auto_masking.main(input_image, output_region_file=region_file, snr=2., min_pix=5, region_type='polygon', sextr_setup=None, galaxy_ellipse=None, offset_size=offset_size, offset_pix=0., verbosity=verbosity)
        output_image_tmp, output_mask_tmp = convert_reg_to_mask.mask(input_image, region_file, output_image='clean.fits', output_mask= region_file.split('.reg')[0]+'.fits', mask_value=1, verbosity=verbosity, mask_DN=None)
        
        if masking_tool=='sextractor':
            if user_interact:
                file_1 = os.stat(region_file).st_mtime
                subprocess.call("ds9 %s -scale histequ -regions %s" % (input_image, region_file), shell=True)  
                file_2 = os.stat(region_file).st_mtime
                if file_1 != file_2:
                    output_image_tmp, output_mask_tmp = convert_reg_to_mask.mask(input_image, region_file, output_image=None, output_mask= region_file.split('.reg')[0]+'.fits', mask_value=1, verbosity=verbosity, mask_DN=None)
            
            shutil.move(region_file.split('.reg')[0]+'.fits', mask_image)
            return float('nan')
        
        else:
            output_image, mean, median, std = determine_sky.sky_subtraction(input_image, output_mask_tmp, polynomial_degree=0, output_image=None, output_sky=None, hdu_inp=0)
            
            gain = mto_func.main(input_image, out='segm_mto.fits', bg_mean=median, bg_variance=std**2, verbosity=verbosity)
        
            if user_interact:
                show_segm_over_image.main(input_image, 'segm_mto.fits', out_region_file=region_file)
                #subprocess.call("ds9 %s -scale histequ -regions %s" % (input_image, region_file), shell=True)  

            output_image_tmp, output_mask_tmp = convert_reg_to_mask.mask(input_image, region_file, output_image=None, output_mask=region_file.split('.reg')[0]+'.fits', mask_value=1, verbosity=verbosity, mask_DN=None)

            arithm_operations.main('segm_mto.fits', 1., 'add', 'tmp_mask.fits')
            shutil.move('tmp_mask.fits', 'segm_mto.fits')
            
            arithm_operations.main(mask_image, 'segm_mto.fits', 'add', 'tmp_mask.fits')
            shutil.move('tmp_mask.fits', mask_image)
            
            os.remove('segm_mto.fits')
            os.remove('segm.fits')
            os.remove('parameters.csv')
            os.remove('mask_tmp.reg')
            
            return gain
            
    if action=='reg_to_mask':
        output_image,output_mask = convert_reg_to_mask.mask(input_image, region_file, output_image=output_image, output_mask=None, mask_value=1, verbosity=verbosity, mask_DN=None)
        
        #if user_interact:
        #    ds9Proc = subprocess.Popen(["ds9", output_image,
        #                                "-scale", "hist"])
        #    ds9Proc.wait()
        #    subprocess.call("ds9 %s -scale histequ -cmap Cool" % (output_image), shell=True)  

    if action=='sky_correction':
        output_image, mean, median, std = determine_sky.sky_subtraction(input_image, mask_image, polynomial_degree=sky_degree, output_image=output_image, output_sky=None, hdu_inp=0, sampling_factor=factor, sigma_smooth=sigma_smooth, verbosity=verbosity)
        
                
        if user_interact:
          answer='no'
          while answer=='no':
            check_background.main(output_image, mask_image, sigma_smooth=10)
            subprocess.call("ds9 %s -scale histequ %s -scale histequ" % (input_image, output_image), shell=True)

            answer = 'yes'
            answer = str(input('Are you happy with this (YES/no)?') or 'yes')
            if answer=='no':
                sky_degree = str(input('Choose a degree for the polynom to fit the background, 0 by default?') or 0)
                output_image, mean, median, std = determine_sky.sky_subtraction(input_image, mask_image, polynomial_degree=int(sky_degree), output_image=output_image, output_sky=None, hdu_inp=0, sampling_factor=factor, sigma_smooth=sigma_smooth, verbosity=verbosity)
                answer='no'
            
        return output_image, mean, median, std, sky_degree

    if action=='sky_galaxy':
        sky,std,bkg_median_std = sky_around_galaxy.main(input_image, output_image, mask_image, manual=user_interact, degree=0, annulus_width=annulus_width, galaxy_ellipse=region_file, Npix=100000, offset_size=offset_size, show_running=verbosity) 
        if user_interact:
          answer='no'
          while answer=='no':
            answer = 'yes'

            check_background.main(output_image, mask_image, sigma_smooth=5)
            
            answer = str(input('Are you happy with this (YES/no)?') or 'yes')
            if answer=='no':
                sky,std,bkg_median_std = sky_around_galaxy.main(input_image, output_image, mask_image, manual=user_interact, degree=0, annulus_width=annulus_width, galaxy_ellipse=region_file, Npix=100000, offset_size=offset_size, show_running=verbosity)       
        
        return sky,std,bkg_median_std

    if action=='determine_galaxy':
        position,a,b,PA = get_galaxy_ellipse.main(input_image, segm_image=mask_image, output_region=region_file,  xc=xc, yc=yc, sky_background=sky_mean, min_radius=smb, verbosity=verbosity)
        if user_interact:
            print('Determine galaxy ellipse')
            subprocess.call("ds9 %s -scale histequ -regions %s" % (input_image, region_file), shell=True)            
        return position[0],position[1],a,b,PA


    if action=='crea_azim_profile':
        output_model, azim_tab = azimProfile.main(input_image, xcen=xc, ycen=yc, ell=1.-smb/sma, posang=PA, sma_min=1., sma_max=2.5*sma, step=5.)
        return output_model, azim_tab

    if action=='plot_azim_profile':
        r_lim = plot_azim_averaged_profile.main(azim_tab,1.,20.,plot_add=True, bin_box=3)
        r_lim = 1.5*r_lim
        get_galaxy_ellipse.create_ellipse_region(xc, yc, r_lim, r_lim*smb/sma, PA, region_file)
        #ds9Proc = subprocess.Popen(["ds9", input_image,
        #                                "-regions", region_file,
        #                                "-scale", "hist"])
        #ds9Proc.wait()
        subprocess.call("ds9 %s -scale histequ -cmap Cool -regions %s" % (input_image, region_file), shell=True)  
        
        return r_lim




    if action=='crop_galaxy':
        output_image = crop_galaxy_image.main(input_image, sma, smb, galaxy_region=region_file, offset_size=offset_size, crop_box=crop_box, output_image=output_image)
        if user_interact:
          answer='no'
          while answer=='no':
            subprocess.call("ds9 %s -scale histequ" % (output_image), shell=True)
            answer = str(input('Are you happy with this (YES/no)?') or 'yes')
            if answer=='no':
                #new_offset_size = str(input('Enter new offset size (was %.1f)' % (offset_size)) or 2.)
                #if new_offset_size is None:
                #    answer='yes'
                #else:
                output_image = crop_galaxy_image.crop_manual(input_image, output_image, offset_size=offset_size, region_file=region_file)
                answer='no'
        return output_image


    if action=='rotate_galaxy':
        if sky_std is not None:
            I_DN_min = 3.* sky_std
        else:
            I_DN_min = None
        if mask_image is None:
            rotate_galaxy.main([input_image], xc, yc, I_DN_min=I_DN_min, I_DN_max=None, output_images = None)
        else:
            rotate_galaxy.main([input_image,mask_image], xc, yc, I_DN_min=I_DN_min, I_DN_max=None, output_images = None)

    if action=='crea_psf':
        crop_psf_star.crop_star(input_image, region_file, output_star_image=psf_image, star_number=number, factor=1, verbosity=verbosity, azim_model=False)
        if user_interact:
          answer='no'
          while answer=='no':
            subprocess.call("ds9 %s -scale histequ" % (psf_image), shell=True)
            answer = 'yes'
            answer = str(input('Are you happy with this (YES/no)?') or 'yes')
            if answer=='no':
                subprocess.call("ds9 %s -scale histequ -regions %s" % (input_image, region_file), shell=True)
                answer = None
                answer = str(input('Choose the number of star to be taken as PSF?') or None)
                if answer is None:
                    print('No number was given. Exiting...')
                    exit()
                else:
                    crop_psf_star.crop_star(input_image, region_file, output_star_image=psf_image, star_number=int(answer), factor=1, azim_model=False)
                    answer='no'
            
    if action == 'final_masking':
        if masking_tool=='sextractor':
            xc,yc,sma,smb,PA = crop_galaxy_image.read_region(region_file)  
            output_region_file = auto_masking.main(input_image, output_region_file='mask_outside.reg', snr=2., min_pix=5, region_type='polygon', sextr_setup=None, galaxy_ellipse=[[xc,yc], sma/4., smb/4., PA], offset_size=1.2, offset_pix=0.)
        
            inner_masking.detect_stars(input_image, psf_image, output_region='mask_inside.reg', median_sky=0., std_sky=sky_std, sigma_numb=5., min_distance=15, offset_size=1.5, min_pix=5, galaxy_ellipse=[[xc,yc], sma, smb, PA])
            
            merge_masks.main(['mask_outside.reg','mask_inside.reg'], 'galaxy_mask.reg')
            if user_interact:
                subprocess.call("ds9 %s -scale histequ -regions galaxy_mask.reg" % (input_image), shell=True)
                answer='no'
                while answer=='no':
                    convert_reg_to_mask.mask(input_image, 'galaxy_mask.reg', output_image=output_image, output_mask=mask_image, mask_value=1, verbosity=verbosity, mask_DN=None)
                    subprocess.call("ds9 %s -scale histequ -regions galaxy_mask.reg" % (input_image.split('.fits')[0] + '_clean.fits'), shell=True)
                    answer = 'yes'
                    answer = str(input('Are you happy with this (YES/no)?') or 'yes')
                    if answer!='no':
                        convert_reg_to_mask.mask(input_image, 'galaxy_mask.reg', output_image=output_image, output_mask=mask_image, mask_value=1, verbosity=verbosity, mask_DN=None)

            else:
                convert_reg_to_mask.mask(input_image, 'galaxy_mask.reg', output_image=output_image, output_mask=mask_image, mask_value=1, verbosity=verbosity, mask_DN=None)
        
        elif masking_tool=='mto':
            #mto_func.main(input_image)
            if verbosity:
                    subprocess.call("python3 %s/mto.py %s -out %s -bg_mean 0. -verbosity 1" % (os.path.join(IMAN_DIR, 'imp/masking/mto-lib'),input_image,'segm_mto.fits'), shell=True)
            else:
                    subprocess.call("python3 %s/mto.py %s -out %s -bg_mean 0. -verbosity 0" % (os.path.join(IMAN_DIR, 'imp/masking/mto-lib'),input_image,'segm_mto.fits'), shell=True)   


            xc,yc,sma,smb,PA = crop_galaxy_image.read_region(region_file)  
            
            coordinates_to_unmask = mark_mto_regions_of_galaxy.main(xc, yc, 'parameters.csv')
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
                
                #convert_reg_to_mask.mask(input_image, 'galaxy_mask.reg', output_image=None, output_mask='tmp_mask.fits', mask_value=1, verbosity=verbosity, mask_DN=None)
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
                coordinates_to_unmask = mark_mto_regions_of_galaxy.main(xc, yc, 'parameters.csv')
                
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


    if action=='fitting':
        if not os.path.exists(out_directory):
            os.mkdir(out_directory)
        else:
            shutil.rmtree(out_directory)
            os.mkdir(out_directory)
        
        if region_file is not None:
                    hdulist = pyfits.open(input_image)
                    data = hdulist[0].data
                    nx, ny = data.shape[1], data.shape[0]

                    xc, yc, sma, smb, PA = crop_galaxy_image.read_region(region_file)
                    [[x_l, y_l], [x_r, y_r]] = crop_galaxy_image.ellipse_borders([xc, yc], 2. * offset_size * sma, 2. * offset_size * smb, PA, nx=nx,
                                                   ny=ny, square=False)
            
            
                    results = sersic_fitting.main(input_image, None, psf_image, mask_image, m0, scale, xc=xc, yc=yc, sampling=1, output_dir = out_directory, xmin=x_l+1, xmax=x_r+1, ymin=y_l+1, ymax=y_r+1, verbose=verbosity)
        else:
                    hdulist = pyfits.open(input_image)
                    data = hdulist[0].data
                    nx, ny = data.shape[1], data.shape[0]
                    xc = nx/2.
                    yc = ny/2.
                    results = sersic_fitting.main(input_image, None, psf_image, mask_image, m0, scale, xc=xc, yc=yc, sampling=1, output_dir = out_directory, verbose=verbosity)
            
        return results

    if action=='ellipse_fitting':
        try:
            iraf_ellipse.main_ell(input_image, xc, yc, ellip=1.-smb/sma, pa=PA, sma0=sma, m0=m0, pix2sec=scale, step=step, minsma=rmin, maxsma=rmax, outp_format='png', ell_file='ellipse.txt', fits_mask=mask_image, fflag=1.0, olthresh=0.0, linear=linear)
            return 0
        except:
          if user_interact:
            maxsma = str(input('Crashed. Enter maximum semimajor-axis?') or '')
            if maxsma!='':
                iraf_ellipse.main_ell(input_image, xc, yc, ellip=1.-smb/sma, pa=PA, sma0=sma, m0=m0, pix2sec=scale, step=step, minsma=rmin, maxsma=float(maxsma), outp_format='png', ell_file='ellipse.txt', fits_mask=mask_image, fflag=1.0, olthresh=0.0, linear=linear,hcenter=hcenter)
                return 0
          else:
              return 1

    if action=='general_pars':
        # wavelength in mum, e.g. 0.661 mum = 661 nm
        l,b = retrieve_l_b(name, RA=RA, DEC=DEC)
        Aext = schlafly_extinction.main([[l,b]], wavelengths=[wavelength], frame='galactic', units='deg', output_file=None)[0]
        #Av = Arenou_model.Arenou([[l,b,100000000000.]])[0][0]
        #Aext = ext_func(Av,3.1,wavelength)

        res = find_iraf_ellipse_magnitude.main(input_image, ellipse_file, m0, scale, rms_sky=sky_std, apert_cor=0., sigma_image=None, sigma_cal=sigma_cal, sigma_apert=0., bin_box=step, SB_lim=SB_lim, rmin=10., sky_med_unc=sky_med_unc, Aext=Aext)
        #print(res)
        return res, Aext
    

    if action=='get_conversion_pars':
        pass
    
    if action=='smooth_picture':
        plot_smoothed_image.main(input_image, mask_image, galaxy_name, m0, scale, SB_bright=24., SB_faint=27., sigma_smooth=sigma_smooth, add_axes=True)



#IMP_main('NGC2683.phot.2.fits', 'image_masking', region_file='general_mask.reg', mask_image='general_mask.fits')
#IMP_main('NGC2683.phot.2.fits', 'sky_correction', region_file='general_mask.reg', mask_image='general_mask.fits', sky_degree=1, output_image='I2_skysub.fits')
#rotate_galaxy.main([input_image,mask_image], xc, yc, I_DN_min=I_DN_min, I_DN_max=None, output_images = None)
#IMP_main('combined.fits', 'rotate_galaxy', mask_image=None, xc=1072,yc=1161,sky_std=0.03)






#IMP_main('combined_someskysub_rot_crop.fits', 'final_masking', region_file=os.path.abspath('galaxy_ellipse.reg'), psf_image='psf_r.fits', output_image=None, mask_image='mask.fits',sky_std=0.00427001, user_interact=True, verbosity=True, masking_tool='combined')










'''


galaxy_name = 'NGC520'
xc,yc = get_galaxy_center_ned.main('new-image.fits', name=galaxy_name, RA=None, DEC=None)

#main('new-image.fits', 'image_masking', user_interact = True)
#main('new-image.fits', 'reg_to_mask', region_file='general_mask.reg', user_interact = True)
#output_image = main('new-image.fits', 'sky_correction', mask_image='general_mask.fits', user_interact = True)
#xc_ell,yc_ell,sma_ell,smb_ell,PA_ell = main('sky_subtr.fits', 'determine_galaxy', xc=xc, yc=yc, mask_image=None, region_file=None, user_interact = True)

#output_model, azim_tab = main('sky_subtr.fits', 'crea_azim_profile', xc=xc, yc=yc, mask_image=None, region_file=None, user_interact = True, sma=sma_ell, smb=smb_ell, PA=PA_ell)
#r_lim = main('sky_subtr.fits', 'plot_azim_profile', xc=xc, yc=yc, azim_tab='azim_model.txt', sma=sma_ell, smb=smb_ell, PA=PA_ell, user_interact = False)
r_lim = 429.
main('sky_subtr.fits', 'sky_galaxy', output_image=None, mask_image='general_mask.fits', region_file='galaxy_ellipse_max.reg', user_interact = True, annulus_width=r_lim*0.414)
exit()
#main('sky_subtr.fits', 'crop_galaxy', xc=xc, yc=yc, sma=100, smb=100, PA=0.)


xc,yc,sma,smb,PA = crop_galaxy_image.read_region('galaxy_ellipse_max.reg')
xc,yc = get_galaxy_center_ned.main('sky_subtr_crop.fits', name=galaxy_name, RA=None, DEC=None)
get_galaxy_ellipse.create_ellipse_region(xc, yc, sma, smb, PA, 'galaxy_ellipse_final.reg')
#main('new-image.fits', 'crea_psf', region_file='best_stars.reg', user_interact = True, number=1)
shutil.copy('sky_subtr_crop.fits', 'galaxy.fits')
main('galaxy.fits', 'final_masking', galaxy_ellipse=[[xc,yc],sma,smb,PA], output_image='galaxy_clean.fits', mask_image='mask.fits')

shutil.copy('galaxy_clean.fits', 'galaxy_clean_galf.fits')



hdulis = pyfits.open('galaxy_clean_galf.fits')
prihdr = hdulis[0].header

pix2sec,note = rebin_image.resolution('galaxy_clean_galf.fits')
pix2sec = float(pix2sec)  

  
  
FWHM = float(prihdr['PSF_FWHM'])
m0 = float(prihdr['ZP'])
RON = 10.
NCOMBINE = 1

if 'EGAIN' in prihdr.keys():
    GAIN = prihdr['EGAIN']
else:
    GAIN = 1.

if 'EXPTIME' in prihdr.keys():
    EXPTIME = prihdr['EXPTIME']
else:
    EXPTIME = 1.

add_keyw_to_header.add_to_header('galaxy_clean_galf.fits',EXPTIME,GAIN,NCOMBINE,RON,m0,pix2sec,FWHM,sky_level=0.,sky_subtr=1,xc=xc,yc=yc)

dec_path = '/home/amosenko/Toshiba_1/CurrentWork/Rich/FINAL_SAMPLE/Decomposition/results/%s' % (galaxy_name)
main('galaxy_clean_galf.fits', 'fitting', output_image=None, psf_image='psf.fits', xc=None, yc=None, mask_image='mask.fits', region_file=None, user_interact = True, sma=None, smb=None, PA=None, azim_tab=None,galaxy_ellipse=None,number=1,m0=20.,scale=1.0,band='R',catalogue='NOMAD', out_directory=dec_path)



#if not os.path.exists(dec_path):
#os.mkdir(dec_path)

#sersic_fitting.main('galaxy_clean_galf.fits', None, 'psf.fits', 'mask.fits', m0, pix2sec, xc=xc, yc=yc, sampling=1, output_dir = dec_path)

'''

