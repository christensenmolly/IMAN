#!/usr/bin/python
# DESCRIPTION:
# Script to prepare galaxy images, which have been initially reduced. One should create a special table which should be given as an input.
# MINIMAL USAGE: python prepare_galaxies.py [input_table] [number_of_item_in_table]

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
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from astropy import wcs
from astroquery import ned
import pickle
import collections
import glob

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
sys.path.append(os.path.join(IMAN_DIR, 'imp'))
sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))
sys.path.append(os.path.join(IMAN_DIR, 'sextractor'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/ds9_regions'))


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
import arithm_operations
import read_data
import convert_coord_regions

from IMP import IMP_main

#*** Colour fonts ***
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

def check_fits_layers(input_image, verbosity):
        # Check that fits-image has more than one layer. If yes, then exit.
        # An image should have just one layer! Otherwise, use parse_cube_fits.py to parse the fits file.
        hdulist = pyfits.open(input_image)
        N_layers = len(hdulist)
        N_sublayers = len(np.shape(hdulist[0].data))
        
        if verbosity: print('N layers = %i' % (N_layers))
        if verbosity: print('N sublayers = %i' % (N_sublayers))
        if N_layers>1 or N_sublayers>2:
            if verbosity: print('An image should have just one layer! Otherwise, use parse_cube_fits.py to parse the fits file befores running this script. Exiting...')
            exit()

def wcs_existence(input_image, verbosity):
        hdulist = pyfits.open(input_image)
        image = hdulist[0].data
        header = hdulist[0].header

        if 'CTYPE1' not in header:
            if verbosity: print('No wcs found! We set the given coordinates ra, dec as the center of the image, with the scale from the table. N is up, E is left')
            return False
        else:
            return True

def add_wcs_to_header(input_image, output_image, x_ref, y_ref, RA_ref, DEC_ref, cdelt1, cdelt2, angle, verbosity):
    if verbosity: print('WARNING: As no WCS was found, a fake wcs was added!')
    hdulist = pyfits.open(input_image)
    referenceheader = pyfits.getheader(input_image, 0)
    outframe = pyfits.getdata(input_image, 0)
    header = referenceheader

    header['EQUINOX'] = 2.000000000000E+03
    header['RADECSYS'] = 'FK5'
    
    header['CTYPE1'] = 'RA---TAN'
    header['CUNIT1'] = 'deg'    
    header['CRVAL1'] = RA_ref
    header['CRPIX1'] = x_ref

    
    header['CTYPE2'] = 'DEC--TAN'
    header['CUNIT2'] = 'deg'  
    header['CRVAL2'] = DEC_ref
    header['CRPIX2'] = y_ref

    header['CDELT1'] = cdelt1
    header['CDELT2'] = cdelt2
    header['CROTA2'] = angle    
    
    header['WCS_NOTE'] = 'This WCS is not real!'   
    
    hdu = pyfits.PrimaryHDU(outframe, header)
    hdu.writeto(output_image, overwrite=True)
    if verbosity: print('New wcs has been added to the image!')



def prep_for_galfit(xc, yc):
    hdulis = pyfits.open('galaxy_clean_galf.fits')
    prihdr = hdulis[0].header

    pix2sec,note = rebin_image.resolution('galaxy_clean_galf.fits')
    pix2sec = float(pix2sec)  

    results = read_pickle()
    
    FWHM = results['PSF_FWHM [pix]']
    m0 = results['M0_AVERAGE [mag]']
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
    hdulis.close()
    return m0,pix2sec


def read_pickle():
    with open('results.pkl','rb') as file_object:
        raw_data = file_object.read()

    results = pickle.loads(raw_data)
    return results  

def write_pickle(results=None):
    output = open('results.pkl', 'wb')
    if results is None:
        results = collections.OrderedDict()
        results['Collected parameters']=''
   
    pickle.dump(results, output)
    output.close()    


def convert_region_ellipse(galaxy_name=None, RA=None, DEC=None, input_image=None, region_file=None):
    #hdulis = pyfits.open(input_image)
    #data = hdulis[0].data
    #ny,nx = np.shape(data)
    
    xc,yc = get_galaxy_center_ned.main(input_image, name=galaxy_name, RA=RA, DEC=DEC)
    # xc = nx/2.
    # yc = ny/2.

    
    f = open(region_file, "r")    
    for line in f:
        if "ellipse" in line:
            params = line.split(",")
            cen = [float(params[0].split('(')[1]),float(params[1])]
            ellA = float(params[2])
            ellB = float(params[3])
            ellPA = float(params[4].split(')')[0])
            if ellA < ellB:
                ellA, ellB = ellB, ellA
                ellPA += 90    
            break
    f.close()

    f_galaxy = open(region_file.split('.reg')[0]+'_rot.reg', 'w')
    f_galaxy.write('%s\n' % ('image') )
    f_galaxy.write('ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # color=red\n' % (xc, yc, ellA, ellB, 0.))
    f_galaxy.close()    



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

def return_band_by_wavelength(wavelength):
    # TODO: Add other filters
    if wavelength>0.600 and wavelength<0.700:
        band = 'R'
    
    return band

def return_wavelength_by_band(band):
    # TODO: Add other filters
    if band=='r':
        wavelength = 0.620
    
    return wavelength

def show_steps(input_image, sky_med=None, add_sky=0., extract_dir = None):
    print(bcolors.OKBLUE + '\nSTEPS:' + bcolors.ENDC)
    if os.path.exists('results.pkl'):
        results = read_pickle()
    
    # -1
    print(bcolors.OKBLUE + '-1  ' + bcolors.ENDC + 'Show input image in saoimage ds9.')
    choose_step = -1
    
    # 0
    if os.path.exists('segm_ini.fits') or os.path.exists('best_stars.reg'):
        print(bcolors.OKBLUE + '0   ' + bcolors.ENDC + 'Run Sextractor and choose good PSF stars. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
        choose_step = 0
    else:
        print(bcolors.OKBLUE + '0   ' + bcolors.ENDC + 'Run Sextractor and choose good PSF stars.')
    
    
    # 1
    if os.path.exists('results.pkl'):
        try:
            if results['M0_AVERAGE [mag]'] is not None and results['M0_AVERAGE [mag]']!='nan':
                print(bcolors.OKBLUE + '1*  ' + bcolors.ENDC + 'Do photometric calibration based on step 0. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
                choose_step = 1
        except:
            print(bcolors.OKBLUE + '1*  ' + bcolors.ENDC + 'Do photometric calibration based on step 0.')
            
            
    # 2
    if os.path.exists('results.pkl'):
        try:
            zz = results['SB_deep [mag/arcsec^2]']
            print(bcolors.OKBLUE + '2   ' + bcolors.ENDC + 'Find some general characteristics of the image: pixel scale, sky background (mean, median and std) and deepness. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
            choose_step = 2
        except:
            print(bcolors.OKBLUE + '2   ' + bcolors.ENDC + 'Find some general characteristics of the image: pixel scale, sky background (mean, median and std) and deepness.')    

    
    # 3
    if os.path.exists('psf.fits'):
        print(bcolors.OKBLUE + '3   ' + bcolors.ENDC + 'Create small PSF image based on the unsaturated stars from step 0. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
        choose_step = 3
    else:
        print(bcolors.OKBLUE + '3   ' + bcolors.ENDC + 'Create small PSF image based on the unsaturated stars from step 0.')



    # 4
    if os.path.exists(input_image.split('.fits')[0] + '_crop.fits'):
        print(bcolors.OKBLUE + '4*  ' + bcolors.ENDC + 'Initial cropping to get square frame 5 times larger the galaxy. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
        choose_step = 4
    else:
        print(bcolors.OKBLUE + '4*  ' + bcolors.ENDC + 'Initial cropping to get square frame 5 times larger the galaxy.')
        
        
    # 5
    if os.path.exists('general_mask.reg'):
        print(bcolors.OKBLUE + '5*  ' + bcolors.ENDC + 'General masking. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC) 
        choose_step = 5
    else:
        print(bcolors.OKBLUE + '5*  ' + bcolors.ENDC + 'General masking.')        
        
        
    # 6
    if sky_med is None:
        if os.path.exists('sky_subtr.fits'):
            print(bcolors.OKBLUE + '6   ' + bcolors.ENDC + 'Sky fitting by a polynomial. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
            choose_step = 6
        else:
            print(bcolors.OKBLUE + '6   ' + bcolors.ENDC + 'Sky fitting by a polynomial.')
    else:
            print(bcolors.OKBLUE + '6   ' + bcolors.ENDC + 'Sky fitting by a polynomial. NOT REQUIRED as sky_med=0.')


    # 7
    if os.path.exists('galaxy_ellipse.reg'):
        print(bcolors.OKBLUE + '7   ' + bcolors.ENDC + 'Determine galaxy ellipse '+ bcolors.WARNING+'(skip if you run 8 next).  '+bcolors.ENDC + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
        #choose_step = 7
    else:
        print(bcolors.OKBLUE + '7   ' + bcolors.ENDC + 'Determine galaxy ellipse '+ bcolors.WARNING+'(Skip if you run 8 next).'+bcolors.ENDC)    
    
    # 8
    if os.path.exists('sky_subtr_galaxy.fits'):
        print(bcolors.OKBLUE + '8   ' + bcolors.ENDC + 'Determine sky within an annulus around the galaxy. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
        choose_step = 8
    else:
        print(bcolors.OKBLUE + '8   ' + bcolors.ENDC + 'Determine sky within an annulus around the galaxy.')


    # 9
    if os.path.exists('sky_subtr_galaxy_rot.fits') or os.path.exists('sky_subtr_rot.fits'):
        print(bcolors.OKBLUE + '9   ' + bcolors.ENDC + 'Galaxy rotation. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
        choose_step = 9
    else:
        print(bcolors.OKBLUE + '9   ' + bcolors.ENDC + 'Galaxy rotation.')
     
    
    # 10
    if os.path.exists('sky_subtr_galaxy_crop.fits') or os.path.exists('sky_subtr_galaxy_rot_crop.fits') or os.path.exists('sky_subtr_crop.fits') or os.path.exists('sky_subtr_rot_crop.fits') or os.path.exists('galaxy.fits'):    
        print(bcolors.OKBLUE + '10  ' + bcolors.ENDC + 'Final cropping. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
        choose_step = 10
    else:    
        print(bcolors.OKBLUE + '10  ' + bcolors.ENDC + 'Final cropping.')
    
    
    # 11
    if add_sky!=0.:
        if os.path.exists('galaxy_fix.fits'):
            print(bcolors.OKBLUE + '11  ' + bcolors.ENDC + 'Add constant sky level to the image. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC) 
            choose_step = 11
        else:
            print(bcolors.OKBLUE + '11  ' + bcolors.ENDC + 'Add constant sky level to the image.') 
    else:
            print(bcolors.OKBLUE + '11  ' + bcolors.ENDC + 'Add constant sky level to the image.' + bcolors.WARNING+' NOT APPLICABLE (change add_sky)'+bcolors.ENDC) 
            
    
    # 12
    if os.path.exists('mask.fits'):
        print(bcolors.OKBLUE + '12  ' + bcolors.ENDC + 'Determine galaxy ellipse for the cropped image and do final masking for it. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC) 
        choose_step = 12
    else:
        print(bcolors.OKBLUE + '12  ' + bcolors.ENDC + 'Determine galaxy ellipse for the cropped image and do final masking for it.') 
    
    
    
    # 13
    if os.path.exists('./sersic_fitting/galfit.01'):
        print(bcolors.OKBLUE + '13  ' + bcolors.ENDC + 'GALFIT fitting to a Sersic function. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
        choose_step = 13
    else:
        print(bcolors.OKBLUE + '13  ' + bcolors.ENDC + 'GALFIT fitting to a Sersic function.')
        
        
    # 14
    if os.path.exists('ellipse.txt'):
        print(bcolors.OKBLUE + '14  ' + bcolors.ENDC + 'Run IRAF/ELLIPSE fitting. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
        choose_step = 14
    else:
        print(bcolors.OKBLUE + '14  ' + bcolors.ENDC + 'Run IRAF/ELLIPSE fitting.')
        
        
    # 15 
    try:
        zz = results['asympt_mag [mag]']
        print(bcolors.OKBLUE + '15  ' + bcolors.ENDC + 'Calculate general galaxy parameters: total magnitude, limiting radius, effective radius etc. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
        choose_step = 15
    except:
        print(bcolors.OKBLUE + '15  ' + bcolors.ENDC + 'Calculate general galaxy parameters: total magnitude, limiting radius, effective radius etc.')
    
    
    # 16
    if os.path.exists('results.pkl'):
        print(bcolors.OKBLUE + '16  ' + bcolors.ENDC + 'Display all collected values from results.pkl')
    else:
        print(bcolors.OKBLUE + '16  ' + bcolors.ENDC + 'Display all collected values from results.pkl. NOT APPLICABLE')


    # 17
    if extract_dir is not None:
        print(bcolors.OKBLUE + '17  ' + bcolors.ENDC + 'Copy all important output to %s' % (extract_dir))
    else:
        print(bcolors.OKBLUE + '17  ' + bcolors.ENDC + 'Copy output.'  + bcolors.WARNING+' NOT APPLICABLE (specify extract_dir)'+bcolors.ENDC)

    # 18
    if os.path.exists('galaxy_resample.fits'):
        print(bcolors.OKBLUE + '18  ' + bcolors.ENDC + 'Resampling to a different frame. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
        choose_step = 18
    else:
        print(bcolors.OKBLUE + '18  ' + bcolors.ENDC + 'Resampling to a different frame.')    
        
    
    # 19
    if os.path.exists('galaxy_rebin.fits'):
        print(bcolors.OKBLUE + '19  ' + bcolors.ENDC + 'Rebinning to a different frame. ' + bcolors.OKGREEN+'DONE'+bcolors.ENDC)
        choose_step = 19
    else:
        print(bcolors.OKBLUE + '19  ' + bcolors.ENDC + 'Rebinning to a different frame.')    
    
    
    # 19
        print(bcolors.OKBLUE + '20  ' + bcolors.ENDC + 'Create deep filtered picture.')        
    
    
    
    
    
    
    
    
    
    print(bcolors.OKBLUE + 'ANY LETTER  ' + bcolors.ENDC + 'EXIT'  + bcolors.WARNING + bcolors.ENDC)
    print('-----------------------------------')
    print('* OPTIONAL')
    print('NOTE: Output magnitudes and SB are corrected for Galactic extinction, except for mtot_sersic.')
    
    if choose_step+1==11 and add_sky==0.:
        return str(choose_step+2)
    else:
        return str(choose_step+1)


def get_galaxy_center(input_image, galaxy_name=None, RA=None, DEC=None):
        # Determine galaxy center:
        if galaxy_name is not None or (RA is not None and DEC is not None):
            xc,yc = get_galaxy_center_ned.main(input_image, name=galaxy_name, RA=RA, DEC=DEC)
            return xc,yc
        else:
            if verbosity: print('Specify galaxy name or RA and DEC! Exiting.')
            exit()

def check_nan(par):
    try:
        if np.isnan(par):
            return None
    except:
        z=1
    
    if par=='nan':
        return None
    else:
        return par


def choose_image_to_do(images):
    for image in images:
        if os.path.exists(image):
            return image
    
    return None


def main_do(directory,
         galaxy_name=None,
         input_image='imp-image.fits',
         sky_degree=1,
         steps=None,
         SB_lim=28.,
         add_sky=0.,
         bin_box=1,
         wavelength=0.620,
         m0 = None, m0_err = None, scale=None,
         RA=None, DEC=None,
         sky_med=None,
         user_interact=True,
         calibr_band = 'r',
         extract_dir = None,
         exptime=1.,
         Ncombine=1,
         gain='nan',
         ron='nan',
         rebin_to_galaxy_directory='nan',
         resample_to_galaxy_directory='nan', 
         verbosity=True
         ):

    # Go to the galaxy preparation directory
    os.chdir(directory)

    if not os.path.exists('results.pkl'):
        # Create pickle file where all output will be saved
        write_pickle(results=None)

    # Check parameters which have been given:
    galaxy_name = check_nan(galaxy_name)
    sky_degree = check_nan(sky_degree)
    SB_lim = check_nan(SB_lim)
    wavelength = check_nan(wavelength)
    calibr_band = check_nan(calibr_band)
    m0 = check_nan(m0)
    m0_err = check_nan(m0_err)
    scale = check_nan(scale)
    sky_med = check_nan(sky_med)
    extract_dir = check_nan(extract_dir)
    exptime = check_nan(exptime)
    Ncombine = check_nan(Ncombine)
    gain = check_nan(gain)
    ron = check_nan(ron)
    

    
    # and write them:
    results = read_pickle()
    #results['scale [arcsec/pix]'] = scale 
    #results['wavelength (um)'] = wavelength
    #results['M0_AVERAGE [mag]'] = m0
    #results['M0_STD [mag]'] = m0_err

    if m0 is not None:
        results['M0_AVERAGE [mag]'] = float(m0)    
    write_pickle(results)

    if steps is not None:
        steps = np.array(np.array(steps, dtype=float), dtype=int)
 
    # Determine band by given wavelength:
    #band = return_band_by_wavelength(wavelength)
    
    
    
    if steps is None:
        choose_step = show_steps(input_image, sky_med=sky_med, add_sky=add_sky, extract_dir=extract_dir)
        step = str(input('\nChoose step(s) from the list above (separated by comma) ' + bcolors.OKBLUE + '[%s] ' % (choose_step) + bcolors.ENDC) or choose_step)
        for sstep in step:
            if sstep.isalpha():
                exit()
        
        steps = np.array(step.split(','), int) 

       


    # If sky_med=0. then no sky subtraction will be done at all and we just copy input file to sky_subtr.fits:
    if sky_med==0.:
        shutil.copy(input_image, 'sky_subtr.fits')
        if verbosity: print('WARNING: Input image is considered as sky-subtracted (sky_med=0.) => sky_subtr.fits')



    if -1 in steps:
        # Show the input image
        subprocess.call("ds9 %s -scale histequ" % (input_image), shell=True)

    
    if 0 in steps:
        # Find good PSF stars by running Sextractor (creates catalog and segmentation map) 
        results = read_pickle()
        fwhm_med, ell_med, theta_med = IMP_main(input_image, 'find_good_stars', band='R', catalogue='NOMAD', user_interact=user_interact, verbosity=verbosity)
        results['PSF_FWHM [pix]'] = fwhm_med
        results['PSF_ell'] = ell_med
        results['PSF_PA [deg,CCW/x]'] = theta_med
        write_pickle(results)        
        
        
        
    if 1 in steps:        
        # Do photometric calibration using step 0: OPTIONAL
        results = read_pickle()
        if not os.path.exists('segm_ini.fits') or not os.path.exists('best_stars.reg'):
            if verbosity: print('Run step 0 first! Exiting.')
            exit()
        res = IMP_main(input_image, 'phot_calibration', region_file='best_stars.reg', mask_image='segm_ini.fits', band=['r','r','r'], catalogue=['I/345/gaia2','V/139','II/349/ps1'], wavelength=wavelength, user_interact=user_interact, verbosity=verbosity)
        
        results['M0_AVERAGE [mag]'] = np.nanmean([res['M0:V/139:r:None'],res['M0:II/349/ps1:r:None'],res['M0:I/345/gaia2:r:None']])
        results['M0_STD [mag]'] = np.nanmean([res['M0_STD:V/139:r:None'],res['M0_STD:II/349/ps1:r:None'],res['M0_STD:I/345/gaia2:r:None']]) 
        results['M0:V/139:r:None [mag]'] = res['M0:V/139:r:None']
        results['M0:II/349/ps1:r:None [mag]'] = res['M0:II/349/ps1:r:None']
        results['M0:I/345/gaia2:r:None [mag]'] = res['M0:I/345/gaia2:r:None']
        results['M0_STD:V/139:r:None [mag]'] = res['M0_STD:V/139:r:None']
        results['M0_STD:II/349/ps1:r:None [mag]'] = res['M0_STD:II/349/ps1:r:None']
        results['M0_STD:I/345/gaia2:r:None [mag]'] = res['M0_STD:I/345/gaia2:r:None']

        if m0 is not None:
            results['M0_AVERAGE [mag]'] = float(m0)
            
        write_pickle(results)  


    if 2 in steps:
        # Get some info about the image
        results = read_pickle()
        
        if not os.path.exists('segm_ini.fits'):
            if verbosity: print('Run step 0 first as no segmentation image found! Exiting.')
            exit()
        try:
            ZP = results['M0_AVERAGE [mag]']
        except:
            if m0 is not None:
                ZP = m0
            else:
                if verbosity: print('No m0 given. Specify m0 or run step 1 first! Exiting.')
                exit()                
            

        limit, median_MEAN, median_MEDIAN, std_MEDIAN, median_STD, scale = IMP_main(input_image, 'ima_stats', m0=ZP, mask_image='segm_ini.fits', user_interact=user_interact, verbosity=verbosity)
        
        results['SB_deep [mag/arcsec^2]'] = limit
        results['Sky_mean [ADU]'] = median_MEAN
        results['Sky_median [ADU]'] = median_MEDIAN
        results['Sky_std_median [ADU]'] = std_MEDIAN
        results['Sky_std [ADU]'] = median_STD 
        results['scale [arcsec/pix]'] = float(scale) 
        results['wavelength (um)'] = float(wavelength)

        
        write_pickle(results)


    if 3 in steps:
        # Create small PSF image based on unsaturated stars from step 0
        # TODO: Add PSFex to find more precise PSF dependent on the position on the frame
        if not os.path.exists('best_stars.reg'):
            if verbosity: print('Run step 0 first! Exiting.')
            exit()
            
        IMP_main(input_image, 'crea_psf', region_file='best_stars.reg', psf_image='psf.fits', user_interact = user_interact, number=None, verbosity=verbosity) # You can choose your own number of PSF star if needed


    if 4 in steps:
        results = read_pickle()

        # Initial cropping: 
        # Determine galaxy center:
        xc,yc = get_galaxy_center_ned.convert_radec_to_image(input_image, RA, DEC) #get_galaxy_center(input_image, galaxy_name=galaxy_name, RA=RA, DEC=DEC)


        if not os.path.exists('segm_ini.fits'):
            if verbosity: print('Run step 0 first as no segmentation image found! Exiting.')
            exit()

        #Determine galaxy ellipse
        try:
            sky_mean = results['Sky_median [ADU]']
        except:
            if sky_med is not None:
                sky_mean = sky_med
            else:
                if verbosity: print('Specify sky_med or run step 2 first')


        # Rebinning to a different frame
        
        if rebin_to_galaxy_directory!='nan': 
            if not os.path.exists(rebin_to_galaxy_directory + '/' + input_image.split('.fits')[0] + '_crop.fits'):
                print('%s does not exist! Exiting...' % (rebin_to_galaxy_directory + '/' + input_image.split('.fits')[0] + '_crop.fits'))
                exit()
            rebin_image.rebin(rebin_to_galaxy_directory + '/' + input_image.split('.fits')[0] + '_crop.fits', input_image, output_image=input_image.split('.fits')[0] + '_crop.fits', hdu_ref=0, hdu_inp=0, preserve_bad_pixels=True)
            shutil.copy(rebin_to_galaxy_directory + '/galaxy_ellipse_wcs.reg', 'galaxy_ellipse_wcs.reg')
        else:        
            xc_ell,yc_ell,sma_ell,smb_ell,PA_ell = IMP_main(input_image, 'determine_galaxy', xc=xc, yc=yc, mask_image='segm_ini.fits', sky_mean=sky_mean, smb=5., user_interact=user_interact, region_file='galaxy_ellipse.reg', verbosity=verbosity)  #### WARNING: If code crashes here, maybe it's because it cannot find object. Change smth!
        
            convert_coord_regions.convert_ellipse_image_to_wcs(input_image, 'galaxy_ellipse.reg', output_region_file='galaxy_ellipse_wcs.reg')

            # Do initial cropping
            IMP_main(input_image, 'crop_galaxy', region_file='galaxy_ellipse.reg', offset_size=5., crop_box=True, user_interact=user_interact, verbosity=verbosity) # Change offset_size if you want to change the size of the output image
        


    if 5 in steps:
        # If cropped image exists, then we work with it
        if os.path.exists(input_image.split('.fits')[0] + '_crop.fits'):
            input_image = input_image.split('.fits')[0] + '_crop.fits'
        
        # Do general masking of the image
        mto_gain = IMP_main(input_image, 'image_masking', region_file = 'general_mask.reg', mask_image='general_mask.fits', user_interact=user_interact, verbosity=verbosity, masking_tool='sextractor', offset_size=1.5)
        
        results = read_pickle()
        if gain is None:
            results['GAIN'] = mto_gain
        else:
            results['GAIN'] = gain
        write_pickle(results)

    if 6 in steps:
        # Do overall frame sky subtraction  
        results = read_pickle()

        # If cropped image exists, then we work with it
        if os.path.exists(input_image.split('.fits')[0] + '_crop.fits'):
            input_image = input_image.split('.fits')[0] + '_crop.fits'
            
        output_image, mean, median, std, sky_degree = IMP_main(input_image, 'sky_correction', mask_image='general_mask.fits', output_image='sky_subtr.fits', sky_degree=sky_degree, user_interact=user_interact, verbosity=verbosity) # NOTE: sky_degree!!!

        results['Sky_subtr_mean [ADU]'] = mean
        results['Sky_subtr_median [ADU]'] = median
        results['Sky_subtr_std [ADU]'] = std
        results['Sky_fit_polynom'] = sky_degree
                
        write_pickle(results)

    if 7 in steps or 8 in steps:
            
        if not os.path.exists('sky_subtr.fits'):
            if verbosity: print('Sky subtraction was not done! Run step 6 first. Exiting.')
            exit()            

        if not os.path.exists('general_mask.fits'):
            if verbosity: print('Run step 5 first. Exiting.')
            exit()   

        # Determine galaxy center:
        xc,yc = get_galaxy_center_ned.convert_radec_to_image('sky_subtr.fits', RA, DEC) #get_galaxy_center('sky_subtr.fits', galaxy_name=galaxy_name, RA=RA, DEC=DEC)

        
        if not os.path.exists('galaxy_ellipse.reg'):
            # Determine galaxy ellipse
            xc_ell,yc_ell,sma_ell,smb_ell,PA_ell = IMP_main('sky_subtr.fits', 'determine_galaxy', xc=xc, yc=yc, mask_image='segm_ini.fits', smb=5., user_interact=user_interact, region_file='galaxy_ellipse.reg', verbosity=verbosity) #### WARNING: If code crashes here, maybe it's because it cannot find object. Change smb!
            convert_coord_regions.convert_ellipse_image_to_wcs('sky_subtr.fits', 'galaxy_ellipse.reg', output_region_file='galaxy_ellipse_wcs.reg')
        else:
            convert_coord_regions.convert_ellipse_wcs_to_image('sky_subtr.fits', 'galaxy_ellipse_wcs.reg', output_region_file='galaxy_ellipse.reg')


    if 8 in steps:
        results = read_pickle()
        
        # Do galaxy sky subtraction (within an annulus) - more accurate:
        sky,std,bkg_median_std = IMP_main('sky_subtr.fits', 'sky_galaxy', output_image='sky_subtr_galaxy.fits', mask_image='general_mask.fits', region_file='galaxy_ellipse.reg', user_interact=user_interact, annulus_width=32, offset_size = 1.5, verbosity=verbosity) #### WARNING: Change annulus_width and offset_size
        results['Sky_subtr_annulus_median [ADU]'] = sky
        results['Sky_subtr_annulus_std [ADU]'] = std
        results['Sky_subtr_annulus_median_std [ADU]'] = bkg_median_std ################ UNCERTAINTY OF THE SKY
        write_pickle(results)


    if 9 in steps:
        results = read_pickle()
        
        do_image = choose_image_to_do(['sky_subtr_galaxy.fits','sky_subtr.fits'])

        if do_image == 'sky_subtr_galaxy.fits':
            sky_std=results['Sky_subtr_annulus_std [ADU]']
        elif do_image == 'sky_subtr.fits':
            sky_std = results['Sky_subtr_std [ADU]']
        elif do_image is None:
            if verbosity: print('Sky subtraction should be done first. Exiting...')
            exit()

        if not os.path.exists('general_mask.fits'):
            if verbosity: print('Run step 5 first. Exiting.')
            exit() 
  
        
        # Determine galaxy center:
        xc,yc = get_galaxy_center_ned.convert_radec_to_image(do_image, RA, DEC) #get_galaxy_center(input_image, galaxy_name=galaxy_name, RA=RA, DEC=DEC)

        
        # Do rotation:
        IMP_main(do_image, 'rotate_galaxy', xc=xc, yc=yc, mask_image='general_mask.fits', sky_std=sky_std, user_interact=user_interact, verbosity=verbosity)
        


    if 10 in steps:
        # Crop galaxy image
        do_image = choose_image_to_do(['sky_subtr_galaxy_rot.fits','sky_subtr_galaxy.fits','sky_subtr_rot.fits','sky_subtr.fits'])
        
        if rebin_to_galaxy_directory!='nan': 
            if not os.path.exists(rebin_to_galaxy_directory + '/galaxy.fits'):
                print('%s does not exist! Exiting...' % ('galaxy.fits'))
                exit()
            rebin_image.rebin(rebin_to_galaxy_directory + '/galaxy.fits', do_image, output_image='galaxy.fits', hdu_ref=0, hdu_inp=0, preserve_bad_pixels=True)
            #convert_coord_regions.convert_ellipse_wcs_to_image(rebin_to_galaxy_directory + '/galaxy.fits', rebin_to_galaxy_directory + '/galaxy_ellipse_wcs.reg', output_region_file='galaxy_ellipse.reg')
        else:
            # Do final cropping
            if os.path.exists(do_image):
                # Convert galaxy_ellipse.reg
                #convert_region_ellipse(galaxy_name, RA, DEC, 'sky_subtr_galaxy_rot.fits', 'galaxy_ellipse.reg')
                convert_coord_regions.convert_ellipse_wcs_to_image(do_image, 'galaxy_ellipse_wcs.reg', output_region_file='galaxy_ellipse.reg')
                
                output_image = IMP_main(do_image, 'crop_galaxy', region_file='galaxy_ellipse.reg', offset_size=2., user_interact=user_interact, verbosity=verbosity)
                shutil.copy(output_image, 'galaxy.fits')            
        

    if 11 in steps and add_sky!=0.:
        # Fix galaxy --- add sky:
        arithm_operations.main('galaxy.fits', add_sky, 'add', 'galaxy_fix.fits')
   
   
    if 12 in steps:
        results = read_pickle()

        if os.path.exists('galaxy_fix.fits') and add_sky!=0.:
            input_image = 'galaxy_fix.fits'
        else:
            input_image = 'galaxy.fits'
        
        # Convert galaxy ellipse to the cropped image
        #if not os.path.exists('galaxy_ellipse.reg'):
        #    if verbosity: print('Galaxy ellipse is not determined. Exiting.')
        #xc,yc,sma,smb,PA = crop_galaxy_image.read_region('galaxy_ellipse.reg')
        #if os.path.exists('galaxy_ellipse_rot.reg'):
        #    PA = 0.
        
        #xc,yc = get_galaxy_center(input_image, galaxy_name=galaxy_name, RA=RA, DEC=DEC)
            
        #get_galaxy_ellipse.create_ellipse_region(input_image, xc, yc, sma, smb, PA, 'galaxy_ellipse_final.reg')



        if os.path.exists('sky_subtr_galaxy.fits'):
            sky_std=results['Sky_subtr_annulus_std [ADU]']
        else:
            sky_std = results['Sky_subtr_std [ADU]']

        if rebin_to_galaxy_directory=='nan':
            convert_coord_regions.convert_ellipse_wcs_to_image(input_image, 'galaxy_ellipse_wcs.reg', output_region_file='galaxy_ellipse.reg')
        else:
            convert_coord_regions.convert_ellipse_wcs_to_image(rebin_to_galaxy_directory + '/galaxy.fits', rebin_to_galaxy_directory + '/galaxy_ellipse_wcs.reg', output_region_file='galaxy_ellipse.reg')
        
        if user_interact:
            subprocess.call("ds9 %s -scale histequ -regions %s" % (input_image, 'galaxy_ellipse.reg'), shell=True)
        
        # Do final masking: inner and outer
        IMP_main(input_image, 'final_masking', region_file='galaxy_ellipse.reg', psf_image='psf.fits', output_image='galaxy_clean.fits', mask_image='mask.fits',sky_std=sky_std, user_interact=user_interact, verbosity=verbosity) # 'galaxy_mask.reg' will be created        
        

    if 13 in steps:
        results = read_pickle()

        if os.path.exists('galaxy_fix.fits') and add_sky!=0.:
            input_image = 'galaxy_fix.fits'
        else:
            input_image = 'galaxy.fits'
        
        #convert_reg_to_mask.mask('galaxy.fits', 'galaxy_mask.reg', output_image='galaxy_clean.fits', output_mask='mask.fits', mask_value=1, show_running=True, mask_DN=None)

        # Prepare image for GALFIT
        xc,yc = get_galaxy_center_ned.convert_radec_to_image(input_image, RA, DEC) #get_galaxy_center(input_image, galaxy_name=galaxy_name, RA=RA, DEC=DEC)
        
        
        shutil.copy(input_image, 'galaxy_clean_galf.fits')
        m0,pix2sec = prep_for_galfit(xc, yc)

        
        # Do GALFIT fitting
        dec_path = './sersic_fitting'
        res = IMP_main('galaxy_clean_galf.fits', 'fitting', psf_image='psf.fits', mask_image='mask.fits', m0=m0, scale=pix2sec, out_directory=dec_path, user_interact=user_interact, verbosity=verbosity)
        xc,yc,mtot,re,n,q,PA,C0 = res
        results['xc_sersic [pix]'] = xc
        results['yc_sersic [pix]'] = yc
        results['mtot_sersic [mag]'] = mtot
        results['re_sersic [pix]'] = re
        results['n_sersic'] = n
        results['q_sersic'] = q
        results['PA_sersic [deg,CCW/y]'] = PA
        results['C0_sersic'] = C0
        write_pickle(results)

    if 14 in steps:
        results = read_pickle()
        if os.path.exists('galaxy_fix.fits') and add_sky!=0.:
            input_image = 'galaxy_fix.fits'
        else:
            input_image = 'galaxy.fits'
        if verbosity: print(RA,DEC)    
        xc,yc = get_galaxy_center_ned.convert_radec_to_image(input_image, RA, DEC) #get_galaxy_center(input_image, galaxy_name=galaxy_name, RA=RA, DEC=DEC)
        
        try:
            ZP = results['M0_AVERAGE [mag]']
        except:
            if m0 is not None:
                ZP = m0
            else:
                if verbosity: print('No m0 given. Specify m0 or run step 1 first! Exiting.')
                exit()     

        try:
            pix2sec = results['scale [arcsec/pix]']
        except:
            if scale is not None:
                pix2sec = scale
            else:
                if verbosity: print('No scale given! Exiting.')
                exit() 
                
        
        try:
            sma = 10.#results['re_sersic [pix]']
            smb = results['q_sersic'] * sma
            PA = results['PA_sersic [deg,CCW/y]']
        except:
            sma = 10.
            smb = 10.*0.8
            PA = 20.
        
        # IRAF/ELLIPSE fitting
        IMP_main(input_image, 'ellipse_fitting', xc=xc, yc=yc, sma=sma, smb=smb, PA=PA, mask_image='mask.fits', m0=ZP, scale=pix2sec, step=1., linear='yes', user_interact=user_interact, verbosity=verbosity) # step=0.03

    if 15 in steps:
        results = read_pickle()
        if os.path.exists('galaxy_fix.fits') and add_sky!=0.:
            input_image = 'galaxy_fix.fits'
        else:
            input_image = 'galaxy.fits'
        
        if galaxy_name is not None:
            result_table = ned.Ned.query_object(galaxy_name)
            try:
                RA,DEC = float(result_table['RA(deg)']),float(result_table['DEC(deg)']) # an astropy.table.Table
            except:
                RA,DEC = float(result_table['RA']),float(result_table['DEC']) # an astropy.table.Table

        if os.path.exists('sky_subtr_galaxy.fits'):
            sky_std = results['Sky_subtr_annulus_std [ADU]']
            sky_med_unc = results['Sky_subtr_annulus_median_std [ADU]']
        else:
            sky_std = results['Sky_subtr_std [ADU]']
            sky_med_unc = results['Sky_std_median [ADU]']
        
        try:
            sigma_cal = results['M0_STD [mag]']
        except:
            if m0_err is not None:
                sigma_cal = m0_err
            else:
                if verbosity: print('No m0_err is given. Specify m0_err or run step 1 first! Exiting.')
                exit()  

        try:
            ZP = results['M0_AVERAGE [mag]']
        except:
            if m0 is not None:
                ZP = m0
            else:
                if verbosity: print('No m0 given. Specify m0 or run step 1 first! Exiting.')
                exit()  

        try:
            pix2sec = results['scale [arcsec/pix]']
        except:
            if scale is not None:
                pix2sec = scale
            else:
                if verbosity: print('No scale given! Exiting.')
                exit() 
        


        res, Aext = IMP_main(input_image, 'general_pars', m0=ZP, scale=pix2sec, sky_std=sky_std, sigma_cal=sigma_cal,SB_lim=SB_lim, ellipse_file='ellipse.txt',step=bin_box, sky_med_unc=sky_med_unc, RA=RA, DEC=DEC, wavelength=return_wavelength_by_band(calibr_band), user_interact=user_interact, verbosity=verbosity)
        
        results['Aext [mag]'] = Aext
        results['min_radius [pix]'] = res['min_radius']
        results['min_SB [mag/arcsec^2]'] = res['min_SB']
        results['lim_radius [pix]'] = res['lim_radius']
        results['lim_radius_err [pix]'] = res['lim_radius_err']
        results['lim_mag [mag]'] = res['lim_mag']
        results['asympt_radius [pix]'] = res['asympt_radius']
        results['asympt_mag [mag]'] = res['asympt_mag']
        results['eff_radius [pix]'] = res['eff_radius']
        results['galaxy_aperture (x0,y0,ell,PA)'] = res['galaxy_aperture']
        try:
            results['mtot_sersic_cor [mag]'] = results['mtot_sersic [mag]'] - Aext
        except:
            zz=1
        write_pickle(results)
        
    if 16 in steps:
        results = read_pickle()
        for key, value in results.items(): 
            if verbosity: print('%s: %s' % (key, str(value))) 


    if 17 in steps:
        # Extract most important output data:
        if verbosity: print('Copying important data...')
        if extract_dir is None:
            if verbosity: print('Error. Please specify extract_dir. Exiting.')
        else:
            files = ['results.pkl','galaxy.fits','galaxy_fix.fits','galaxy_clean.fits','galaxy_mask.reg','psf.fits','iraf.png','azim_aver_28.0.png','sersic_fitting/galfit.01']
            
            if not os.path.exists(extract_dir):
                os.mkdir(extract_dir)

            if not os.path.exists(extract_dir + '/%s' % (directory.split('/')[-1])):
                os.mkdir(extract_dir + '/%s' % (directory.split('/')[-1]))
            
            for file in files:
                if os.path.exists(file):
                    shutil.copy(file, extract_dir + '/%s/%s' % (directory.split('/')[-1], file.split('/')[-1]))
            
    if 18 in steps:            
        # Resampling to a different frame    
        if os.path.exists(resample_to_galaxy_directory + '/galaxy.fits'):
            #TODO!
            zz=1
            
    if 19 in steps:
        # Rebinning to a different frame
        if os.path.exists(rebin_to_galaxy_directory + '/galaxy.fits'):
            rebin_image.rebin(rebin_to_galaxy_directory + '/galaxy.fits', 'galaxy.fits', output_image='galaxy_rebin.fits', hdu_ref=0, hdu_inp=0, preserve_bad_pixels=True)            
            
    if 20 in steps:
        # Create a deep smooth image
        
        results = read_pickle()

        try:
            ZP = results['M0_AVERAGE [mag]']
        except:
            if m0 is not None:
                ZP = m0
            else:
                if verbosity: print('No m0 given. Specify m0 or run step 1 first! Exiting.')
                exit()  

        try:
            pix2sec = results['scale [arcsec/pix]']
        except:
            if scale is not None:
                pix2sec = scale
            else:
                if verbosity: print('No scale given! Exiting.')
                exit() 

        IMP_main('galaxy.fits', 'smooth_picture', mask_image='mask.fits', m0=ZP, scale=pix2sec, sigma_smooth = 1., verbosity=verbosity)    
            
            
            
            
            
            

def get_input_output_directories(sample_file, verbosity):
        f = open(sample_file, 'r')
        lines = f.readlines()
        
        if 'directory_from' in lines[0]:
            ss = lines[0].replace(" ", "")
            directory_from = ss.split('directory_from=')[-1].split('\n')[0].strip()     
            if directory_from != '':
                if not os.path.exists(directory_from):
                    if verbosity: print('\033[91m ERROR! Input directory %s is not found! Exiting... \033[0m ' % (directory_from))
                    exit()
        else:
            directory_from = ''

        if 'directory_to' in lines[1]:
            ss = lines[1].replace(" ", "")
            directory_to = ss.split('directory_to=')[-1].split('\n')[0].strip()
        else:
            directory_to = os.getcwd() + '/data_preparation'


        if 'directory_res' in lines[2]:
            ss = lines[2].replace(" ", "")
            directory_res = ss.split('directory_res=')[-1].split('\n')[0].strip()
        else:
            directory_res = None
            
        f.close()
        
        if verbosity: print('\n\033[92m --------------------- \033[0m ')
        if verbosity: print('\033[92m Input directory:     \033[0m %s ' % (directory_from))
        if verbosity: print('\033[92m Output directory:    \033[0m %s ' % (directory_to))
        if verbosity: print('\033[92m Resultant directory: \033[0m %s ' % (directory_res))
        if verbosity: print('\033[92m --------------------- \033[0m \n\n')
        
        return directory_from, directory_to, directory_res



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

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")




def main(number, sample_file, verbosity=True):
    # Greeting:
    if verbosity: print('\n\033[94m ****  GALAXY IMAGE PREPARATION   **** \033[0m')
    if verbosity: print('\033[94m **** Aleksandr Mosenkov, 2019 (C)**** \033[0m')
    
    # Get current directory from where we run the script:
    current_dir = os.getcwd()
    
    # Read in the input file with the sample:
    directory_from, directory_to, directory_res = get_input_output_directories(sample_file, verbosity)


    # Create output directory:
    if not os.path.exists(directory_to):
        os.mkdir(directory_to)
    
    # Create res directory:
    if directory_res is not None:
        if not os.path.exists(directory_res):
            os.mkdir(directory_res)


    # Reading the sample table:
    data,units = read_data.main(sample_file, delimiter='\t', header_line=2, units_line=None, skip_lines = [], comment='#')

    # Mandatory parameters:
    Number = np.array(data['Number'], dtype=int)
    galaxy_names = np.array(data['Obj_name'], dtype=str)
    galaxy_frames = np.array(data['File'], dtype=str)
 

    # Optional parameters:
    try:
        rebin_to = np.array(data['Rebin_to'], dtype=str)
    except:
        rebin_to = np.array(len(Number)*['nan'], dtype=str)
    

    try:
        resample_to = np.array(data['Resample_to'], dtype=str)
    except:
        resample_to = np.array(len(Number)*['nan'], dtype=str)

    try:
        X = np.array(data['X[pix]'], dtype=str)
    except:
        X = np.array(len(Number)*['nan'], dtype=str)


    try:
        Y = np.array(data['Y[pix]'], dtype=str)
    except:
        Y = np.array(len(Number)*['nan'], dtype=str)

    try:
        RA = np.array(data['Ra[deg]'], dtype=str)
    except:
        RA = np.array(len(Number)*['nan'], dtype=str)



    try:
        DEC = np.array(data['Dec[deg]'], dtype=str)
    except:
        DEC = np.array(len(Number)*['nan'], dtype=str)
        

    try:
        sky_degrees = np.array(data['Sky_polynomial'], dtype=int)
    except:
        sky_degrees = np.array(len(Number)*[1], dtype=int)
    
    try:
        stepss = np.array(data['Steps'], dtype=str)
    except:
        stepss = np.array(len(Number)*['nan'], dtype=str)        
    
    steps = []
    for j in range(len(stepss)):
        if stepss[j]!='nan':
            steps.append(stepss[j].split(','))
        else:
            steps.append(None)

    steps = np.array(steps)

    
    try:
        user_inter = np.array(data['User_interact'])
        user_interact = []
        for ii in range(len(user_inter)):
            user_interact.append(str2bool(user_inter[ii]))
        
        user_interact = np.array(user_interact, dtype=bool)
    except:
        user_interact = np.array(len(Number)*[True], dtype=bool)

    try:
        exptime = np.array(data['Exptime[sec]'], dtype=float)
    except:
        exptime = np.array(len(Number)*[1.], dtype=float)
        
    try:
        Ncombine = np.array(data['Ncombine'], dtype=float)
    except:
        Ncombine = np.array(len(Number)*[1.], dtype=float)
    
    try:
        Sblim = np.array(data['Sblim[mag/arcsec2]'], dtype=float)
    except:
        Sblim = np.array(len(Number)*[28.], dtype=float)
    
    try:
        m0 = np.array(data['Zero-point[mag]'], dtype=str)
    except:
        m0 = np.array(len(Number)*['nan'], dtype=str)
    
    try:
        pix2sec = np.array(data['Scale[arcsec/pix]'], dtype=str)
    except:
        pix2sec = np.array(len(Number)*['nan'], dtype=str)
        
    try:
        gain = np.array(data['Gain[e-/ADU]'], dtype=str)
    except:
        gain = np.array(len(Number)*['nan'], dtype=str)

    try:
        ron = np.array(data['Ron[e-]'], dtype=str)
    except:
        ron = np.array(len(Number)*['nan'], dtype=str)
        
    
    try:
        wavelengths = np.array(data['Wavelength[um]'], dtype=str)
    except:
        wavelengths = np.array(len(Number)*['nan'], dtype=str)
        

    try:        
        calibr_band = np.array(data['Calibr_band'], dtype=str)
    except:
        calibr_band = np.array(len(Number)*['nan'], dtype=str)


    # Work with each galaxy in the sample:
    if number in Number:
            # Now we work with a specific item (number) from the input table:
            k = number - 1
            
            galaxy_name = galaxy_names[k]
            
            
            print('\n\033[92m ********GALAXY #%i: %s (%s)******** \033[0m' % (number, galaxy_name, galaxy_frames[k]) )

            # We need to specify wcs galaxy coordinates
            if RA[k] == 'nan' or DEC[k] == 'nan':
                ra,dec = get_galaxy_center_ned.get_coords_by_name(galaxy_name)
                if np.isnan(ra) or np.isnan(dec):
                    if verbosity: print('\033[91m ERROR! Ra and Dec are nan. Exiting... \033[0m')
                    exit()
            else:
                ra = float(RA[k])
                dec = float(DEC[k])
            


            if os.path.isdir(galaxy_frames[k]) and os.path.exists(galaxy_frames[k]) and os.path.exists(galaxy_frames[k]+'/results.pkl'):
                # THIS IS FOR CREATUNG SUBFRAMES (E.G. COMPACT GROUPS)
                # If a frame has several objects under study, this option can be used.
                # First you have to do some basic steps to prepare the frame (e.g. calibrarion, sky subtraction, intial cropping)
                # Then you can add next line to the sample file: galaxy_name path_to_dir_with_done_preparation
                # In this case new directory for each object will be created within the parent directory and all files
                # will be copied in there. The script will process each galaxy individually.
                galaxy_directory = galaxy_frames[k] + '/%s' % (galaxy_name)
                if not os.path.exists(galaxy_directory):
                    # Create special directory for each galaxy frame:
                    all_files = glob.glob(galaxy_frames[k] + '/*')
                    os.mkdir(galaxy_directory) 
                    # Copy all files from the parent directory
                    for file in all_files:
                        try:
                            check_fits_layers(file, verbosity)
                            shutil.copy(file, galaxy_directory + '/%s' % (file.split('/')[-1]))
                        except:
                            zz=1

                if os.path.exists(galaxy_directory + '/%s' % ('sky_subtr_galaxy.fits')):
                    input_image = 'sky_subtr_galaxy.fits'
                else:
                    input_image = 'sky_subtr.fits'             
            else:
                # THIS IS FOR REGULAR WORKING 
                galaxy_directory = directory_to+'/%s_%s' % (galaxy_name, galaxy_frames[k].split('.fi')[0].split('/')[-1])
                
                if not os.path.exists(galaxy_directory):
                    # Create special directory for each galaxy frame:
                    os.mkdir(galaxy_directory)
                    
                if not os.path.exists(galaxy_directory+'/imp-image.fits'):
                    if True:
                        # Check that the input image has one layer, otherwise the error occures.
                        check_fits_layers(directory_from+'/%s' % (galaxy_frames[k]), verbosity)
                        
                        # Copy the input file to the galaxy directory
                        shutil.copy(directory_from+'/%s' % (galaxy_frames[k]), galaxy_directory+'/imp-image.fits')
                        
                        # Go into the galaxy directory
                        os.chdir(galaxy_directory)
                        
                        # Check if wcs is found in the image:
                        if not wcs_existence('imp-image.fits', verbosity):
                            if np.isnan(float(X[k])) or np.isnan(float(Y[k])) or np.isnan(float(pix2sec[k])):
                                if verbosity: print('\033[91m ERROR! No WCS found. In this case you should specify world (Ra[deg],Dec[deg]), image (X[pix],Y[pix]) and Scale[arcsec/pix]  columns! Exiting... \033[0m ')
                                exit()
                            
                            # Add wcs using the coordinates of the target and input pixel scale:
                            add_wcs_to_header('imp-image.fits', 'imp-image.fits', float(X[k]), float(Y[k]), ra, dec, -float(pix2sec[k])/3600., float(pix2sec[k])/3600., 0., verbosity)
                        
                        # SWarping the image to reduce all image to a similar wcs:
                        if verbosity: print('SWarping to fix the image...')
                        swarpName = swarp_name(verbosity)
                        callSt = "%s -BACK_TYPE MANUAL -COMBINE_TYPE AVERAGE -VERBOSE_TYPE QUIET -BACK_DEFAULT 0.0 -INTERPOLATE N -RESAMPLE N " % (swarpName)
                        callSt += " ".join(["%s[0]" % (s) for s in ['imp-image.fits','imp-image.fits']])
                        subprocess.call(callSt, shell="True")
                        if verbosity: print('Done!')
                        
                        # Rename and remove SWarp tmp files 
                        shutil.move('coadd.fits', 'imp-image.fits')
                        os.remove('coadd.weight.fits')
                        os.remove('swarp.xml')
                        
                        imageHDU = pyfits.open('imp-image.fits', mode='update')
                        header = imageHDU[0].header
                        if 'GAIN' in header:
                            if header['GAIN'] == 0.:
                                header['GAIN'] = 10000.
                                
                        if 'EXPTIME' in header:
                            header['EXPTIME'] = float(header['EXPTIME'])/2.
                        
                        imageHDU.flush()
                        
                        # Go back to the initial directory
                        os.chdir(current_dir)
                    else:
                        zz=1
                        # This is a directory

            

            if ',' in galaxy_name:
                # FOR CREATING SUBFRAMES:
                try:
                    [ra,dec] = galaxy_name.split(',')
                    galaxy_name = None
                    ra = float(ra)
                    dec = float(dec)
                except:
                    if verbosity: print('\033[91m ERROR! In the galaxy name you can use comma only in this case: RA,DEC. Exiting... \033[0m ')


            
            # Repeat the script each time to show the menu:
            repeat_menu = True

            if rebin_to[k]!='nan':
                rebin_to_galaxy_directory = directory_to+'/%s_%s' % (galaxy_names[int(float(rebin_to[k]))-1], galaxy_frames[int(float(rebin_to[k]))-1].split('.fi')[0].split('/')[-1])
            else:
                rebin_to_galaxy_directory = 'nan'
                
            if resample_to[k]!='nan':
                resample_to_galaxy_directory = directory_to+'/%s_%s' % (galaxy_names[int(float(resample_to[k]))-1], galaxy_frames[int(float(resample_to[k]))-1].split('.fi')[0].split('/')[-1])
            else:
                resample_to_galaxy_directory = 'nan'                
                

            while repeat_menu==True:
                main_do(galaxy_directory, galaxy_name, input_image='imp-image.fits', sky_degree=sky_degrees[k], steps=steps[k], SB_lim=Sblim[k], wavelength=wavelengths[k], m0=m0[k], scale=pix2sec[k], user_interact=user_interact[k], RA=ra, DEC=dec, calibr_band=calibr_band[k], exptime=exptime[k], Ncombine=Ncombine[k], gain=gain[k], ron=ron[k], rebin_to_galaxy_directory=rebin_to_galaxy_directory, resample_to_galaxy_directory=resample_to_galaxy_directory,verbosity=verbosity)  
                
                
                if user_interact[k]:
                    repeat_menu = True
                else:
                    repeat_menu = False




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Image preparation: pipeline")
    parser.add_argument("input_sample", help="tsv file with a sample")
    parser.add_argument("sample_number", help="Item number to be done, starting from 1", type=int)
    parser.add_argument("-v", "--v", action="store_true", help="Optional: Verbosity", default=False) 
    args = parser.parse_args()

    input_sample = args.input_sample
    sample_number = args.sample_number
    verbosity = args.v
    
    
    main(sample_number, input_sample, verbosity=verbosity)



#sample_file = '/home/amosenko/MyWork/Test_image_preparation/sample_pavel.tsv'  # Specify the input file with a list of galaxies (galaxy_name\tgalaxy_frame) 
#number = 1
#main(number, sample_file, verbosity=True)






