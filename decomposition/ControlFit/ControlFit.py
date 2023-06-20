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

LOCAL_DIR = "/decomposition/ControlFit"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'decomposition/make_model'))
sys.path.append(os.path.join(IMAN_DIR, 'plotting/2dprofile'))
sys.path.append(os.path.join(IMAN_DIR, 'plotting/1dprofile'))
sys.path.append(os.path.join(IMAN_DIR, 'deca/deca_tk_lib'))
import make_model_ima_galfit
import make_model_ima_imfit
import plot_2d_profile
import plot_profile
import radial_profile
import tex_creator
import plot_oned_profile

tmp_out = sys.stdout
FNULL = open(os.devnull, 'w')

galfit_path = ''


def get_input_images(input_file, code, m0, pix2sec):
    f = open(input_file, 'r')

    lines = f.readlines()
    if code == 'galfit':
        for line in lines:
            if 'A)' in line and '# Input data image (FITS file)' in line:
                input_image = line.split()[1]
            if 'C)' in line and '# Sigma image name' in line:
                sigma_image = line.split()[1]
            if 'D)' in line and '# Input PSF image and (optional) diffusion kernel' in line:
                psf_image = line.split()[1]
            if 'F)' in line and '# Bad pixel mask (FITS image or ASCII coord list)' in line:
                mask_image = line.split()[1]
            if 'J)' in line and '# Magnitude photometric zeropoint' in line:
                MagZP = float(line.split()[1])
                hdulist = pyfits.open(input_image)
                header = hdulist[0].header
                if 'EXPTIME' in header:
                    exptime = float(header['EXPTIME'])
                    m0 = MagZP + 2.5 * math.log10(exptime)
                else:
                    m0 = MagZP
            if 'K)' in line and '# Plate scale (dx dy)   [arcsec per pixel]' in line:
                pix2sec = float(line.split()[1])
    else:
        for line in lines:
            if '#' in line and 'imfit' in line:
                input_image = line.split()[2].split('[')[0]
                if '--psf=' in line:
                    psf_image = line.split('--psf=')[1].split()[0].split('[')[0]
                else:
                    psf_image = 'none'
                if '--mask=' in line:
                    mask_image = line.split('--mask=')[1].split()[0].split('[')[0]
                else:
                    mask_image = 'none'

                if '--noise=' in line:
                    sigma_image = line.split('--noise=')[1].split()[0].split('[')[0]
                else:
                    sigma_image = 'none'

                # TODO: Add m0 and pix2sec from the imfit input file as I did for HEROES

    input_image = os.path.abspath(input_image)

    if psf_image == 'none':
        psf_image = None
    else:
        psf_image = os.path.abspath(psf_image)

    if sigma_image == 'none':
        sigma_image = None
    else:
        sigma_image = os.path.abspath(sigma_image)

    if mask_image == 'none':
        mask_image = None
    else:
        mask_image = os.path.abspath(mask_image)

    f.close()

    return input_image, sigma_image, psf_image, mask_image, pix2sec, m0


def get_pars(input_file, code, full_image=False):
    f = open(input_file, 'r')

    lines = f.readlines()
    xc = [];
    yc = []
    if code == 'galfit':
        for line in lines:
            if '# Image region to fit (xmin xmax ymin ymax)' in line:
                xmin = int(line.split()[1])
                xmax = int(line.split()[2])
                ymin = int(line.split()[3])
                ymax = int(line.split()[4])

            if '1)' in line and '#  Position x, y' in line:
                if not full_image:
                    xc.append(float(line.split()[1]) - xmin)
                    yc.append(float(line.split()[2]) - ymin)
                else:
                    xc.append(float(line.split()[1]))
                    yc.append(float(line.split()[2]))

    else:
        for line in lines:
            if 'imfit' in line and '-c' in line and '[' in line and ']' in line:
                coords = line.split('[')[1].split(']')[0]
                xmin = int(coords.split(':')[0])
                xmax = int(coords.split(':')[1].split(',')[0])
                ymin = int(coords.split(',')[1].split(':')[0])
                ymax = int(coords.split(':')[2])
            else:
                xmin = 0.
                ymin = 0.

            if 'X0' in line:
                xc.append(float(line.split()[1]) - xmin)
            if 'Y0' in line:
                yc.append(float(line.split()[1]) - ymin)

    f.close()

    XC = np.mean(xc)
    YC = np.mean(yc)

    return XC, YC


def recognize_code(input_file):
    f = open(input_file, 'r')

    code = None

    lines = f.readlines()
    for line in lines:
        if '#' in line and 'imfit' in line:
            code_run = line.split('#')[1]
            code = 'imfit'
            output_file = 'bestfit_parameters_imfit.dat'
        if 'GALFIT' in line:
            code_run = '%sgalfit %s' % (galfit_path, input_file.replace(" ", "\ "))
            code = 'galfit'
            output_file = 'galfit.01'

    if code is None:
        print('Error! Input file is not recognized as galfit or imfit. Exiting!')
        exit()

    f.close()
    return code, code_run, output_file


def determine_new_dir(code, output_directory):
    dirs = glob.glob('%s/%s_*' % (output_directory, code))

    dir_numbers = []
    for di in dirs:
        dir_numbers.append(int(di.split('_')[-1]))
    if dir_numbers == []:
        new_dir = '%s/%s_1' % (output_directory, code)
    else:
        new_dir = '%s/%s_%i' % (output_directory, code, max(dir_numbers) + 1)

    return os.path.abspath(new_dir)


def main(input_file, output_directory='./trial_fits', crea_model=True, geom_units='arcsec', lum_units='lum',
         SB_units='mag/arcsec2', verbose=False, min_level=27., m0=None, pix2sec=None, new_dir=None, comp_names=[], sigma=None, show_negative=True, oned=False):

    current_dir = os.getcwd()

    input_file = os.path.abspath(input_file)

    code, code_run, output_file = recognize_code(input_file)
    input_image, sigma_image, psf_image, mask_image, pix2sec, m0 = get_input_images(input_file, code, m0, pix2sec)
    
    if new_dir is None:
        new_dir = determine_new_dir(code, output_directory)

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    if not os.path.exists(new_dir):
        os.mkdir(new_dir)

    os.chdir(new_dir)
    
    if not os.path.exists(output_file):
        for file in [input_image, sigma_image, psf_image, mask_image, input_file]:
            if file != None and file != 'None' and file!= 'none':
                shutil.copy(file, file.split('/')[-1])

        # Start fitting
        if verbose:
            subprocess.call(code_run, shell=True, stdout=FNULL)
        else:
            subprocess.call(code_run, shell=True)

        if code == 'galfit':
            if os.path.exists('galfit.01'):
                if crea_model:
                    make_model_ima_galfit.main(input_image, 'galfit.01', composed_model_file='composed_model.fits',
                                            subtract_sky=True, galfitPath=galfit_path)
                os.remove('fit.log')
                #output_file = 'galfit.01'
            else:
                print('Galfit failed! Exiting.')
                exit()
        elif code == 'imfit':
            if os.path.exists('bestfit_parameters_imfit.dat'):
                if crea_model:
                    make_model_ima_imfit.main(input_image, 'bestfit_parameters_imfit.dat', psf_image=psf_image,
                                            composed_model_file='composed_model.fits',
                                            comp_names=comp_names, imfitPath='', oned=oned) #'disc', 'AGN', 'bulge', 'ring', 'ring', 'bar'
                #output_file = 'bestfit_parameters_imfit.dat'
            else:
                print('Imfit failed! Exiting.')
                exit()

    xc, yc = get_pars(output_file, code)

    if not oned:
        # Plot 2d 
        plot_2d_profile.main('composed_model.fits', min_level, pix2sec, m0, mask_file='mask_cropped.fits',
                            borders='0,0,0,0', output_file='plot_2d.png', view='line', color='nipy_spectral_r',
                            scale_bar=None, grid=None, show_bar_label='yes', region=None, text=None, show_negative=show_negative,
                            sigma=sigma)

        # Plot 1d (azimuthally averaged)
        # if mask_image is not None:
        #    hdu_mask = pyfits.open(mask_image)
        #    mask = hdu_mask[0].data        
        #    mask_astropy = np.ones_like(np.array(mask,dtype=float),dtype=bool)
        # else:
        #    mask_astropy = None

        plot_profile.main('composed_model.fits', m0, pix2sec, mask_image='mask_cropped.fits', profile='azim', xc=xc, yc=yc, PA=0., Rmin=0., Rmax=0., step=1., zmin=0., zmax=0., output_file='azim.png', AX=None, geom_units='arcsec', SB_units='mag/arcsec2', Scale=0.1, legend_size=10, interp=False, FWHM=3., max_SB=None, min_SB=None, do_not_show_full_model=False, plot_symbs='o', text=None)

        # bin_centers, radial_prof = radial_profile.azimuthalAverage(data, center=[xc,yc], returnradii=True, binsize=step, weights=None, interpnan=False, mask=mask_astropy )

        # Create table with the output data

        '''
        pictures = ['plot_2d.png', 'plot_prof_azim.png']
        # Create tex-file with all these pictures and the table
        tex_creator.main(object_name,output_file,code,luminosities,pictures,m0,pix2sec,Distance,Scale,Filter,Aext,Kcorr,geom_units,lum_units,SB_units, continue_file=False, last_object == False)
                
        # TODO:        
        object_name
        luminosities
        Distance
        Scale
        Filter
        Aext
        Kcorr
        '''
    else:
        plot_oned_profile.main('composed_model.fits', m0, pix2sec, output_file='azim_profile.png', geom_units='arcsec', SB_units='mag/arcsec2', x_scale='linear')
    os.chdir(current_dir)


# main('imfit.inp', m0=20.472, pix2sec=0.75)
# main('galfit.inp', m0=29.55288, pix2sec=0.83)
# main('galfit.inp', m0=29.45153, pix2sec=0.833)
# main('galfit.inp', m0=20.472, pix2sec=0.75)
# main('galfit.inp', m0=8.9, pix2sec=0.6)

# xc,yc = get_pars('galfit.01', 'galfit')
# plot_profile.main('composed_model.fits', 20.472, 0.75, mask_image='mask_cropped.fits', profile = 'azim',xc=xc,yc=yc,PA=0.,Rmin=0.,Rmax=0.,step=1.,zmin=0.,zmax=0.,output_file='azim.png',AX=None, geom_units='arcsec',SB_units='mag/arcsec2',Scale=0.1, legend_size=6, interp=False, FWHM=3.,max_SB=26.,min_SB=None, do_not_show_full_model=False, plot_symbs='o',text=None)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Script to control fitting using an existing galfit file")
    parser.add_argument("input_file", help="Input galfit file (all images should be in the same directory)")
    parser.add_argument("m0", help="Zero-point")
    parser.add_argument("pix2sec", help="Scale")
    parser.add_argument("--replot", help="Optional: Path to directory for replotting the output image", type=str, default=None)
    parser.add_argument("--comps", help="Optional: Names of components in the order", type=str, default=None)
    parser.add_argument("--sigma", help="Optional: Smoothing of the output image", type=float, default=None)
    parser.add_argument("--SBfaint", help="Optional: Faintest SB level to show", type=float, default=25.)
    parser.add_argument("--neg", action="store_true", default=False,
                        help="Show negative values in the output 2D image")      
    parser.add_argument("--oned", action="store_true", default=False,
                        help="Fitting of 1D profile")   
    
    args = parser.parse_args()

    input_file = args.input_file
    m0 = float(args.m0)
    pix2sec = float(args.pix2sec)
    new_dir = args.replot
    if new_dir is not None:
        new_dir = os.path.abspath(new_dir)
    
    comp_names = args.comps
    
    if comp_names is None:
        comp_names = []
    else:
        comp_names = comp_names.split(',')
    
    sigma = args.sigma
    SBfaint = args.SBfaint
    show_negative = args.neg
    oned = args.oned

    
    main(input_file, m0=m0, pix2sec=pix2sec, new_dir=new_dir, comp_names=comp_names, sigma=sigma, min_level=SBfaint, show_negative=show_negative, oned=oned)
