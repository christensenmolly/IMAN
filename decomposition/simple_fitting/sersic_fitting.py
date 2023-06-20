#!/usr/bin/python
# DESCRIPTION:
# GALFIT fitting with a simple sersic function.
# MINIMAL USAGE: python sersic_fitting.py [input_image]
# EXAMPLE: python sersic_fitting.py galaxy.fits

import sys
import math
import numpy as np
from scipy import stats
import scipy as sp
from numpy import *
import os
import shutil
import subprocess
from astropy.io import fits as pyfits
import argparse

galfit_path = ''

LOCAL_DIR = "/decomposition/simple_fitting"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))
sys.path.append(os.path.join(IMAN_DIR, 'decomposition/make_model'))
sys.path.append(os.path.join(IMAN_DIR, 'plotting/2dprofile'))
sys.path.append(os.path.join(IMAN_DIR, 'plotting/1dprofile'))

import make_model_ima_galfit
import make_model_ima_imfit
import plot_2d_profile
import plot_profile
import get_galaxy_polygon


tmp_out = sys.stdout
FNULL = open(os.devnull, 'w')


def read_galfit_file(galfit_file):
    # Function to parse galfit file
    fff = open(galfit_file, 'r')
    lines = fff.readlines()
    C0 = float('nan')
    for line in lines:
        #print line
        if '1)' in line and '#  Position x, y' in line:
            xc = float(line.split()[1])
            yc = float(line.split()[2])
        if '3)' in line and '#  Integrated magnitude' in line:
            mtot = float(line.split()[1])
        if '4)' in line and '#  R_e (effective radius)   [pix]' in line:
            re = float(line.split()[1])
        if '5)' in line and '#  Sersic index n (de Vaucouleurs n=4) ' in line:
            n = float(line.split()[1])
        if '9)' in line and '#  Axis ratio (b/a)' in line:
            q = float(line.split()[1])
        if '10)' in line and '#  Position angle (PA) [deg: Up=0, Left=90]' in line:
            PA = float(line.split()[1])            
        if 'C0)' in line and '#  Diskyness(-)/Boxyness(+)' in line:
            C0 = float(line.split()[1])
    
    return xc,yc,mtot,re,n,q,PA,C0


def Bulge(component,xc,yc,mb,reb,n,q,C0=None,PA=90.,fix_center=False):
        # Function to define the simple sersic component
        print("\n# Component %i: General" % (component))
        print(" 0) sersic             # Object type")
        if fix_center:
            print(" 1) %.1f  %.1f  0 0    #  position x, y        [pixel]" % (xc,yc))
        else:
            print(" 1) %.1f  %.1f  1 1    #  position x, y        [pixel]" % (xc,yc))
        print(" 3) %.2f       1       #  total magnitude"  % (mb))   
        print(" 4) %.1f       1       #    R_e              [Pixels]" % (reb))
        print(" 5) %.1f       1       #  Sersic exponent (deVauc=4, expdisk=1)" % (n))
        print(" 9) %.1f        1       # axis ratio (b/a)" % (q))
        print("10) %.1f       1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA))
        print(" Z) 0                  #  Skip this model in output image?  (yes=1, no=0)")
        if C0 is not None:
            print(" C0) %f       1          #  Diskyness(-)/Boxyness(+)" % (C0))

 
def header(file_Image, file_Sigma, file_Psf, file_Mask, magZP, generalScaleValue=1., sampling=1, xmin=None, xmax=None, ymin=None, ymax=None, nx_psf=None, ny_psf=None):
        if file_Sigma is None:
            file_Sigma = 'none'

        #if file_Psf is None:
        #    file_Psf = 'none'

        if file_Mask is None:
            file_Mask = 'none'
            
        # Function to create header for galfit file
        hdulist = pyfits.open(file_Image)
        image = hdulist[0].data    
        ySize, xSize = image.shape
        if file_Psf is not None:
            hdulist = pyfits.open(file_Psf)
            image = hdulist[0].data    
            ySize_psf, xSize_psf = image.shape
        
            #psf_box = min([min([ySize_psf, xSize_psf]),ySize])
            psf_box = min([ySize_psf, xSize_psf])
        else:
            psf_box = min([ySize, xSize])
            #file_Psf
        
        if xmin is None:
            xmin = 1
        if xmax is None:
            xmax = xSize

        if ymin is None:
            ymin = 1
            
        if ymax is None:
            ymax = ySize
            
        if nx_psf is None:
            nx_psf = xSize
            
        if ny_psf is None:
            ny_psf = ySize            

        if file_Psf is None:
            nx_psf = 1
            ny_psf = 1
        
        print("\n===============================================================================")
        print("# IMAGE and GALFIT CONTROL PARAMETERS")
        print("A) %s                  # Input data image (FITS file)" % (file_Image))
        print("B) %s                  # Output data image block" % ('model.fits'))
        print("C) %s                # Sigma image name (made from data if blank or none)" % (str(file_Sigma)))
        print("D) %s                # Input PSF image and (optional) diffusion kernel" % (str(file_Psf)))
        print("E) %i                   # PSF fine sampling factor relative to data" % (sampling))
        print("F) %s                # Bad pixel mask (FITS image or ASCII coord list)" % (str(file_Mask)))
        print("G) %s                # File with parameter constraints (ASCII file)" % ('none'))
        print("H) %i    %i   %i    %i   # Image region to fit (xmin xmax ymin ymax)" % (xmin, xmax, ymin, ymax))
        print("I) %i    %i              # Size of the convolution box (x y)" % (psf_box,psf_box))
        print("J) %.3f              # Magnitude photometric zeropoint" % (magZP))
        print("K) %.3f  %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (generalScaleValue,generalScaleValue))
        print("O) regular             # Display type (regular, curses, both)")
        print("P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n")

        print("# INITIAL FITTING PARAMETERS")
        print("#")
        print("#   For object type, the allowed functions are:")
        print("#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat,")
        print("#       ferrer, powsersic, sky, and isophote.") 
        print("#")  
        print("#   Hidden parameters will only appear when they're specified:")
        print("#       C0 (diskyness/boxyness),") 
        print("#       Fn (n=integer, Azimuthal Fourier Modes),")
        print("#       R0-R10 (PA rotation, for creating spiral structures).")
        print("#") 
        print("# -----------------------------------------------------------------------------")
        print("#   par)    par value(s)    fit toggle(s)    # parameter description") 
        print("# -----------------------------------------------------------------------------\n")

# -----------------------------------------------------------------
# FUNCTION FOR SKY
def Sky(component,sky_level,fit_sky=1):
        # Function to define sky
        if fit_sky==1:
          fit_sky = 0
        else:
          fit_sky = 1
        print("\n# Component %i: SKY" % (component))
        print(" 0) sky")
        print(" 1) %.5f       %i       # sky background       [ADU counts]" % (sky_level, fit_sky))
        print(" 2) 0.000      %i       # dsky/dx (sky gradient in x)" % (fit_sky))
        print(" 3) 0.000      %i       # dsky/dy (sky gradient in y)" % (fit_sky)) 
        print(" Z) 0                  #  Skip this model in output image?  (yes=1, no=0)")
# -----------------------------------------------------------------  


def convert_segm_to_boolean(mask):
    return np.ma.make_mask(mask, shrink=False)


def main(input_image, sigma_image, psf_image, mask_image, ZP, pixel_scale, xc=None, yc=None, sampling=1, output_dir = './fitting', C0=None, initial_pars=None, xmin=None, xmax=None, ymin=None, ymax=None, verbosity=True, plot_graphs=True, nx_psf=None, ny_psf=None, fix_center=False, make_model=True, I_min=None):
   dirpath = os.getcwd()

   if output_dir!='.':
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.mkdir(output_dir)
   
        # Copy files:
        for file in [input_image, sigma_image, psf_image, mask_image]:
            if file is not None:
                shutil.copy(file, output_dir + '/' + file)
   
        # Change directory
        os.chdir(output_dir)
  
   hdulist = pyfits.open(input_image)
   prihdr = hdulist[0].header
   data = hdulist[0].data
   nx = int(prihdr['NAXIS1'])        
   ny = int(prihdr['NAXIS2'])

   hdulist_mask = pyfits.open(mask_image)
   mask = hdulist_mask[0].data
   
   mask_astropy = convert_segm_to_boolean(mask)
   
   
   try:
       EXPTIME = float(prihdr['EXPTIME'])
   except:
       EXPTIME = 1.
       
   if EXPTIME==0.:
       EXPTIME = 1.

   if yc is None:
        yc=ny/2.
   if xc is None:
        xc=nx/2.



   radius = min([xc,yc,nx-xc,ny-yc])



   if initial_pars is None:
        try:
            x0,y0,sma,smb,PA,Flux,mask_astropy = get_galaxy_polygon.main(input_image, mask_image, output_region='galaxy_polygon.reg', xc=xc, yc=yc, I_min=I_min, min_radius=10., verbosity=verbosity)
            PA = PA-90.
            if not fix_center: 
                xc = x0
                yc = y0
            q = smb/sma
            re = sma/2.
            n=2.
            mtot = ZP - 2.5*log10(Flux)
        except:
            re = radius/4.
            Flux = np.nansum(data * ~mask_astropy)
            n = 2.
            q = float(np.min([nx,ny])) / float(np.max([nx,ny]))
            PA = 0.
            mtot = ZP - 2.5*log10(Flux)
   else:
        [xc,yc,mtot,re,n,q,PA] = initial_pars

   f = open('galfit.inp', "w") 
   sys.stdout = f
   header(input_image, sigma_image, psf_image, mask_image, ZP-2.5*log10(EXPTIME), pixel_scale, sampling, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, nx_psf=nx_psf, ny_psf=ny_psf)    
   
   Bulge(1,xc,yc,mtot,re,n,q,C0=C0,PA=PA, fix_center=fix_center)
   Sky(2,0.,fit_sky=1)

 
   sys.stdout = tmp_out
   f.close()
   
   if verbosity: print('Simple Sersic Galfit fitting...')
   if os.path.exists('galfit.01'):
        os.remove('galfit.01')
   if os.path.exists('fit.log'):
        os.remove('fit.log')
   
   if verbosity:
         subprocess.call("%sgalfit galfit.inp" % (galfit_path), shell=True) 
   else:
        subprocess.call("%sgalfit galfit.inp" % (galfit_path), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
   
   for file in ['fit.log','model.fits','galaxy_polygon.reg','galaxy_polygon_ell.reg']:
        if os.path.exists(file):
                os.remove(file)
   
   if make_model:
        make_model_ima_galfit.main(input_image, 'galfit.01', composed_model_file='composed_model.fits',
                                                    subtract_sky=True, galfitPath=galfit_path)
   if plot_graphs:
        plot_2d_profile.main('composed_model.fits', None, pixel_scale, ZP, mask_file='mask_cropped.fits',
                                    borders='0,0,0,0', output_file='plot_2d.png', view='line', color='nipy_spectral_r',
                                    scale_bar=None, grid=None, show_bar_label='yes', region=None, text=None, show_negative=True,
                                    sigma=5.)
        
        plot_profile.main('composed_model.fits', ZP, pixel_scale, mask_image='mask_cropped.fits', profile='azim', xc=xc, yc=yc, PA=0., Rmin=0., Rmax=0., step=1., zmin=0., zmax=0., output_file='azim.png', AX=None, geom_units='arcsec', SB_units='mag/arcsec2', Scale=0.1, legend_size=6, interp=False, FWHM=3., max_SB=None, min_SB=None, do_not_show_full_model=False, plot_symbs='o', text=None)
   
   
   if os.path.exists('galfit.01'):
        xc,yc,mtot,re,n,q,PA,C0 = read_galfit_file('galfit.01')
        if verbosity: print('Done!')
        os.chdir(dirpath)
        return xc,yc,mtot,re,n,q,PA,C0
   else:
        if verbosity: print('Failed!')
        os.chdir(dirpath)
        return float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')
     
    
  
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simple Sersic fitting")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("--sigma_image", nargs='?', help="Optional: Input sigma image",type=str, default=None)     
    parser.add_argument("--psf_image", nargs='?', help="Optional: Input psf image",type=str, default=None)   
    parser.add_argument("--mask_image", nargs='?', help="Optional: Input mask image",type=str, default=None)       
    parser.add_argument("--output_dir", nargs='?', help="Optional: Output directory",type=str, default='./fitting')       
    
    parser.add_argument("--xcen", nargs='?', help="Optional: Galaxy center x",type=float,default=None)  
    parser.add_argument("--ycen", nargs='?', help="Optional: Galaxy center y",type=float,default=None)     

    parser.add_argument("--m0", nargs='?', help="Optional: Zero-point",type=float,default=20.0)  
    parser.add_argument("--scale", nargs='?', help="Optional: Pixel scale",type=float,default=1.0) 
    parser.add_argument("--sampling", nargs='?',  help="Optional: Sampling",type=int,default=1) 
    
    args = parser.parse_args()

    input_image = args.input_image
    sigma_image = args.sigma_image
    psf_image = args.psf_image
    mask_image = args.mask_image
    output_dir = args.output_dir
    xcen = args.xcen
    ycen = args.ycen
    m0 = args.m0
    scale = args.scale
    sampling = args.sampling

    main(input_image, sigma_image, psf_image, mask_image, m0, scale, xc=xcen, yc=ycen, sampling=sampling, output_dir = output_dir)  
  
  
  
