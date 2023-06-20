#! /usr/bin/env python
#  python3 /home/amosenko/MyGit/IMAN/deca_lite/deca_lite.py galaxy.fits --mask_image mask.fits --psf psf.fits --comps sersic,edgedisk --distance 30.0 --Scale 0.1
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

LOCAL_DIR = "/deca_lite"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'decomposition/make_model'))
sys.path.append(os.path.join(IMAN_DIR, 'decomposition/simple_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'plotting/2dprofile'))
sys.path.append(os.path.join(IMAN_DIR, 'plotting/1dprofile'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))
sys.path.append(os.path.join(IMAN_DIR, 'deca/deca_tk_lib'))
import make_model_ima_galfit
import make_model_ima_imfit
import plot_2d_profile
import plot_profile
import radial_profile
import plot_oned_profile
import sersic_fitting
import get_galaxy_ellipse
import merge_masks
import tex_creator
tmp_out = sys.stdout
FNULL = open(os.devnull, 'w')

galfit_path = ''



def Bulge(component,xc,yc,mb,reb,n,q,C0=None,PA=90., fix_PA=1):
        # Function to define the simple sersic component
        print("\n# Component %i: General" % (component))
        print(" 0) sersic             # Object type")
        print(" 1) %.1f  %.1f  1 1    #  position x, y        [pixel]" % (xc,yc))
        print(" 3) %.2f       1       #  total magnitude"  % (mb))   
        print(" 4) %.1f       1       #    R_e              [Pixels]" % (reb))
        print(" 5) %.1f       1       #  Sersic exponent (deVauc=4, expdisk=1)" % (n))
        print(" 9) %.1f        1       # axis ratio (b/a)" % (q))
        print("10) %.1f       %i       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA, fix_PA))
        print(" Z) 0                  #  Skip this model in output image?  (yes=1, no=0)")
        if C0 is not None:
            print(" C0) %f       1          #  Diskyness(-)/Boxyness(+)" % (C0))


def read_galfit_edge_on(galfit_file):
    # Function to parse galfit file
    fff = open(galfit_file, 'r')
    lines = fff.readlines()

    for k in range(len(lines)):
        if "0) edgedisk" in lines[k]:
            start_eon_ind = k+1

    xc = float(lines[start_eon_ind].split()[1])
    yc = float(lines[start_eon_ind].split()[2])
    mu0d = float(lines[start_eon_ind+1].split()[1])
    z0 = float(lines[start_eon_ind+2].split()[1])
    h = float(lines[start_eon_ind+3].split()[1])
    PA = float(lines[start_eon_ind+8].split()[1])            

    return xc,yc,mu0d,z0,h,PA


def read_galfit_sersic(galfit_file):
    # Function to parse galfit file
    fff = open(galfit_file, 'r')
    lines = fff.readlines()
    
    for k in range(len(lines)):
        if "0) sersic" in lines[k]:
            start_ser_ind = k+1
            end_ser_ind = k+12
    
    C0 = float('nan')
    for k in range(start_ser_ind,end_ser_ind):
        line = lines[k]
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


def Edge_on_disk(component, xc, yc, mu0d, z0, h, PA=90., fix_PA=1):
        # Function to define the edge-on disk
        print("\n# Component %i: DISK" % (component))
        print("  0) edgedisk               #  Component type")
        print("  1) %.1f  %.1f 1 1  #  Position x, y"  % (xc, yc))
        print("  3) %.2f     1          #     Mu(0)   [mag/arcsec^2]" % (mu0d))
        print("  4) %.1f      1          #  h_s (disk scale-height)   [pix]" % (z0))
        print("  5) %.1f     1          #  R_s (disk scale-length)   [pix]" % (h))
        print(" 10) %.1f     %i          #  Position angle (PA) [deg: Up=0, Left=90]" % (PA, fix_PA))
        print("  Z) 0                      #  Skip this model in output image?  (yes=1, no=0)")



def run_fitting(xc, yc, input_image, ZP, pixel_scale, dirpath, verbosity=True, plot_graphs=True):
   if verbosity: print('DECA fitting...')
   if os.path.exists('galfit.01'):
        os.remove('galfit.01')
   if os.path.exists('fit.log'):
        os.remove('fit.log')
   
   if verbosity:
        subprocess.call("%sgalfit galfit.inp" % (galfit_path), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)   
   else:
        subprocess.call("%sgalfit galfit.inp" % (galfit_path), shell=True)
   
   for file in ['fit.log','model.fits','galaxy_polygon.reg','galaxy_polygon_ell.reg']:
        if os.path.exists(file):
                os.remove(file)
   

   
   
   if os.path.exists('galfit.01'):
        xc_d,yc_d,mu0d,z0,h,PA_d = read_galfit_edge_on('galfit.01')
        try:
            xc_b,yc_b,mtot,re,n,q,PA_b,C0 = read_galfit_sersic('galfit.01')
        except:
            xc_b,yc_b,mtot,re,n,q,PA_b,C0 = float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')
            #if verbosity: print('Failed!')
            #os.chdir(dirpath)            
            #return [float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')],[float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')]
   else:
        if verbosity: print('Failed!')
        #os.chdir(dirpath)
        return [float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')],[float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')]

   if plot_graphs:
        make_model_ima_galfit.main(input_image, 'galfit.01', composed_model_file='composed_model.fits',
                                                    subtract_sky=True, galfitPath=galfit_path)

        plot_2d_profile.main('composed_model.fits', None, pixel_scale, ZP, mask_file='mask_cropped.fits',
                                    borders='0,0,0,0', output_file='plot_2d.png', view='line', color='nipy_spectral_r',
                                    scale_bar=None, grid=None, show_bar_label='yes', region=None, text=None, show_negative=True,
                                    sigma=5.)#, rotate_PA=90.+PA_d)
        
        plot_profile.main('composed_model.fits', ZP, pixel_scale, mask_image='mask_cropped.fits', profile='cut', xc=xc_d, yc=yc_d, PA=90.+PA_d, Rmin=0., Rmax=6.*h, step=1., zmin=0., zmax=0., output_file='major_axis_cut.png', AX=None, geom_units='arcsec', SB_units='mag/arcsec2', Scale=0.1, legend_size=6, interp=False, FWHM=3., max_SB=None, min_SB=None, do_not_show_full_model=False, plot_symbs='o', text=None)

   if verbosity: print('Done!')
   #os.chdir(dirpath)
   return [xc_d,yc_d,mu0d,z0,h,PA_d],[xc_b,yc_b,mtot,re,n,q,PA_b,C0]     


# -----------------------------------------------------------------
# EDGE-ON DISC SB
def m0_disc_edge_on_f(mag_disc,h,z0):
	return +2.5*math.log10(2.*math.pi) + mag_disc + 2.5*math.log10(z0*h)
# -----------------------------------------------------------------


def create_and_change_dir(output_dir, input_image, sigma_image, psf_image, mask_image):
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


def main(input_image, sigma_image, psf_image, mask_image,
         ZP, pixel_scale,
         xc=None, yc=None,
         comps = ['sersic'],
         comp_names = ['disk'],
         sampling=1,
         C0 = None,
         initial_pars=None,
         xmin=None, xmax=None, ymin=None, ymax=None,
         verbosity=True,
         plot_graphs=True,
         mask_cen_for_pure_disks=False,
         Distance=None,
         Scale=None,
         Aext=0.,
         Kcorr=0.,
         object_name=None,
         run_sersic=True):

   start = time.time()
   
   if run_sersic:
        # A necessary fit
        xc_ser,yc_ser,mtot_ser,re_ser,n_ser,q_ser,PA_ser,C0_ser = sersic_fitting.main(input_image, sigma_image, psf_image, mask_image, ZP, pixel_scale, xc=xc, yc=yc, sampling=sampling, output_dir = './deca_sersic', C0=C0, initial_pars=initial_pars, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, verbosity=verbosity, plot_graphs=plot_graphs)
   else:
       [xc_ser,yc_ser,mtot_ser,re_ser,n_ser,q_ser,PA_ser] = initial_pars
       C0_ser = 0.



   dirpath = os.getcwd()





   hdulist = pyfits.open(input_image)
   prihdr = hdulist[0].header
   data = hdulist[0].data
   nx = int(prihdr['NAXIS1'])        
   ny = int(prihdr['NAXIS2'])
   
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


   
   if 'edgedisk' in comps:
        if True:    
            create_and_change_dir('deca_edgedisk', input_image, sigma_image, psf_image, mask_image)

            if mask_cen_for_pure_disks:
                get_galaxy_ellipse.create_ellipse_region(input_image, xc_ser, yc_ser, re_ser, re_ser, 0., 'mask_center.reg')
                convert_reg_to_mask.mask(input_image, 'mask_center.reg', output_image=None, output_mask='mask_center.fits', mask_value=1, mask_DN=None, verbosity=verbosity)
                if mask_image is not None:
                        merge_masks.main(['mask_center.fits',mask_image], 'deca_mask.fits', verbosity=True)
                else:
                        shutil.copy('mask_center.fits','deca_mask.fits')
                
                mask_image = 'deca_mask.fits'           

       
        f = open('galfit.inp', "w") 
        sys.stdout = f
        sersic_fitting.header(input_image, sigma_image, psf_image, mask_image, ZP-2.5*math.log10(EXPTIME), pixel_scale, sampling, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)    
        
        #Bulge(1,xc,yc,mtot,re,n,q,C0=C0,PA=PA)
        
        mu0d = m0_disc_edge_on_f(mtot_ser, re_ser/1.68, q_ser*re_ser/(1.68))
        
        Edge_on_disk(1, xc_ser, yc_ser, mu0d, q_ser*re_ser/(1.68), re_ser/1.68, PA=PA_ser)
        sersic_fitting.Sky(2,0.,fit_sky=1)
  
        sys.stdout = tmp_out
        f.close()
        
        [xc_d,yc_d,mu0d,z0,h,PA_d],[xc_b,yc_b,mtot,re,n,q,PA_b,C0_b] = run_fitting(xc_ser, yc_ser,input_image, ZP, pixel_scale, dirpath, verbosity=verbosity)
        
        if Distance is not None and Scale is not None and os.path.exists('galfit.01'):
            tex_creator.main(object_name,'galfit.01','composed_model.fits','GALFIT',['plot_2d.png','major_axis_cut.png'],ZP,pixel_scale,Distance,Scale,float('nan'),Aext,Kcorr,'kpc','mag','mag/arcsec2')
        os.chdir(dirpath)

 
        if 'sersic' in comps:
            create_and_change_dir('deca_sersic_edgedisk', input_image, sigma_image, psf_image, mask_image)
       

       
        f = open('galfit.inp', "w") 
        sys.stdout = f
        sersic_fitting.header(input_image, sigma_image, psf_image, mask_image, ZP-2.5*math.log10(EXPTIME), pixel_scale, sampling, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)    
        
        Bulge(1, xc_ser, yc_ser, mtot_ser+0.75, re_ser/5., 2., 0.8, C0=C0, PA=PA_d, fix_PA=0)        
        Edge_on_disk(2, xc_d, yc_d, mu0d, z0, h, PA=PA_d, fix_PA=0)
        sersic_fitting.Sky(3, 0., fit_sky=1)
  
        sys.stdout = tmp_out
        f.close()

        [xc_d,yc_d,mu0d,z0,h,PA_d],[xc_b,yc_b,mtot,re,n,q,PA_b,C0_b] = run_fitting(xc_ser, yc_ser,input_image, ZP, pixel_scale, dirpath, verbosity=verbosity)

        if Distance is not None and Scale is not None and os.path.exists('galfit.01'):
            tex_creator.main(object_name,'galfit.01','composed_model.fits','GALFIT',['plot_2d.png','major_axis_cut.png'],ZP,pixel_scale,Distance,Scale,float('nan'),Aext,Kcorr,'kpc','mag','mag/arcsec2')
        
        os.chdir(dirpath)
   print('Time: %g s.' % (time.time() - start))
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Simple Sersic fitting")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("--sigma_image", nargs='?', help="Optional: Input sigma image",type=str, default=None)     
    parser.add_argument("--psf_image", nargs='?', help="Optional: Input psf image",type=str, default=None)   
    parser.add_argument("--mask_image", nargs='?', help="Optional: Input mask image",type=str, default=None)       
    parser.add_argument("--mask_cen", nargs='?', help="Optional: Mask center while fitting edge-on disk",type=bool, default=False)       
    
    parser.add_argument("--xcen", nargs='?', help="Optional: Galaxy center x",type=float,default=None)  
    parser.add_argument("--ycen", nargs='?', help="Optional: Galaxy center y",type=float,default=None)     

    parser.add_argument("--m0", nargs='?', help="Optional: Zero-point",type=float,default=20.0)  
    parser.add_argument("--scale", nargs='?', help="Optional: Pixel scale",type=float,default=1.0) 
    parser.add_argument("--sampling", nargs='?', help="Optional: Sampling",type=int,default=1) 
    parser.add_argument("--comps", help="Optional: GALFIT components in the order (Default: sersic). Other availables components: edgedisk", type=str, default='sersic')
    parser.add_argument("--comp_names", help="Optional: Names of components in the order (Default: model). Other availables components: disk, bulge", type=str, default='model')
    
    parser.add_argument("--distance", nargs='?', help="Optional: Distance (Mpc)", type=float, default=None) 
    parser.add_argument("--Scale", nargs='?', help="Optional: Scale (kpc/arcsec)", type=float, default=None)
    parser.add_argument("--Aext", nargs='?', help="Optional: Galactic dust attenuation", type=float, default=0.)
    parser.add_argument("--Kcorr", nargs='?', help="Optional: K-correction", type=float, default=0.)
    parser.add_argument("--name", nargs='?', help="Optional: Object name",type=str, default=None) 
    
    args = parser.parse_args()

    input_image = args.input_image
    sigma_image = args.sigma_image
    psf_image = args.psf_image
    mask_image = args.mask_image
    mask_cen_for_pure_disks = args.mask_cen
    xcen = args.xcen
    ycen = args.ycen
    m0 = args.m0
    scale = args.scale
    sampling = args.sampling
    comps = args.comps
    comp_names = args.comp_names
    distance = args.distance
    Scale = args.Scale
    Aext = args.Aext
    Kcorr = args.Kcorr
    object_name = args.name
    
    comps = comps.split(',')
    comp_names = comp_names.split(',')

  
    main(input_image, sigma_image, psf_image, mask_image,
         m0, scale,
         xc=xcen, yc=ycen,
         comps = comps,
         comp_names = comp_names,
         sampling=sampling,
         C0 = None,
         initial_pars=None,
         xmin=None, xmax=None, ymin=None, ymax=None,
         verbosity=True,
         plot_graphs=True,
         mask_cen_for_pure_disks=mask_cen_for_pure_disks,
         Distance=distance,
         Scale=Scale,
         Aext=Aext,
         Kcorr=Kcorr,
         object_name=object_name
         )  
          
    
    
    
    
