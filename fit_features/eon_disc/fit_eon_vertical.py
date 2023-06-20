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
import warnings
import random
warnings.filterwarnings("ignore")

FNULL = open(os.devnull, 'w')
sys.path.append( '/home/amosenko/MEGA/MyPrograms/IMAN/imp/masking')
sys.path.append( '/home/amosenko/MEGA/MyPrograms/IMAN/decomposition/make_model')
import make_model_ima_galfit
import merge_masks

path = ''#'/Users/mosenkov/MEGA/HERON/ESO240-G011/'
imfit_path = path#'/home/amosenko/MEGA/MyPrograms/imfit-1.6.1/'
galfit_path = path
fitscopy = path

####
n_times = 10 # Max number of randomly selected cuts from a given interval
code = 'galfit'
####



def crea_galfit_input(input_image, sigma_image, psf_image, mask_image, xc, yc, mu0d1, z01=None, h1=None,  mu0d2=None, z02=None, h2=None):
    if sigma_image is None:
        sigma_image = 'none'

    if mask_image is None:
        mask_image = 'none'

    hdu = pyfits.open(input_image)
    image = hdu[0].data
    ny,nx =  np.shape(image)



    if psf_image is None:
        psf_image = 'none'
        n_psf = 1
    else:
        n_psf = ny

    
    if mu0d2 is None:
            xc_ch = 0
            yc_ch = 1 #### TODO 1 
    else:
            xc_ch = 0
            yc_ch = 0        
    
    h1_ch = 0
    h2_ch = 0
    z01_ch = 1 
    z02_ch = 1    
    
    
    #############
    #h1 = 40.
    #############
        
    s1="""

================================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) %s         # Input data image (FITS file)
B) model.fits          # Output data image block
C) %s          # Sigma image name (made from data if blank or "none") 
D) %s            # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data 
F) %s           # Bad pixel mask (FITS image or ASCII coord list)
G) none                # File with parameter constraints (ASCII file) 
H) 1  %i  1  %i  # Image region to fit (xmin xmax ymin ymax)
I) %i    %i          # Size of the convolution box (x y)
J) 0              # Magnitude photometric zeropoint 
K) 1.0  1.0        # Plate scale (dx dy)   [arcsec per pixel]
O) regular             # Display type (regular, curses, both)
P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

# INITIAL FITTING PARAMETERS
#
#   For component type, the allowed functions are: 
#       sersic, expdisk, edgedisk, devauc, king, nuker, psf, 
#       gaussian, moffat, ferrer, and sky. 
#  
#   Hidden parameters will only appear when they're specified:
#       Bn (n=integer, Bending Modes).
#       C0 (diskyness/boxyness), 
#       Fn (n=integer, Azimuthal Fourier Modes).
#       R0-R10 (coordinate rotation, for creating spiral structures).
#       To, Ti, T0-T10 (truncation function).
# 
# ------------------------------------------------------------------------------
#   par)    par value(s)    fit toggle(s)    # parameter description 
# ------------------------------------------------------------------------------
""" % (input_image, sigma_image, psf_image, mask_image, nx, ny, 3*n_psf, 3*n_psf)

    s2="""
# Component number: 1
 0) edgedisk               #  Component type
 1) %f %f %i %i  #  Position x, y
 3) %f     1          #     Mu(0)   [mag/arcsec^2]
 4) %f      %i          #  h_s (disk scale-height)   [pix]
 5) %f      %i          #  R_s (disk scale-length)   [pix]
10) 90.0     0          #  Position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)
""" % (xc, yc, xc_ch, yc_ch, mu0d1, z01, z01_ch, h1, h1_ch)

    if mu0d2 is not None:
        s3="""
# Component number: 2
 0) edgedisk               #  Component type
 1) %f %f %i %i  #  Position x, y
 3) %f     1          #     Mu(0)   [mag/arcsec^2]
 4) %f      %i          #  h_s (disk scale-height)   [pix]
 5) %f      %i          #  R_s (disk scale-length)   [pix]
10) 90.0     0          #  Position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)
================================================================================
""" % (xc, yc, xc_ch, yc_ch, mu0d2, z02, z02_ch, h2, h2_ch)
    else:
        s3 = "\n================================================================================"
        
    s = s1 + s2 + s3
    
    f = open("input_galf.txt", 'w')
    f.write(s)
    f.close()



def crea_galfit_input_tripple(input_image, sigma_image, psf_image, mask_image, xc, yc, mu0d1, z01=None, h1=None,  mu0d2=None, z02=None, h2=None, mu0d3=None, z03=None, h3=None):
    if sigma_image is None:
        sigma_image = 'none'

    if psf_image is None:
        psf_image = 'none'
        n_psf=1
    else:
        hdu = pyfits.open(psf_image)
        image_psf = hdu[0].data 
        ny_psf,nx_psf =  np.shape(image_psf)
        n_psf = max([ny_psf,nx_psf])

    if mask_image is None:
        mask_image = 'none'

    hdu = pyfits.open(input_image)
    image = hdu[0].data
    ny,nx =  np.shape(image)
    

    xc_ch = 0
    yc_ch = 0        
    
    
    h1_ch = 0
    h2_ch = 0
    h3_ch = 0
    z01_ch = 1 
    z02_ch = 1
    z03_ch = 1   
        
        
        
    s1="""

================================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) %s         # Input data image (FITS file)
B) model.fits          # Output data image block
C) %s          # Sigma image name (made from data if blank or "none") 
D) %s            # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data 
F) %s           # Bad pixel mask (FITS image or ASCII coord list)
G) none                # File with parameter constraints (ASCII file) 
H) 1  %i  1  %i  # Image region to fit (xmin xmax ymin ymax)
I) %i    %i          # Size of the convolution box (x y)
J) 0              # Magnitude photometric zeropoint 
K) 1.0  1.0        # Plate scale (dx dy)   [arcsec per pixel]
O) regular             # Display type (regular, curses, both)
P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

# INITIAL FITTING PARAMETERS
#
#   For component type, the allowed functions are: 
#       sersic, expdisk, edgedisk, devauc, king, nuker, psf, 
#       gaussian, moffat, ferrer, and sky. 
#  
#   Hidden parameters will only appear when they're specified:
#       Bn (n=integer, Bending Modes).
#       C0 (diskyness/boxyness), 
#       Fn (n=integer, Azimuthal Fourier Modes).
#       R0-R10 (coordinate rotation, for creating spiral structures).
#       To, Ti, T0-T10 (truncation function).
# 
# ------------------------------------------------------------------------------
#   par)    par value(s)    fit toggle(s)    # parameter description 
# ------------------------------------------------------------------------------
""" % (input_image, sigma_image, psf_image, mask_image, nx, ny, n_psf, n_psf)

    s2="""
# Component number: 1
 0) edgedisk               #  Component type
 1) %f %f %i %i  #  Position x, y
 3) %f     1          #     Mu(0)   [mag/arcsec^2]
 4) %f      %i          #  h_s (disk scale-height)   [pix]
 5) %f      %i          #  R_s (disk scale-length)   [pix]
10) 90.0     0          #  Position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)
""" % (xc, yc, xc_ch, yc_ch, mu0d1, z01, z01_ch, h1, h1_ch)


    s3="""
# Component number: 2
 0) edgedisk               #  Component type
 1) %f %f %i %i  #  Position x, y
 3) %f     1          #     Mu(0)   [mag/arcsec^2]
 4) %f      %i          #  h_s (disk scale-height)   [pix]
 5) %f      %i          #  R_s (disk scale-length)   [pix]
10) 90.0     0          #  Position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)
================================================================================
""" % (xc, yc, xc_ch, yc_ch, mu0d2, z02, z02_ch, h2, h2_ch)

    s4="""
# Component number: 3
 0) edgedisk               #  Component type
 1) %f %f %i %i  #  Position x, y
 3) %f     1          #     Mu(0)   [mag/arcsec^2]
 4) %f      %i          #  h_s (disk scale-height)   [pix]
 5) %f      %i          #  R_s (disk scale-length)   [pix]
10) 90.0     0          #  Position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)
================================================================================
""" % (xc, yc, xc_ch, yc_ch, mu0d3, z03, z03_ch, h3, h3_ch)
        
    s = s1 + s2 + s3 + s4
    
    f = open("input_galf.txt", 'w')
    f.write(s)
    f.close()



def read_galf(galf_file):
  f = open(galf_file,'r')
  L0 = 99999.
  chi2 = 1.
  for line in f:
    if "#  Chi^2/nu" in line:
      chi2 = float(line.split('#  Chi^2/nu =')[-1].split(',')[0])
    if "#  Position x, y" in line:
        x0 = float(line.split()[1])
        y0 = float(line.split()[2])
    if '#     Mu(0)   [mag/arcsec^2]' in line:
      L0 = float(line.split()[1])
    if '#  h_s (disk scale-height)   [pix]' in line and L0!=99999.:
      z0 = float(line.split()[1])
    if '#  R_s (disk scale-length)   [pix]' in line:
      h = float(line.split()[1])
  f.close()
  
  L0 = 10**(0.4*(-L0))/(2.*h)
  return x0,y0,L0,z0/2.,h,chi2


def read_galf_double(galf_file):
  f = open(galf_file,'r')
  func = 0
  z01 = 99999.
  for line in f:
    if "#  Chi^2/nu" in line:
      chi2 = float(line.split('#  Chi^2/nu =')[-1].split(',')[0])
    if "0) edgedisk               #  Component type" in line and func==0:
      func=1
    if '#     Mu(0)   [mag/arcsec^2]' in line and func==1:
      L01 = float(line.split()[1])      
    if '#  h_s (disk scale-height)   [pix]' in line and func==1:
      z01 = float(line.split()[1])
    if '#  R_s (disk scale-length)   [pix]' in line and func==1:
      h1 = float(line.split()[1])
    
    if "0) edgedisk               #  Component type" in line and func==1 and z01!=99999.:
      func=2
    if '#     Mu(0)   [mag/arcsec^2]' in line and func==2:
      L02 = float(line.split()[1])   
    if '#  h_s (disk scale-height)   [pix]' in line and func==2:
      z02 = float(line.split()[1])
    if '#  R_s (disk scale-length)   [pix]' in line and func==2:
      h2 = float(line.split()[1])    
  f.close()
  
  L01 = 10**(0.4*(-L01))/(2.*h1)
  L02 = 10**(0.4*(-L02))/(2.*h2)
  if z02>z01:
      return L01,h1,z01/2.,L02,h2,z02/2.,chi2
  else:
      return L02,h2,z02/2.,L01,h1,z01/2.,chi2      


def read_galf_tripple(galf_file):
  f = open(galf_file,'r')
  func = 0
  z01 = 99999.
  z02 = 99999.
  for line in f:
    if "#  Chi^2/nu" in line:
      chi2 = float(line.split('#  Chi^2/nu =')[-1].split(',')[0])
    if "0) edgedisk               #  Component type" in line and func==0:
      func=1
    if '#     Mu(0)   [mag/arcsec^2]' in line and func==1:
      L01 = float(line.split()[1])      
    if '#  h_s (disk scale-height)   [pix]' in line and func==1:
      z01 = float(line.split()[1])
    if '#  R_s (disk scale-length)   [pix]' in line and func==1:
      h1 = float(line.split()[1])
    
    if "0) edgedisk               #  Component type" in line and func==1 and z01!=99999.:
      func=2
    if '#     Mu(0)   [mag/arcsec^2]' in line and func==2:
      L02 = float(line.split()[1])   
    if '#  h_s (disk scale-height)   [pix]' in line and func==2:
      z02 = float(line.split()[1])
    if '#  R_s (disk scale-length)   [pix]' in line and func==2:
      h2 = float(line.split()[1])    

    if "0) edgedisk               #  Component type" in line and func==2 and z02!=99999.:
      func=3
    if '#     Mu(0)   [mag/arcsec^2]' in line and func==3:
      L03 = float(line.split()[1])   
    if '#  h_s (disk scale-height)   [pix]' in line and func==3:
      z03 = float(line.split()[1])
    if '#  R_s (disk scale-length)   [pix]' in line and func==3:
      h3 = float(line.split()[1])   
  f.close()
  
  L01 = 10**(0.4*(-L01))/(2.*h1)
  L02 = 10**(0.4*(-L02))/(2.*h2)
  L03 = 10**(0.4*(-L03))/(2.*h3)


  RES= [[L01,h1,z01/2.],[L02,h2,z02/2.],[L03,h3,z03/2.]]
  zzz = [z01,z02,z03]
  inds_sorted = sorted(range(len(zzz)), key=lambda k: zzz[k])
  
  return RES[inds_sorted[0]],RES[inds_sorted[1]],RES[inds_sorted[2]],chi2 



def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

def chi2_func(ref,mod):
    return sum( (ref-mod)**2)

def double_disc_func(z, I0d1, z01, I0d2, z02 ):
    #*** For edge-on disk SB in mag/arcsec^2 (along z-axis). n - index of this law. ***
    return I0d1 * (1.0 / np.cosh(1.*np.fabs(z)/(2.*z01)) ) **(2./1.) + I0d2 * (1.0 / np.cosh(1.*np.fabs(z)/(2.*z02)) ) **(2./1.)

def find_double_disc_guess(zz, II, y0, h1, h2):#, z0, I0):
      zz = zz - y0
      ZZ = zz
      z = []; I = []
      for k in range(len(zz)):
          if zz[k]>=0.:
            try:
              I.append(math.log10(II[k]))
              z.append(zz[k])
            except:
                zzz=1
      z = np.array(z); I = np.array(I)  
      #plt.plot(z, I)

      Res_chi2 = []
      Res_p = []
      Res_e = []
      for zz in z:
          try:
            ind = list(z).index(zz)
            p , e = curve_fit(piecewise_linear, z, I,[zz,I[ind],1.,1.])
            chi2 = chi2_func(I,piecewise_linear(z, *p))
            Res_chi2.append(chi2)
            Res_p.append(p)
            Res_e.append(np.sqrt(np.diag(e)))
          except:
              ppp =1
      try:
        results_pos = Res_p[Res_chi2.index(min(Res_chi2))]
        results_err = Res_e[Res_chi2.index(min(Res_chi2))]
        
        z01 = np.fabs(-1.086/results_pos[2])
        z02 = np.fabs(-1.086/results_pos[3])
        if z01>z02:
            z0_thick=z01
            z0_thin=z02    
            I0_thin = 10**(results_pos[1]-results_pos[3]*results_pos[0])
            I0_thick = 10**(results_pos[1]-results_pos[2]*results_pos[0])
            z0_thin_err = results_err[3]
            z0_thick_err = results_err[2]
        else:
            z0_thick=z02
            z0_thin=z01    
            I0_thin = 10**(results_pos[1]-results_pos[2]*results_pos[0])
            I0_thick = 10**(results_pos[1]-results_pos[3]*results_pos[0])
            z0_thin_err = results_err[2]
            z0_thick_err = results_err[3]
      except:
          I0_thin = np.nanmax(II)
          z0_thin = len(zz)/40.
          I0_thick = I0_thin/10.
          z0_thick = len(zz)/20.
      #print()
      #print(y0)
      
      #plt.plot(z, piecewise_linear(z, results_pos[0], results_pos[1], results_pos[2], results_pos[3]))
      

      p,e = curve_fit(double_disc_func, ZZ, II, [I0_thin, z0_thin, I0_thick, z0_thick] )

      I0_thin = p[0]
      z0_thin = np.fabs(p[1])

      I0_thick = p[2]
      z0_thick = np.fabs(p[3])
      #print(I0_thin,I0_thick,z0_thin,z0_thick)
      #plt.plot(ZZ, -2.5*np.log10(II))
      #plt.plot(ZZ, -2.5*np.log10(double_disc_func(ZZ, I0_thin, z0_thin, I0_thick, z0_thick)))
      #plt.show()
      #exit()
      if z0_thin<z0_thick:
            return I0_thin/(2.*h1),I0_thick/(2.*h1),z0_thin,z0_thick,
      else:
            return I0_thick/(2.*h2),I0_thin/(2.*h2),z0_thick,z0_thin         

def disk_exp(r, I0d, h, x_c):
    return I0d * np.exp(-np.fabs(r-x_c)/h)

def crea_hor_average_cut(image_data, mask_data, sigma_data, z_min_ind, z_max_ind, x_min_ind, x_max_ind):
    I_image = []; I_mask = []; I_sigma = []; rr = []; r_no_mask = []; I_no_mask = []
    ny,nx = np.shape(image_data)
    
    for i in range(x_min_ind, x_max_ind):  
        image_array = list(image_data[z_min_ind:z_max_ind,i])
        mask_array = list(mask_data[z_min_ind:z_max_ind,i])
        sigma_array = list(sigma_data[z_min_ind:z_max_ind,i])
        
        a_image = np.ma.array(image_array, mask = mask_array)
        a_sigma = np.ma.array(sigma_array, mask = mask_array)
        try:
            ii = float(np.ma.average(a_image, weights=1./(np.array(sigma_array)**2))) #WARNING
            if np.isnan(ii) or np.isinf(ii):
                I_image.append(0.)
                I_mask.append(True)
                I_sigma.append(0.)
            else:
                I_image.append(ii)
                I_mask.append(False)
                I_sigma.append(float(np.ma.average(a_sigma, weights=1./(np.array(sigma_array)**2)))) #WARNING
                r_no_mask.append(i+1)
                I_no_mask.append(ii)
        except:
            I_image.append(0.)
            I_mask.append(True)
            I_sigma.append(0.)
    return np.array(r_no_mask), np.array(I_no_mask)

def crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind, y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up):
    #print(y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up)

    I_image = []; I_mask = []; I_sigma = []; rr = []; r_no_mask = []; I_no_mask = []
    ny,nx = np.shape(image_data)
    xc = nx/2.
    yc = ny/2.
    
    '''
    for i in range(y_min_ind, y_max_ind):  
        image_array = list(image_data[i,x_min_ind:x_max_ind])
        mask_array = list(mask_data[i,x_min_ind:x_max_ind])
        sigma_array = list(sigma_data[i,x_min_ind:x_max_ind])
    '''

    for i in range(0, ny):
      rr.append(i+1)
      if (i >=y_min_ind_dn and i<y_max_ind_dn) or (i>=y_min_ind_up and i<y_max_ind_up):
        image_array = list(image_data[i,x_min_ind:x_max_ind])
        mask_array = list(mask_data[i,x_min_ind:x_max_ind])
        sigma_array = list(sigma_data[i,x_min_ind:x_max_ind])
        a_image = np.ma.array(image_array, mask = mask_array)
        a_sigma = np.ma.array(sigma_array, mask = mask_array)

        try:
            ii = float(np.ma.average(a_image, weights=1./(np.array(sigma_array)**2))) #WARNING
            if np.isnan(ii) or np.isinf(ii):
                I_image.append(0.)
                I_mask.append(True)
                I_sigma.append(0.)
                
            else:
                I_image.append(ii)
                I_mask.append(False)
                I_sigma.append(float(np.ma.average(a_sigma, weights=1./(np.array(sigma_array)**2)))) #WARNING
                r_no_mask.append(i+1)
                I_no_mask.append(ii)
        except:
            I_image.append(0.)
            I_mask.append(True)
            I_sigma.append(0.)
      else:
            I_image.append(0.)
            I_mask.append(True)
            I_sigma.append(0.) 
    
    x_max_ind = nx ##### WARNING!!!
    #x_min_ind = ds9_to_np(xc)
    
    r_plane,I_plane = crea_hor_average_cut(image_data, mask_data, sigma_data, ds9_to_np(ny/2.), ds9_to_np(ny/2.+ny/6.), x_min_ind, x_max_ind)
    popt_soph, pcov_soph = curve_fit(lambda r, I0, h: disk_exp(r, I0, h, xc), r_plane, I_plane, p0=(max(I_plane), nx/10.) )
    h1 = popt_soph[1] 
    #plt.plot(r_plane, -2.5*np.log10(I_plane))
    #plt.plot(r_plane, -2.5*np.log10(disk_exp(r_plane, popt_soph[0], 40., xc)))
    #plt.show()
    #print(hd1)
    #exit()
    
    r_plane,I_plane = crea_hor_average_cut(image_data, mask_data, sigma_data, ds9_to_np(ny/2.+ny/6.), ds9_to_np(ny/2.+ny/3.), x_min_ind, x_max_ind)
    popt_soph, pcov_soph = curve_fit(lambda r, I0, h: disk_exp(r, I0, h, xc), r_plane, I_plane, p0=(max(I_plane), nx/10.) )
    h2 = popt_soph[1]     

    r_plane,I_plane = crea_hor_average_cut(image_data, mask_data, sigma_data, ds9_to_np(ny/2.+ny/3.), ds9_to_np(ny), x_min_ind, x_max_ind)
    popt_soph, pcov_soph = curve_fit(lambda r, I0, h: disk_exp(r, I0, h, xc), r_plane, I_plane, p0=(max(I_plane), nx/10.) )
    h3 = popt_soph[1]    

    return np.array(rr), np.array(I_image), np.array(I_mask), np.array(I_sigma), np.array(r_no_mask), np.array(I_no_mask), h1, h2, h3
        



def psf_slice(input_psf, axis='yaxis'):
    imageHDU = pyfits.open(input_psf)[0]
    image = imageHDU.data
    ny,nx = image.shape
    if axis=='yaxis': 
        if nx%2:
            x0 = int(nx/2.+1.)
        else:
            x0 = int(nx/2.+0.5)
        array = image[:,x0-1]
        if os.path.exists('psf_yline.fits'):
            os.remove('psf_yline.fits')
        subprocess.call('%sfitscopy \'%s[%i:%i,%i:%i]\' psf_yline.fits' % (fitscopy,input_psf,x0,x0,1,ny), shell=True)

    else: 
        if ny%2:
            y0 = int(ny/2.+1.)
        else:
            y0 = int(ny/2.+0.5)
        array = image[y0-1,:]
        if os.path.exists('psf_xline.fits'):
            os.remove('psf_xline.fits')
        subprocess.call('%sfitscopy \'%s[%i:%i,%i:%i]\' psf_xline.fits' % (fitscopy,input_psf,1,nx,y0,y0), shell=True)            

    return array/sum(array)


def disk_edge_soph(z, I0d, z0, z_c):
    #*** For edge-on disk SB in mag/arcsec^2 (along z-axis). n - index of this law. ***
    return I0d * (1.0 / np.cosh(1.*np.fabs(z-z_c)/(2.*z0)) ) **(2./1.)




def crea_imf(L0, n, z0, y0, n_fixed='fixed', ron=None, gain=None):
    f = open('input_imf.txt','w')
    if ron is not None and gain is not None:
        print >>f,  "GAIN	%.2f" % (gain)
        print >>f,  "READNOISE	%.2f" % (ron)
    '''
    print >>f,  "X0	0.5	fixed"
    print >>f,  "Y0	%.2f	fixed" % (y0)

    print >>f,  "FUNCTION Z_eon"
    print >>f,  "I_0	%.10f" % (I0)
    print >>f,  "n	%.2f	0.1,100" % (n)
    print >>f,  "z_0	%.2f" % (z0)
    '''
    s = """
X0		%.2f	fixed
Y0		%.2f   
FUNCTION EdgeOnDisk
PA		%f	fixed
L_0		%f
h		%f	fixed
n		%f	%s
z_0		%f    
""" % (0.5,y0,90.,L0,10.,n,n_fixed,z0)  #### n is fixed!!!!!!!
    print >>f, s
    f.close()
    #exit()


def crea_imf_double(L01, n1, z01, y0, L02, n2, z02, n_fixed='fixed', ron=None, gain=None):
    f = open('input_imf.txt','w')
    if ron is not None and gain is not None:
        print >>f,  "GAIN	%.2f" % (gain)
        print >>f,  "READNOISE	%.2f" % (ron)
    s = """
X0		%.2f	fixed
Y0		%.2f	fixed    
FUNCTION EdgeOnDisk
PA		%f	fixed
L_0		%f
h		%f	fixed
n		%f	%s
z_0		%f

FUNCTION EdgeOnDisk
PA		%f	fixed
L_0		%f
h		%f	fixed
n		%f	%s
z_0		%f
""" % (0.5,y0,90.,L01,10.,n1,n_fixed,z01,90.,L02,10.,n1,n_fixed,z02)
    print >>f, s
    f.close()


def crea_imf_tripple(L01, n1, z01, y0, L02, n2, z02, L03, n3, z03, n_fixed='fixed', ron=None, gain=None):
    f = open('input_imf.txt','w')
    if ron is not None and gain is not None:
        print >>f,  "GAIN	%.2f" % (gain)
        print >>f,  "READNOISE	%.2f" % (ron)
    s = """
X0		%.2f	fixed
Y0		%.2f	fixed    
FUNCTION EdgeOnDisk
PA		%f	fixed
L_0		%f
h		%f	fixed
n		%f	%s
z_0		%f

FUNCTION EdgeOnDisk
PA		%f	fixed
L_0		%f
h		%f	fixed
n		%f	%s
z_0		%f

FUNCTION EdgeOnDisk
PA		%f	fixed
L_0		%f
h		%f	fixed
n		%f	%s
z_0		%f
""" % (0.5,y0,90.,L01,10.,n1,n_fixed,z01,90.,L02,10.,n2,n_fixed,z02,90.,L03,10.,n3,n_fixed,z03)
    print >>f, s
    f.close()

def read_imf(imf_file):
  f = open(imf_file,'r')
  L0 = 99999.
  for line in f:
    if "#   Reduced value:" in line:
      chi2 = float(line.split()[3])
    if "Y0\t" in line:
        y0 = float(line.split()[1])
    if 'L_0' in line:
      L0 = float(line.split()[1])
    if 'n' in line and L0!=99999.:
      n = float(line.split()[1])
    if 'z_0' in line:
      z0 = float(line.split()[1])
  f.close()
  return y0,L0,n,z0,chi2


def read_imf_double(imf_file):
  f = open(imf_file,'r')
  func = 0
  n1 = 99999.
  for line in f:
    if "#   Reduced value:" in line:
      chi2 = float(line.split()[3])
    if "FUNCTION EdgeOnDisk" in line and func==0:
      func=1
    if 'L_0' in line and func==1:
      L01 = float(line.split()[1])      
    if 'n\t' in line and func==1:
      n1 = float(line.split()[1])
    if 'z_0' in line and func==1:
      z01 = float(line.split()[1])
    if "FUNCTION EdgeOnDisk" in line and func==1 and n1!=99999.:
      func=2
    if 'L_0' in line and func==2:
      L02 = float(line.split()[1])   
    if 'n\t' in line and func==2:
      n2 = float(line.split()[1])
    if 'z_0' in line and func==2:
      z02 = float(line.split()[1])      
  f.close()
  if z01<z02:
      return L01,n1,z01,L02,n2,z02,chi2
  else:
      return L02,n2,z02,L01,n1,z01,chi2        

def read_imf_tripple(imf_file):
  f = open(imf_file,'r')
  func = 0
  n1 = 99999.
  n2 = 99999.
  for line in f:
    if "#   Reduced value:" in line:
      chi2 = float(line.split()[3])
    if "FUNCTION EdgeOnDisk" in line and func==0:
      func=1
    if 'L_0' in line and func==1:
      L01 = float(line.split()[1])      
    if 'n\t' in line and func==1:
      n1 = float(line.split()[1])
    if 'z_0' in line and func==1:
      z01 = float(line.split()[1])
      
    if "FUNCTION EdgeOnDisk" in line and func==1 and n1!=99999.:
      func=2
    if 'L_0' in line and func==2:
      L02 = float(line.split()[1])   
    if 'n\t' in line and func==2:
      n2 = float(line.split()[1])
    if 'z_0' in line and func==2:
      z02 = float(line.split()[1])    
      
    if "FUNCTION EdgeOnDisk" in line and func==2 and n2!=99999.:
      func=3
    if 'L_0' in line and func==3:
      L03 = float(line.split()[1])   
    if 'n\t' in line and func==3:
      n3 = float(line.split()[1])
    if 'z_0' in line and func==3:
      z03 = float(line.split()[1])   

  f.close()
  RES= [[L01,n1,z01],[L02,n2,z02],[L03,n3,z03]]
  zzz = [z01,z02,z03]
  inds_sorted = sorted(range(len(zzz)), key=lambda k: zzz[k])
  
  return RES[inds_sorted[0]],RES[inds_sorted[1]],RES[inds_sorted[2]],chi2 


def run_galf(R, L0, n, z0, y0, h, sigma_image, psf_image, mask_image):
  try:
    crea_galfit_input('vertical_line.fits', sigma_image, psf_image, mask_image, -R, y0, -2.5*math.log10(L0*(2.*h)), z01=2.*z0, h1=h,  mu0d2=None, z02=None, h2=None)
    
    if os.path.exists('galfit.01'):
        os.remove('galfit.01')
        
    subprocess.call('%sgalfit %s' % (galfit_path, 'input_galf.txt'), shell=True, stdout=FNULL)
    if os.path.exists('galfit.01'):
        status = 0
        shutil.copy('galfit.01','galfit_lm.01')
        os.remove('fit.log')
        os.remove('model.fits')
    else:
        status = 1
    #exit()
    return status    
  except:
      return 1

def run_galf_double(R, L01, n1, z01, y0, h1, L02, n2, z02, h2, sigma_image, psf_image, mask_image, ron, gain):
  try:
    crea_galfit_input('vertical_line.fits', sigma_image, psf_image, mask_image, -R, y0, -2.5*math.log10(L01*(2.*h1)), z01=2.*z01, h1=h1,  mu0d2=-2.5*math.log10(L02*(2.*h2)), z02=2.*z02, h2=h2)
    #exit()

    if os.path.exists('galfit.01'):
        os.remove('galfit.01')
    #exit()    
    subprocess.call('%sgalfit %s' % (galfit_path, 'input_galf.txt'), shell=True)#, stdout=FNULL)
    if os.path.exists('galfit.01'):
        status = 0
        shutil.copy('galfit.01','galfit_lm.01')
        os.remove('fit.log')
        os.remove('model.fits')
    else:
        status = 1
    #exit()
    return status 
  except:
      return 1

def run_galf_tripple(L01, n1, z01, y0, h1, L02, n2, z02, h2, L03, n3, z03, h3, sigma_image, psf_image, mask_image, ron, gain):
  try:
    crea_galfit_input_tripple(input_image, sigma_image, psf_image, mask_image, 0.5, y0, -2.5*math.log10(L01*(2.*h)), z01=2.*z01,  mu0d2=-2.5*math.log10(L02*(2.*h)), z02=2.*z02, mu0d3=-2.5*math.log10(L03*(2.*h)), z03=2.*z03)  


    if os.path.exists('galfit.01'):
        os.remove('galfit.01')
        
    subprocess.call('%sgalfit %s' % (galfit_path, 'input_galf.txt'), shell=True, stdout=FNULL)
    if os.path.exists('galfit.01'):
        status = 0
        shutil.copy('galfit.01','galfit_lm.01')
        os.remove('fit.log')
        os.remove('model.fits')
    else:
        status = 1
    #exit()
    return status 
  except:
      return 1

def run_imf(L0, n, z0, y0, n_fixed, sigma_image, psf_image, mask_image, ron, gain):
  try:
    crea_imf(L0, n, z0, y0, n_fixed=n_fixed, ron=ron, gain=gain)

    imfit_line = '%simfit vertical_line.fits -c input_imf.txt --ftol 0.00001 --max-threads 1' % (imfit_path)
    if sigma_image!=None:
        imfit_line += ' --noise sigma_line.fits '

    if psf_image!=None:
        imfit_line += ' --psf %s ' % (psf_image)       

    if mask_image!=None:
        imfit_line += ' --mask mask_line.fits '   


    subprocess.call(imfit_line, shell=True, stdout=FNULL)
    try:    
        shutil.copy('bestfit_parameters_imfit.dat','bestfit_parameters_imfit_lm.dat')
        status = 0
    except:
        status = 1
    return status
  except:
      return 1

def run_imf_double(L01, n1, z01, y0, L02, n2, z02, n_fixed, sigma_image, psf_image, mask_image, ron, gain):
  try:
    #print('Input pars:',L01, n1, z01, y0, L02, n2, z02)
    crea_imf_double(L01, n1, z01, y0, L02, n2, z02, n_fixed=n_fixed, ron=ron, gain=gain)

    imfit_line = '%simfit vertical_line.fits -c input_imf.txt --ftol 0.00001 --max-threads 1' % (imfit_path)
    if sigma_image!=None:
        imfit_line += ' --noise sigma_line.fits '

    if psf_image!=None:
        imfit_line += ' --psf psf_yline.fits '        

    if mask_image!=None:
        imfit_line += ' --mask mask_line.fits '   


    subprocess.call(imfit_line, shell=True, stdout=FNULL)
    try:
        shutil.copy('bestfit_parameters_imfit.dat','bestfit_parameters_imfit_lm.dat')
        status = 0
    except:
        status = 1
    #exit()
    return status
  except:
      return 1


def run_imf_tripple(L01, n1, z01, y0, L02, n2, z02, L03, n3, z03, n_fixed, sigma_image, psf_image, mask_image, ron, gain):
  try:
    #print('Input pars:',L01, n1, z01, y0, L02, n2, z02)
    crea_imf_tripple(L01, n1, z01, y0, L02, n2, z02, L03, n3, z03, n_fixed=n_fixed, ron=ron, gain=gain)

    imfit_line = '%simfit vertical_line.fits -c input_imf.txt --ftol 0.00001 --max-threads 1' % (imfit_path)
    if sigma_image!=None:
        imfit_line += ' --noise sigma_line.fits '

    if psf_image!=None:
        imfit_line += ' --psf psf_yline.fits '        

    if mask_image!=None:
        imfit_line += ' --mask mask_line.fits '   


    subprocess.call(imfit_line, shell=True, stdout=FNULL)
    try:
        shutil.copy('bestfit_parameters_imfit.dat','bestfit_parameters_imfit_lm.dat')
        status = 0
    except:
        status = 1
    #exit()
    return status
  except:
      return 1    


def single_disc_fit(R, L0, n, z0, y0, h, n_fixed, input_image, sigma_image, psf_image, mask_image, ron, gain):
    if code=='imfit':
        # Levenberg-Marquardt:
        status = run_imf(L0, n, z0, y0, n_fixed, sigma_image, psf_image, mask_image, ron, gain)

        if status == 0:
            Y0,L00,N,Z0,chi2 = read_imf('bestfit_parameters_imfit_lm.dat')
        else:
            L00=float('nan'); Y0=float('nan'); N=float('nan'); Z0=float('nan')
    elif code=='galfit':
        # Levenberg-Marquardt:
        status = run_galf(R, L0, n, z0, y0, h, sigma_image, psf_image, mask_image)

        if status == 0:
            X0,Y0,L00,Z0,H,chi2 = read_galf('galfit_lm.01')
            N = 1.
            make_model_ima_galfit.main(input_image, 'galfit_lm.01', composed_model_file = 'composed_model.fits', subtract_sky=True, galfitPath=galfit_path)
        else:
            L00=float('nan'); Y0=float('nan'); N=float('nan'); Z0=float('nan')        
        
    return L00,Y0,N,Z0



def double_disc_fit(R, L01, n1, z01, y0, h1, L02, n2, z02, h2, n_fixed, input_image, sigma_image, psf_image, mask_image, ron, gain):
    z01 = np.fabs(z01)
    z02 = np.fabs(z02)
    if code=='imfit':
        # Levenberg-Marquardt:
        status = run_imf_double(L01, n1, z01, y0, L02, n2, z02, n_fixed, sigma_image, psf_image, mask_image, ron, gain)

        if status == 0:
            L01,N1,Z01,L02,N2,Z02,chi2 = read_imf_double('bestfit_parameters_imfit_lm.dat')
        else:
            L01=float('nan'); N1=float('nan'); Z01=float('nan'); L02=float('nan'); N2=float('nan'); Z02=float('nan'); chi2=float('nan')
    elif code=='galfit':
        status = run_galf_double(R, L01, n1, z01, y0, h1, L02, n2, z02, h2, sigma_image, psf_image, mask_image, ron, gain)
        if status == 0:
            L01,h1,Z01,L02,h2,Z02,chi2 = read_galf_double('galfit_lm.01')
            N1 = 1
            N2 = 1
            make_model_ima_galfit.main(input_image, 'galfit_lm.01', composed_model_file = 'composed_model.fits', subtract_sky=True, galfitPath=galfit_path)
        else:
            L01=float('nan'); N1=float('nan'); Z01=float('nan'); L02=float('nan'); N2=float('nan'); Z02=float('nan'); chi2=float('nan')    
    
    return L01,N1,Z01,L02,N2,Z02
    

def tripple_disc_fit(L01, n1, z01, y0, h1, L02, n2, z02, h2, L03, n3, z03, h3, n_fixed, sigma_image, psf_image, mask_image, ron, gain):
    z01 = np.fabs(z01)
    z02 = np.fabs(z02)
    z03 = np.fabs(z03)
    if code=='imfit':
        # Levenberg-Marquardt:
        status = run_imf_tripple(L01, n1, z01, y0, L02, n2, z02, L03, n3, z03, n_fixed, sigma_image, psf_image, mask_image, ron, gain)

        if status == 0:
            [L01,N1,Z01],[L02,N2,Z02],[L03,N3,Z03],chi2 = read_imf_tripple('bestfit_parameters_imfit_lm.dat')
        else:
            L01=float('nan'); N1=float('nan'); Z01=float('nan'); L02=float('nan'); N2=float('nan'); Z02=float('nan'); L03=float('nan'); N3=float('nan'); Z03=float('nan'); chi2=float('nan')
    elif code=='galfit':
        status = run_galf_tripple(L01, n1, z01, y0, h1, L02, n2, z02, h2, L03, n3, z03, h3, sigma_image, psf_image, mask_image, ron, gain)
        if status == 0:
            [L01,h1,Z01],[L02,h2,Z02],[L03,h3,Z03],chi2 = read_galf_tripple('galfit_lm.01')
            N1 = 1
            N2 = 1
            N3 = 1
        else:
            L01=float('nan'); N1=float('nan'); Z01=float('nan'); L02=float('nan'); N2=float('nan'); Z02=float('nan'); L03=float('nan'); N3=float('nan'); Z03=float('nan'); chi2=float('nan')    
    
    return L01,N1,Z01,L02,N2,Z02,L03,N3,Z03

def find_single_disc_guess(rr, I_image, I0, y_min_ind):
        #print(rr, I_image, I0, y_min_ind)
        
        # Find first guess on the parameters:
        try:
            y0 = float(np.where(I_image==I0)[0])+y_min_ind
        except:
            y0 = float(np.where(I_image==I0)[0][0])+y_min_ind
        try:
            #plt.plot(rr, np.log10(I_image))
            #plt.show()
            #exit()
            popt_soph, pcov_soph = curve_fit( disk_edge_soph, rr, I_image, p0=(I0, len(I_image)/4., y0))
            L0 = popt_soph[0]
            z0 = np.fabs(popt_soph[1])
            y0 = popt_soph[2]
            #print(math.log10(L0),z0,y0)
            N = 1
            #plt.plot(rr, np.log10(disk_edge_soph(rr, L0, z0, y0)))
            #plt.show()
            #exit()
            
        except:
            L0 = I0
            z0 = len(I_image)/4.
            y0 = y0
            N = 1

        if z0>len(I_image):
            L0 = I0
            z0 = len(I_image)/4.
            y0 = y0
            N = 1
        return L0, z0, y0, N


def write_profile(R, rr, I_image, disc_comps, components=[1]):
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        hdu_model = pyfits.open('composed_model.fits')
        
        model_data = []
        for comp in components:
            model_data.append(hdu_model[comp].data) 
        
        Model_data = []
        for kk in range(len(rr)):
            for ii in range(len(model_data)):
                if ii==0:
                    ss = '%.10f' % (model_data[ii][kk])
                else:
                    ss = ss + '\t%.10f' % (model_data[ii][kk])
            Model_data.append(ss)
        
        
        
        # Write profile to file:
        ff = open('vertical_profiles_%s_%s.txt' % ('1', disc_comps),'a') # TODO
        ff.write('#RADIUS:%f\n' % (R))

        for kk in range(len(rr)):
            ff.write('%i\t%.10f\t%s\n' % (rr[kk],I_image[kk],Model_data[kk]))
        ff.write('#END\n')
        ff.close()
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
 


def fit_each_vert_slice(R, image_data, mask_data, sigma_data, psf_image, x_min_ind, x_max_ind,  y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up, ron, gain, n_fixed='fixed', L0=None, N=None, z0=None, y0=None, h=None, L02=None, N2=None, z02=None, h2=None, L03=None, N3=None, z03=None, h3=None,  disc_comps = 'single', h1_given=None, h2_given=None, z01_given=None, z02_given=None):


    rr, I_image, I_mask, I_sigma, r_no_mask, I_no_mask, H1, H2, H3 = crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind, y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up)     
    
    
    if h is None:
        h = H1
    
    if h2 is None:
        h2 = H2
    
    if h3 is None:
        h3 = H3

        
    if h1_given is not None:
        h = h1_given

    if h2_given is not None:
        h2 = h2_given

    if h1_given is not None:
        h3 = h2_given        
    
    I0 = max(I_image)

    
    input_image, sigma_image, mask_image = prepare_slice_images(I_image, I_sigma, I_mask)
    

    if L0 is None and z0 is None and N is None:
        # Find first guess on the parameters:
        L0, z0, y0, N = find_single_disc_guess(r_no_mask, I_no_mask, I0, 0.)

    if z01_given is not None:
        z0 = z01_given

    # Single disc fit:
    if disc_comps=='single':
        L0,Y0,N,Z0 = single_disc_fit(R, L0, N, z0, y0, h, n_fixed, input_image, sigma_image, psf_image, mask_image, ron, gain)
        write_profile(R, rr, I_image, disc_comps, components=[1])
 
        return L0,Y0,N,Z0,h,h2,h3


    elif disc_comps=='double1' or disc_comps=='tripple':
        '''
        if L02 is None and N2 is None and z02 is None: # WARNING: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            L02 = L0/30.
            N2 = N
            z02 = 3.*z0
        '''
        if L02 is None:
            L01,L02,z01,z02 = find_double_disc_guess(r_no_mask, I_no_mask, y0, h, h2)  #### WARNING: DOES NOT WORK PROPERLY!!!
        else:
            L01 = L0
            z01 = z0
        
        if z01>len(I_image) or z02>len(I_image):
            #####################
            # Double disc fit:
            L01=L0
            z01=z0
            L02=L0/10.
            z02=z0*3.
            #####################
        N=1
        
        if z01_given is not None:
            z01 = z01_given
        if z02_given is not None:
            z02 = z02_given

        L01,N1,Z01,L02,N2,Z02 = double_disc_fit(R, L01, N, z01, y0, h, L02, N, z02, h2, n_fixed, input_image, sigma_image, psf_image, mask_image, ron, gain)
        if not np.isnan(L01):
            write_profile(R, rr, I_image, disc_comps, components=[1,4,5]) # TODO: CHECK
        
        if disc_comps=='double1':
            return L01,N1,Z01,L02,N2,Z02,h,h2,h3

    if disc_comps=='tripple':
            I_1 = disk_edge_soph(rr, 2.*h*L01, Z01, y0)
            I_2 = disk_edge_soph(rr, 2.*h2*L02, Z02, y0)
            I_sum = I_1 + I_2
            I_3 = I_image - I_sum
            L_3 = np.trapz(I_3, x=rr)
            Ltot = np.trapz(I_image, x=rr)

            if L_3/Ltot>0.01:
                if L03 is None:
                        # Find first guess on the parameters:
                        L03, z03, y03, N3 = find_single_disc_guess(rr, I_3, np.max(I_3), max(rr)/2.)
                        #print(L03, z03, y03, N3)
                
                L01,N1,Z01,L02,N2,Z02,L03,N3,Z03 = tripple_disc_fit(L01, N1, Z01, y0, L02, N2, Z02, L03, N3, z03, n_fixed, sigma_image, psf_image, mask_image, ron, gain)
                #print(L01,N1,Z01,L02,N2,Z02,L03,N3,Z03)
                write_profile(R, rr, I_image, disc_comps, components=[1,4,5,6]) # TODO: CHECK
                return L01,N1,Z01,L02,N2,Z02,L03,N3,Z03
            else:
                return L01,N1,Z01,L02,N2,Z02,float('nan'),float('nan'),float('nan')
            

def prepare_slice_images(I_image, I_sigma, I_mask):
    #### Create slice images: ###
    # Data image:
    fict_data = np.random.random((len(I_image),2))
    fict_data[:,0] = I_image
    fict_data[:,1] = np.zeros(len(I_image))
    hdu = pyfits.PrimaryHDU(fict_data)
    hdu.writeto('slice.fits',clobber=True)

    if os.path.exists('vertical_line.fits'):
        os.remove('vertical_line.fits')
    subprocess.call('%sfitscopy \'%s[%i:%i,%i:%i]\' vertical_line.fits' % (fitscopy,'slice.fits',1,1,1,len(I_image)), shell=True)
    os.remove('slice.fits')
    input_image = 'vertical_line.fits'
    
    
    # Sigma image:
    fict_data = np.random.random((len(I_sigma),2))
    fict_data[:,0] = I_sigma
    fict_data[:,1] = np.zeros(len(I_sigma))
    hdu = pyfits.PrimaryHDU(fict_data)
    hdu.writeto('slice_sigma.fits',clobber=True)      

    if os.path.exists('sigma_line.fits'):
        os.remove('sigma_line.fits')
    subprocess.call('%sfitscopy \'%s[%i:%i,%i:%i]\' sigma_line.fits' % (fitscopy,'slice_sigma.fits',1,1,1,len(I_image)), shell=True)
    os.remove('slice_sigma.fits')
    sigma_image = 'sigma_line.fits'
    

    # Mask image:
    I_mask = np.multiply(I_mask,1)
    fict_data = np.random.random((len(I_mask),2))
    fict_data[:,0] = I_mask
    fict_data[:,1] = np.zeros(len(I_mask))
    hdu = pyfits.PrimaryHDU(fict_data)
    hdu.writeto('slice_mask.fits',clobber=True)      

    if os.path.exists('mask_line.fits'):
        os.remove('mask_line.fits')
    subprocess.call('%sfitscopy \'%s[%i:%i,%i:%i]\' mask_line.fits' % (fitscopy,'slice_mask.fits',1,1,1,len(I_image)), shell=True)
    os.remove('slice_mask.fits')
    mask_image = 'mask_line.fits'
    #### Create slice images END ###   
    return input_image,sigma_image,mask_image



def fit_one_disc(R, r_no_mask, I_no_mask, I_image, I_sigma, I_mask, h1, n_fixed, psf_image, ron, gain, z01=None, Y0=None):   
    I0 = max(I_no_mask)

    input_image,sigma_image,mask_image = prepare_slice_images(I_image, I_sigma, I_mask)

    L0, z0, y0, N = find_single_disc_guess(r_no_mask, I_no_mask, I0, 0.)
    
    if Y0 is not None:
        if abs(y0-Y0)/Y0>0.2:
            y0=Y0
    
    if z01 is not None:
        z0 = z01
        
    L0,Y0,N,Z0 = single_disc_fit(R, L0, N, z0, y0, h1, n_fixed, input_image, sigma_image, psf_image, mask_image, ron, gain)
    return L0,Y0,N,Z0


def write_profile_func(R, rr, I_image, total_model, model_thin, model_thick):
    # Write profile to file:
    ff = open('vertical_profiles_%s_%s.txt' % ('1', 'double'),'a') # TODO
    ff.write('#RADIUS:%f\n' % (R))

    for kk in range(len(rr)):
            ff.write('%i\t%.10f\t%.10f\t%.10f\t%.10f\n' % (rr[kk],I_image[kk],total_model[kk],model_thin[kk],model_thick[kk]))
    ff.write('#END\n')
    ff.close()    

def fit_each_vert_slice_double(R, image_data, mask_data, sigma_data, psf_image, x_min_ind, x_max_ind,  y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up, ron, gain, n_fixed='fixed', error_est = False, h1_given=None, h2_given=None, z01_given=None, z02_given=None):
    #print(y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up)
    #exit()
    L01 = float('nan')
    L01_err = float('nan')
    N1 = float('nan')
    N1_err = float('nan')
    Z01 = float('nan')
    Z01_err = float('nan')
    L02 = float('nan')
    L02_err = float('nan')
    N2 = float('nan')
    N2_err = float('nan')
    Z02 = float('nan')
    Z02_err = float('nan')


    # Thick disc:
    rr_2, I_image_2, I_mask_2, I_sigma_2, r_no_mask_2, I_no_mask_2, h1, h2, h3 = crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind, y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up)     
    
   
    if h2_given is not None:
        h2 = h2_given

    
    L02,Y02,N2,Z02 = fit_one_disc(R, r_no_mask_2, I_no_mask_2, I_image_2, I_sigma_2, I_mask_2, h2, n_fixed, psf_image, ron, gain, z01=z02_given)


    #plt.plot(rr_2, -np.log10(I_image_2))
    #plt.plot(rr_2, -np.log10(disk_edge_soph(rr_2, 2.*h2*L02, Z02, Y02)))
    #plt.show()
    #exit()
    
    
    
    
    
    try:
        hdu_model = pyfits.open('composed_model.fits')
        model_thick = hdu_model[1].data[:,0] 
    except:
        return L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err
    
 
    # Thin disc
    rr, I_image, I_mask, I_sigma, r_no_mask, I_no_mask, h1, h2, h3 = crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind, y_min_ind_dn, y_max_ind_up, y_min_ind_up, y_max_ind_up)     #### TODO:CHECK!!!   
    #print(y_min_ind_dn, y_max_ind_up, y_min_ind_up, y_max_ind_up)

    
    Ltot = np.trapz(I_image, x=rr)
    I_image_1 = I_image - disk_edge_soph(rr, 2.*h2*L02, Z02, Y02) # TODO: Y02 should be real y0!!!!!!!!!!!!!!!!!!!!!!
    

    L_1 = np.trapz(I_image_1, x=rr)
    #plt.plot(rr_2, -np.log10(I_image_2))
    #plt.plot(rr, -np.log10(I_image_1))
    #plt.plot(rr, -np.log10(I_image))
    #plt.plot(rr, -np.log10(disk_edge_soph(rr, 2.*h2*L02, Z02, Y02)))
    #plt.show()
    #exit()


    if abs(L_1/Ltot)<0.01: #### TODO: CHECK
        # Just one disc
        total_model = model_thick
        model_thin = np.zeros(len(rr))

        if h1_given is not None:
            h1 = h1_given

        L01,Y01,N1,Z01 = fit_one_disc(R, r_no_mask, I_no_mask, I_image, I_sigma, I_mask, h1, n_fixed, psf_image, ron, gain, z01=z01_given)
        
        write_profile_func(R, rr, I_image, total_model, model_thin, model_thick)
        
        if error_est:
            L001 = []; Y001 = []; NN1 = []; Z001 = []   
            if abs(x_max_ind - x_min_ind)<=n_times:
                slice_list = range(x_min_ind, x_max_ind-1)
            else:
                slice_list = random.sample(range(x_min_ind, x_max_ind-1),n_times)
            
            for k in slice_list:
                # Thin disc
                rr, I_image, I_mask, I_sigma, r_no_mask, I_no_mask, h1, h2, h3 = crea_vert_average_cut(image_data, mask_data, sigma_data, k, k+1, y_min_ind_dn, y_min_ind_dn+2, y_min_ind_dn+1, y_max_ind_up)     #### TODO:CHECK!!!   

                if h1_given is not None:
                    h1 = h1_given
                
                p1,p2,p3,p4 = fit_one_disc(R, r_no_mask, I_no_mask, I_image, I_sigma, I_mask, h1, n_fixed, psf_image, ron, gain, z01=z01_given)
                L001.append(p1)
                Y001.append(p2)
                NN1.append(p3)
                Z001.append(p4)                   

            L01_err = np.std(L001)
            Y01_err = np.std(Y001)
            N1_err = np.std(NN1)
            Z01_err = np.std(Z001)    
            
    
    
    else:
        if h1_given is not None:
            h1 = h1_given
            
        L01,Y01,N1,Z01 = fit_one_disc(R, r_no_mask, I_no_mask, I_image_1, I_sigma, I_mask, h1, n_fixed, psf_image, ron, gain, z01=z01_given)
        hdu_model = pyfits.open('composed_model.fits')
        model_thin = hdu_model[1].data[:,0]    
        total_model = model_thin + model_thick
        
        write_profile_func(R, rr, I_image, total_model, model_thin, model_thick)
        
        if error_est:
            L001 = []; Y001 = []; NN1 = []; Z001 = [] 
            L002 = []; Y002 = []; NN2 = []; Z002 = []      
            if abs(x_max_ind - x_min_ind)<=n_times:
                slice_list = range(x_min_ind, x_max_ind-1)
            else:
                slice_list = random.sample(range(x_min_ind, x_max_ind-1),n_times)
            
            for k in slice_list:
                # Thick disc:
                rr_2, I_image_2, I_mask_2, I_sigma_2, r_no_mask_2, I_no_mask_2, h1, h2, h3 = crea_vert_average_cut(image_data, mask_data, sigma_data, k, k+1, y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up)     

                if h2_given is not None:
                    h2 = h2_given
                
                p1,p2,p3,p4 = fit_one_disc(R, r_no_mask_2, I_no_mask_2, I_image_2, I_sigma_2, I_mask_2, h2, n_fixed, psf_image, ron, gain, z01=z02_given)
                L002.append(p1)
                Y002.append(p2)
                NN2.append(p3)
                Z002.append(p4)                
                

                # Thin disc
                rr, I_image, I_mask, I_sigma, r_no_mask, I_no_mask, h1, h2, h3 = crea_vert_average_cut(image_data, mask_data, sigma_data, k, k+1, y_min_ind_dn, y_min_ind_dn+2, y_min_ind_dn+1, y_max_ind_up)     #### TODO:CHECK!!!   

                I_image_1 = I_image - disk_edge_soph(rr, 2.*h2*p1, p4, p2) # TODO: Y02 should be real y0!!!!!!!!!!!!!!!!!!!!!!    

                if h1_given is not None:
                    h1 = h1_given
                
                p1,p2,p3,p4 = fit_one_disc(R, r_no_mask, I_no_mask, I_image_1, I_sigma, I_mask, h1, n_fixed, psf_image, ron, gain, z01=z01_given)
                L001.append(p1)
                Y001.append(p2)
                NN1.append(p3)
                Z001.append(p4)                   

            L01_err = np.std(L001)
            Y01_err = np.std(Y001)
            N1_err = np.std(NN1)
            Z01_err = np.std(Z001)    

            L02_err = np.std(L002)
            Y02_err = np.std(Y002)
            N2_err = np.std(NN2)
            Z02_err = np.std(Z002)    
    



    return L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err



def fit_each_vert_slice_double1(R, image_data, mask_data, sigma_data, psf_image, x_min_ind, x_max_ind,  y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up, ron, gain, n_fixed='fixed', error_est = False, h1_given=None, h2_given=None, z01_given=None, z02_given=None, yc=None):

    L01 = float('nan')
    L01_err = float('nan')
    N1 = float('nan')
    N1_err = float('nan')
    Z01 = float('nan')
    Z01_err = float('nan')
    L02 = float('nan')
    L02_err = float('nan')
    N2 = float('nan')
    N2_err = float('nan')
    Z02 = float('nan')
    Z02_err = float('nan')


    # Thin disc:
    rr_1, I_image_1, I_mask_1, I_sigma_1, r_no_mask_1, I_no_mask_1, h1, h2, h3 = crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind, y_max_ind_dn, y_max_ind_dn+10, y_max_ind_dn+10, y_min_ind_up)     
    
   
    if h1_given is not None:
        h1 = h1_given

    
    L01,Y01,N1,Z01 = fit_one_disc(R, r_no_mask_1, I_no_mask_1, I_image_1, I_sigma_1, I_mask_1, h1, n_fixed, psf_image, ron, gain, z01=z01_given, Y0=yc)


    #plt.plot(rr_1, -np.log10(I_image_1))
    #plt.plot(rr_1, -np.log10(disk_edge_soph(rr_1, 2.*h1*L01, Z01, Y01)))
    #plt.show()
    #exit()
    
    
    
    
    

    hdu_model = pyfits.open('composed_model.fits')
    model_thin = hdu_model[1].data[:,0]     
    
 
    # Thick disc
    rr, I_image, I_mask, I_sigma, r_no_mask, I_no_mask, h1, h2, h3 = crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind, y_min_ind_dn, y_max_ind_up, y_min_ind_up, y_max_ind_up)     #### TODO:CHECK!!!   
    

    
    
    Ltot = np.trapz(I_image, x=rr)
    I_image_2 = I_image - disk_edge_soph(rr, 2.*h1*L01, Z01, Y01) # TODO: Y02 should be real y0!!!!!!!!!!!!!!!!!!!!!!
    
    

    L_2 = np.trapz(I_image_2, x=rr)
    #plt.plot(rr_1, -np.log10(I_image_1))
    #plt.plot(rr, -np.log10(I_image_2))
    #plt.plot(rr, -np.log10(I_image))
    #plt.plot(rr, -np.log10(disk_edge_soph(rr, 2.*h2*L02, Z02, Y02)))
    #plt.show()
    #exit()
    print('\t',L_2/Ltot)
    #exit()

    if abs(L_2/Ltot)<0.01 or abs(L_2/Ltot)>100.: #### TODO: CHECK
        print('HERE')
        # Just one disc
        total_model = model_thin
        model_thick = np.zeros(len(rr))

        if h1_given is not None:
            h1 = h1_given

        L01,Y01,N1,Z01 = fit_one_disc(R, r_no_mask, I_no_mask, I_image, I_sigma, I_mask, h1, n_fixed, psf_image, ron, gain, z01=z01_given, Y0=yc)
        
        write_profile_func(R, rr, I_image, total_model, model_thin, model_thick)
        
        if error_est:
            L001 = []; Y001 = []; NN1 = []; Z001 = []   
            if abs(x_max_ind - x_min_ind)<=n_times:
                slice_list = range(x_min_ind, x_max_ind-1)
            else:
                slice_list = random.sample(range(x_min_ind, x_max_ind-1),n_times)
            
            for k in slice_list:
                # Thin disc
                rr, I_image, I_mask, I_sigma, r_no_mask, I_no_mask, h1, h2, h3 = crea_vert_average_cut(image_data, mask_data, sigma_data, k, k+1, y_min_ind_dn, y_min_ind_dn+2, y_min_ind_dn+1, y_max_ind_up)     #### TODO:CHECK!!!   

                if h1_given is not None:
                    h1 = h1_given
                
                p1,p2,p3,p4 = fit_one_disc(R, r_no_mask, I_no_mask, I_image, I_sigma, I_mask, h1, n_fixed, psf_image, ron, gain, z01=z01_given, Y0=yc)
                L001.append(p1)
                Y001.append(p2)
                NN1.append(p3)
                Z001.append(p4)                   

            L01_err = np.std(L001)
            Y01_err = np.std(Y001)
            N1_err = np.std(NN1)
            Z01_err = np.std(Z001)    
            
    
    
    else:
        if h2_given is not None:
            h2 = h2_given
        
        rr_2, I_image_22, I_mask_2, I_sigma_2, r_no_mask_2, I_no_mask_2, h1, h2, h3 = crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind, y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up)  
        
        L02,Y02,N2,Z02 = fit_one_disc(R, r_no_mask_2, I_no_mask_2, I_image_2, I_sigma_2, I_mask_2, h2, n_fixed, psf_image, ron, gain, z01=z02_given, Y0=yc)
        hdu_model = pyfits.open('composed_model.fits')
        model_thick = hdu_model[1].data[:,0]    
        total_model = model_thin + model_thick
        
        write_profile_func(R, rr, I_image, total_model, model_thin, model_thick)
        
        if error_est:
            L001 = []; Y001 = []; NN1 = []; Z001 = [] 
            L002 = []; Y002 = []; NN2 = []; Z002 = []      
            if abs(x_max_ind - x_min_ind)<=n_times:
                slice_list = range(x_min_ind, x_max_ind-1)
            else:
                slice_list = random.sample(range(x_min_ind, x_max_ind-1),n_times)
            
            for k in slice_list:
                # Thin disc:               
                rr_1, I_image_1, I_mask_1, I_sigma_1, r_no_mask_1, I_no_mask_1, h1, h2, h3 = crea_vert_average_cut(image_data, mask_data, sigma_data, k, k+1, y_max_ind_dn, y_max_ind_dn+10, y_max_ind_dn+10, y_min_ind_up) 
                

                if h1_given is not None:
                    h1 = h1_given
                
                p1,p2,p3,p4 = fit_one_disc(R, r_no_mask_1, I_no_mask_1, I_image_1, I_sigma_1, I_mask_1, h1, n_fixed, psf_image, ron, gain, z01=z01_given, Y0=yc)
                L001.append(p1)
                Y001.append(p2)
                NN1.append(p3)
                Z001.append(p4)                
                

                # Thick disc
                rr, I_image, I_mask, I_sigma, r_no_mask, I_no_mask, h1, h2, h3 = crea_vert_average_cut(image_data, mask_data, sigma_data, k, k+1, y_min_ind_dn, y_max_ind_up, y_min_ind_up, y_max_ind_up)

                I_image_2 = I_image - disk_edge_soph(rr, 2.*h1*p1, p4, p2) # TODO: Y02 should be real y0!!!!!!!!!!!!!!!!!!!!!!    

                if h2_given is not None:
                    h2 = h2_given
                
                rr_2, I_image_22, I_mask_2, I_sigma_2, r_no_mask_2, I_no_mask_2, h1, h2, h3 = crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind, y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up) 
                
                p1,p2,p3,p4 = fit_one_disc(R, r_no_mask_2, I_no_mask_2, I_image_2, I_sigma_2, I_mask_2, h2, n_fixed, psf_image, ron, gain, z01=z02_given, Y0=yc)
                L002.append(p1)
                Y002.append(p2)
                NN2.append(p3)
                Z002.append(p4)                   

            L01_err = np.std(L001)
            Y01_err = np.std(Y001)
            N1_err = np.std(NN1)
            Z01_err = np.std(Z001)    

            L02_err = np.std(L002)
            Y02_err = np.std(Y002)
            N2_err = np.std(NN2)
            Z02_err = np.std(Z002)    
    



    return L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err



def do_each_bin_x(R, image_data, sigma_data, mask_data, psf_image, x_min_ind, x_max_ind, y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up, disc_comps='single', error_est=True, ron=None, gain=None, n_fixed='fixed', side=1, h1_given=None, h2_given=None, z01_given=None, z02_given=None, h_given=None, yc=None):

    print('SLICE X=%i:%i (R=%.1f pix)' % (x_min_ind, x_max_ind, R) )

    L0 = float('nan')
    L0_err = float('nan')
    Y0 = float('nan') 
    Y0_err = float('nan') 
    N = float('nan')
    N_err = float('nan')
    Z0 = float('nan')
    Z0_err = float('nan')

    L01 = float('nan')
    L01_err = float('nan')
    N1 = float('nan')
    N1_err = float('nan')
    Z01 = float('nan')
    Z01_err = float('nan')
    L02 = float('nan')
    L02_err = float('nan')
    N2 = float('nan')
    N2_err = float('nan')
    Z02 = float('nan')
    Z02_err = float('nan')
    L03 = float('nan')
    L03_err = float('nan')
    N3 = float('nan')
    N3_err = float('nan')
    Z03 = float('nan')
    Z03_err = float('nan')  


    # Make average slices of image and sigma image:
    if x_min_ind>x_max_ind:
      x_min_ind,x_max_ind = x_max_ind,x_min_ind  
    
    if disc_comps!='double': 
        L0,Y0,N,Z0,H,H2,H3 = fit_each_vert_slice(R, image_data, mask_data, sigma_data, psf_image, x_min_ind, x_max_ind, y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up, ron, gain, n_fixed, h=h_given, h1_given=h1_given, h2_given=h2_given, z01_given=z01_given, z02_given=z02_given)


        if error_est:
            L00 = []; Y00 = []; NN = []; Z00 = []
            if abs(x_max_ind - x_min_ind)<=n_times:
                slice_list = range(x_min_ind, x_max_ind-1)
            else:
                slice_list = random.sample(range(x_min_ind, x_max_ind-1),n_times)

            for k in slice_list:
                p1,p2,p3,p4,p5,p6,p7 = fit_each_vert_slice(R, image_data, mask_data, sigma_data, psf_image, k, k+1,  y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up, ron, gain, n_fixed, L0, N, Z0, Y0, H, h1_given=h1_given, h2_given=h2_given, z01_given=z01_given, z02_given=z02_given)
                L00.append(p1)
                Y00.append(p2)
                NN.append(p3)
                Z00.append(p4)
            L0_err = np.std(L00)
            Y0_err = np.std(Y00)
            N_err = np.std(NN)
            Z0_err = np.std(Z00)    
        else:
            L0_err = float('nan')
            Y0_err = float('nan')
            N_err = float('nan')
            Z0_err = float('nan')           


    if disc_comps=='double1':
        # Double disc fit:
        L01,N1,Z01,L02,N2,Z02,H,H2,H3 = fit_each_vert_slice(R, image_data, mask_data, sigma_data, psf_image, x_min_ind, x_max_ind,  y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up, ron, gain, n_fixed, L0, N, Z0, Y0, H, disc_comps='double1', h1_given=h2_given, h2_given=h2_given, z01_given=z01_given, z02_given=z02_given)

        if error_est:
            L001 = []; Y001 = []; NN1 = []; Z001 = []
            L002 = []; Y002 = []; NN2 = []; Z002 = []            
            if abs(x_max_ind - x_min_ind)<=n_times:
                slice_list = range(x_min_ind, x_max_ind-1)
            else:
                slice_list = random.sample(range(x_min_ind, x_max_ind-1),n_times)

            for k in slice_list:
                p1,p2,p3,p4,p5,p6,p7,p8,p9 = fit_each_vert_slice(R, image_data, mask_data, sigma_data, psf_image, k, k+1,  y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up, ron, gain, n_fixed, L01, N1, Z01, Y0, H, L02, N2, Z02, H2, disc_comps='double1', h1_given=h2_given, h2_given=h2_given, z01_given=z01_given, z02_given=z02_given)
                #exit()
                L001.append(p1)
                NN1.append(p2)
                Z001.append(p3)
                L002.append(p4)
                NN2.append(p5)
                Z002.append(p6)
            L01_err = np.std(L001)
            N1_err = np.std(NN1)
            Z01_err = np.std(Z001)    
            L02_err = np.std(L002)
            N2_err = np.std(NN2)
            Z02_err = np.std(Z002)  

    if disc_comps=='double': 
        L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err = fit_each_vert_slice_double(R, image_data, mask_data, sigma_data, psf_image, x_min_ind, x_max_ind,  y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up, ron, gain, n_fixed='fixed', error_est =error_est, h1_given=h1_given, h2_given=h2_given, z01_given=z01_given, z02_given=z02_given)


    if disc_comps=='doublee':         
        L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err = fit_each_vert_slice_double1(R, image_data, mask_data, sigma_data, psf_image, x_min_ind, x_max_ind,  y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up, ron, gain, n_fixed='fixed', error_est =error_est, h1_given=h1_given, h2_given=h2_given, z01_given=z01_given, z02_given=z02_given, yc=yc)        
        
 

    if disc_comps=='tripple':
        # Tripple disc fit:
        L01,N1,Z01,L02,N2,Z02,L03,N3,Z03 = fit_each_vert_slice(R, image_data, mask_data, sigma_data, psf_image, x_min_ind, x_max_ind,  y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up, ron, gain, n_fixed, L0, N, Z0, Y0, disc_comps='tripple', h1_given=h2_given, h2_given=h2_given, z01_given=z01_given, z02_given=z02_given)


        if error_est:
            L001 = []; Y001 = []; NN1 = []; Z001 = []
            L002 = []; Y002 = []; NN2 = []; Z002 = []  
            L003 = []; Y003 = []; NN3 = []; Z003 = []  
            if abs(x_max_ind - x_min_ind)<=n_times:
                slice_list = range(x_min_ind, x_max_ind-1)
            else:
                slice_list = random.sample(range(x_min_ind, x_max_ind-1),n_times)

            for k in slice_list:
                p1,p2,p3,p4,p5,p6,p7,p8,p9 = fit_each_vert_slice(R, image_data, mask_data, sigma_data, psf_image, k, k+1,  y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up, ron, gain, n_fixed, L01, N1, Z01, Y0, L02, N2, Z02, L03, N3, Z03, tripple_disc=tripple_disc, h1_given=h2_given, h2_given=h2_given, z01_given=z01_given, z02_given=z02_given)
                #exit()
                L001.append(p1)
                NN1.append(p2)
                Z001.append(p3)
                L002.append(p4)
                NN2.append(p5)
                Z002.append(p6)
                L003.append(p7)
                NN3.append(p8)
                Z003.append(p9)                
            L01_err = np.std(L001)
            N1_err = np.std(NN1)
            Z01_err = np.std(Z001)    
            L02_err = np.std(L002)
            N2_err = np.std(NN2)
            Z02_err = np.std(Z002)  
            L03_err = np.std(L003)
            N3_err = np.std(NN3)
            Z03_err = np.std(Z003)  

    

    for file in ['bestfit_parameters_imfit.dat',
                 'modelimage.fits',
                 'vertical_line.fits',
                 'bestfit_parameters_imfit_lm.dat',
                 'sigma_line.fits',
                 'galfit.01',
                 'galfit_lm.01',
                 'composed_model.fits']:
      if os.path.exists(file):
        os.remove(file)
    #exit()
    return [L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err],[L01,L01_err,N1,N1_err,Z01,Z01_err],[L02,L02_err,N2,N2_err,Z02,Z02_err],[L03,L03_err,N3,N3_err,Z03,Z03_err] 



























def define_xbins(rmin, rmax, Rmax, bin_type='log'):
    if rmax is None:
        rmax = Rmax

    if bin_type=='log':
        return Rmax*np.power(10, np.arange(math.log10(1./Rmax), math.log10(rmax/Rmax), 0.1))
    if bin_type=='linear':
        return Rmax*np.arange(rmin/Rmax, rmax/Rmax, 0.1)        
    if bin_type=='pix':
        return np.arange(rmin, rmax, 10)     








def ds9_to_np(x):
    '''
    Function to convert ds9 pixel coordinates to numpy format, i.e
    0.5 pix is 0, 1.5 pix is 1
    '''
    return int(round(x)) - 1


def create_add_mask(nx, ny, R_in_l, R_out_l, R_in_r, R_out_r, z_in_dn, z_out_dn, z_in_up, z_out_up):
    xc = nx/2.
    yc = ny/2.
    
    # define numpy indices:
    if R_out_l is None:
        x_min_ind_l = 0
    else:
        x_min_ind_l = ds9_to_np(xc - R_out_l)
        
    x_max_ind_l = ds9_to_np(xc - R_in_l)        


    if R_out_r is None:
        x_max_ind_r = nx
    else:
        x_max_ind_r = ds9_to_np(xc + R_out_r)
        
    x_min_ind_r = ds9_to_np(xc + R_in_r)   


    if z_out_dn is None:
        z_min_ind_dn = 0
    else:
        z_min_ind_dn = ds9_to_np(yc - z_out_dn)

    z_max_ind_dn = ds9_to_np(yc - z_in_dn)        

    #print(yc,z_in_dn,z_max_ind_dn)
    #exit()

    if z_out_up is None:
        z_max_ind_up = ny
    else:
        z_max_ind_up = ds9_to_np(yc + z_out_up)
        
    z_min_ind_up = ds9_to_np(yc + z_in_up)
    
    add_mask_data = np.ones((ny,nx))
    
    for k in range(ny):
        for i in range(nx):
            #if (k >= z_min_ind_dn and k<=z_max_ind_dn) or (k>=z_min_ind_up and k<=z_max_ind_up):
            if k >= z_min_ind_dn and k<=z_max_ind_up:
                if (i>=x_min_ind_l and i<=x_max_ind_l) or (i>=x_min_ind_r and i<=x_max_ind_r):
                    add_mask_data[k,i]=0
    
    hdu = pyfits.PrimaryHDU(add_mask_data)
    hdu.writeto('add_mask_tmp.fits', clobber=True)
    return x_min_ind_l,x_max_ind_l,x_min_ind_r,x_max_ind_r,z_min_ind_dn,z_max_ind_dn,z_min_ind_up,z_max_ind_up
    





def main(input_image, mask_image=None, weight_image=None, psf_image=None, R_in_l=0., R_out_l=None, R_in_r=0., R_out_r=None, z_in_dn = 0., z_out_dn=None, z_in_up=0., z_out_up=None, bin_type='linear', n_fixed='fixed', side='1', error_est=True, disc_comps='single', h1_given=None, h2_given=None, z01_given=None, z02_given=None, h_given=None, R_opt=None):
    # Output z0 is, in fact, hz, i.e. real_z0 = 2*z0
    
    # Sides:
    # 1 - right
    # 2 - left
    

    # The input images should be rotated, cropped and centered so that the major axis is horizontal
    hdu = pyfits.open(input_image)
    image_data = hdu[0].data 
    ny,nx =  np.shape(image_data)
    
    xc = nx/2.
    yc = ny/2.

    if R_opt is None:
        R_opt = nx/2.


    ##### MASK ########
    # Create additional mask:
    x_min_ind_l,x_max_ind_l,x_min_ind_r,x_max_ind_r,z_min_ind_dn,z_max_ind_dn,z_min_ind_up,z_max_ind_up = create_add_mask(nx, ny, R_in_l, R_out_l, R_in_r, R_out_r, z_in_dn, z_out_dn, z_in_up, z_out_up)

    #z_min_ind_dn = z_min_ind_dn #### WARNING
    #z_max_ind_dn = z_max_ind_up #### WARNING

    #print(z_min_ind_dn, z_max_ind_dn, z_min_ind_up, z_max_ind_up)
    #exit()

    if mask_image is None:
        shutil.copy('add_mask_tmp.fits', 'mask_tmp.fits')
    else:
        # Merge all masks
        #merge_masks.main([mask_image, 'add_mask_tmp.fits'], 'mask_tmp.fits') #### WARNING: UNCOMMENT THIS!!!!!!
        shutil.copy(mask_image, 'mask_tmp.fits')
    hdu = pyfits.open('mask_tmp.fits')
    mask_data = hdu[0].data
    mask_data = (mask_data>0)
    ##### MASK END ########        




    ##### PSF ########    
    # Create slice of PSF:
    if psf_image is not None:
            psf_slice(psf_image, axis='yaxis')
            psf_image_cut = 'psf_yline.fits'
     
    else:
        psf_image_cut = None
    ##### PSF END ########     


    psf_image_cut = 'psf.fits' ######################## WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! It's better than psf_yline.fits !!!!!!!!!!!!!!!!


    ##### WEIGHT IMAGE ######## 
    if weight_image is not None:
        hdu = pyfits.open(weight_image)
        sigma_data = hdu[0].data
    else: 
        sigma_data = np.ones((ny,nx))
    ##### WEIGHT IMAGE END ######## 

    
    ind_x_bins = []
    if side=='2':
        if R_out_l is None:
            R_out_l = nx/2.
            
        r_bins = define_xbins(R_in_l, R_out_l, R_opt, bin_type)
        for k in range(len(r_bins)):
            ind_x_bins.append( ds9_to_np(xc-r_bins[k]) )
    else:
        if R_out_r is None:
            R_out_r = nx/2.
            
        r_bins = define_xbins(R_in_r, R_out_r, R_opt, bin_type)
        for k in range(len(r_bins)):
            ind_x_bins.append( ds9_to_np(xc+r_bins[k]) )


    print('Fit code is: %s' % (code))
    ff = open('vertical_profiles_%s_%s.txt' % (side, disc_comps),'w')
    ff.close()


    print('Fitting averaged vertical cuts in the %s side:' % (side))    
    f = open('vertical_fits_%s_%s.dat' % (side, disc_comps), 'w')
    f.write('R\tL0\tL0_err\tY0\tY0_err\tN\tN_err\tZ0\tZ0_err\tL01\tL01_err\tN1\tN1_err\tZ01\tZ01_err\tL02\tL02_err\tN2\tN2_err\tZ02\tZ02_err\tL03\tL03_err\tN3\tN3_err\tZ03\tZ03_err\n')
    for k in range(0, len(ind_x_bins)-1):
            R = (r_bins[k]+r_bins[k+1])/2.
            [L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err],[L01,L01_err,N1,N1_err,Z01,Z01_err],[L02,L02_err,N2,N2_err,Z02,Z02_err],[L03,L03_err,N3,N3_err,Z03,Z03_err] = do_each_bin_x(R, image_data, sigma_data, mask_data, psf_image_cut, ind_x_bins[k], ind_x_bins[k+1], z_min_ind_dn, z_max_ind_dn, z_min_ind_up, z_max_ind_up, disc_comps=disc_comps, error_est=error_est, n_fixed=n_fixed, side=side, h1_given=h1_given, h2_given=h2_given, z01_given=z01_given, z02_given=z02_given, h_given=h_given, yc=yc)

            f.write('%f\t%.10f\t%.10f\t%f\t%f\t%f\t%f\t%f\t%f\t%.10f\t%.10f\t%f\t%f\t%f\t%f\t%.10f\t%.10f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n'   % (R,L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err,L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err,L03,L03_err,N3,N3_err,Z03,Z03_err ) )
    f.close()    

    # Remove all tmp files:
    for file in ['add_mask_tmp.fits',
                 'input_galf.txt',
                 'input_imf.txt',
                 'mask_line.fits',
                 'mask_tmp.fits',
                 'psf_xline.fits',
                 'psf_yline.fits',
                 'model.fits']:
        if os.path.exists(file):
            os.remove(file)

'''
# Test_eon/one_disc_galaxy
if __name__ == "__main__":    
    
    input_image = 'model_gal.fits'
    mask_image = None
    weight_image = None
    R_in_l = 10.
    R_out_l = None
    R_in_r = 10.
    R_out_r = None
    z_in_dn = 0.
    z_out_dn = None
    z_in_up = 0.
    z_out_up = None
    n_fixed = 'fixed'
  
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, bin_type='linear', n_fixed=n_fixed, side='1', error_est=False, disc_comps='single')
'''
'''
# Test_eon/two_disc_galaxy
if __name__ == "__main__":    
    
    input_image = 'model_gal.fits'
    mask_image = None
    weight_image = None
    R_in_l = 10.
    R_out_l = None
    R_in_r = 10.
    R_out_r = None
    z_in_dn = 0.
    z_out_dn = None
    z_in_up = 0.
    z_out_up = None
    n_fixed = 'fixed'
  
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, bin_type='linear', n_fixed=n_fixed, side='1', error_est=False, disc_comps='double1')
'''

'''
# Test_eon/one_disc_galaxy
if __name__ == "__main__":    
    
    input_image = 'model_gal.fits'
    mask_image = None
    weight_image = None
    R_in_l = 10.
    R_out_l = None
    R_in_r = 10.
    R_out_r = None
    z_in_dn = 0.
    z_out_dn = None
    z_in_up = 0.
    z_out_up = None
    n_fixed = 'fixed'
  
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, bin_type='linear', n_fixed=n_fixed, side='1', error_est=False, disc_comps='single')


    
    input_image = 'model_gal.fits'
    mask_image = None
    weight_image = None
    R_in_l = 10.
    R_out_l = None
    R_in_r = 10.
    R_out_r = None
    z_in_dn = 0.
    z_out_dn = None
    z_in_up = 0.
    z_out_up = None
    n_fixed = 'fixed'
  
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, bin_type='linear', n_fixed=n_fixed, side='1', error_est=False, disc_comps='double1')
'''    
'''   
# NGC4302_s4g_model
if __name__ == "__main__":       
    input_image = 'model_rot_crop.fits'
    mask_image = 'mask_rot_crop.fits'
    weight_image = 'sigma_rot_crop.fits'
    R_in_l = 10.
    R_out_l = None
    R_in_r = 10.
    R_out_r = None
    z_in_dn = 0.
    z_out_dn = None
    z_in_up = 0.
    z_out_up = None
    n_fixed = 'fixed'
  
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, bin_type='linear', n_fixed=n_fixed, side='1', error_est=False, disc_comps='single')
'''
'''
# NGC4302_s4g_model
if __name__ == "__main__":       
    input_image = 'model_rot_crop.fits'
    mask_image = 'mask_rot_crop.fits'
    weight_image = 'sigma_rot_crop.fits'
    R_in_l = 10.
    R_out_l = None
    R_in_r = 10.
    R_out_r = None
    z_in_dn = 0.
    z_out_dn = None
    z_in_up = 0.
    z_out_up = None
    n_fixed = 'fixed'
  
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, bin_type='linear', n_fixed=n_fixed, side='1', error_est=False, disc_comps='double1')
'''

'''
# NGC4302
if __name__ == "__main__":       
    input_image = 'NGC4302.phot.1_nonan_rot_crop_aver.fits'
    mask_image = 'NGC4302.1.finmask_nonan_rot_crop_aver.fits'
    weight_image = 'NGC4302_sigma2014_rot_crop_aver.fits'
    R_in_l = 10.
    R_out_l = None
    R_in_r = 10.
    R_out_r = None
    z_in_dn = 67.
    z_out_dn = None
    z_in_up = 67.
    z_out_up = None
    n_fixed = 'fixed'

    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, bin_type='linear', n_fixed=n_fixed, side='1', error_est=False, disc_comps='double')
'''