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
    
    if mu0d2 is None:
        if xc==0.5:
            xc_ch = 0
            yc_ch = 1
        if yc==0.5:
            xc_ch = 0 # FIXED!!!
            yc_ch = 0
    else:
            xc_ch = 0
            yc_ch = 0        
    
    if z01 is None:
        z01_ch = 0
        z02_ch = 0
        z01 = 5.
        z02 = 5.
        h1_ch = 1
        h2_ch = 1
    
    if h1 is None:
        h1_ch = 0
        h2_ch = 0
        h1 = 10.
        h2 = 10.
        z01_ch = 1 
        z02_ch = 1    
        
        
        
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




def read_galf(galf_file):
  f = open(galf_file,'r')
  L0 = 99999.
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
  if z01<z02:
      return L01,h1,z01/2.,L02,h2,z02/2.,chi2
  else:
      return L02,h2,z02/2.,L01,h1,z01/2.,chi2      




def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

def chi2_func(ref,mod):
    return sum( (ref-mod)**2)

def double_disc_func(z, I0d1, z01, I0d2, z02 ):
    #*** For edge-on disk SB in mag/arcsec^2 (along z-axis). n - index of this law. ***
    return I0d1 * (1.0 / np.cosh(1.*np.fabs(z)/(2.*z01)) ) **(2./1.) + I0d2 * (1.0 / np.cosh(1.*np.fabs(z)/(2.*z02)) ) **(2./1.)

def find_double_disc_guess(zz, II, y0):#, z0, I0):
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

      #print(I0_thin,I0_thick,z0_thin,z0_thick,z0_thin_err,z0_thick_err)
      
      #plt.plot(z, piecewise_linear(z, results_pos[0], results_pos[1], results_pos[2], results_pos[3]))
      

      p,e = curve_fit(double_disc_func, ZZ, II, [I0_thin, z0_thin, I0_thick, z0_thick] )

      I0_thin = p[0]
      z0_thin = np.fabs(p[1])

      I0_thick = p[2]
      z0_thick = np.fabs(p[3])
      #print(I0_thin,I0_thick,z0_thin,z0_thick)
      #plt.plot(ZZ, II)
      #plt.plot(ZZ, double_disc_func(ZZ, I0_thin, z0_thin, I0_thick, z0_thick))
      #plt.show()
      #exit()
      if z0_thin<z0_thick:
            return I0_thin/(2.*10.),I0_thick/(2.*10.),z0_thin,z0_thick,
      else:
            return I0_thick/(2.*10.),I0_thin/(2.*10.),z0_thick,z0_thin         


def crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind, y_min_ind, y_max_ind):
    I_image = []; I_mask = []; I_sigma = []
    ny,nx = np.shape(image_data)
    
    for i in range(y_min_ind, y_max_ind):  
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
        except:
            I_image.append(0.)
            I_mask.append(True)
            I_sigma.append(0.)
    return np.array(I_image), np.array(I_mask), np.array(I_sigma)
        



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
""" % (0.5,y0,90.,L01,10.,n1,n_fixed,z01,90.,L02,10.,n2,n_fixed,z02)
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


def run_galf(L0, n, z0, y0, sigma_image, psf_image, mask_image):
  try:
    crea_galfit_input('vertical_line.fits', sigma_image, psf_image, mask_image, 0.5, y0, -2.5*math.log10(L0*(2.*10.)), z01=2.*z0, h1=None,  mu0d2=None, z02=None, h2=None)
    
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

def run_galf_double(L01, n1, z01, y0, L02, n2, z02, sigma_image, psf_image, mask_image, ron, gain):
  try:
    crea_galfit_input('vertical_line.fits', sigma_image, psf_image, mask_image, 0.5, y0, -2.5*math.log10(L01*(2.*10.)), z01=2.*z01, h1=None,  mu0d2=-2.5*math.log10(L02*(2.*10.)), z02=2.*z02, h2=None)

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
    


def single_disc_fit(L0, n, z0, y0, n_fixed, sigma_image, psf_image, mask_image, ron, gain):
    if code=='imfit':
        # Levenberg-Marquardt:
        status = run_imf(L0, n, z0, y0, n_fixed, sigma_image, psf_image, mask_image, ron, gain)

        if status == 0:
            Y0,L00,N,Z0,chi2 = read_imf('bestfit_parameters_imfit_lm.dat')
        else:
            L00=float('nan'); Y0=float('nan'); N=float('nan'); Z0=float('nan')
    elif code=='galfit':
        # Levenberg-Marquardt:
        status = run_galf(L0, n, z0, y0, sigma_image, psf_image, mask_image)

        if status == 0:
            X0,Y0,L00,Z0,H,chi2 = read_galf('galfit_lm.01')
            N = 1.
        else:
            L00=float('nan'); Y0=float('nan'); N=float('nan'); Z0=float('nan')        
        
    return L00,Y0,N,Z0



def double_disc_fit(L01, n1, z01, y0, L02, n2, z02, n_fixed, sigma_image, psf_image, mask_image, ron, gain):
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
        status = run_galf_double(L01, n1, z01, y0, L02, n2, z02, sigma_image, psf_image, mask_image, ron, gain)
        if status == 0:
            L01,h1,Z01,L02,h2,Z02,chi2 = read_galf_double('galfit_lm.01')
            N1 = 1
            N2 = 1
        else:
            L01=float('nan'); N1=float('nan'); Z01=float('nan'); L02=float('nan'); N2=float('nan'); Z02=float('nan'); chi2=float('nan')    
    
    return L01,N1,Z01,L02,N2,Z02
    


def fit_each_vert_slice(image_data, mask_data, sigma_data, psf_image, x_min_ind, x_max_ind, y_min_ind, y_max_ind, ron, gain, n_fixed='fixed', L0=None, N=None, z0=None, y0=None, L02=None, N2=None, z02=None,  double_disc=False):
    I_image, I_mask, I_sigma = crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind, y_min_ind, y_max_ind)     

    I0 = max(I_image)
    #y0 = float(np.where(I_image==I0)[0])+y_min_ind
    rr = np.arange(1, len(I_image)+1, 1)+y_min_ind
    #plt.plot(rr, I_image)
    #plt.show()
    #exit()
    
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

    if L0 is None and z0 is None and N is None:
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

    #exit()
    # Single disc fit:
    if not double_disc:
        L0,Y0,N,Z0 = single_disc_fit(L0, N, z0, y0, n_fixed, sigma_image, psf_image, mask_image, ron, gain)
        return L0,Y0,N,Z0,I_image,rr
    else:
        '''
        if L02 is None and N2 is None and z02 is None: # WARNING: !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            L02 = L0/30.
            N2 = N
            z02 = 3.*z0
        '''
        if L02 is None:
            L01,L02,z01,z02 = find_double_disc_guess(rr, I_image, y0)  #### WARNING: DOES MOT WORK PROPERLY!!!
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
        
        L01,N1,Z01,L02,N2,Z02 = double_disc_fit(L01, N, z01, y0, L02, N, z02, n_fixed, sigma_image, psf_image, mask_image, ron, gain)
        #exit()
        return L01,N1,Z01,L02,N2,Z02






def do_each_bin_x(R, image_data, sigma_data, mask_data, psf_image, x_min_ind, x_max_ind, y_min_ind, y_max_ind, double_disc=True, error_est=True, ron=None, gain=None, n_fixed='fixed', side=1):
    print('SLICE X=%i:%i (R=%.1f pix)' % (x_min_ind, x_max_ind, R) )
    # Make average slices of image and sigma image:
    if x_min_ind>x_max_ind:
      x_min_ind,x_max_ind = x_max_ind,x_min_ind  
    
    
    L0,Y0,N,Z0,I_image,rr = fit_each_vert_slice(image_data, mask_data, sigma_data, psf_image, x_min_ind, x_max_ind, y_min_ind, y_max_ind, ron, gain, n_fixed)

    # Write profile to file:
    ff = open('vertical_profiles_%s.txt' % (side),'a')
    ff.write('#RADIUS:%f\n' % (R))

    for kk in range(len(rr)):
        ff.write('%i\t%.10f\n' % (rr[kk],I_image[kk]))
    ff.write('#END\n')
    ff.close()



    if error_est:
        L00 = []; Y00 = []; NN = []; Z00 = []
        if abs(x_max_ind - x_min_ind)<=n_times:
            slice_list = range(x_min_ind, x_max_ind-1)
        else:
            slice_list = random.sample(range(x_min_ind, x_max_ind-1),n_times)

        for k in slice_list:
            p1,p2,p3,p4,p5,p6 = fit_each_vert_slice(image_data, mask_data, sigma_data, psf_image, k, k+1, y_min_ind, y_max_ind, ron, gain, n_fixed, L0, N, Z0, Y0)
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

    #print(L0,Z0,Y0)
    #exit()
    if double_disc:
        # Double disc fit:
        L01,N1,Z01,L02,N2,Z02 = fit_each_vert_slice(image_data, mask_data, sigma_data, psf_image, x_min_ind, x_max_ind, y_min_ind, y_max_ind, ron, gain, n_fixed, L0, N, Z0, Y0, double_disc=double_disc)

        #exit()
        if error_est:
            L001 = []; Y001 = []; NN1 = []; Z001 = []
            L002 = []; Y002 = []; NN2 = []; Z002 = []            
            if abs(x_max_ind - x_min_ind)<=n_times:
                slice_list = range(x_min_ind, x_max_ind-1)
            else:
                slice_list = random.sample(range(x_min_ind, x_max_ind-1),n_times)

            for k in slice_list:
                p1,p2,p3,p4,p5,p6 = fit_each_vert_slice(image_data, mask_data, sigma_data, psf_image, k, k+1, y_min_ind, y_max_ind, ron, gain, n_fixed, L01, N1, Z01, Y0, L02, N2, Z02, double_disc=double_disc)
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
        else: 
            L01_err = float('nan') 
            N1_err = float('nan') 
            Z01_err = float('nan')    
            L02_err = float('nan') 
            N2_err = float('nan') 
            Z02_err = float('nan')   






    else:
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
      

    for file in ['bestfit_parameters_imfit.dat',
                 'modelimage.fits',
                 'vertical_line.fits',
                 'bestfit_parameters_imfit_lm.dat',
                 'sigma_line.fits',
                 'galfit.01',
                 'galfit_lm.01']:
      if os.path.exists(file):
        os.remove(file)
    #exit()
    return [L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err],[L01,L01_err,N1,N1_err,Z01,Z01_err],[L02,L02_err,N2,N2_err,Z02,Z02_err] 



























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


    if z_out_up is None:
        z_max_ind_up = ny
    else:
        z_max_ind_up = ds9_to_np(yc + z_out_up)
        
    z_min_ind_up = ds9_to_np(yc + z_in_up)
    
    add_mask_data = np.ones((ny,nx))
    
    for k in range(ny):
        for i in range(nx):
            if (k >= z_min_ind_dn and k<=z_max_ind_dn) or (k>=z_min_ind_up and k<=z_max_ind_up):
                if (i>=x_min_ind_l and i<=x_max_ind_l) or (i>=x_min_ind_r and i<=x_max_ind_r):
                    add_mask_data[k,i]=0
    
    hdu = pyfits.PrimaryHDU(add_mask_data)
    hdu.writeto('add_mask_tmp.fits', clobber=True)
    return x_min_ind_l,x_max_ind_l,x_min_ind_r,x_max_ind_r,z_min_ind_dn,z_max_ind_dn,z_min_ind_up,z_max_ind_up
    





def main(input_image, mask_image=None, weight_image=None, psf_image=None, R_in_l=0., R_out_l=None, R_in_r=0., R_out_r=None, z_in_dn = 0., z_out_dn=None, z_in_up=0., z_out_up=None, bin_type='linear', n_fixed='fixed', side='1', error_est=True, double_disc=True):
    # Sides:
    # 1 - right
    # 2 - left

    

    # The input images should be rotated, cropped and centered so that the major axis is horizontal
    hdu = pyfits.open(input_image)
    image_data = hdu[0].data 
    ny,nx =  np.shape(image_data)
    xc = nx/2.
    yc = ny/2.



    ##### MASK ########
    # Create additional mask:
    x_min_ind_l,x_max_ind_l,x_min_ind_r,x_max_ind_r,z_min_ind_dn,z_max_ind_dn,z_min_ind_up,z_max_ind_up = create_add_mask(nx, ny, R_in_l, R_out_l, R_in_r, R_out_r, z_in_dn, z_out_dn, z_in_up, z_out_up)

    z_min_ind_dn = z_min_ind_dn #### WARNING
    z_max_ind_dn = z_max_ind_up #### WARNING

    #print(z_min_ind_dn, z_max_ind_dn)
    #exit()

    if mask_image is None:
        shutil.copy('add_mask_tmp.fits', 'mask_tmp.fits')
    else:
        # Merge all masks
        merge_masks.main([mask_image, 'add_mask_tmp.fits'], 'mask_tmp.fits') #### WARNING: UNCOMMENT THIS!!!!!!
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





    ##### WEIGHT IMAGE ######## 
    if weight_image is not None:
        hdu = pyfits.open(weight_image)
        sigma_data = hdu[0].data
    else: 
        sigma_data = np.ones((ny,nx))
    ##### WEIGHT IMAGE END ######## 

    
    ind_x_bins = []
    if side=='2':
        r_bins = define_xbins(R_in_l, R_out_l, nx/2., bin_type)
        for k in range(len(r_bins)):
            ind_x_bins.append( ds9_to_np(xc-r_bins[k]) )
    else:
        r_bins = define_xbins(R_in_r, R_out_r, nx/2., bin_type)
        for k in range(len(r_bins)):
            ind_x_bins.append( ds9_to_np(xc+r_bins[k]) )

    #print(r_bins, ind_x_bins)
    #exit()
    print('Fit code is: %s' % (code))
    ff = open('vertical_profiles_%s.txt' % (side),'w')
    ff.close()


    print('Fitting averaged vertical cuts in the %s side:' % (side))    
    f = open('vertical_fits_%s.dat' % (side), 'w')
    f.write('R\tL0\tL0_err\tY0\tY0_err\tN\tN_err\tZ0\tZ0_err\tL01\tL01_err\tN1\tN1_err\tZ01\tZ01_err\tL02\tL02_err\tN2\tN2_err\tZ02\tZ02_err\n')
    for k in range(0, len(ind_x_bins)-1):
            R = (r_bins[k]+r_bins[k+1])/2.
            [L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err],[L01,L01_err,N1,N1_err,Z01,Z01_err],[L02,L02_err,N2,N2_err,Z02,Z02_err] = do_each_bin_x(R, image_data, sigma_data, mask_data, psf_image_cut, ind_x_bins[k], ind_x_bins[k+1], z_min_ind_dn, z_max_ind_dn, double_disc=double_disc, error_est=error_est, n_fixed=n_fixed, side=side)

            f.write('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n'   % (R,L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err,L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err ) )
    f.close()    

    # Remove all tmp files:
    for file in ['add_mask_tmp.fits',
                 'input_galf.txt',
                 'input_imf.txt',
                 'mask_line.fits',
                 'mask_tmp.fits',
                 'psf_xline.fits',
                 'psf_yline.fits']:
        if os.path.exists(file):
            os.remove(file)




  
if __name__ == "__main__":    
    input_image = 'NGC3628.phot.1_nonan_rot_crop_aver.fits'
    mask_image = 'new_mask_rot_crop_aver.fits'
    weight_image = 'NGC3628_sigma2014_rot_crop_aver.fits'
    R_in_l = 107.
    R_out_l = None
    R_in_r = 107.
    R_out_r = None
    z_in_dn = 0.
    z_out_dn=None
    z_in_up=0.
    z_out_up = None
    n_fixed = 'fixed'
    #main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, bin_type='linear', n_fixed=n_fixed, side='1', error_est=True, double_disc=True)   
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, bin_type='linear', n_fixed=n_fixed, side='1', error_est=False, double_disc=True)    

'''
if __name__ == "__main__":    
    input_image = 'NGC4302.phot.1_nonan_rot_crop_aver.fits'
    mask_image = 'NGC4302.1.finmask_nonan_rot_crop_aver.fits'
    weight_image = 'NGC4302_sigma2014_rot_crop_aver.fits'
    R_in_l = 0.
    R_out_l = None
    R_in_r = 0.
    R_out_r = None
    z_in_dn = 0.
    z_out_dn=None
    z_in_up=0.
    z_out_up = None
    n_fixed = 'fixed'
  
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, bin_type='linear', n_fixed=n_fixed, side='1', error_est=False, double_disc=True)   
'''
'''
if __name__ == "__main__":    
    input_image = 'NGC891_coadd_rot_crop_aver.fits'
    mask_image = 'NGC891_coadd_mask_rot_crop_aver.fits'
    weight_image = 'NGC891_coadd_sigma_rot_crop_aver.fits'
    R_in_l = 130.
    R_out_l = None
    R_in_r = 130.
    R_out_r = None
    z_in_dn = 0.
    z_out_dn=None
    z_in_up=0.
    z_out_up = None
    n_fixed = 'fixed'
  
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, bin_type='linear', n_fixed=n_fixed, side='1', error_est=False, double_disc=True)    
'''













