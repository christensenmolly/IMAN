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
import collections

warnings.filterwarnings("ignore")

FNULL = open(os.devnull, 'w')

import fit_eon_vertical

sys.path.append( '/home/amosenko/MEGA/MyPrograms/IMAN/imp/masking')
sys.path.append( '/home/amosenko/MEGA/MyPrograms/IMAN/decomposition/make_model')
import make_model_ima_galfit
import merge_masks


path = ''#'/Users/mosenkov/MEGA/HERON/ESO240-G011/'
imfit_path = path#'/home/amosenko/MEGA/MyPrograms/imfit-1.6.1/'
galfit_path = path
fitscopy = path


####
n_times = 5 # Max number of randomly selected cuts from a given interval
code = 'galfit' # galfit check with and without psf
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
        n_psf=1
    else:
        #hdu = pyfits.open(psf_image)
        #image_psf = hdu[0].data 
        #ny_psf,nx_psf =  np.shape(image_psf)
        if xc==0.5:
            n_psf = ny
        if yc==0.5:
            n_psf = nx
    
    if mu0d2 is None:
        if xc==0.5:
            xc_ch = 0
            yc_ch = 1 #### TODO 1 
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
""" % (input_image, sigma_image, psf_image, mask_image, nx, ny, 3*ny, 3*ny)

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


def crea_hor_average_cut(image_data, mask_data, sigma_data, z_min_ind, z_max_ind, x_min_ind, x_max_ind):
    I_image = []; I_mask = []; I_sigma = []
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


def disk_exp(r, I0d, h, x_c):
    return I0d * np.exp(-np.fabs(r-x_c)/h)




def crea_imf(L0, h, x0, ron=None, gain=None):
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
h		%f
n		%f	fixed
z_0		%f	fixed    
""" % (x0,0.5,90.,L0,h,1.,5.)  
    print >>f, s
    f.close()
    #exit()



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
    if 'h\t' in line and L0!=99999.:
      h = float(line.split()[1])
    if 'z_0' in line:
      z0 = float(line.split()[1])
  f.close()
  return L0,h,chi2

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
  return L0,h,chi2


def run_imf(L0, h, x0, sigma_image, psf_image, mask_image, ron, gain):
    crea_imf(L0, h, x0, ron, gain)

    imfit_line = '%simfit radial_line.fits -c input_imf.txt --ftol 0.00001 --max-threads 1' % (imfit_path)
    if sigma_image!=None:
        imfit_line += ' --noise sigma_line.fits '

    if psf_image!=None:
        imfit_line += ' --psf %s ' % (psf_image)        

    if mask_image!=None:
        imfit_line += ' --mask mask_line.fits '   


    subprocess.call(imfit_line, shell=True, stdout=FNULL)

    try:
        shutil.copy('bestfit_parameters_imfit.dat','bestfit_parameters_imfit_lm.dat')
        return 0
    except:
        return 1


def run_galf(L0, h, x0, sigma_image, psf_image, mask_image):
  try:
    crea_galfit_input('radial_line.fits', sigma_image, psf_image, mask_image, x0, 0.5, -2.5*math.log10(abs(L0*2.*h)), z01=None, h1=h,  mu0d2=None, z02=None, h2=None) #### TODO: ADD module!
    
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

def single_disc_fit(L0, h, x0, input_image, sigma_image, psf_image, mask_image, ron, gain):
    if code=='imfit':
        # Levenberg-Marquardt:
        status = run_imf(L0, h, x0, sigma_image, psf_image, mask_image, ron, gain)
        if status==0:
            L00,H,chi2 = read_imf('bestfit_parameters_imfit_lm.dat')
        else:
            L00 = float('nan'); H = float('nan')
    if code=='galfit':
        status = run_galf(L0, h, x0, sigma_image, psf_image, mask_image)
        if status==0:
            L00,H,chi2 = read_galf('galfit_lm.01')
            make_model_ima_galfit.main(input_image, 'galfit_lm.01', composed_model_file = 'composed_model.fits', subtract_sky=True, galfitPath=galfit_path)
        else:
            L00 = float('nan'); H = float('nan')    
    return L00,H


def fit_each_hor_slice(Z, image_data, mask_data, sigma_data, psf_image, z_min_ind, z_max_ind, x_min_ind, x_max_ind, xc, ron, gain, L0=None, h=None, z0=None, y0=None, L02=None, N2=None, z02=None,  double_disc=False, side='1'):
    I_image, I_mask, I_sigma = crea_hor_average_cut(image_data, mask_data, sigma_data, z_min_ind, z_max_ind, x_min_ind, x_max_ind)     
    #print(z_min_ind, z_max_ind, x_min_ind, x_max_ind)
    #exit()
    I0 = max(I_image)
    #exit()
    if side=='2' or side=='3':
        rr = np.arange(1, len(I_image)+1, 1) + x_min_ind
    else:
        rr = np.arange(1, len(I_image)+1, 1) + x_min_ind 
    

    #### Create slice images: ###
    # Data image:
    fict_data = np.random.random((2,len(I_image)))
    fict_data[0,:] = I_image
    fict_data[1,:] = np.zeros(len(I_image))
    hdu = pyfits.PrimaryHDU(fict_data)
    hdu.writeto('slice.fits',clobber=True)
    #exit()
    if os.path.exists('radial_line.fits'):
        os.remove('radial_line.fits')
    subprocess.call('%sfitscopy \'%s[%i:%i,%i:%i]\' radial_line.fits' % (fitscopy,'slice.fits',1,len(I_image),1,1), shell=True)
    os.remove('slice.fits')
    input_image = 'radial_line.fits'
    #exit()
    
    # Sigma image:
    fict_data = np.random.random((2,len(I_image)))
    fict_data[0,:] = I_sigma
    fict_data[1,:] = np.zeros(len(I_sigma))
    hdu = pyfits.PrimaryHDU(fict_data)
    hdu.writeto('slice_sigma.fits',clobber=True)      

    if os.path.exists('sigma_line.fits'):
        os.remove('sigma_line.fits')
    subprocess.call('%sfitscopy \'%s[%i:%i,%i:%i]\' sigma_line.fits' % (fitscopy,'slice_sigma.fits',1,len(I_image),1,1), shell=True)
    os.remove('slice_sigma.fits')
    sigma_image = 'sigma_line.fits'
    

    # Mask image:
    fict_data = np.random.random((2,len(I_image)))
    fict_data[0,:] = I_mask
    fict_data[1,:] = np.zeros(len(I_mask))
    hdu = pyfits.PrimaryHDU(fict_data)
    hdu.writeto('slice_mask.fits',clobber=True)      

    if os.path.exists('mask_line.fits'):
        os.remove('mask_line.fits')
    subprocess.call('%sfitscopy \'%s[%i:%i,%i:%i]\' mask_line.fits' % (fitscopy,'slice_mask.fits',1,len(I_image),1,1), shell=True)
    os.remove('slice_mask.fits')
    mask_image = 'mask_line.fits'
    #### Create slice images END ###

    #print(I0)
    #exit()
    #plt.plot(rr, I_image)
    #plt.show()
    #exit()
    #print('here',I_image, h)
    if L0 is None and h is None:
      try:
        # Find first guess on the parameters:
        if side=='2' or side=='3':
            popt_soph, pcov_soph = curve_fit(lambda r, I0, h: disk_exp(r, I0, h, xc), rr, I_image, p0=(I0, len(I_image)/4.) )
        else:
            popt_soph, pcov_soph = curve_fit(lambda r, I0, h: disk_exp(r, I0, h, xc), rr, I_image, p0=(I0, len(I_image)/4.) )
        h = popt_soph[1]
        L0 = popt_soph[0]/(2.*h)
      except:
        h = len(I_image)/4.
        L0 = I0/(2.*h)
    #print(L0,h)
    '''
    print(h,L0)

    plt.plot(rr, I_image, 'o')
    if side=='2' or side=='3':
        plt.plot(rr, disk_exp(rr, L0, h, len(I_image)), '-')
    else:
        plt.plot(rr, disk_exp(rr, L0, h, 0.), '-')
    plt.show()


    exit()
    '''
    if side=='2' or side=='3':
        L0,h = single_disc_fit(L0, h, xc-x_min_ind, input_image, sigma_image, psf_image, mask_image, ron, gain)
    else:
        L0,h = single_disc_fit(L0, h, xc-x_min_ind, input_image, sigma_image, psf_image, mask_image, ron, gain)
    write_profile(Z, rr, I_image, 'single', components=[1])
    return L0,h,I_image,rr

def write_profile(Z, rr, I_image, disc_comps, components=[1]):
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        hdu_model = pyfits.open('composed_model.fits')
        
        model_data = []
        for comp in components:
            #print(hdu_model[comp].data[0,:])
            #exit()
            model_data.append(hdu_model[comp].data[0,:]) 
        
        Model_data = []
        for kk in range(len(rr)):
            for ii in range(len(model_data)):
                if ii==0:
                    #print(model_data[ii][kk])
                    ss = '%.10f' % (model_data[ii][kk])
                else:
                    ss = ss + '\t%.10f' % (model_data[ii][kk])
            Model_data.append(ss)
        
        
        
        # Write profile to file:
        ff = open('radial_profiles_%s_%s.txt' % ('1', disc_comps),'a') # TODO
        ff.write('#Z:%f\n' % (Z))

        for kk in range(len(rr)):
            ff.write('%i\t%.10f\t%s\n' % (rr[kk],I_image[kk],Model_data[kk]))
        ff.write('#END\n')
        ff.close()
        # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

def do_each_bin_z(Z, image_data, sigma_data, mask_data, psf_image, z_min_ind, z_max_ind, x_min_ind, x_max_ind, xc, break_radii=None, error_est=True, ron=None, gain=None, side='1'):
    print('SLICE Y=%i:%i (Z=%.1f pix)' % (z_min_ind, z_max_ind, Z) )
    # Make average slices of image and sigma image:
    if z_min_ind>z_max_ind:
      z_min_ind,z_max_ind = z_max_ind,z_min_ind  
    

    L0,h,I_image,zz = fit_each_hor_slice(Z, image_data, mask_data, sigma_data, psf_image, z_min_ind, z_max_ind, x_min_ind, x_max_ind, xc, ron, gain, side=side)



    if error_est:
        L00 = []; H00 = []
        if abs(z_max_ind - z_min_ind)<=n_times:
            slice_list = range(z_min_ind, z_max_ind-1)
        else:
            slice_list = random.sample(range(z_min_ind, z_max_ind-1),n_times)

        for k in slice_list:
            p1,p2,p3,p4 = fit_each_hor_slice(Z, image_data, mask_data, sigma_data, psf_image, k, k+1, x_min_ind, x_max_ind, xc, ron, gain, L0, h, side=side)
            L00.append(p1)
            H00.append(p2)

        L0_err = np.std(L00)
        h_err = np.std(H00)
  
    else:
        L0_err = float('nan')
        h_err = float('nan')

    Breaks_res = collections.OrderedDict()
    if break_radii is not None:
        #print('Breaks', x_min_ind, x_max_ind)
        if side=='2' or side=='3':
            x_nodes = [x_min_ind] + list(np.array(xc-break_radii, int)) + [x_max_ind]
        else:
            break_radii = np.array(sorted(break_radii))
            x_nodes = [x_min_ind] + list(np.array(xc+break_radii, int)) + [x_max_ind]
    
        for k in range(len(x_nodes)-1):
            #print(x_nodes[k], x_nodes[k+1])
            L000,h00,I_image,zz = fit_each_hor_slice(Z, image_data, mask_data, sigma_data, psf_image, z_min_ind, z_max_ind, x_nodes[k], x_nodes[k+1], xc, ron, gain, side=side)
            #print('\t',L000,h00)
            Breaks_res['L0%i' % (k+1)] = L000
            Breaks_res['h%i' % (k+1)] = h00
            if k<len(break_radii):
                Breaks_res['Rbr%i' % (k+1)] = break_radii[k]
            
            if error_est:
                L00 = []; H00 = []
                if abs(z_max_ind - z_min_ind)<=n_times:
                    slice_list = range(z_min_ind, z_max_ind-1)
                else:
                    slice_list = random.sample(range(z_min_ind, z_max_ind-1),n_times)

                for kk in slice_list:
                    p1,p2,p3,p4 = fit_each_hor_slice(Z, image_data, mask_data, sigma_data, psf_image, kk, kk+1, x_nodes[k], x_nodes[k+1], xc, ron, gain, L0, h, side=side)
                    L00.append(p1)
                    H00.append(p2)

                Breaks_res['L0%i_err' % (k+1)] = np.std(L00)
                Breaks_res['h%i_err' % (k+1)] = np.std(H00)
        
            else:
                Breaks_res['L0%i_err' % (k+1)] = float('nan')
                Breaks_res['h%i_err' % (k+1)] = float('nan')


    else:
        for k in range(3):
            Breaks_res['L0%i' % (k+1)] = float('nan')
            Breaks_res['h%i' % (k+1)] = float('nan')
            Breaks_res['L0%i_err' % (k+1)] = float('nan')
            Breaks_res['h%i_err' % (k+1)] = float('nan')
            Breaks_res['Rbr%i' % (k+1)] = float('nan')
    
     

    for file in ['bestfit_parameters_imfit.dat',
                 'modelimage.fits',
                 'radial_line.fits',
                 'bestfit_parameters_imfit_lm.dat',
                 'sigma_line.fits',
                 'galfit.01',
                 'galfit_lm.01',
                 'fit.log']:
        if os.path.exists(file):
            os.remove(file)
    #exit()
    return [L0,L0_err,h,h_err],[Breaks_res['L01'],Breaks_res['L01_err'],Breaks_res['h1'],Breaks_res['h1_err'],Breaks_res['Rbr1']],[Breaks_res['L02'],Breaks_res['L02_err'],Breaks_res['h2'],Breaks_res['h2_err'],Breaks_res['Rbr2']],[Breaks_res['L03'],Breaks_res['L03_err'],Breaks_res['h3'],Breaks_res['h3_err']] 




def define_zbins(z_min, z_max, Zmax, bin_type='log'):
    if z_max is None:
        z_max = Zmax

    if bin_type=='log':
        return Zmax*np.power(10, np.arange(math.log10(1./Zmax), math.log10(z_max/Zmax), 0.1))
    if bin_type=='linear':
        return Zmax*np.arange(z_min/Zmax, z_max/Zmax, 0.1)        
    if bin_type=='pix':
        return np.arange(z_min, z_max, 10) 


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
    





def main(input_image, mask_image=None, weight_image=None, psf_image=None, R_in_l=0., R_out_l=None, R_in_r=0., R_out_r=None, z_in_dn = 0., z_out_dn=None, z_in_up=0., z_out_up=None, break_radii=None, bin_type='linear', side='1', error_est=False):
    # Sides:
    # 1 - up, right
    # 2 - up, left
    # 3 - bottom, left
    # 4 - bottom, right
    
    
    if break_radii is not None:
        break_radii = np.array(sorted(break_radii, reverse=True))


    # The input images should be rotated, cropped and centered so that the major axis is horizontal.
    hdu = pyfits.open(input_image)
    image_data = hdu[0].data 
    ny,nx =  np.shape(image_data)
    xc = nx/2.
    yc = ny/2.



    ##### MASK ########
    # Create additional mask:
    x_min_ind_l,x_max_ind_l,x_min_ind_r,x_max_ind_r,z_min_ind_dn,z_max_ind_dn,z_min_ind_up,z_max_ind_up = create_add_mask(nx, ny, R_in_l, R_out_l, R_in_r, R_out_r, z_in_dn, z_out_dn, z_in_up, z_out_up)
    #print(x_min_ind_l,x_max_ind_l)
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
            psf_slice(psf_image, axis='xaxis')
            psf_image_cut = 'psf_xline.fits'            
    else:
        psf_image_cut = None
    ##### PSF END ########     
    if code=='galfit':
        psf_image_cut = psf_image # For the radial profile GALFIT crashes with psf_xline.fits. Why???




    ##### WEIGHT IMAGE ######## 
    if weight_image is not None:
        hdu = pyfits.open(weight_image)
        sigma_data = hdu[0].data
    else: 
        sigma_data = np.ones((ny,nx))
    ##### WEIGHT IMAGE END ######## 

    if side=='1' or side=='2':
        z_bins = define_zbins(z_in_up, z_out_up, ny/2., bin_type=bin_type)
        ind_z_bins = []
        for k in range(len(z_bins)):
            ind_z_bins.append( ds9_to_np(yc+z_bins[k]) )   
    if side=='3' or side=='4':
        z_bins = define_zbins(z_in_dn, z_out_dn, ny/2., bin_type=bin_type)
        ind_z_bins = []
        for k in range(len(z_bins)):
            ind_z_bins.append( ds9_to_np(yc-z_bins[k]) )
        ind_z_bins = sorted(ind_z_bins)
        z_bins = sorted(z_bins, reverse=True)

    if side=='1' or side=='4':
        x_min_ind = x_min_ind_r
        x_max_ind = x_max_ind_r
    else:
        x_min_ind = x_min_ind_l
        x_max_ind = x_max_ind_l    
    
    #print(z_bins)
    #exit()
    

    print('Fit code is: %s' % (code))    
    ff = open('radial_profiles_%s_single.txt' % (side),'w')
    ff.close()
  
    
    print('Fitting averaged radial cuts in the %s side:' % (side))
    #print(z_bins, ind_z_bins)
    #exit()
    f = open('radial_fits_%s_single.dat' % (side), 'w')
    f.write('Z\tL0\tL0_err\th\th_err\tL01\tL01_err\th1\th1_err\tRbr1\tL02\tL02_err\th2\th2_err\tRbr2\tL03\tL03_err\th3\th3_err\n')
    for k in range(0, len(ind_z_bins)-1):
        Z = (z_bins[k]+z_bins[k+1])/2.
        [L0,L0_err,h,h_err],[L01,L01_err,h1,h1_err,Rbr1],[L02,L02_err,h2,h2_err,Rbr2],[L03,L03_err,h3,h3_err] = do_each_bin_z(Z, image_data, sigma_data, mask_data, psf_image_cut, ind_z_bins[k], ind_z_bins[k+1], x_min_ind, x_max_ind, ds9_to_np(xc), break_radii=break_radii, error_est=error_est, side=side)

        f.write('%f\t%.20f\t%.20f\t%f\t%f\t%.20f\t%.20f\t%f\t%f\t%f\t%.20f\t%.20f\t%f\t%f\t%f\t%.20f\t%.20f\t%f\t%f\n'   % (Z,L0,L0_err,h,h_err,L01,L01_err,h1,h1_err,Rbr1,L02,L02_err,h2,h2_err,Rbr2,L03,L03_err,h3,h3_err) )
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
















'''    
if __name__ == "__main__":      
    input_image = 'galaxy.fits'
    mask_image = None#'mask.fits'
    weight_image = None#'sigma.fits'
    R_in_l = 20.
    R_out_l = 170.
    R_in_r = 20.
    R_out_r = 170.
    z_in_dn = 0.
    z_out_dn=None
    z_in_up=0.
    z_out_up = None
    break_radii = None
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, break_radii=[50,100], bin_type='linear', side='1', error_est=True)    

'''

'''
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
    break_radii = None
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, break_radii=None, bin_type='linear', side='1', error_est=True)  
'''

'''
if __name__ == "__main__":      
    input_image = 'NGC4302.phot.1_nonan_rot_crop_aver.fits'
    mask_image = 'NGC4302.1.finmask_nonan_rot_crop_aver.fits'
    weight_image = 'NGC4302_sigma2014_rot_crop_aver.fits'
    R_in_l = 10.
    R_out_l = None
    R_in_r = 10.
    R_out_r = None
    z_in_dn = 0.
    z_out_dn=None
    z_in_up=0.
    z_out_up = None
    break_radii = None
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, break_radii=None, bin_type='linear', side='1', error_est=True)  
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
    break_radii = None
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, break_radii=None, bin_type='linear', side='1', error_est=False)  
'''

'''
if __name__ == "__main__":      
    input_image = 'model_rot_crop.fits'
    mask_image = 'mask_rot_crop.fits'
    weight_image = 'sigma_rot_crop.fits'
    R_in_l = 10.
    R_out_l = None
    R_in_r = 10.
    R_out_r = None
    z_in_dn = 0.
    z_out_dn=None
    z_in_up=0.
    z_out_up = None
    break_radii = None
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, break_radii=None, bin_type='linear', side='1', error_est=False)  
'''

if __name__ == "__main__":      
    input_image = 'model_gal.fits'
    mask_image = None
    weight_image = None
    R_in_l = 10.
    R_out_l = None
    R_in_r = 10.
    R_out_r = None
    z_in_dn = 0.
    z_out_dn=None
    z_in_up=0.
    z_out_up = None
    break_radii = None
    main(input_image, mask_image=mask_image, weight_image=weight_image, psf_image='psf.fits', R_in_l=R_in_l, R_out_l=R_out_l, R_in_r=R_in_r, R_out_r=R_out_r, z_in_dn = z_in_dn, z_out_dn=z_out_dn, z_in_up=z_in_up, z_out_up=z_out_up, break_radii=None, bin_type='linear', side='1', error_est=False)  
