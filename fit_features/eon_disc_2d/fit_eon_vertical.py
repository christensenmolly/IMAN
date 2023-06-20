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


def crea_galfit_input(input_image, sigma_image, psf_image, mask_image, xc, yc, mu0d, z0, h):
    if sigma_image is None:
        sigma_image = 'none'

    if mask_image is None:
        mask_image = 'none'

    hdu = pyfits.open(input_image)
    image = hdu[0].data
    ny,nx =  np.shape(image)
   
        
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
""" % (input_image, sigma_image, psf_image, mask_image, nx, ny, nx, ny)

    s2="""
# Component number: 1
 0) edgedisk               #  Component type
 1) %f %f %i %i  #  Position x, y
 3) %f     1          #     Mu(0)   [mag/arcsec^2]
 4) %f      1          #  h_s (disk scale-height)   [pix]
 5) %f      1          #  R_s (disk scale-length)   [pix]
10) 90.0     0          #  Position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)
 ================================================================================
""" % (xc, yc, 0, 1, mu0d, z0, z01, h1)

 
    s = s1 + s2
    
    f = open("input_galf.txt", 'w')
    f.write(s)
    f.close()


def crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind, y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up):
    I_image = []; I_mask = []; I_sigma = []; rr = []; r_no_mask = []; I_no_mask = []
    ny,nx = np.shape(image_data)

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
            

    return np.array(rr), np.array(I_image), np.array(I_mask), np.array(I_sigma), np.array(r_no_mask), np.array(I_no_mask)

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

def ds9_to_np(x):
    '''
    Function to convert ds9 pixel coordinates to numpy format, i.e
    0.5 pix is 0, 1.5 pix is 1
    '''
    return int(round(x)) - 1


def crea_border_mask(nx, ny, R_in=0, R_out=None, z_in=0, z_out=None):
    '''
    Create border mask:
    '''
    xc = nx/2.
    yc = ny/2.
    
    if R_out is None:
        R_out = nx/2.
    
    if z_out is None:
        z_out = ny/2.
    
    mask = np.zeros((ny,nx))
    
    for k in range(ny):
        for i in range(nx):
            if k<ds9_to_np(yc-z_out) or k>ds9_to_np(yc+z_out) or (k>ds9_to_np(yc-z_in) and k<ds9_to_np(yc+z_in)): 
                mask[k,i]=1
            if i<ds9_to_np(xc-R_out) or i>ds9_to_np(xc+R_out) or (i>ds9_to_np(xc-R_in) and i<ds9_to_np(xc+R_in)): 
                mask[k,i]=1
    
    #hdu = pyfits.PrimaryHDU(mask)
    #hdu.writeto('border_mask.fits', clobber=True)    
    #return 'border_mask.fits'
    return mask
    

def define_xbins(rmin, rmax, Rmax, bin_type='log'):
    if rmax is None:
        rmax = Rmax

    if bin_type=='log':
        return Rmax*np.power(10, np.arange(math.log10(1./Rmax), math.log10(rmax/Rmax), 0.1))
    if bin_type=='linear':
        return Rmax*np.arange(rmin/Rmax, rmax/Rmax, 0.1)        
    if bin_type=='pix':
        return np.arange(rmin, rmax, 10)   

def disk_exp(r, I0d, h, x_c):
    return I0d * np.exp(-np.fabs(r-x_c)/h)

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


def do_each_bin_x(R, input_image, sigma_image, mask_data, psf_image, x_min_ind, x_max_ind, y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up, disc_comps='single', error_est=True, ron=None, gain=None, n_fixed='fixed', side=1):
    print('SLICE X=%i:%i (R=%.1f pix)' % (x_min_ind, x_max_ind, R) )

    hdu = pyfits.open(input_image)
    image_data = hdu[0].data 

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
    
    # Mask all data except slice and add mask_data
    ny,nx = np.shape(mask_data)
    xc = nx/2.
    yc = ny/2.
    except_scice_mask = np.zeros((ny,nx))
    for k in range(ny):
        for i in range(nx):
            if i<x_min_ind or i>x_max_ind:
                except_scice_mask[k,i] = 1
    mask_data = mask_data + except_scice_mask
    hdu = pyfits.PrimaryHDU(mask_data)
    hdu.writeto('total_mask_tmp.fits', clobber=True)      
    
    # Estimate L0, h:
    '''
    I_plane = image_data[ds9_to_np(yc), x_min_ind:x_max_ind]
    r_plane = np.arange(x_min_ind, x_max_ind, 1) - ds9_to_np(xc) + 1.
    popt_soph, pcov_soph = curve_fit(lambda r, I0, h: disk_exp(r, I0, h, 0.), r_plane, I_plane, p0=(np.max(image_data), nx/10.) )
    h = popt_soph[1]
    L0 = popt_soph[0]/(2.*h)    
    '''
    r_plane,I_plane = crea_hor_average_cut(image_data, mask_data, sigma_data, y_min_ind_up, y_max_ind_up, x_min_ind, x_max_ind)
    popt_soph, pcov_soph = curve_fit(lambda r, I0, h: disk_exp(r, I0, h, xc), r_plane, I_plane, p0=(max(I_plane), nx/10.) )
    h = popt_soph[1]
    L0 = popt_soph[0]/(2.*h)  
    
    
    # Estimate z0
    rr, I_image, I_mask, I_sigma, r_no_mask, I_no_mask = crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind, y_min_ind_dn, y_max_ind_dn, y_min_ind_up, y_max_ind_up)
    L0_vert, z0, y0, N = find_single_disc_guess(r_no_mask, I_no_mask, np.max(I_no_mask), 0.)
    
    
    
    crea_galfit_input(input_image, sigma_image, psf_image, 'total_mask_tmp.fits', xc, yc, -2.5*math.log10(2.*L0*h), 2.*z0, h)
    
    
    
    
    
    
    
    
    
    
    
    
    






def main(input_image, mask_image=None, sigma_image=None, psf_image=None, R_in=0., R_out=None, z_in = 0., z_out=, bin_type='linear', n_fixed='fixed', disc_comps='single'):
    # The input images should be rotated, cropped and centered so that
    # the major axis is horizontal and the center of teh galaxy coincides with the frame center
    side = '1'
    
    print('Fit code is: %s' % (code))
    ff = open('vertical_profiles_%s.txt' % (side),'w')
    ff.close()
    
    
    # Open input image:
    hdu = pyfits.open(input_image)
    image_data = hdu[0].data 
    ny,nx =  np.shape(image_data)
    
    # Define galaxy center:
    xc = nx/2. 
    yc = ny/2.



    # Create border mask:
    border_mask_data = crea_border_mask(nx, ny, R_in=0, R_out=None, z_in=0, z_out=None)

    # Merge border mask with initial mask:
    if mask_image is not None:
        ini_hdu = pyfits.open(mask_image)
        ini_mask_data = ini_hdu[0].data
        mask_data = ini_mask_data + border_mask_data
    else:
        mask_data = border_mask_data



    # Define bins:
    r_bins = define_xbins(R_in, R_out, nx/2., bin_type)
    for k in range(len(r_bins)):
        ind_x_bins.append( ds9_to_np(xc+r_bins[k]) )




    print('Fitting 2D vertical slices in the %s side:' % (side))   
    
    f = open('vertical_fits_%s.dat' % (side), 'w')
    f.write('R\tL0\tL0_err\tY0\tY0_err\tN\tN_err\tZ0\tZ0_err\tL01\tL01_err\tN1\tN1_err\tZ01\tZ01_err\tL02\tL02_err\tN2\tN2_err\tZ02\tZ02_err\tL03\tL03_err\tN3\tN3_err\tZ03\tZ03_err\n')
    
    for k in range(0, len(ind_x_bins)-1):
            R = (r_bins[k]+r_bins[k+1])/2.
            [L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err],[L01,L01_err,N1,N1_err,Z01,Z01_err],[L02,L02_err,N2,N2_err,Z02,Z02_err],[L03,L03_err,N3,N3_err,Z03,Z03_err] = do_each_bin_x(R, input_image, sigma_image, mask_data, psf_image, ind_x_bins[k], ind_x_bins[k+1], disc_comps=disc_comps, error_est=error_est, ron=None, gain=None, n_fixed=n_fixed, side=side)
            

            f.write('%f\t%.10f\t%.10f\t%f\t%f\t%f\t%f\t%f\t%f\t%.10f\t%.10f\t%f\t%f\t%f\t%f\t%.10f\t%.10f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n'   % (R,L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err,L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err,L03,L03_err,N3,N3_err,Z03,Z03_err ) )
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