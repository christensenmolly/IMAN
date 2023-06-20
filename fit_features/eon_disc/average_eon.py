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

sys.path.append( '/home/amosenko/MEGA/MyPrograms/IMAN/imp/misc')
import flip_image



def main(input_image, mask_image=None, weight_image=None):
    flip_image.main(input_image, output_image='ima1_tmp.fits', flip='lr')
    flip_image.main(input_image, output_image='ima2_tmp.fits', flip='ud')
    flip_image.main('ima2_tmp.fits', output_image='ima3_tmp.fits', flip='lr')
    
    if mask_image is not None:
        flip_image.main(mask_image, output_image='ma1_tmp.fits', flip='lr')
        flip_image.main(mask_image, output_image='ma2_tmp.fits', flip='ud')
        flip_image.main('ma2_tmp.fits', output_image='ma3_tmp.fits', flip='lr')        
    
    if weight_image is not None:
        flip_image.main(weight_image, output_image='wt1_tmp.fits', flip='lr')
        flip_image.main(weight_image, output_image='wt2_tmp.fits', flip='ud')
        flip_image.main('wt2_tmp.fits', output_image='wt3_tmp.fits', flip='lr')

        wt0 = pyfits.open(weight_image)[0].data
        header_wt0 = pyfits.open(weight_image)[0].header
        wt1 = pyfits.open('wt1_tmp.fits')[0].data
        wt2 = pyfits.open('wt2_tmp.fits')[0].data
        wt3 = pyfits.open('wt3_tmp.fits')[0].data 
        wt_aver = (wt0 + wt1 + wt2 + wt3)/4.
    

    ima0 = pyfits.open(input_image)[0].data
    header = pyfits.open(input_image)[0].header
    ima1 = pyfits.open('ima1_tmp.fits')[0].data
    ima2 = pyfits.open('ima2_tmp.fits')[0].data
    ima3 = pyfits.open('ima3_tmp.fits')[0].data 
    ima_aver = (ima0 + ima1 + ima2 + ima3)/4.
    


    
    
    if mask_image is not None:
        ny,nx =  np.shape(ima0) 
        
        ma0 = pyfits.open(mask_image)[0].data
        ma1 = pyfits.open('ma1_tmp.fits')[0].data
        ma2 = pyfits.open('ma2_tmp.fits')[0].data
        ma3 = pyfits.open('ma3_tmp.fits')[0].data
        
        ma_sum = ma0 + ma1 + ma2 + ma3
        
        if weight_image is None:
            for k in range(ny):
                for i in range(nx):
                    #if ma_sum[k,i]<4:
                    if ma0[k,i]==0 or ma1[k,i]==0 or ma2[k,i]==0 or ma3[k,i]==0:
                            sum_val = 0.; N_val = 0
                            mask_vals = [ma0[k,i],ma1[k,i],ma2[k,i],ma3[k,i]]
                            vals = [ima0[k,i],ima1[k,i],ima2[k,i],ima3[k,i]]
                            for l in range(len(vals)):
                                if mask_vals[l]==0:
                                    ma_sum[k,i]=0
                                    sum_val = sum_val + vals[l]
                                    N_val = N_val + 1
                            ima_aver[k,i] = sum_val/N_val
        else:
            for k in range(ny):
                for i in range(nx):
                    #if ma_sum[k,i]<4:
                    if ma0[k,i]==0 or ma1[k,i]==0 or ma2[k,i]==0 or ma3[k,i]==0:
                            sum_val = 0.; N_val = 0; sum_wt = 0.
                            mask_vals = [ma0[k,i],ma1[k,i],ma2[k,i],ma3[k,i]]
                            vals = [ima0[k,i],ima1[k,i],ima2[k,i],ima3[k,i]]
                            vals_wt = [wt0[k,i],wt1[k,i],wt2[k,i],wt3[k,i]]
                            for l in range(len(vals)):
                                if mask_vals[l]==0:
                                    ma_sum[k,i]=0
                                    sum_val = sum_val + vals[l]
                                    sum_wt = sum_wt + vals_wt[l]
                                    N_val = N_val + 1
                            ima_aver[k,i] = sum_val/N_val
                            wt_aver[k,i] = sum_wt/N_val
                                                
    hdu = pyfits.PrimaryHDU(ima_aver, header)
    hdu.writeto(input_image.split('.fits')[0]+'_aver.fits', clobber=True)           
    
    if mask_image is not None:
        hdu = pyfits.PrimaryHDU(ma_sum, header)
        hdu.writeto(mask_image.split('.fits')[0]+'_aver.fits', clobber=True)        
    
    if weight_image is not None:
        hdu = pyfits.PrimaryHDU(wt_aver, header)
        hdu.writeto(weight_image.split('.fits')[0]+'_aver.fits', clobber=True)   
        
'''
input_image = 'NGC3628.phot.1_nonan_rot_crop.fits'
#mask_image = 'NGC3628.1.finmask_nonan_rot_crop.fits'
mask_image = 'new_mask_rot_crop.fits'
weight_image = 'NGC3628_sigma2014_rot_crop.fits'
'''
'''
input_image = 'NGC4302.phot.1_nonan_rot_crop.fits'
mask_image = 'NGC4302.1.finmask_nonan_rot_crop.fits'
weight_image = 'NGC4302_sigma2014_rot_crop.fits'
'''
input_image = 'NGC891_coadd_rot_crop.fits'
mask_image = 'NGC891_coadd_mask_rot_crop.fits'
weight_image = 'NGC891_coadd_sigma_rot_crop.fits'
main(input_image, mask_image=mask_image, weight_image=weight_image)