#!/usr/bin/python
# -*- coding:  cp1251 -*-

import sys
import math
import numpy as np
from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import NullFormatter
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
import pyfits

radius_find = 10
FIND_MAJOR_OBJECT = 0


def read_sextr_find(RA,DEC,file_out,coords=1):
        # FUNCTION TO FIND THE OBJECT AND RETURN ITS PARAMETERS FROM THE SEXTRACTOR CATALOGUE
        number,x_image,y_image,x_world,y_world,XPEAK_IMAGE,YPEAK_IMAGE,XPEAK_WORLD,YPEAK_WORLD,flux_radius25,flux_radius50,flux_radius75,flux_radius99,mag_auto,a_image,b_image,theta_image,ellipticity,kron_radius,backgr,class_star = loadtxt(file_out, usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], unpack=True, skiprows = 18)


        #if type(x_image)!='numpy.float64':        #dtype for 1 object!
        if type(x_image)==np.ndarray:        #dtype for 1 object!
                #print x_world


                if fabs(max(x_world)-fabs(min(x_world)))>1. and fabs(RA-fabs(max(x_world)))>1.:                #### !!!!!!!!!!!!!!!! EQUAL to 1 DEGREE !!!
                        radius_find1 = radius_find
                        print 'I am here coords=1'
                elif fabs(max(x_world)-fabs(min(x_world)))<1. and fabs(RA-fabs(max(x_world)))<1.:        
                        radius_find1 = radius_find/3600.        
                        print 'I am here coords=2'
                elif fabs(max(x_world)-fabs(min(x_world)))<1. and fabs(RA-fabs(max(x_world)))>1.:        
                        radius_find1 = radius_find
                        x_world = x_image
                        y_world = y_image
                        print 'I am here coords=3'



                n_iter = len(x_image)

                number_found = 0
                PA1 = [0.,0.]; xc1 = [0.,0.]; yc1 = [0.,0.]; RA_c1 = [0.,0.]; DEC_c1 = [0.,0.]; a1 = [0.,0.]; b1 = [0.,0.]
                flu_rad1 = [0.,0.]; ma_auto1 = [0.,0.]; ell1 = [0.,0.]; kron_r1 = [0.,0.]; C311 = [0.,0.]; rmax1 = [0.,0.]


                
                if FIND_MAJOR_OBJECT == 0:
                                a_image = list(a_image)
                                iii = a_image.index(max(a_image))

                                print '**0**'
                                PA = - theta_image[iii]
                                xc = x_image[iii]
                                yc = y_image[iii]
                                RA_c = x_world[iii]
                                DEC_c = y_world[iii]
                                a = float(a_image[iii])
                                b = float(b_image[iii])
                                flu_rad = flux_radius50[iii]
                                ma_auto = mag_auto[iii]
                                ell = ellipticity[iii]
                                kron_r = kron_radius[iii]
                                C31 = flux_radius75[iii]/flux_radius25[iii]
                                rmax = flux_radius99[iii]                         



                if FIND_MAJOR_OBJECT == 1:
                        for k in range(n_iter):
                                if fabs(RA-x_world[k])< radius_find1 and fabs(DEC-y_world[k])< radius_find1 and a_image[k]*kron_radius[k] > semimajor_find and class_star[k]<0.1:
                                        PA1.append(- theta_image[k])
                                        xc1.append(x_image[k])
                                        yc1.append(y_image[k])
                                        RA_c1.append(x_world[k])
                                        DEC_c1.append(y_world[k])
                                        a1.append(float(a_image[k]))
                                        b1.append(float(b_image[k]))
                                        flu_rad1.append(flux_radius50[k])
                                        ma_auto1.append(mag_auto[k])
                                        ell1.append(ellipticity[k])
                                        kron_r1.append(kron_radius[k])
                                        C311.append(flux_radius75[k]/flux_radius25[k]) 
                                        rmax1.append(flux_radius99[k])
                                        number_found = number_found + 1
                                        print '**1**'



                        if number_found==0:
                                        PA = - theta_image[0]
                                        xc = x_image[0]
                                        yc = y_image[0]
                                        RA_c = x_world[0]
                                        DEC_c = y_world[0]
                                        a = float(a_image[0])
                                        b = float(b_image[0])
                                        flu_rad = flux_radius50[0]
                                        ma_auto = mag_auto[0]
                                        ell = ellipticity[0]
                                        kron_r = kron_radius[0]
                                        C31 = flux_radius75[0]/flux_radius25[0]
                                        rmax = flux_radius99[0] 
                                        print '**2**'

                        elif number_found==1:
                                        PA = PA1[2]
                                        xc = xc1[2]
                                        yc = yc1[2]
                                        RA_c = RA_c1[2]
                                        DEC_c = DEC_c1[2]
                                        a = a1[2]
                                        b = b1[2]
                                        flu_rad = flu_rad1[2]
                                        ma_auto = ma_auto1[2]
                                        ell = ell1[2]
                                        kron_r = kron_r1[2]
                                        C31 = C311[2]
                                        rmax = rmax1[2]
                                        print '**3**'

                        elif number_found>1:
                                        iii = a1.index(max(a1))
                                        PA = PA1[iii]
                                        xc = xc1[iii]
                                        yc = yc1[iii]
                                        RA_c = RA_c1[iii]
                                        DEC_c = DEC_c1[iii]
                                        a = a1[iii]
                                        b = b1[iii]
                                        flu_rad = flu_rad1[iii]
                                        ma_auto = ma_auto1[iii]
                                        ell = ell1[iii]
                                        kron_r = kron_r1[iii]
                                        C31 = C311[iii]
                                        rmax = rmax1[iii]
                                        print '**4**'
        else:
                                        PA = - theta_image
                                        xc = x_image
                                        yc = y_image
                                        RA_c = x_world
                                        DEC_c = y_world
                                        a = float(a_image)
                                        b = float(b_image)
                                        flu_rad = flux_radius50
                                        ma_auto = mag_auto
                                        ell = ellipticity
                                        kron_r = kron_radius
                                        C31 = flux_radius75/flux_radius25
                                        rmax = flux_radius99
                                        print '**5**'


                                
        return xc,yc,RA_c,DEC_c,PA,a,b,flu_rad,ma_auto,ell,kron_r,C31,rmax
