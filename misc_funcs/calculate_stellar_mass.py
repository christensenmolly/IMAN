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
import re
import glob
from scipy.interpolate import splev, splrep

from scipy.optimize import curve_fit
from scipy import optimize
from scipy import special
from uncertainties import ufloat
from uncertainties.umath import *


def calc_Mag(mag, Dist, Aext):
    return mag - 5.*np.log10(Dist) - Aext - 25.
    
def calc_scale(Dist):
    return Dist*1000./206265. # kpc per arcsec

def Jansky_to_AB(flux):
    return -2.5*np.log10(flux) + 8.90

def AB_to_Jansky(AB_mag):    
    return 3631.*10**(-0.4*AB_mag)

def nuLnu_Mstar_coef(Type='Late'):
  # Table 1 from Wen et al. 2013 (Mon Not R Astron Soc (2013) 433 (4): 2946-2957), Eq. (1) http://adsabs.harvard.edu/abs/2013MNRAS.433.2946W
  if Type=='Total':
    A,B = [-0.040,1.12]
  if Type=='HII':
    A,B = [0.779,1.019]    
  if Type=='Composites':
    A,B = [1.081,1.003]
  if Type=='AGNs':
    A,B = [1.132,1.000]
  if Type=='Absorptions':
    A,B = [0.774,1.041]
  if Type=='Emissions':
    A,B = [0.685,1.051]
  if Type=='Early':
    A,B = [0.761,1.044]
  if Type=='Late':
    A,B = [0.679,1.033]
  return A,B

def find_Lsun_Msun(mag_3p4_AB, Dist, Type):
    #http://adsabs.harvard.edu/abs/2013MNRAS.433.2946W
    Jy = 1E-6 * 10.0**( ( 23.9 - mag_3p4_AB )  /2.5 )

    Wt = Jy * 10**(-26.)
    Lnu = 4.*pi*(Dist*10**6*(3.0857*10**16))**2*Wt
    nuLnu =  88235294117647.06*Lnu

    Lsun = nuLnu / (3.828*10**26)
    log_Mstar_Msun = []

    if Type<=0.:
            TYPE = 'Early'
    else:
            TYPE = 'Late'
        
    A,B = nuLnu_Mstar_coef(Type=TYPE)    
    log_Mstar_Msun = A + B*log10(Lsun)
    
    return log10(Lsun),log_Mstar_Msun

def find_Lsun_Msun_Meidt(mag_3p6_AB, Dist):
    #http://faculty.ucr.edu/~gillianw/cal.html
    #https://arxiv.org/pdf/1505.03534.pdf
    # https://www.researchgate.net/figure/Empirical-mass-to-light-relation-M-L-derived-from-WISE-34-m-m-in-band-luminosity_fig14_232244453
    mag_3p6_Vega = mag_3p6_AB - 2.79
    Mag_3p6_Vega = calc_Mag(mag_3p6_Vega, Dist, 0.)

    Lsun = 10**(0.4*(3.24-Mag_3p6_Vega))
    M_L = 0.65#0.53 # New value taken from https://arxiv.org/pdf/1709.08316.pdf
    log_Mstar_Msun_Meidt = log10( Lsun * M_L)
    
    return log_Mstar_Msun_Meidt


def find_Lsun_Msun_Eskew(F_3p6, F_4p5, Dist):
    #https://arxiv.org/pdf/1204.0524.pdf
    Msun = 10**(5.65) * F_3p6**(2.85) * F_4p5**(-1.85) * (Dist/0.05)**2    
    return log10(Msun)

def calc_incl(logr25, logr0):
    return degrees(math.asin(math.sqrt(  (1.-10**(-2.*logr25)) / (1.-10**(-2.*logr0))    )))
















def main(mag_ap, Dist):
    
    mag_ap_AB = float(mag_ap) + 2.699  # Convert: Vega -> AB
    mag_ap_AB = mag_ap_AB + 0.034 # Correction for non-PSF source


    Aext_W1 = 0.

    
    mag_ap_AB = mag_ap_AB - Aext_W1
    mag_3p6_AB = mag_ap_AB + 2.5*log10(1.027)
    

    
    log_Mstar_Msun_Meidt = find_Lsun_Msun_Meidt(mag_3p6_AB, Dist)
    print(log_Mstar_Msun_Meidt)
    return log_Mstar_Msun_Meidt

main(10.46, 32.89)
main(9.04, 16.8)
main(7.77, 16.75)