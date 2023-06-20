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
from distutils.dir_util import copy_tree
from collections import OrderedDict 
import random
from astropy.stats import sigma_clipped_stats
import warnings
warnings.filterwarnings("ignore")


def main(input_image, imfit_output_files, m0):
    Components = []
    Fluxes = []
    Magnitudes = []
    Fractions = []
    for imfit_output_file in imfit_output_files:
        print(imfit_output_file)
        subprocess.call('makeimage -print-fluxes %s --refimage %s --zero-point %s > err_model.txt' % (imfit_output_file, input_image, str(m0)), shell=True)
        
        components = []
        fluxes = []
        magnitudes = []
        fractions = []
        f = open('err_model.txt', 'r')
        lines = f.readlines()
        for k in range(len(lines)):
            if 'Component                 Flux        Magnitude  Fraction' in lines[k]:
                k_start = k + 1
            if 'Total' in lines[k]:
                k_end = k - 1
        
        for k in range(k_start, k_end, 1):
            s = lines[k].split()
            components.append(s[0])
            fluxes.append(float(s[1]))
            try:
                magnitudes.append(float(s[2]))
            except:
                magnitudes.append(float('nan'))
            fractions.append(float(s[3].split('\n')[0]))
            
        
        Components.append(components)
        Fluxes.append(fluxes)
        Magnitudes.append(magnitudes)
        Fractions.append(fractions)
    
    all_components = components
    Components = np.array(Components)
    Fluxes = np.array(Fluxes)
    Magnitudes = np.array(Magnitudes)
    Fractions = np.array(Fractions)
    #print(Components)
    #print(Fluxes)
    #print(Magnitudes)
    #print(Fractions)
    
    for k in range(len(all_components)):
        mean_flux, median_flux, stddev_flux = sigma_clipped_stats(Fluxes.T[k])
        mean_mag, median_mag, stddev_mag = sigma_clipped_stats(Magnitudes.T[k])
        mean_frac, median_frac, stddev_frac = sigma_clipped_stats(Fractions.T[k])
        print('%s:' % (all_components[k]))
        print('\tFlux:%.5f+/-%.5f' % (mean_flux, stddev_flux))
        print('\tMag:%.2f+/-%.2f' % (mean_mag, stddev_mag))
        print('\tFrac:%.3f+/-%.3f\n' % (mean_frac, stddev_frac))

'''
#NGC4452
m0 = 28.098
input_image = 'combined_reconv_res_rot_crop.fits'
'''

#'''
#NGC4469
m0 = 28.294
input_image = 'combined.fits'
#'''

#NGC509
#m0 = 28.15
#input_image = 'combined_rot.fits'

#main(input_image, ['bestfit_parameters_imfit.dat'], m0)
#exit()

imfit_output_files = []
for k in range(1,11):
    imfit_output_files.append('%i_bestfit_parameters_imfit.dat' % (k))
main(input_image, imfit_output_files, m0)