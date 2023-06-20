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

import random

LOCAL_DIR = "/decomposition/ControlFit"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'decomposition/make_model'))
sys.path.append(os.path.join(IMAN_DIR, 'plotting/2dprofile'))
sys.path.append(os.path.join(IMAN_DIR, 'plotting/1dprofile'))
sys.path.append(os.path.join(IMAN_DIR, 'deca/deca_tk_lib'))
import make_model_ima_galfit
import make_model_ima_imfit
import plot_2d_profile
import plot_profile
import radial_profile
import tex_creator

tmp_out = sys.stdout
FNULL = open(os.devnull, 'w')

galfit_path = ''

import ControlFit

def main(input_dir, sky_rms, n_trials, input_file = 'imfit.inp', verbose=False):
    
    
    os.chdir(input_dir)
    if os.path.exists(input_dir+'/error_trials'):
        shutil.rmtree(input_dir+'/error_trials')
    
    copy_tree(input_dir, input_dir+'/error_trials')

    os.chdir(input_dir+'/error_trials')
    
    code, code_run, output_file = ControlFit.recognize_code(input_file)
    input_image, sigma_image, psf_image, mask_image, pix2sec, m0 = ControlFit.get_input_images(input_file, code, 1., 1.)
  
    shutil.copy(input_image, 'ini_image.fits')
    galaxyHDU = pyfits.open('ini_image.fits')
    galaxyData = galaxyHDU[0].data        
    galaxyHeader = galaxyHDU[0].header 

      
    for k in range(n_trials):
        print('Trial #%i' % (k+1))
        
        sky_error = random.gauss(0., sky_rms)
        print(sky_error)
        
        newData = np.copy(galaxyData)
        
        newData = newData + sky_error

        hdu = pyfits.PrimaryHDU(newData, galaxyHeader)
        hdu.writeto(input_image, overwrite=True)

        # Start fitting
        if verbose:
            subprocess.call(code_run, shell=True, stdout=FNULL)
        else:
            subprocess.call(code_run, shell=True)
            
        os.rename('bestfit_parameters_imfit.dat', '%i_bestfit_parameters_imfit.dat' % (k+1))
    
    os.rename('ini_image.fits', input_image)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Script to control fitting using an existing galfit file")
    parser.add_argument("input_dir", help="Input directory with imfit file and all data")

    parser.add_argument("--sky_rms", help="Optional: RMS of the sky background", type=float, default=0.1)
    parser.add_argument("--n_trials", help="Optional: Number of trials", type=int, default=10)
   
    
    args = parser.parse_args()

    input_dir = args.input_dir
    sky_rms = float(args.sky_rms)
    n_trials = int(args.n_trials)


    
    main(input_dir, sky_rms, n_trials)
