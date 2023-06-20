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

def main(input_dir, verbose=False):
    os.chdir(input_dir)
    if os.path.exists('mean_bestfit_parameters_imfit.dat'):
        os.remove('mean_bestfit_parameters_imfit.dat')
    
    files = glob.glob('*_bestfit_parameters_imfit.dat')
    #print(files)
    #exit()
    res = OrderedDict() 
    
    for i in range(len(files)):
        file = files[i]
        print(file)
        f = open(file, 'r')
        lines = f.readlines()
        for k in range(len(lines)):
            if lines[k][0]!='#' and 'FUNCTION' not in lines[k] and lines[k]!='\n':
                print(lines[k])
                pars = lines[k].split()
                key = pars[0]
                value = float(pars[1])
                if i==0:
                    res['%i' % (k)] = []
                    res['%i' % (k)].append(value) 
                else:
                    res['%i' % (k)].append(value) 
        f.close()

    
    ff = open('mean_bestfit_parameters_imfit.dat', 'w')
    f = open(files[0], 'r')
    lines = f.readlines()
    for k in range(len(lines)):
        try:
            value = np.mean(res['%i' % (k)])
            value_std = np.std(res['%i' % (k)])
            pars = lines[k].split()
            key = pars[0]
                
            ff.write('%s\t%.5f # +/- %.5f\n' % (key, value, value_std))
        except:
            ff.write('%s\n' % (lines[k]))
    ff.close()
    f.close()
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Script to calculate mean and error of decomposition")
    parser.add_argument("input_dir", help="Input directory with imfit file and all data")

   
    
    args = parser.parse_args()

    input_dir = args.input_dir

    
    main(input_dir, verbose=False)
