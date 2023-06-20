#!/usr/bin/python
# -*- coding:  cp1251 -*-
import random as random_number
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
#import subprocess
import signal
import pyfits



tmp_out = sys.stdout

import subprocess as sub


def crea_input_file(IMAGES,INFO,ADD_INFO,nx,ny):
  [input_image,mask_image,sigma_image,psf_image] = IMAGES
  [scale,m0,gain,read_out_noise,ncombine,exptime,fwhm] = ADD_INFO
  [number,FILTER,name,sampling,distance,fitting_proc] = INFO
  
  f = open('input.txt', "w") 
  sys.stdout = f  
  header(input_image,sigma_image,psf_image,mask_image,'none','model.fits',1,nx,1,ny,m0,scale,sampling)
  