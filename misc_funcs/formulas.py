# python /home/amosenko/CurrentWork/ImaPrep/imfit_to_skirt_conv/crea_ski.py imf_output.dat --crea galaxy.fits,3.55,psf_rebin.fits
#! /usr/bin/env python
# ./SkirtMakeUp 
import pylab
import random as random_number
import sys
import os
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
tmp_out = sys.stdout
from scipy import stats
import re
from scipy.optimize import fsolve
import pyfits
import subprocess



G = 4.302e-3	# pc Msun^(-1) (km/s)^2





def find_grav_mass(V,hR,dist,pix2sec):
   # hR in [pix]
   # V in [km/s]
   # dist in Mpc
   # pix2sec in [arcsec]
   scale_arc = math.tan(radians(1./3600.)) * float(dist) * 1000000. # 1"=scale_arc [pc]

   scale_pix = scale_arc*pix2sec	# 1pix=scale_pix [pc]
   hR = hR*scale_pix # Now in [pc]
   print "Total mass in [Msun]: %.2E" % (4.*V**2*hR/G)
   return 4.*V**2*hR/G
 
 
print '%.2E' % (find_grav_mass(80.8,51.75,17.4,0.75)/10**8.947)