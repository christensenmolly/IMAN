# Import standard modules
import pylab
import sys
import os
import shutil
import math
import numpy as np
import scipy as sp
from numpy import *
from pylab import *
import subprocess
import argparse
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import optimize
from scipy.odr.odrpack import *
from scipy import special




# -----------------------------------------------------------------
# SERSIC EFFECTIVE SB
def meb_bulge_f(reb,mag_bulge,n,q):
	nu = 1.9987*n-0.3267
	An = 2.5*log10(2.*math.pi*n/(nu**(2.*n)) * special.gamma(2.0*n)) + 1.0857*nu
	return mag_bulge+5.*log10(reb) + An + 2.5*log10(q)
# -----------------------------------------------------------------

reb = 64
mag_bulge = 9.92
n = 0.17
q = 0.73



pix2sec = 0.75
m0 = 20.472

reb = reb*pix2sec

meb = meb_bulge_f(reb,mag_bulge,n,q)
Ieb = (pix2sec**2) * 10**(0.4*(m0-meb)) 
print(Ieb.round(2))