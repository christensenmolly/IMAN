#! /usr/bin/env python
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

def incl_dust_lane(a,z_cen,delta_z):
  '''
  Estimation of the inclination angle for a near edge-on galaxy
  '''
  # a is the major axis of the outer isophote
  # z_cen is the center of the dust lane
  # delta_z is the width of the dust lane
  incl = degrees(math.acos(2.*z_cen/a))
  #incl_plus = degrees(math.acos(2.*(z_cen-delta_z)/a)) - incl  
  #incl_minus = degrees(math.acos(2.*(z_cen+delta_z)/a)) - incl  
  d_incl = fabs(-1./sqrt(1.-4.*z_cen*z_cen/(a*a)) * 2./a * delta_z)
  return incl,degrees(d_incl)


def incl_dust_torus(a,b,d1,d2):
  '''
  Estimation of the inclination angle for a torus
  '''
  # a is the major axis of the torus
  # b is the minor axis of the torus
  # d1 is the viewed heighth of the down or up side of the torus
  # d2 is the viewed length of the left of right side of the torus
  alfa = b - 2.*d1
  def f(x):
    return a/2. - d2 - (a - 2.*(d1+a*x-b)/x)*sqrt(1.-((b-a*x)/(2.*(b-a*x+alfa)))**2)
  x_initial_guess = 0.15
  x_solution = fsolve(f, x_initial_guess)
  incl1 = degrees(math.acos(x_solution))
  D = (d1-b+a*cos(radians(incl1)))/cos(radians(incl1))
  h = (b-a*cos(radians(incl1)))/sin(radians(incl1))    
  return incl1,D,h
  

def Hubble_formula(q,q0=0.25,q0_sigma=0.05):
  '''
  Estimation of the inclination angle for an inclined disc
  '''  
  # q is the apparent axis ratio
  # q0 is the true axis ratio
  # q0_sigma is the dispersion of the true axis ratio
  AA = (q*q-q0*q0)/(1.-q0*q0)
  incl = degrees(math.acos(sqrt(AA)))
  d_incl = -1./sqrt(1.-AA) * 1./(2.*sqrt(AA)) * ( 2.*q*q*q0/((1.-q0*q0)**2) - 2.*q0/(1.-q0*q0) - 2.*q0*q0*q0/((1.-q0*q0)**2)  ) * q0_sigma
  return degrees(math.acos(sqrt((q*q-q0*q0)/(1.- q0*q0)))),degrees(d_incl)


def Hubble_formula_eon_ring(a,b,h):
  '''
  Simple estimation of the inclination of the near edge-on ring
  '''
  # a is the major axis
  # b is the minor axis
  # h is the viewed heigth of the near eon ring 
  q = b/a
  q0 = h/a
  incl = Hubble_formula(q,q0)
  # Correction of the heigth of the ring for the inclination
  h = h*sin(radians(incl))
  incl = Hubble_formula(q,h/a)
  return incl

def main(sma, delta_z):
    incl = degrees(math.acos(delta_z/sma))
    print('Inclination (deg) %.1f' % (round(incl)))
    return incl


if __name__ == '__main__':


    sma = float(sys.argv[1])
    delta_z = float(sys.argv[2]) 
    
    main(sma, delta_z)
