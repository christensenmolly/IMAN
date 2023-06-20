#!/usr/bin/python
# -*- coding:  cp1251 -*-

# Import standard modules
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
import warnings
from scipy import special
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy import optimize
import heapq
import pyfits
from scipy import interpolate
import argparse

#from scipy.optimize import curve_fit
from scipy import optimize


sys.path.append('/home/amosenko/MyPrograms/Stats')
import regression_analysis

def line_intersect(k1, b1, k2, b2):
    if k1 == k2:
        print ("These lines are parallel!!!")
        return None
    # y = kx + b
    # Set both lines equal to find the intersection point in the x direction
    # k1 * x + b1 = k2 * x + b2
    # k1 * x - k2 * x = b2 - b1
    # x * (k1 - k2) = b2 - b1
    # x = (b2 - b1) / (k1 - k2)
    x = (b2 - b1) / (k1 - k2)
    # Now solve for y -- use either line, because they are equal here
    # y = mx + b
    y = k1 * x + b1
    return x,y



# Function to find intersection between two lines (with errors)
def find_intersection_of_two_lines(y1, y2, y1_err=None, y2_err=None, Number_of_trials=1000):
     [k1,b1] = y1
     [k2,b2] = y2
     
     if y1_err!=None:
        [k1_err,b1_err] = y1_err
      
     if y2_err!=None:
        [k2_err,b2_err] = y2_err
     
     # Find true intersection:
     x_int,y_int = line_intersect(k1, b1, k2, b2)
     
     
     # Find std:
     if y1_err!=None or y2_err!=None:
        X_int = []; Y_int = []
        for N in range(Number_of_trials):
            if y1_err!=None:
                K1 = random.gauss(k1, k1_err)
                B1 = random.gauss(b1, b1_err)
            else:
                K1 = k1
                B1 = b1

            if y2_err!=None:
                K2 = random.gauss(k2, k2_err)
                B2 = random.gauss(b2, b2_err)
            else:
                K2 = k2
                B2 = b2
                
            xx_int,yy_int = line_intersect(K1, B1, K2, B2)
            X_int.append(xx_int)
            Y_int.append(yy_int)

        x_int_err = np.std(X_int)
        y_int_err = np.std(Y_int)
        
        return [x_int,y_int],[x_int_err,y_int_err]
     else:
        return [x_int,y_int],[None,None]







'''
def piecewise_linear(x, x_break, y_break, k1, k2):
    return np.piecewise(x, [x < x_break], [lambda x:k1*x + y_break-k1*x_break, lambda x:k2*x + y_break-k2*x_break])


def main(x, y, x_break=None, y_break=None, method_to_fit='from_center'):
    if method_to_fit=='from_center':
                        # Taken from bps_disc_com_regions:
			if inter_mode==True and repeat==2:
			  x_break_left = R_break_left_man
			  y_break_left = Mag_break_left_man
			  def piecewise_linear_fix2(x, k1, k2):
			      return np.piecewise(x, [x < x_break_left], [lambda x:k1*x + y_break_left-k1* x_break_left, lambda x:k2*x +y_break_left-k2* x_break_left])
			  results_neg1, e = optimize.curve_fit(piecewise_linear_fix2, R, Mag,[1,1])     
'''


def main(x, y, y_err=None, x_break=None, method_to_fit='from_center'):
    
    if x_break!=None:
        x_left = []
        y_left = []
        if y_err!=None:
            y_left_err = []
        else:
            y_left_err = None

        x_right = []
        y_right = []
        if y_err!=None:
            y_right_err = []
        else:
            y_right_err = None

        
        for k in range(len(x)):
            if x[k]<x_break:
                x_left.append(x[k])
                y_left.append(y[k])
                
                if y_err!=None:
                    y_left_err.append(y_err[k])

            else:
                x_right.append(x[k])
                y_right.append(y[k])
                
                if y_err!=None:
                    y_right_err.append(y_err[k])                    

    # Fit both ranges:
    [k_left,b_left],[k_left_err,b_left_err],r_left = regression_analysis.polyfit_func(x_left, y_left, y_err=y_left_err, degree=1)
    [k_right,b_right],[k_right_err,b_right_err],r_right = regression_analysis.polyfit_func(x_right, y_right, y_err=y_right_err, degree=1)
    
    [x_break_new,y_break_new],[x_break_new_err,y_break_new_err] = find_intersection_of_two_lines([k_left,b_left], [k_right,b_right], y1_err=[k_left_err,b_left_err], y2_err=[k_right_err,b_right_err], Number_of_trials=1000)
    
    return [x_break_new,y_break_new],[x_break_new_err,y_break_new_err]



    
    