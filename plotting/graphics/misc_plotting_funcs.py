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
import matplotlib.gridspec as gridspec
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
import re
import glob

import collections
from scipy import special
import matplotlib
from astropy.stats import sigma_clip
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['pdf.fonttype'] = 42                                                                                                              
matplotlib.rcParams['ps.fonttype'] = 42   



def find_regression(xx, yy):
        k1, b1, r_value, p_value, std_err = stats.linregress(xx,yy)
        k2, b2, r_value, p_value, std_err = stats.linregress(yy,xx)

        t =  (1.0-k1*k2)/(k1+k2)
        t2 =  t/(1.0+sqrt(1.0+t*t))
        k4 = (t2+k1)/(1.0-t2*k1)
        xc = ((b2+b1*k2)/(1.0-k1*k2))
        yc = k1*xc + b1
        b4 = yc - k4*xc
        return k1,b1,k2,b2,k4,b4,r_value,std_err
    
def RegrLine(xxx, yyy, xx1, xx2, yy1, yy2, i, sigma_outliers=3.):
        xx = []
        yy = []
        if xx1>xx2:
            xx1,xx2 = xx2,xx1
        if yy1>yy2:
            yy1,yy2 = yy2,yy1


        for k in range(0,len(xxx)):
            if xxx[k]>xx1 and xxx[k]<xx2 and yyy[k]>yy1 and yyy[k]<yy2:
                xx.append(xxx[k])
                yy.append(yyy[k])

        k1,b1,k2,b2,k4,b4,r_value,std_err = find_regression(xx,yy)

        if sigma_outliers is not None:
            XX = []; YY = []; Dist = []
            for k in range(len(xx)):
                Dist.append( abs(k4*xx[k]-yy[k]+b4) / ( sqrt(k4*k4+1.) ) )
                
            std_dist = np.std(Dist)    
            for k in range(len(xx)):                
                if Dist[k]<=sigma_outliers*std_dist:
                    XX.append(xx[k])
                    YY.append(yy[k])
            k1,b1,k2,b2,k4,b4,r_value,std_err = find_regression(XX,YY)    

        if i==1:
                return k1,b1,r_value
        if i==2:
                return 1./k2,-b2/k2,r_value
        if i==3:
                return k4,b4,r_value










def plot_regr_line(x, y, sel_objects=None, minx=None, miny=None, maxx=None, maxy=None, ax=None, line_type='-', lw=2., color='black', zorder=1):
    x = np.array(x)
    y = np.array(y)
    
    if sel_objects is not None:
        x = x[sel_objects]
        y = y[sel_objects]
    
    if minx is None:
        minx = np.min(x)

    if maxx is None:
        maxx = np.max(x)

    if miny is None:
        miny = np.min(y)        

    if maxy is None:
        maxy = np.max(y)

        
    k,b,r_value = RegrLine(x, y, minx, maxx, miny, maxy, 3)
    
    
    if b>=0.:
        s = '$y = %.3f x + %.3f\,(r=%.2f)$' % (k,b,r_value)
        print('\t Regression line: Y = %.3f * X + %.3f (r=%.2f)\n' % (k,b,r_value))
    else:
        print('\t Regression line: Y = %.3f * X - %.3f (r=%.2f)\n' % (k,abs(b),r_value))
        s = '$y = %.3f x - %.3f\,(r=%.2f)$' % (k,abs(b),r_value)
        
    x_regr = np.arange(minx, maxx, (maxx-minx)/100.)
    y_regr = k*x_regr + b
    
    if ax is None:
        plt.plot(x_regr, y_regr, line_type, lw=lw, color=color, zorder=zorder)
    else:
        ax.plot(x_regr, y_regr, line_type, lw=lw, color=color, zorder=zorder)
    return s    