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
from matplotlib.patches import Polygon
import pyfits
from scipy import interpolate
import argparse

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

warnings.filterwarnings("ignore")

#***************************
scale_bar_length = 10
offset_x_factor = 0.98
offset_y_factor = 0.1
format_of_output_file = 'eps'
fsize = 15
#***************************

def crea_scale_bar(x0,x1,y0,y1,pix2sec,orientation):
  if orientation=='N':
    xc = fabs(x1)-scale_bar_length/pix2sec - (1.-offset_x_factor)*(x1-x0)
    yc = fabs(y0) + (y1-y0)* offset_y_factor
  else:
    xc = fabs(x0)+scale_bar_length/pix2sec + (1.-offset_x_factor)*(x1-x0)
    yc = fabs(y1) - (y1-y0)* offset_y_factor   
  plt.errorbar(xc, yc, xerr=scale_bar_length/pix2sec,color='lime',capsize=2,c='black')

def main(input_image,iso_reg,lines_reg,points_reg,name='',orientation='N',x0=0,x1=0,y0=0,y1=0):
    plt.clf()
    plt.close()  
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data 
    nx = data.shape[1]
    ny = data.shape[0]
    pix2sec = 0.395
    AX = plt.subplot()
  
    reg = open(iso_reg,'r')

    for Line in reg:
      x = []
      y = []
      P = []
      if 'polygon' in Line:
	    coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]
	    for k in range(0,len(coords)-1,2):
	      x.append(float(coords[k]))
	      y.append(float(coords[k+1]))
      else:
	continue
      plt.plot(x,y,color='black')
      plt.plot([x[0],x[-1]],[y[0],y[-1]],color='black')
      plt.gca().set_aspect('equal', adjustable='box')

    
    # Read points:
    reg_points = open(points_reg,'r')

    x_points = []; y_points = []
    for Line in reg_points:
	  if 'circle(' in Line:
	      x_points.append(float(Line.split('(')[1].split(',')[0]))
	      y_points.append(float(Line.split(',')[1]))
	  if 'point(' in Line:
	  	  coords = Line.split('(')[1].split(')')[0].split(',')
	  	  x_points.append(float(coords[0]))
	  	  y_points.append(float(coords[1]))
    x_points = np.array(x_points)
    y_points = np.array(y_points)    

    ZZ = np.polyfit(x_points, y_points, 3.)

    f_cube = np.poly1d(ZZ)    
    #plt.plot(x_points,y_points,'x',color='red')
    #plt.plot(x_points,f_cube(x_points),'x',color='red')
    
    # Read lines:
    reg = open(lines_reg,'r')

    Line = []
    for LL in reg:
      if 'line' in LL:
        Line.append(LL)

    print lines_reg
    x0_lines_left = float(Line[0].split('(')[1].split(',')[0])
    y0_lines_left = float(Line[0].split(',')[1])
    x1_lines_left = float(Line[0].split(',')[2])
    y1_lines_left = float(Line[0].split(',')[3].split(')')[0])

    x0_lines_mid_left = float(Line[1].split('(')[1].split(',')[0])
    y0_lines_mid_left = float(Line[1].split(',')[1])
    x1_lines_mid_left = float(Line[1].split(',')[2])
    y1_lines_mid_left = float(Line[1].split(',')[3].split(')')[0])

    x0_lines_mid_right = float(Line[2].split('(')[1].split(',')[0])
    y0_lines_mid_right = float(Line[2].split(',')[1])
    x1_lines_mid_right = float(Line[2].split(',')[2])
    y1_lines_mid_right = float(Line[2].split(',')[3].split(')')[0])

    x0_lines_right = float(Line[3].split('(')[1].split(',')[0])
    y0_lines_right = float(Line[3].split(',')[1])
    x1_lines_right = float(Line[3].split(',')[2])
    y1_lines_right = float(Line[3].split(',')[3].split(')')[0])    

    
    plt.plot([x0_lines_left,x1_lines_left],[y0_lines_left,y1_lines_left],color='r',lw=4)
    plt.plot([x1_lines_left,x1_lines_mid_left],[y1_lines_left,y1_lines_mid_left],color='r',lw=4)
    plt.plot([x0_lines_mid_right,x1_lines_mid_right],[y0_lines_mid_right,y1_lines_mid_right],color='r',lw=4)
    plt.plot([x0_lines_right,x1_lines_right],[y0_lines_right,y1_lines_right],color='r',lw=4)

    if orientation=='N':
      if x0==0 and x1==0 and y0==0 and y1==0:
	crea_scale_bar(min(x)-5,max(x)+5,min(y)-5,max(y)+5,pix2sec,orientation)
      else:
	crea_scale_bar(x0,x1,y0,y1,pix2sec,orientation) 
    else:
      if x0==0 and x1==0 and y0==0 and y1==0:
	crea_scale_bar(min(x)-5,max(x)+5,min(y)-5,max(y)+5,pix2sec,orientation)
      else:
	crea_scale_bar(x0,x1,y0,y1,pix2sec,orientation)       
      
      
      
    if x0==0 and x1==0 and y0==0 and y1==0:
      xlim(min(x)-5,max(x)+5)
      ylim(min(y)-6,max(y)+6)    
    else:
      xlim(x0,x1)
      ylim(y0,y1)      
      
    if orientation=='S':
      #print orientation
      plt.gca().invert_xaxis()
      plt.gca().invert_yaxis()
    #plt.axis('off')
    plt.gca().axes.get_xaxis().set_ticks([])
    plt.gca().axes.get_yaxis().set_ticks([])
    #xlim()
    
    #if ',' in name:
    #  name = name.split(',')[0] + ', ' + name.split(',')[1]
    
    plt.text(0.05, 0.95,name, transform=AX.transAxes,fontsize=fsize+2, fontweight='bold', va='top')
    #plt.show()

    plt.savefig(str(input_image.split('.fits')[0])+'_iso_map.png', bbox_inches='tight', pad_inches=0.01, dpi = 300)
    plt.clf()
    plt.close()   
    return str(input_image.split('.fits')[0])+'_iso_map.png',min(x)-5,max(x)+5,min(y)-5,max(y)+5
  
#main('IC 194','cropped_i.fits','isophotes.reg','cropped_i_warp_lines.reg','cropped_i_warp_points.reg','N')  


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Creating plot with isophotes")
    parser.add_argument("input_image", help="Input image")
    parser.add_argument("iso_reg", help="Input region file with isophotes")
    parser.add_argument("lines_reg", help="Input region file with lines describing the warps")
    parser.add_argument("points_reg", help="Input region file with points describing the warps")
    parser.add_argument("name", nargs='?', const=0., help="Optional: Input name of the galaxy",type=str,default='')
    parser.add_argument("Or", nargs='?', const=0., help="Optional: Orientation of the image",type=str,default='N')
    parser.add_argument("x0", nargs='?', const=0., help="Optional: Input x-coordinate of the left warp",type=float,default=0.)
    parser.add_argument("y0", nargs='?', const=0., help="Optional: Input y-coordinate of the left warp",type=float,default=0.)    
    parser.add_argument("x1", nargs='?', const=0., help="Optional: Input x-coordinate of the right warp",type=float,default=0.)
    parser.add_argument("y1", nargs='?', const=0., help="Optional: Input y-coordinate of the right warp",type=float,default=0.)     
    
    args = parser.parse_args()
    input_image = args.input_image
    iso_reg = args.iso_reg
    lines_reg = args.lines_reg
    points_reg = args.points_reg
    name = args.name
    Or = args.Or
    x0 = float(args.x0)
    y0 = float(args.y0)
    x1 = float(args.x1)
    y1 = float(args.y1)
    
    main(input_image,iso_reg,lines_reg,points_reg,name,orientation='N',x0=x0,x1=x1,y0=y0,y1=y1)