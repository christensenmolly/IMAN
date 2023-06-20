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
import matplotlib.ticker as ticker
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

import pyfits
from scipy import interpolate
import matplotlib.colors as colors


def crea_hist_image(input_image,output_image,sigma,m0,pix2sec):
	#im = pylab.imread('lena.png').sum(axis=2) # make grayscale
	levels_step = 0.3
	hdulist = pyfits.open(input_image)
	im = hdulist[0].data
	ny,nx = np.shape(im)

	contour(im, origin='lower',cmap='gray_r',linewidths=2)#,vmin=sigma, vmax=3.*sigma)
	axis('equal')
	axis('off')
	#plt.show()
	plt.savefig(output_image, bbox_inches='tight', pad_inches=0)
	plt.clf()
	plt.close()
	'''
	exit()

	
	data = sp.ndimage.zoom(im, 3)


	f = plt.figure(0)
	ax1 = subplot(111)
	x = []
	y = []
	I = []
	I_crit = sigma

	for k in range(ny):
		for i in range(nx):
			if im[k,i]>I_crit:
				x.append(i)
				y.append(k)
				I.append(im[k,i])
	x = np.array(x)
	y = np.array(y)
	I = np.array(I)

	m = -2.5*log10(I/(pix2sec**2)) + m0 
	m = np.around(m,decimals=1)
	m_crit= m0 - 2.5*log10(I_crit/(pix2sec**2))

	con_levels = np.linspace(m_crit-0.5,m_crit,9)

	xi = np.linspace(min(x),max(x),100)
	yi = np.linspace(min(y),max(y),100)

	Ii = griddata(x,y,m,xi,yi)

	#ax1.contour(xi,yi,Ii,con_levels,linewidths=1,colors='black')


	'''


def crea_histogram_image(input_image,output_image,sigma,cmap="grey",scale="histequ", invert=True, mode="zscale",nlevels=9,smooth=4):
          #outp_format = output_image.split('.')[-1]
          lower_limit = str(0.)#0.7
          upper_limit = str(sigma/1000.)
          
          if invert==False:
            p = subprocess.Popen(["ds9",input_image,"-scale","histequ","-scale","mode",mode,"-cmap",cmap,"-regions","system","image","-contour","yes","-contour","scale","log","-contour","limits",lower_limit,upper_limit,"-contour","smooth",str(smooth),"-contour","nlevels",str(nlevels),"-contour","convert","-regions","save","con1.reg","-export",output_image,"-exit"])
	  else:
            p = subprocess.Popen(["ds9",input_image,"-scale","histequ","-scale","limits","0.0",upper_limit,"-cmap",cmap,"-cmap","invert","yes","-regions","system","image","-contour","yes","-contour","scale","log","-contour","limits",lower_limit,upper_limit,"-contour","smooth",str(smooth),"-contour","nlevels",str(nlevels),"-contour","convert","-regions","save","con1.reg","-export",output_image,"-exit"])
	  p.wait()    
 
 	  #subprocess.call("ds9 %s -scale log -cmap b -contour yes -contour scale log -contour limits %.1f %.1f -contour smooth %i -contour nlevels %i -contour convert -regions save con.reg -exit" % (input_image,outer_level,inner_level,smooth,nLevels), shell=True)
 
def main(input_image_g,input_image_r,input_image_i,mask='no',show=False,iso='no'):
	# ds9
	hdulist = pyfits.open(input_image_g)
	inframe = hdulist[0].data
	nx, ny =inframe.shape[1], inframe.shape[0]
	max_g = inframe[int(ny/2),int(nx/2)]
	
	hdulist = pyfits.open(input_image_r)
	inframe = hdulist[0].data
	max_r = inframe[int(ny/2),int(nx/2)]
	
	hdulist = pyfits.open(input_image_i)
	inframe = hdulist[0].data
	max_i = inframe[int(ny/2),int(nx/2)]
	
	line_hor_x1 = 0
	line_hor_x2 = nx
	line_hor_y = ny/2.
	line_vert_y1 = 0
	line_vert_y2 = ny
	line_vert_x = nx/2.
	
	#lines_ds9 = open('cent_lines.reg','w')
	#print >>lines_ds9, "line(%.1f,%.1f,%.1f,%.1f) # line=0 0 color=green width=2" % (line_hor_x1,line_hor_y,line_hor_x2,line_hor_y)
	#print >>lines_ds9, "line(%.1f,%.1f,%.1f,%.1f) # line=0 0 color=green width=2" % (line_vert_x,line_vert_y1,line_vert_x,line_vert_y2)
	#lines_ds9.close()

	zoom = 1.
	command1 = "ds9 -geometry %ix%i -zoom %.2f " % (1.01*nx*zoom,260+ny*zoom,zoom)
	#command1 = command1 + " -rgb -red %s -scale log -scale limits 0 %.1f -invert -green %s -scale log -scale limits 0 %.1f -invert -blue %s -scale log -scale limits 0 %.1f -regions load cent_lines.reg -regions load general_mask.reg " % (input_image_i,max_i,input_image_r,max_r,input_image_g,max_g)
	command1 = command1 + " -rgb -red %s -scale log -scale limits 0 %.1f -invert -green %s -scale log -scale limits 0 %.1f -invert -blue %s -scale log -scale limits 0 %.1f " % (input_image_i,max_i,input_image_r,max_r,input_image_g,max_g)
	if mask!='no':
	  command1 = command1 + "-region load %s " % (mask)
	if iso!='no':
	  command1 = command1 + "-contour yes -contour color red -contour width 2 -contour limits %.1f %.1f -contour smooth 7 -contour nlevels 1 " % (float(iso),float(iso))
	if show==True:
	  command1 = command1 + "-colorbar no -saveimage png %s" % ('rgb.png')
	else:
	  command1 = command1 + "-colorbar no -saveimage png %s -exit" % ('rgb.png')
	subprocess.call(command1, shell=True)
	#subprocess.call("ds9 -rgb -red %s -scale log -scale limits 0 %.1f -invert -green %s -scale log -scale limits 0 %.1f -invert -blue %s -scale log -scale limits 0 %.1f -regions load cent_lines.reg -regions load general_mask.reg" % (input_image_i,max_i,input_image_r,max_r,input_image_g,max_g), shell=True)
	


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Create RGB image from three images")
    parser.add_argument("input_image_g", help="Input B-image")
    parser.add_argument("input_image_r", help="Input G-image")  
    parser.add_argument("input_image_i", help="Input R-image")
    parser.add_argument("--mask", action="store_true", default='no',
                        help="Add region with masks")
    parser.add_argument("--show", action="store_true", default=False,
                        help="Show ds9 window without exit")    
    parser.add_argument("--iso", action="store_true", default='no',
                        help="Plot isophote of given level in [DN]")  
    args = parser.parse_args()
    
    input_image_g = args.input_image_g
    input_image_r = args.input_image_r
    input_image_i = args.input_image_i
    mask = args.mask
    show = args.show
    iso = args.iso
    main(input_image_g,input_image_r,input_image_i,mask=mask,show=show,iso=iso)
