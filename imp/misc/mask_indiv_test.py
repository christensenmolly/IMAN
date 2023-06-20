#!/usr/bin/python
# -*- coding:  cp1251 -*-
# Import the necessary modules
import pyfits
import numpy as np
import math
import itertools
import matplotlib.pyplot as plt
from scipy import ndimage
import sys
from matplotlib.path import Path
from itertools import product
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import imp_setup
import os
import subprocess

try:
    import find_objects
    import calibration
except:
    z=1

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]



Cov_coeff = imp_setup.Cov_coeff

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


def convert_to_pl(mask_image,pl_output_file):
  from pyraf import iraf

  hdulist = pyfits.open(mask_image)
  mask = hdulist[0].data
  ny,nx = np.shape(mask)

  f = open('tmp.txt','w')
  for i in range(ny):
    for k in range(nx):
      if mask[i,k]>0:
	f.write(str(k+1)+' '+str(i+1)+'\n')
  f.close()
  
  iraf.proto(_doprint=0)
  iraf.text2mask.unlearn()
  iraf.text2mask(text='tmp.txt',mask=pl_output_file,ncols=1024,nlines=1024)
  print("Binary .pl file of Bad pixels created in "+pl_output_file)
  print("To view mask in IRAF : display "+mask_image+" 1 overlay="+pl_output_file)  
  os.remove('tmp.txt')











def rot_point(p, orig, angle):
    """ performs rotation of a point (x,y) around point (x0,y0)
        by angle radians"""
    x1 = cos(angle)*(p.x-orig.x) - sin(angle)*(p.y-orig.y) + orig.x
    y1 = sin(angle)*(p.x-orig.x) + cos(angle)*(p.y-orig.y) + orig.y
    return Point(x1, y1)

def listRightIndex(alist, value):
    return len(alist) - alist[-1::-1].index(value) -1

def show_complete(i,N):
  percent = 100. * float(i) / float(N)
  sys.stdout.write("\r%2d%%" % percent)
  sys.stdout.flush()

def bad_pix_mask(input_image, file_segm, xc, yc, output_image=None, output_mask=None):
		hdulist_segm = pyfits.open(file_segm, do_not_scale_image_data=True)
		img_segm = hdulist_segm[0].data


		(dimy,dimx) = img_segm.shape

		shutil.copy(file_segm,output_mask) 
		hdulist3 = pyfits.open(output_mask, do_not_scale_image_data=True, mode='update')
		img3 = hdulist3[0].data

		shutil.copy(input_image,output_image) 
		hdulist4 = pyfits.open(output_image, do_not_scale_image_data=True, mode='update')
		img4 = hdulist4[0].data

		for i in range(dimy):
		   for k in range(dimx):
		      if img_segm[i,k]>0 and img_segm[i,k]!=img_segm[int(yc),int(xc)]:
			img3[i,k] = 1.
			img4[i,k] = 0.
		      else:
			img3[i,k] = 0.  
		hdulist3.flush()
		hdulist4.flush()




def bad_pix_mask_txt(input_image, file_segm, mask_txt, xc, yc, output_image=None, output_mask=None):
		x,y = np.loadtxt(mask_txt, usecols=[0,1],dtype=float, unpack=True)
		
		hdulist_segm = pyfits.open(file_segm, do_not_scale_image_data=True)
		img_segm = hdulist_segm[0].data


		(dimy,dimx) = img_segm.shape

		shutil.copy(file_segm,output_mask) 
		hdulist3 = pyfits.open(output_mask, do_not_scale_image_data=True, mode='update')
		img3 = hdulist3[0].data

		shutil.copy(input_image,output_image) 
		hdulist4 = pyfits.open(output_image, do_not_scale_image_data=True, mode='update')
		img4 = hdulist4[0].data

		for i in range(len(x)):
		      if img_segm[y[i],x[i]]!=img_segm[int(yc),int(xc)]:
			img3[y[i],x[i]] = 1.
			img4[y[i],x[i]] = 0.
  
		hdulist3.flush()
		hdulist4.flush()


# Remove object from region file which lie within the ellipse whihc describes the outer isophote of the galaxy under study:
def remove_galaxy_masks(input_image,FWHM,reg_file, cen, ellA, ellB, ellPA, xSize, ySize, mask_type='circles',reg_file_new = 'non_cov_gal_mask.reg'):
	f = open(reg_file, "r")

        cospa = cos(radians(ellPA))
        sinpa = sin(radians(ellPA))
        # Check if whole ellipse is inside of the image
        # and obtain size of the ellipse in xy plane
        xMax = 0.0
        yMax = 0.0
        xMin = 1e10
        yMin = 1e10
        for e in linspace(0, 4*pi, 1000):
            cose = cos(e)
            sine = sin(e)
            x = cen.x + ellA * cose * cospa - ellB * sine * sinpa
            y = cen.y + ellB * sine * cospa + ellA * cose * sinpa
            if x > xMax:
                xMax = x
            if y > yMax:
                yMax = y
            if x < xMin:
                xMin = x
            if y < yMin:
                yMin = y
        xMin = max(0, int(round(xMin)))
        xMax = min(xSize, int(round(xMax)))
        yMin = max(0, int(round(yMin)))
        yMax = min(ySize, int(round(yMax)))
        focusR = (ellA ** 2.0 - ellB ** 2.0) ** 0.5
        focus10 = Point(cen.x + focusR, cen.y)  # Unrotated
        focus20 = Point(cen.x - focusR, cen.y)  #
        focus1 = rot_point(focus10, cen, radians(ellPA))
        focus2 = rot_point(focus20, cen, radians(ellPA))
        # Find pixels inside of the ellipse
        dEll = 2 * ellA
        
	I = np.zeros(shape=(ySize,xSize))
	for x in xrange(xMin, xMax+1):
	    for y in xrange(yMin, yMax+1):
		  dFocus1 = hypot(x-focus1.x, y-focus1.y)
		  dFocus2 = hypot(x-focus2.x, y-focus2.y)
		  dPoint = dFocus1 + dFocus2
		  if dPoint < dEll:
		    try:
		      I[y,x] = 1
		    except:
		      zz = 1

	
	f_new = open(reg_file_new,'w')
	f_new.write('image\n')
	
	if mask_type=='circles':
	  for line in f:
	    if 'circle(' in line:
	      xc = float(line.split('(')[-1].split(',')[0])
	      yc = float(line.split('(')[-1].split(',')[1])
	      try:
		if I[int(math.ceil(yc)),int(math.ceil(xc))]==0:
		  f_new.write(line)
	      except:
		z=1
	    else:
		f_new.write(line)
	  #f_new.close()
	  
	if mask_type=='all' or mask_type=='clean' or mask_type=='se':
	  for line in f:
	    if 'circle(' in line or 'ellipse(' in line:
	      xc = float(line.split('(')[-1].split(',')[0])
	      yc = float(line.split('(')[-1].split(',')[1])
	      try:
		if I[int(math.ceil(yc)),int(math.ceil(xc))]==0:
		  f_new.write(line)
	      except:
		z=1
	    else:
		f_new.write(line)
	  #f_new.close()
	
	if mask_type=='se':
	  [[xmin,ymin],[xmax,ymax]] = calibration.ellipse_borders([cen.x,cen.y],ellA, ellB, ellPA)
	  '''
	  additional_mask(input_image,xc=None,yc=None,FWHM=FWHM,output_mask='tmp.reg',sigma_numb=8.0,xmin=xmin,ymin=ymin,xmax=xmax,ymax=ymax)  
	  f_tmp = open('tmp.reg', 'r')
	  for line in f_tmp:
	    f_new.write(line)

	  
	  f_tmp.close()
	  os.remove('tmp.reg')
	  '''  
	  f.close()
	  
	  # crop the image with sub-region
	  import imp_center
	  import backEst
	  if xmin<0:	xmin=0
	  if ymin<0:	ymin=0
	  if xmax>xSize:	xmax = xSize
	  if ymax>ySize:	ymax = ySize
	  imp_center.crop(input_image, 'tmp.fits',int(xmin),int(ymin),int(xmax),int(ymax))
	  
	  # Launch Sextractor for this cut out image
	  if os.path.exists('field.cat'):
	    os.rename('field.cat', 'field1.cat')

	  if os.path.exists('segm.fits'):
	    os.rename('segm.fits', 'segm1.fits')
  
	  backEst.call_SE('tmp.fits', both=False, sextr_setup='overlap.sex')
	  #exit()
	  X_IMAGE,Y_IMAGE,A_IMAGE,B_IMAGE,KRON_RADIUS,THETA_IMAGE,CLASS_STAR,FWHM_WORLD = find_objects.find_sex_column('field.cat', ['X_IMAGE','Y_IMAGE','A_IMAGE','B_IMAGE','KRON_RADIUS','THETA_IMAGE','CLASS_STAR','FWHM_WORLD'], float)


	  for k in range(len(X_IMAGE)):
	    if A_IMAGE[k]*KRON_RADIUS[k]<3.*FWHM:# and CLASS_STAR[k]<0.5:# and FWHM_WORLD[k]<1.5*FWHM:# B_IMAGE[k]/A_IMAGE[k]>0.8:#  CLASS_STAR[k]<0.2:
		#ellA = A_IMAGE[k] * KRON_RADIUS[k] * 1.3 + 2.
		#ellB = B_IMAGE[k] * KRON_RADIUS[k] * 1.3 + 2.
		#ellA = 1.5*FWHM
		#ellB = 1.5*FWHM
		ellA = min([A_IMAGE[k],B_IMAGE[k]]) * KRON_RADIUS[k] * 1.3 + 5.
		ellB = ellA

		f_new.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f)\n" % (int(xmin)+X_IMAGE[k], int(ymin)+Y_IMAGE[k],
									    ellA, ellB, THETA_IMAGE[k]))

	  f_new.close()


	  if os.path.exists('field1.cat'):
	    os.rename('field1.cat', 'field.cat')

	  if os.path.exists('segm1.fits'):
	    os.rename('segm1.fits', 'segm.fits')
	  os.remove('tmp.fits')

	  
	  
	  
	  
	  
	else:
	  f_new.close()
	  
	  
# Main function, this part actually runs when routine is called
def ellipse_mask(cen, ellA, ellB, ellPA, inframe, xSize, ySize, img_inp=None,img_mask=None):
        cospa = cos(radians(ellPA))
        sinpa = sin(radians(ellPA))
        # Check if whole ellipse is inside of the image
        # and obtain size of the ellipse in xy plane
        xMax = 0.0
        yMax = 0.0
        xMin = 1e10
        yMin = 1e10
        for e in linspace(0, 4*pi, 1000):
            cose = cos(e)
            sine = sin(e)
            x = cen.x + ellA * cose * cospa - ellB * sine * sinpa
            y = cen.y + ellB * sine * cospa + ellA * cose * sinpa
            if x > xMax:
                xMax = x
            if y > yMax:
                yMax = y
            if x < xMin:
                xMin = x
            if y < yMin:
                yMin = y
        xMin = max(0, int(round(xMin)))
        xMax = min(xSize, int(round(xMax)))
        yMin = max(0, int(round(yMin)))
        yMax = min(ySize, int(round(yMax)))
        focusR = (ellA ** 2.0 - ellB ** 2.0) ** 0.5
        focus10 = Point(cen.x + focusR, cen.y)  # Unrotated
        focus20 = Point(cen.x - focusR, cen.y)  #
        focus1 = rot_point(focus10, cen, radians(ellPA))
        focus2 = rot_point(focus20, cen, radians(ellPA))
        # Find pixels inside of the ellipse
        dEll = 2 * ellA
        
        if img_inp==None:
	  for x in xrange(xMin, xMax+1):
	      for y in xrange(yMin, yMax+1):
		  dFocus1 = hypot(x-focus1.x, y-focus1.y)
		  dFocus2 = hypot(x-focus2.x, y-focus2.y)
		  dPoint = dFocus1 + dFocus2
		  if dPoint < dEll:
		    try:
			  inframe[y-1,x-1] = 0.0
			  img_mask[y-1,x-1] = 1
		    except:
		      zz=1
        else:
	  for x in xrange(xMin, xMax+1):
	      for y in xrange(yMin, yMax+1):
		  dFocus1 = hypot(x-focus1.x, y-focus1.y)
		  dFocus2 = hypot(x-focus2.x, y-focus2.y)
		  dPoint = dFocus1 + dFocus2
		  if dPoint < dEll:
		    try:
			  inframe[y-1,x-1] = img_inp[y-1,x-1]
		    except:
		      zz=1



def mask_array(data, data_fill, reg_file):
		(dimy,dimx) = data.shape
		mask = np.zeros(shape=(dimy,dimx))
		
		paths=[]
                borders = []
		f = open(reg_file, "r")
                factor = 1
                obj_ellipse = None
		if factor==1:
		    for line in f:
			    # Detect circle regions and fill them with the mask value
			    if 'circle(' in line:
				      params = line.split("(")[-1].split(",")[0]
				      cen = Point(float(line.split("(")[-1].split(",")[0]),
						  float(line.split("(")[-1].split(",")[1]))
				      ellA = float(line.split("(")[-1].split(",")[2].split(')')[0])
				      ellB = ellA
				      ellPA = 0.
				      ellipse_mask(cen, ellA, ellB, ellPA, mask, dimx, dimy)
    
			    # Detect polygon regions and add them to the path list
			    elif 'polygon(' in line:
				      param = map(float, ((line.split('polygon(')[1]).split(')')[0]).split(','))
				      param2 = [None]*(len(param)/2)
				      for i in xrange(0,len(param2)):	param2[i] = (int(param[2*i]),int(param[2*i+1])) 
				      param2.append(param2[0])
				      codes = []
				      codes.append(Path.MOVETO)
				      for i in xrange(1,len(param2)-1): codes.append(Path.LINETO)
				      codes.append(Path.CLOSEPOLY)
				      path = Path(param2, codes)
				      paths.append(path)
				      coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
				      X = []; Y = []
				      for k in range(0,len(coords)-1,2):
					      X.append(int(float(coords[k])))
					      Y.append(int(float(coords[k+1])))

				      borders.append([min(X),min(Y),max(X),max(Y)])    
			    elif "ellipse" in line:
				      params = line[8:-2].split(",")
				      cen = Point(float(params[0]),
						  float(params[1]))
				      ellA = float(params[2])
				      ellB = float(params[3])
				      try:
					ellPA = float(params[4])
				      except:
					ellPA = 0.
				      if ellA < ellB:
					  ellA, ellB = ellB, ellA
					  ellPA += 90
				      ellipse_mask(cen, ellA, ellB, ellPA, mask, dimx, dimy)


		if obj_ellipse==None:
		      # Loop over the image and all polygons and fill them with the mask value	
 
		      for kk in range(len(paths)):
			path = paths[kk]
			border = borders[kk]
			xmin = border[0]
			ymin = border[1]
			xmax = border[2]
			ymax = border[3]
			for i, j in product(range(xmin,xmax), range(ymin,ymax)):
			  try:
				mask[j-1][i-1]=1 if path.contains_point((i,j)) else mask[j-1][i-1]
			  except:
			    zz=1


		      for i in range(dimy):
                            for k in range(dimx):
				  if mask[i,k] == 1:
						data[i,k] = data_fill[i,k]
		      return data,mask




















def mask(input_image, reg_file, factor=1, output_image=None, output_mask=None, obj_ellipse=None, mask_value=1):
		print 'Masking contaminants...'
		if output_image==None:
		  output_image = input_image.split('.fits')[0]+'_clean.fits'  

		shutil.copy(input_image,output_image)
		
		if output_mask!=None:
		  shutil.copy(input_image,output_mask)
		  hdulist2 = pyfits.open(output_mask, do_not_scale_image_data=True, mode='update')
		  img2 = hdulist2[0].data
  
		hdulist3 = pyfits.open(output_image, do_not_scale_image_data=True, mode='update')
		img3 = hdulist3[0].data


		(dimy,dimx) = img3.shape
		value = 0.
		paths=[]
		borders=[]
		f = open(reg_file, "r")
		all_lines = f.readlines()
		Number_of_lines = len(all_lines)
		f.close()
		f = open(reg_file, "r")

		coun = 0
		if factor==1:
		    for line in f:
			    #print line
			    # Detect circle regions and fill them with the mask value
			    if 'circle(' in line:
				      '''
				      param = ((line.split('circle(')[1]).split(')')[0]).split(',')
				      a, b ,r  = int(float(param[0])), int(float(param[1])), int(float(param[2])) 
				      y,x = np.ogrid[-b:img3.shape[0]-b, -a:img3.shape[1]-a]
				      mask = x*x + y*y <= r*r
				      img3[mask] = value
				      '''
				      params = line.split("(")[-1].split(",")[0]
				      cen = Point(float(line.split("(")[-1].split(",")[0]),
						  float(line.split("(")[-1].split(",")[1]))
				      ellA = float(line.split("(")[-1].split(",")[2].split(')')[0])
				      ellB = ellA
				      ellPA = 0.
				      ellipse_mask(cen, ellA, ellB, ellPA, img3, dimx, dimy)
    
			    # Detect polygon regions and add them to the path list
			    elif 'polygon(' in line:
				      param = map(float, ((line.split('polygon(')[1]).split(')')[0]).split(','))
				      param2 = [None]*(len(param)/2)
				      for i in xrange(0,len(param2)):	param2[i] = (int(param[2*i]),int(param[2*i+1])) 
				      param2.append(param2[0])
				      codes = []
				      codes.append(Path.MOVETO)
				      for i in xrange(1,len(param2)-1): codes.append(Path.LINETO)
				      codes.append(Path.CLOSEPOLY)
				      path = Path(param2, codes)
				      paths.append(path)
				      coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
				      X = []; Y = []
				      for k in range(0,len(coords)-1,2):
					      X.append(int(float(coords[k])))
					      Y.append(int(float(coords[k+1])))

				      borders.append([min(X),min(Y),max(X),max(Y)])
				    
			    elif "ellipse" in line:
				      params = line[8:-2].split(",")
				      cen = Point(float(params[0]),
						  float(params[1]))
				      ellA = float(params[2])
				      ellB = float(params[3])
				      try:
					ellPA = float(params[4])
				      except:
					ellPA = 0.
				      if ellA < ellB:
					  ellA, ellB = ellB, ellA
					  ellPA += 90
				      ellipse_mask(cen, ellA, ellB, ellPA, img3, dimx, dimy)
			    coun = coun + 1
			    #print ' Done %.1f perc.' % (float(coun)*100./float(Number_of_lines)), line
			    show_complete(coun,Number_of_lines)
		else:
		    for line in f:
			    #print line
			    # Detect circle regions and fill them with the mask value
			    if 'circle(' in line:
				      '''
				      param = ((line.split('circle(')[1]).split(')')[0]).split(',')
				      a, b ,r  = int(float(param[0])/factor), int(float(param[1])/factor), int(float(param[2])/factor) 
				      y,x = np.ogrid[-b:img3.shape[0]-b, -a:img3.shape[1]-a]
				      mask = x*x + y*y <= r*r
				      img3[mask] = value
				      '''
				      params = line.split("(")[-1].split(",")[0]
				      cen = Point(float(line.split("(")[-1].split(",")[0])/factor,
						  float(line.split("(")[-1].split(",")[1])/factor)
				      ellA = float(line.split("(")[-1].split(",")[2].split(')')[0])
				      ellB = ellA
				      ellPA = 0.
				      ellipse_mask(cen, ellA, ellB, ellPA, img3, dimx, dimy)
			    # Detect polygon regions and add them to the path list
			    elif 'polygon(' in line:
				      param = map(float, ((line.split('polygon(')[1]).split(')')[0]).split(','))
				      param2 = [None]*(len(param)/2)
				      for i in xrange(0,len(param2)):	param2[i] = (int(param[2*i]),int(param[2*i+1])) 
				      param2.append(param2[0]/factor)
				      codes = []
				      codes.append(Path.MOVETO)
				      for i in xrange(1,len(param2)-1): codes.append(Path.LINETO)
				      codes.append(Path.CLOSEPOLY)
				      path = Path(param2, codes)
				      paths.append(path)
				      coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
				      X = []; Y = []
				      for k in range(0,len(coords)-1,2):
					      X.append(int(float(coords[k])))
					      Y.append(int(float(coords[k+1])))

				      borders.append([min(X),min(Y),max(X),max(Y)])    
			    elif "ellipse" in line:
				      params = line[8:-2].split(",")
				      cen = Point(float(params[0])/factor,
						  float(params[1])/factor)
				      ellA = float(params[2])/factor
				      ellB = float(params[3])/factor
				      try:
					ellPA = float(params[4])
				      except:
					ellPA = 0.
				      if ellA < ellB:
					  ellA, ellB = ellB, ellA
					  ellPA += 90
				      ellipse_mask(cen, ellA, ellB, ellPA, img3, dimx, dimy)
			    coun = coun + 1
			    #print ' Done %.1f perc.' % (float(coun)*100./float(Number_of_lines)), line  
			    show_complete(coun,Number_of_lines)
		
		if obj_ellipse==None:
		      # Loop over the image and all polygons and fill them with the mask value	
		      nx, ny =img3.shape[1], img3.shape[0]
		      '''
		      for i, j, path in product(range(0,nx), range(0,ny), paths):
			      #print path
			      #exit()
			      img3[j][i]=value if path.contains_point((i,j)) else img3[j][i]
		      hdulist3.flush()
		      '''      
		      for kk in range(len(paths)):
			path = paths[kk]
			border = borders[kk]
			xmin = border[0]
			ymin = border[1]
			xmax = border[2]
			ymax = border[3]
			for i, j in product(range(xmin,xmax), range(ymin,ymax)):
			  try:
				img3[j-1][i-1]=value if path.contains_point((i,j)) else img3[j-1][i-1]
			  except:
			    zz=1
		      hdulist3.flush()

		      if output_mask!=None:
			for i in range(dimy):
				for k in range(dimx):
				  if img3[i,k] == value:
						img2[i,k] = mask_value
				  else:
						img2[i,k] = 0



			hdulist2.flush()
		else:
		      hdulist_inp = pyfits.open(input_image)
		      img_inp = hdulist_inp[0].data
		      
		      [xc,yc,sma,smb,PA] = obj_ellipse
		      cen = Point(xc/factor,yc/factor)

		      ellipse_mask(cen, sma, smb, PA, img3, dimx, dimy, img_inp)

		      # Loop over the image and all polygons and fill them with the mask value	
		      nx, ny =img3.shape[1], img3.shape[0]
		      for kk in range(len(paths)):
			path = paths[kk]
			border = borders[kk]
			xmin = border[0]
			ymin = border[1]
			xmax = border[2]
			ymax = border[3]
			for i, j in product(range(xmin,xmax), range(ymin,ymax)):
			  try:
				img3[j-1][i-1]=value if path.contains_point((i,j)) else img3[j-1][i-1]
			  except:
			    zz=1
		      hdulist3.flush()

		      if output_mask!=None:
			for i in range(dimy):
				for k in range(dimx):
				  if img3[i,k] == value:
						img2[i,k] = mask_value
				  else:
						img2[i,k] = 0



			hdulist2.flush()


def additional_mask(input_image,xc=None,yc=None,FWHM=3.,output_mask='add_mask.reg',sigma_numb=10.0,xmin=None,ymin=None,xmax=None,ymax=None,box_size=5):
    from astropy.stats import sigma_clipped_stats
    from photutils import find_peaks

    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    
    if xmin!=None and ymin!=None and xmax!=None and ymax!=None:
      mean, median, std = sigma_clipped_stats(data[int(math.ceil(ymin)):int(math.ceil(ymax)),int(math.ceil(xmin)):int(math.ceil(xmax))], sigma=3.0)
    else:
      mean, median, std = sigma_clipped_stats(data, sigma=3.0)
      
    threshold = median + (sigma_numb * std)
    
    if xmin==None and ymin==None and xmax==None and ymax==None:
	tbl = find_peaks(data, threshold, box_size=5, subpixel=False)
	xmin = 0; ymin = 0
    else:
	#print 'here'
	tbl = find_peaks(data[int(math.ceil(ymin)):int(math.ceil(ymax)),int(math.ceil(xmin)):int(math.ceil(xmax))], threshold, box_size=box_size, subpixel=False)      
    
    
    f = open(output_mask,'w')
    if xc==None or yc==None:
      for k in range(len(tbl)):
	x = tbl['x_peak'][k]
	y = tbl['y_peak'][k]
	f.write('circle(%f,%f,%f)\n' % (x+1+xmin,y+1+ymin,FWHM))
    else:
      for k in range(len(tbl)):
	x = tbl['x_peak'][k]
	y = tbl['y_peak'][k]
	if math.sqrt((xc-(x+xmin))**2+(yc-(y+ymin))**2)>=FWHM:
	  f.write('circle(%f,%f,%f) # color=red\n' % (x+1+xmin,y+1+ymin,FWHM))      
    f.close()







def final_mask(input_image, segm_file, sex_catalog, output_region_file, xc,yc, output_image=None, output_mask=None, mode='a',rescale=None,rescale_value=None,check_mask=True,create_masks=True, FWHM=3.,sigma_numb=10.0, type_of_masking='non_poly'):

    # Read SExtractor catalogue
    ellA = []
    for line in open(sex_catalog):
        if line.startswith("#"):
            continue
        params = line.split()
        n = int(params[0])
        xCen = float(params[1])
        yCen = float(params[2])
        kron = float(params[8])
        ellA.append(kron * float(params[4]) * Cov_coeff)
        ellB = kron * float(params[5]) * Cov_coeff
        PA = float(params[6])

    if rescale=='offset':
      offset = [rescale_value]
    elif rescale=='factor':
      offset = []
      for k in range(len(ellA)):
	offset.append(ellA[k]*abs(rescale_value-1))
    else:
      offset = None
      

    if os.path.exists('final_mask.txt'):
      os.remove('final_mask.txt')


    if type_of_masking=='poly':
	    '''
	    # 1. Open segmentation file
	    hdulist = pyfits.open(segm_file)
	    data = hdulist[0].data
	    ny,nx = np.shape(data)
	    
	    from skimage import measure
	    from skimage.measure import find_contours, approximate_polygon, \
	    subdivide_polygon
	    import pyclipper
	  
	    # Plot contours for the segmentation map
	    contours = measure.find_contours(data, 0)


	    # Save it to region file
	    f_tmp = open('tmp.reg','w')

	    for n, contour_old in enumerate(contours):
	      contour = approximate_polygon(contour_old, tolerance=.9) #### Maybe, this step can be omitted
	      x = contour[:, 1]
	      y = contour[:, 0]
	      if x[0]>x[-1]:
		x = list(reversed(x))
		y = list(reversed(y))

	      for k in range(len(x)):
		if k==0:
		  f_tmp.write('polygon(%.1f,%.1f,' % (x[k]+1.,y[k]+1.) )
		elif k>0 and k<len(x)-1:
		  f_tmp.write('%.1f,%.1f,' % (x[k]+1.,y[k]+1.) )
		else:
		  f_tmp.write('%.1f,%.1f)\n' % (x[k]+1.,y[k]+1.) )
	    f_tmp.close()
	    
	    # If offseting should be done:
	    if offset!=None:
	      f_tmp = open('tmp.reg','r')
	      f_reg = open(output_region_file,'w')
	      count = -1
	      for Line in f_tmp:
			count = count + 1
			if len(offset)==1:
			  Offset = offset[0]
			else:
			  Offset = offset[count]
			X=[]; Y=[]; COORDS=[]; XX=[]; YY=[]
			if 'polygon' in Line:
			  coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]

			  for k in range(0,len(coords)-1,2):
				  COORDS.append([int(float(coords[k])), int(float(coords[k+1]))])
				  X.append(int(float(coords[k])))
				  Y.append(int(float(coords[k+1])))
			  pco = pyclipper.PyclipperOffset()
			  pco.AddPath(COORDS, pyclipper.JT_ROUND, pyclipper.ET_CLOSEDPOLYGON)
			  solution = pco.Execute(Offset)
			  polygon = []
			  try:    
			    for i in range(len(solution[0])):
				      xx = solution[0][i][0]
				      yy = solution[0][i][1]
				      XX.append(xx)
				      YY.append(yy)
				      polygon.append( str(int(round(xx))) + ',' + str(int(round(yy))))
			    #print xc,yc,min(XX)
			    if xc>=min(XX) and xc<=max(XX) and yc>=min(YY) and yc<=max(YY):
			      #print 'here'
			      continue
			  except:
			      z=1
			      
			  for k in range(len(polygon)):
			      if k==0:
				f_reg.write('polygon(%s,' % (polygon[k]) )
			      elif k>0 and k<len(polygon)-1:
				f_reg.write('%s,' % (polygon[k]) )
			      else:
				f_reg.write('%s)\n' % (polygon[k]) )    
	      f_reg.close()
	      f_tmp.close()
	    else:
	      os.rename('tmp.reg',output_region_file)
	    '''
	    import convert_segm_to_region
	    convert_segm_to_region.main(segm_file, output_region_file, output_mask_image=None, fits_slice = 0, offset=offset, xc=xc, yc=yc)
    else:
	    X_IMAGE,Y_IMAGE,A_IMAGE,B_IMAGE,KRON_RADIUS,THETA_IMAGE,CLASS_STAR = find_objects.find_sex_column(sex_catalog, ['X_IMAGE','Y_IMAGE','A_IMAGE','B_IMAGE','KRON_RADIUS','THETA_IMAGE','CLASS_STAR'], float)

	    f_reg = open(output_region_file,'w')
	    #print output_region_file
	    for k in range(len(X_IMAGE)):
		ellA = A_IMAGE[k] * KRON_RADIUS[k] * Cov_coeff + rescale_value
		ellB = B_IMAGE[k] * KRON_RADIUS[k] * Cov_coeff + rescale_value
		f_reg.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f)\n" % (X_IMAGE[k], Y_IMAGE[k],
									    ellA, ellB, THETA_IMAGE[k]))
		#print X_IMAGE[k], Y_IMAGE[k], ellA, ellB, THETA_IMAGE[k]
		#exit()
	    f_reg.close()

    additional_mask(input_image,xc=xc,yc=yc,FWHM=FWHM,output_mask='add_mask.reg',sigma_numb=sigma_numb)
    

    # Merge masks:
    filenames = [output_region_file, 'add_mask.reg']
    with open('tmp1.reg', 'w') as outfile:
	for fname in filenames:
	    with open(fname) as infile:
		outfile.write(infile.read())
    os.rename('tmp1.reg',output_region_file)
    os.remove('add_mask.reg')


    
    if check_mask==True:
	      p = subprocess.Popen(["ds9",input_image,"-scale","histequ","-invert","-regions","load",output_region_file])
	      p.wait()    

    if create_masks==True:
      mask(input_image, output_region_file, factor=1, output_image=output_image, output_mask=output_mask)








if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Isophote analysis")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("reg_file", help="Input region file")
    parser.add_argument("output_image", nargs='?', const=1, help="Optional: Ouput image name",type=str, default=None) 
    parser.add_argument("output_mask", nargs='?', const=1, help="Optional: Ouput mask name",type=str, default=None) 
    args = parser.parse_args()

    input_image = args.inputImage
    reg_file = args.regFile
    output_image = args.output_image
    output_mask = args.output_mask    

    mask(input_image, reg_file, output_image=None, output_mask=None)
