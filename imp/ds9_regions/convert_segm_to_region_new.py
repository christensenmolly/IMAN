#!/usr/bin/python
# Script to convert a segmentation map from Sextractor or PTS to the DS9 region file consisting of polygons 
# EXAMPLE: python ~/MEGA/MyPrograms/IMAN/IMP_NEW/convert_segm_to_region.py segments.fits mask.reg --fits_slice=0 --xc=120 --yc=123 --offset=3

# Import the necessary modules
import pyfits
import numpy as np
import math
from scipy import ndimage
import sys
import shutil
import argparse
import os
from skimage import measure
from skimage.measure import find_contours, approximate_polygon, subdivide_polygon
import pyclipper
from astropy import wcs
from skimage.draw import ellipse
from skimage.measure import label, regionprops
from skimage.morphology import watershed
from skimage.feature import peak_local_max
import matplotlib.pyplot as plt

def unmask_galaxy(data, xc, yc):
      '''
      Function to remove from the mask the area which overlaps the galaxy with the coordinates (xc,yc)
      '''
      xc = float(xc); yc = float(yc)
      object_number = data[int(yc), int(xc)]

      ny,nx = np.shape(data)
      for k in range(ny):
	for i in range(nx):
	  if data[k,i]==object_number:
	    data[k,i] = 0
      return data


def unmask_galaxy_test(Data, xc, yc):
      #hdulist = pyfits.open('segm.fits') # open FITS file
      #Data = hdulist[0].data
      
      ny,nx = np.shape(Data)
      image = np.copy(Data)
      for k in range(ny):
	for i in range(nx):
	  if Data[k,i]>0:
	    image[k,i] = True
	  else:
	    image[k,i] = False 
      
      #from photutils import data_properties, properties_table
      #props = data_properties(Data)
      #columns = ['id', 'xcentroid', 'ycentroid', 'semimajor_axis_sigma','semiminor_axis_sigma', 'orientation']
      #tbl = properties_table(props, columns=columns)
      #print float(tbl['xcentroid'][0])
      #exit()
      # Now we want to separate the two objects in image
      # Generate the markers as local maxima of the distance to the background
      distance = ndimage.distance_transform_edt(image)
      local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((3, 3)),
				  labels=image)
      markers = ndimage.label(local_maxi)[0]
      labels = watershed(-distance, markers, mask=image)

      regions = regionprops(labels)

      #fig, ax = plt.subplots()
      #ax.imshow(labels, cmap=plt.cm.spectral, interpolation='nearest')
      #plt.show()
      #exit()

      object_labels = []
      for kk in range(len(regions)):
	  props = regions[kk]
	  y0, x0 = props.centroid
	  orientation = props.orientation
	  x1 = x0 + math.cos(orientation) * 0.5 * props.major_axis_length
	  y1 = y0 - math.sin(orientation) * 0.5 * props.major_axis_length
	  x2 = x0 - math.sin(orientation) * 0.5 * props.minor_axis_length
	  y2 = y0 - math.cos(orientation) * 0.5 * props.minor_axis_length

	  minr, minc, maxr, maxc = props.bbox
	  bx = (minc, maxc, maxc, minc, minc)
	  by = (minr, minr, maxr, maxr, minr)
	  if xc>=minc and xc<=maxc and yc>=minr and yc<=maxr:
	    #ax.plot(bx, by, '-b', linewidth=2.5)
	    #ax.plot((x0, x1), (y0, y1), '-r', linewidth=2.5)
	    #ax.plot((x0, x2), (y0, y2), '-r', linewidth=2.5)
	    #ax.plot(x0, y0, '.g', markersize=15)
	    object_labels.append(labels[int(y0),int(x0)])
	    
      #ax.axis((0, 600, 600, 0))
      #plt.show()
      for k in range(ny):
	for i in range(nx):
	  if labels[k,i] in object_labels:
	    labels[k,i]=0

      #exit()


      #hdu = pyfits.PrimaryHDU(labels)
      #hdu.writeto('labels.fits',clobber=True)
      #exit()
      return labels


def do_offseting(input_region_file, output_region_file, offset, xc=None, yc=None, system='image'):
	      f_tmp = open(input_region_file,'r')
	      f_reg = open(output_region_file,'w')
	      f_reg.write('%s\n' % (system) )
	      for Line in f_tmp:
			X=[]; Y=[]; COORDS=[]; XX=[]; YY=[]
			if 'polygon' in Line:
			  coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]

			  for k in range(0,len(coords)-1,2):
				  COORDS.append([int(float(coords[k])), int(float(coords[k+1]))])
				  X.append(int(float(coords[k])))
				  Y.append(int(float(coords[k+1])))
			  if len(offset)==1:
			    Offset = float(offset[0])
			  else:
			    Offset = max([abs(max(X)-min(X)), abs(max(Y)-min(Y))]) * abs(float(offset[0])-1)
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

			    if xc!=None and yc!=None:
			      if float(xc)>=min(XX) and float(xc)<=max(XX) and float(yc)>=min(YY) and float(yc)<=max(YY):
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
	      

    
def main(segm_file, output_region_file, output_mask_image=None, fits_slice = 0, offset=None, xc=None, yc=None, system='image', ignore_value=None):
	    # Open the segmentation file
	    hdulist = pyfits.open(segm_file)
	    data = hdulist[0].data
	    if ignore_value!=None:
	      print 'Ignoring masks with values larger than',ignore_value,'!'
	      np.putmask(data, data>=float(ignore_value), 0.)
	      
	    if len(np.shape(data))==2 and xc!=None and yc!=None:
	      data = unmask_galaxy(data, xc, yc)
	      #data = unmask_galaxy_test(data, xc, yc)
	    if fits_slice!='all':
	      if len(np.shape(data))==3:
		data = data[int(fits_slice)]
		if xc!=None and yc!=None:
		  data = unmask_galaxy(data, xc, yc)
		  #data = unmask_galaxy_test(data, xc, yc)
	    else:
	      if len(np.shape(data))==3:
		(n_slices,ySize,xSize) = np.shape(data)
		Data = np.zeros(shape=(ySize,xSize))
		for k in range(n_slices):
		  if xc!=None and yc!=None:
		    data_unmask = unmask_galaxy(data[k], xc, yc)
		    #data_unmask = unmask_galaxy_test(data[k], xc, yc)
		    #print data_unmask[int(yc),int(xc)]
		    data[k] = data_unmask
		  Data = Data + data[k]
		data = Data
	    if output_mask_image!=None:
	      hdu = pyfits.PrimaryHDU(data)
	      hdu.writeto(output_mask_image, clobber=True)

	    ny,nx = np.shape(data)
	    
	    '''
	    # Add a frame
	    data[0,0:nx]=0.
	    data[ny-1,0:nx]=0.

	    data[0:ny,0]=0.
	    data[0:ny,nx-1]=0.    
	    '''
	    # Make contours for the segmentation map:
	    contours = measure.find_contours(data, 0)
	    
	    '''
	    # Plot contours:
	    import matplotlib.pyplot as plt
	    fig, ax = plt.subplots()
	    ax.imshow(data, interpolation='nearest', cmap=plt.cm.gray)

	    for n, contour in enumerate(contours):
		ax.plot(contour[:, 1], contour[:, 0], linewidth=2, color='r')

	    ax.axis('image')
	    ax.set_xticks([])
	    ax.set_yticks([])
	    plt.show()
	    exit()
	    '''
	    # Save them in the region file
	    f_tmp = open('tmp.reg','w')
	    f_tmp.write('%s\n' % (system) )   
	    for n, contour_old in enumerate(contours):
	      contour = contour_old#approximate_polygon(contour_old, tolerance=.9) #### Maybe, this step can be omitted
	      x = contour[:, 1]
	      y = contour[:, 0]

	      if x[0]>x[-1]:
		x = list(reversed(x))
		y = list(reversed(y))
	      else:
		x=list(x)
		y=list(y)
	      
	      # Right top corner
	      if (x[0]>=nx-1 and y[-1]>=ny-1) or (x[-1]>=nx-1 and y[0]>=ny-1):
		x.append(float(nx-1))
		y.append(float(ny-1))

	      # Right bottom corner
	      if (x[0]>=nx-1 and y[-1]<=0) or (x[-1]>=nx-1 and y[0]<=0):
		x.append(float(nx-1))
		y.append(float(0))

	      # Left top corner
	      if (x[0]<=0 and y[-1]>=ny-1) or (x[-1]<=0 and y[0]>=ny-1):
		x.append(float(0))
		y.append(float(ny-1))

	      # Left bottom corner
	      if (x[0]<=0 and y[-1]<=0) or (x[-1]<=0 and y[0]<=0):
		x.append(float(0))
		y.append(float(0))

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
              do_offseting('tmp.reg', output_region_file, offset, xc=xc, yc=yc)
              os.remove('tmp.reg')
              ''' 
	      f_tmp = open('tmp.reg','r')
	      f_reg = open(output_region_file,'w')
	      f_reg.write('%s\n' % (system) )
	      for Line in f_tmp:
			X=[]; Y=[]; COORDS=[]; XX=[]; YY=[]
			if 'polygon' in Line:
			  coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]

			  for k in range(0,len(coords)-1,2):
				  COORDS.append([int(float(coords[k])), int(float(coords[k+1]))])
				  X.append(int(float(coords[k])))
				  Y.append(int(float(coords[k+1])))
			  if len(offset)==1:
			    Offset = float(offset[0])
			  else:
			    Offset = max([abs(max(X)-min(X)), abs(max(Y)-min(Y))]) * abs(float(offset[0])-1)
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

			    if xc!=None and yc!=None:
			      if float(xc)>=min(XX) and float(xc)<=max(XX) and float(yc)>=min(YY) and float(yc)<=max(YY):
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
	      os.remove('tmp.reg')
              '''
	    else:
	      os.rename('tmp.reg',output_region_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convertion from segmentation map to region file")
    parser.add_argument("segm_file", help="Input segmentation map (fits)")
    parser.add_argument("output_region_file", help="Output region file")
    parser.add_argument("--output_mask_image", nargs='?', const=1, help="Optional: Output mask image. Default None.",type=str, default=None) 
    parser.add_argument("--fits_slice", nargs='?', const=1, help="Optional: Fits slice in the segmentation map to be used. Default 0.",type=str, default='0') 
    parser.add_argument("--offset", nargs='?', const=1, help="Optional: Offset to make the regions larger or smaller, in pixels. Default None.",type=str, default=None) 
    parser.add_argument("--xc", nargs='?', const=1, help="Optional: If you don't want to mask the galaxy with the given x-coordinate of the center",type=str, default=None) 
    parser.add_argument("--yc", nargs='?', const=1, help="Optional: If you don't want to mask the galaxy with the given y-coordinate of the center",type=str, default=None)
   
    
    args = parser.parse_args()

    segm_file = args.segm_file
    output_region_file = args.output_region_file
    output_mask_image = args.output_mask_image
    fits_slice = args.fits_slice
    offset = args.offset
    xc = args.xc
    yc = args.yc
    

    main(segm_file, output_region_file, output_mask_image=output_mask_image, fits_slice = fits_slice, offset=offset, xc=xc, yc=yc)

