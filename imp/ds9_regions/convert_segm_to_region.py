#!/usr/bin/python
# Script to convert a segmentation map from Sextractor or PTS to the DS9 region file consisting of polygons 
# EXAMPLE: python convert_segm_to_region.py NGC0150_mask.fits test.reg --ignore_value 60000 --scale 2

# Import the necessary modules:
import pyfits
import numpy as np
import math
import argparse
import os
from skimage import measure
from skimage.measure import find_contours, approximate_polygon, subdivide_polygon
import pyclipper
import matplotlib.pyplot as plt
from shapely.geometry.polygon import Polygon
from shapely.affinity import scale as Scale


def scale_region(input_region_file, output_region_file, scale_value, max_radius=None):
	      f_inp = open(input_region_file,'r')
	      f_out = open(output_region_file,'w')
	      
	      f_out.write('image\n')
	      for Line in f_inp:
			if 'polygon' in Line:
			  coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]
			  pol = []
			  for k in range(0,len(coords)-1,2):
			    pol.append((float(coords[k]),float(coords[k+1])))
			  try:
			      polygon = Polygon(pol)  
			      polygon_a = Scale(polygon, xfact=scale_value, yfact=scale_value)
			      x,y = polygon_a.exterior.xy
				  
			      for k in range(len(x)):
				  if k==0:
				    f_out.write('polygon(%.1f,%.1f,' % (x[k],y[k]) )
				  elif k>0 and k<len(x)-1:
				    f_out.write('%.1f,%.1f,' % (x[k],y[k]) )
				  else:
				    f_out.write('%.1f,%.1f)\n' % (x[k],y[k]) )    
			  except:
			    z=1
			if 'circle(' in Line:
				      params = Line.split("(")[-1].split(",")[0]
				      cen = [float(Line.split("(")[-1].split(",")[0]),float(Line.split("(")[-1].split(",")[1])]
				      ellA = float(Line.split("(")[-1].split(",")[2].split(')')[0])
				      ellB = ellA
				      ellPA = 0.
				      if max_radius==None:
					f_out.write('circle(%f,%f,%f)\n' % (cen[0],cen[1],ellA*scale_value))
				      else:
					if ellA>max_radius:
					  f_out.write('circle(%f,%f,%f)\n' % (cen[0],cen[1],ellA))
					else:
					  f_out.write('circle(%f,%f,%f)\n' % (cen[0],cen[1],ellA*scale_value))
			if 'ellipse(' in Line:
				      params = Line.split(",")
				      cen = [float(params[0].split('(')[1]),float(params[1])]
				      ellA = float(params[2])
				      ellB = float(params[3])
				      ellPA = float(params[4].split(')')[0])
				      if ellA < ellB:
					    ellA, ellB = ellB, ellA
					    ellPA += 90
				      if max_radius==None:
					f_out.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f)\n" % (cen[0],cen[1],
									    ellA*scale_value, ellB*scale_value, ellPA))
				      else:
					if ellA>max_radius:
					  f_out.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f)\n" % (cen[0],cen[1],
									    ellA, ellB, ellPA))
					else:
					  f_out.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f)\n" % (cen[0],cen[1],
									    ellA*scale_value, ellB*scale_value, ellPA))
	      f_out.close()
	      f_inp.close()  
  
    
def main(segm_file, output_region_file, fits_slice = 0, scale=None, offset=None, ignore_value=None):
	    # Open the segmentation file
	    hdulist = pyfits.open(segm_file)
	    Data = hdulist[fits_slice].data
	    
	    # Check it there are some pixels with values larger than ignore_value. They will be skipped.
	    if ignore_value!=None:
	      print('Ignoring masks with values larger than %s!' % (ignore_value))
	      np.putmask(Data, Data>=float(ignore_value), 0.)
	      
	    ny,nx = np.shape(Data)
	    
	    # Create a bit larger image to handle polygons on the edges:
	    data = np.zeros(shape=(ny+2,nx+2))
	    data[1:ny+1,1:nx+1] = Data
	    
	    # Draw contours for the segmentation map:
	    contours = measure.find_contours(data, 0)
	    
	    '''
	    # Plot contours:
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
	    
	    # Save them in a region file
	    f_tmp = open('tmp.reg','w')
	    f_tmp.write('%s\n' % ('image') )   
	    for n, contour in enumerate(contours):
	      x = contour[:, 1] -1.
	      y = contour[:, 0] -1.

	      if x[0]>x[-1]:
		x = list(reversed(x))
		y = list(reversed(y))
	      else:
		x=list(x)
		y=list(y)

	      for k in range(len(x)):
		if k==0:
		  f_tmp.write('polygon(%.1f,%.1f,' % (x[k]+1.,y[k]+1.) )
		elif k>0 and k<len(x)-1:
		  f_tmp.write('%.1f,%.1f,' % (x[k]+1.,y[k]+1.) )
		else:
		  f_tmp.write('%.1f,%.1f)\n' % (x[k]+1.,y[k]+1.) )
	    f_tmp.close()

	    
	    if scale!=None:
		  # If scaling should be done:
		  f_tmp = open('tmp.reg','r')
		  f_reg = open(output_region_file,'w')
		  f_reg.write('%s\n' % ('image') )
		  for Line in f_tmp:
			    X=[]; Y=[]; COORDS=[]; XX=[]; YY=[]
			    if 'polygon' in Line:
			      coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]
			      pol = []
			      for k in range(0,len(coords)-1,2):
				pol.append((int(float(coords[k])),int(float(coords[k+1]))))
			      try:
				  polygon = Polygon(pol)  
				  polygon_a = Scale(polygon, xfact=float(scale), yfact=float(scale))
				  x,y = polygon_a.exterior.xy
				      
				  for k in range(len(x)):
				      if k==0:
					f_reg.write('polygon(%.1f,%.1f,' % (x[k],y[k]) )
				      elif k>0 and k<len(x)-1:
					f_reg.write('%.1f,%.1f,' % (x[k],y[k]) )
				      else:
					f_reg.write('%.1f,%.1f)\n' % (x[k],y[k]) )    
			      except:
				z=1
		  f_reg.close()
		  f_tmp.close()
		  os.remove('tmp.reg')
	    
	    elif offset!=None:
		  # If inward/outward polygon offseting should be done:
		  f_tmp = open('tmp.reg','r')
		  f_reg = open(output_region_file,'w')
		  f_reg.write('%s\n' % ('image') )
		  for Line in f_tmp:
			    X=[]; Y=[]; COORDS=[]; XX=[]; YY=[]
			    if 'polygon' in Line:
			      coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]

			      for k in range(0,len(coords)-1,2):
				      COORDS.append([int(float(coords[k])), int(float(coords[k+1]))])
				      X.append(int(float(coords[k])))
				      Y.append(int(float(coords[k+1])))
			      pco = pyclipper.PyclipperOffset()
			      pco.AddPath(COORDS, pyclipper.JT_ROUND, pyclipper.ET_CLOSEDPOLYGON)
			      solution = pco.Execute(float(offset))
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
	      
	    else:      
	      os.rename('tmp.reg',output_region_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convertion from segmentation map to region file")
    parser.add_argument("segm_file", help="Input segmentation map (fits)")
    parser.add_argument("output_region_file", help="Output region file")
    parser.add_argument("--fits_slice", nargs='?', const=0, help="Optional: Fits slice in the segmentation map to be used. Default 0.",type=int, default=0) 
    parser.add_argument("--scale", nargs='?', const=1, help="Optional: Scaling factor to make the regions larger or smaller. Default None.",type=str, default=None) 
    parser.add_argument("--offset", nargs='?', const=1, help="Optional: Offset to make the regions larger or smaller, in pixels. Default None.",type=str, default=None)
    parser.add_argument("--ignore_value", nargs='?', const=1, help="Optional: Ignore pixels with values larger than this.",type=str, default=None)
    args = parser.parse_args()

    segm_file = args.segm_file
    output_region_file = args.output_region_file
    fits_slice = args.fits_slice
    scale = args.scale
    offset = args.offset
    ignore_value = args.ignore_value

    main(segm_file, output_region_file, fits_slice = fits_slice, scale=scale, offset=offset, ignore_value=ignore_value)

#scale_region('inner_mask_proto.reg', 'inner_mask_100.reg', 3.)