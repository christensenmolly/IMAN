#!/usr/bin/python
# DESCRIPTION:
# Script to mask stars inside of a galaxy image using a given PSF image.
# The target galaxy can be unmasked if needed.
# The masked areas can be presented as ellipses or polygons.
# The masked areas can be enlarged (by multiplication or subtraction).
# MINIMAL USAGE: python inner_masking.py [input_image] [psf_image]
        
from astropy.io import fits as pyfits
from scipy import ndimage as ndi
import matplotlib.pyplot as plt
from skimage.feature import peak_local_max
from skimage import data, img_as_float
from skimage import measure
import numpy as np
import math
from astropy.stats import sigma_clipped_stats
import os
import sys
import shutil
import shapely
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import warnings
warnings.filterwarnings("ignore")
#PATH_TO_MODULE = os.path.dirname(__file__)
#sys.path.append(PATH_TO_MODULE.split('/IMAN')[0]+'/IMAN/imp/create_psf')
#sys.path.append(PATH_TO_MODULE.split('/IMAN')[0]+'/IMAN/imp/masking')
#sys.path.append(PATH_TO_MODULE.split('/IMAN')[0]+'/IMAN/imp/cropping')
import convert_segm_to_region
import argparse



def read_region(reg_file):
    f = open(reg_file, "r")    
    for line in f:
        if "ellipse" in line:
            params = line.split(",")
            cen = [float(params[0].split('(')[1]),float(params[1])]
            ellA = float(params[2])
            ellB = float(params[3])
            ellPA = float(params[4].split(')')[0])
            if ellA < ellB:
                ellA, ellB = ellB, ellA
                ellPA += 90    
            break
    f.close()
    return cen[0],cen[1],ellA,ellB,ellPA

def prepare_psf(psf_image):
    # Upload the PSF image:
    hdulist_psf = pyfits.open(psf_image)
    data_psf = hdulist_psf[0].data
    Imax_psf = np.max(data_psf)
    (ycen_psf, xcen_psf) = np.where(data_psf==Imax_psf)

    return data_psf,Imax_psf,ycen_psf,xcen_psf

def detect_stars(input_image, psf_image, output_region='test_mask.reg', median_sky=None, std_sky=None, sigma_numb=5., min_distance=5, offset_size=1.,offset_pix=0., min_pix=5, galaxy_ellipse=None):
  print('Masking by searching for local maxima...')
  if galaxy_ellipse is not None:
        try:
            if galaxy_ellipse.geom_type == 'Polygon':
                cont = True
                ellipse = galaxy_ellipse
                SMB = None; PA = None
                [(xc_gal,yc_gal)] = galaxy_ellipse.centroid.coords
            else:
                cont = False
                print('ERROR: Some error with the galaxy polygon provided! Exiting.')
                exit()
        except:
            try:
                xc_gal,yc_gal,SMA,SMB,PA = read_region(galaxy_ellipse)
            except:
                [[xc_gal,yc_gal],SMA,SMB,PA] = galaxy_ellipse
                
            circle = Point(xc_gal,yc_gal).buffer(1)
            if SMA is not None:
                if SMB is not None and PA is not None:
                    ellipse = shapely.affinity.scale(circle,SMA,SMB)
                    rot_ellipse = shapely.affinity.rotate(ellipse, PA, origin='center', use_radians=False)
                else:
                    ellipse = shapely.affinity.scale(circle,SMA,SMA)
  else:
      xc_gal=None; yc_gal=None
      
  hdu = pyfits.open(input_image)
  data = hdu[0].data
  header = hdu[0].header
  ny, nx = data.shape

  if xc_gal is None or yc_gal is None:
        xc_gal = nx/2. # was 1000000
        yc_gal = ny/2. # was 1000000

  if median_sky is None:
      median_sky = header['SKY_MED']
      
  if std_sky is None:
      std_sky = header['SKY_STD']    
      
  f_tmp = open(output_region,'w')
  f_tmp.write('%s\n' % ('image') )  
  
  data_psf, Imax_psf, ycen_psf, xcen_psf = prepare_psf(psf_image)
  
  # Detect stars as peaks
  threshold = median_sky + (sigma_numb * std_sky)
  coordinates = peak_local_max(data, min_distance=min_distance)

  for k in range(len(coordinates)):
      yc,xc = coordinates[k][0],coordinates[k][1]
      if galaxy_ellipse is not None:
        if not rot_ellipse.contains(Point(xc,yc)):
            continue
      Imax = data[yc,xc]
      if Imax>threshold:
            if math.sqrt((xc_gal-(xc+1.))**2+(yc_gal-(yc+1.))**2)<=10.:	#### WAS 5.
              # Don't mask the centre of the galaxy!
              continue
      else:
          continue


      # Find background around the star:
      Rin = int(math.ceil(max([xcen_psf,ycen_psf])))
      Rout = int(math.ceil(max([xcen_psf,ycen_psf])))+5
      backgr = []
      X = []; Y = []; I = []
      for y in range(int(yc)-Rout,int(yc)+Rout,1):
        for x in range(int(xc)-Rout,int(xc)+Rout,1):
            try:
                I.append(data[y,x])
                X.append(x)
                Y.append(y)
                if (x-int(xc))**2+(y-int(yc))**2>Rin**2 and (x-int(xc))**2+(y-int(yc))**2<=Rout**2 and x>=0 and y>=0 and x<nx and y<ny:
                    backgr.append(data[y,x])
            except:
                zz = 1
      mean_bckg, median_bckg, std_bckg = sigma_clipped_stats(backgr, sigma=3.0, iters=5)
      Imax = Imax - median_bckg
      psf = data_psf * Imax / Imax_psf

      outer_contour(psf, threshold, xcen_psf, ycen_psf, xc, yc, f_tmp, min_pix)          
  
  f_tmp.close()
  if offset_size!=1. or offset_pix!=0.:
    convert_segm_to_region.do_offseting(output_region, 'tmp.reg', offset_size=offset_size, offset_pix=offset_pix, xc=None, yc=None, system='image')
    shutil.move('tmp.reg', output_region)
  print('Done!')

def outer_contour(old_image_data, level, xc_psf, yc_psf, xc, yc, f_tmp, min_pix): 
  nx, ny =old_image_data.shape[1], old_image_data.shape[0]
  image_data = np.zeros(shape=(ny,nx))
  np.putmask(image_data, old_image_data>=level, 1.)

  contours = measure.find_contours(image_data, 0)


  for n, contour_old in enumerate(contours):
              contour = contour_old
              x = contour[:, 1]
              y = contour[:, 0]

              XX = x - xc_psf

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

              if max(x)-min(x)>min_pix and max(y)-min(y)>min_pix:
                  for k in range(len(x)):
                    if k==0:
                      f_tmp.write('polygon(%.1f,%.1f,' % (x[k]+1.-xc_psf+xc,y[k]+1.-yc_psf+yc) )
                    elif k>0 and k<len(x)-1:
                      f_tmp.write('%.1f,%.1f,' % (x[k]+1.-xc_psf+xc,y[k]+1.-yc_psf+yc) )
                    else:
                      f_tmp.write('%.1f,%.1f)\n' % (x[k]+1-xc_psf+xc,y[k]+1.-yc_psf+yc) )




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Masking by searching for local maxima")
    parser.add_argument("inputImage", help="Input fits image")
    parser.add_argument("psf_image", help="Input PSF fits image")
    parser.add_argument("--std_sky", nargs='?', const=1, help="Optional: Background rms (should be given)",type=float,default=None)  

    parser.add_argument("--outputRegion", help="Optional: Output region file", type=str, default='mask_lm.reg') 
    parser.add_argument("--med_sky", nargs='?', const=1, help="Optional: Median background level",type=float,default=0.)
    parser.add_argument("--snr", nargs='?', const=1, help="Optional: Input signal-to-noise ratio of pixels to mask",type=float,default=3.)
    parser.add_argument("--min_pix", nargs='?', const=1, help="Optional: Number of joint pixels along each x- and y-axis",type=int,default=2)
    parser.add_argument("--min_dist", nargs='?', const=1, help="Optional: Minimal distance in pixels between local maxima",type=int,default=5)
    parser.add_argument("--offset_pix", nargs='?', const=1, help="Optional: Offset to make the regions larger or smaller, in pixels. Default 0..",type=float, default=0.) 
    parser.add_argument("--offset_size", nargs='?', const=1, help="Optional: Offset to make the regions larger or smaller, in units of the region size (multiplication). Default 1.",type=float, default=1.) 
    parser.add_argument("--galaxy_ellipse", help="Optional: Galaxy ellipse (list of values or region)", type=str, default=None) 
    
    args = parser.parse_args()

    input_image = args.inputImage
    psf_image = args.psf_image
    output_region = args.outputRegion
    median_sky = args.med_sky
    std_sky = args.std_sky
    sigma_numb = args.snr
    min_pix = args.min_pix
    min_distance = args.min_dist
    offset_pix = args.offset_pix
    offset_size = args.offset_size
    galaxy_ellipse = args.galaxy_ellipse
    
    detect_stars(input_image, psf_image, output_region=output_region, median_sky=median_sky, std_sky=std_sky, sigma_numb=sigma_numb, min_distance=min_distance, offset_size=offset_size,offset_pix=offset_pix, min_pix=min_pix, galaxy_ellipse=galaxy_ellipse)
    
