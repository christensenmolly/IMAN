#!/usr/bin/python        
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
#import shapely
#from shapely.geometry import Point
#from shapely.geometry.polygon import Polygon
import warnings
warnings.filterwarnings("ignore")
from scipy.signal import argrelextrema
from scipy.signal import find_peaks
import subprocess
import argparse


def reject_outliers(data, m=2):
    return np.where(abs(data - np.mean(data)) < m * np.std(data))[0]

def detect_lane(input_image, output_region_file, SBlevel, m0, pix2sec, mask_image=None, user_interact=True):
  # Convert SBlevel in mag/arcsec^2 to ADU/pix^2. This is done by th Poghson formula.
  # This SBlevel basically defines the borders by Inetnsity where to find dust lane. 
  SBlevel = (pix2sec**2)*10**(0.4*(m0-SBlevel))
  
  # Open and read in fits file
  hdu = pyfits.open(input_image)
  data = hdu[0].data
  header = hdu[0].header
  ny, nx = data.shape
  
  # Open and read in mask image (if any)
  if mask_image is None:
      # In this case mask array will be 0s
      mask = np.zeros((ny,nx))
  else:
      hdu_mask = pyfits.open(mask_image)
      mask = hdu_mask[0].data      


  
  x_min = []; y_min = []
  for x in range(0,nx,1):
    # We make photometric cuts perpendicular to the galaxy plane (which should be horizontal!)
    I = data[0:ny,x] # Just a perpenducilar slice at the x-coordinate
    I_mask = mask[0:ny,x]
    z = np.arange(0,ny,1)
    II = []; zz = []
    for i in range(len(I)):
        if I_mask[i]==0:
            # Take into account non-masked pixels
            II.append(I[i])
            zz.append(z[i])
    II = np.array(II)
    zz = np.array(zz)
    #plt.plot(z,I)
    
    #inds = argrelextrema(I, np.less, mode='wrap')[0]
    inds, _ = find_peaks(-II, distance=5, width=2) # This scipy function is usefull here. Some additional optional parameters might help to trace the dust lane better.
    
    
    for ind in inds:
        if II[ind]>SBlevel: # We consider pixels which are brighter than SBlevel
            x_min.append(x)
            y_min.append(zz[ind])
  y_min = np.array(y_min)        
  
  inds = reject_outliers(y_min, m=2) # Remove outliers

  #plt.plot(z[inds],I[inds],'x')
  
  #plt.show()
  
  f = open(output_region_file,'w')
  for ind in inds:
      f.write('point(%f,%f) # point=x\n' % (x_min[ind]+1,y_min[ind]+1)) # Create ds9 region file to see the result
  f.close()

  if user_interact:
            # Show the image and the overlayed region with the traced dust lane
            ds9Proc = subprocess.Popen(["ds9", input_image,
                                        "-regions", output_region_file,
                                        "-scale", "log"])
            ds9Proc.wait()        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to trace the centerline of a dust lane")
    parser.add_argument("inputImage", help="Input fits image with the centered, horizontally alligned object")
    parser.add_argument("--outputRegion", help="Output region file with centers of the dust lane",default='dust_lane.reg')
    parser.add_argument("--mask", help="Input mask",default=None)
    parser.add_argument("--SBlevel", help="Level for the outer isophote in mag/arcsec^2", type=float, default=23.5)
    parser.add_argument("--scale", help="Pixel scale in arcsec/pix", type=float, default=1.)
    parser.add_argument("--m0", help="Zero-point", type=float, default=28.)
    
    args = parser.parse_args()

    input_image = args.inputImage
    output_region_file = args.outputRegion
    SBlevel = args.SBlevel
    m0 = args.m0
    pix2sec = args.scale
    mask_image = args.mask

    detect_lane(input_image, output_region_file, SBlevel, m0, pix2sec, mask_image=mask_image)
    
