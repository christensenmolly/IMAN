#!/usr/bin/python
# DESCRIPTION:
# Script to fit sky background along the lines (rows) inclined by PA. A mask is required.
# MINIMAL USAGE: python subtract_lines.py [input_image] [input_mask]

import numpy as np
from astropy.io import fits as pyfits
from astropy.modeling import models, fitting
import warnings
warnings.filterwarnings("ignore")
from astropy.stats import sigma_clipped_stats
import argparse
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import subprocess
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from itertools import product
from matplotlib.path import Path
import os
import sys
import shutil

LOCAL_DIR = "/imp/sky_fitting"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/masking'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/rotate'))
sys.path.append(os.path.join(IMAN_DIR, 'plotting/2dprofile'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))

import merge_masks
import plot_smoothed_image
import rotate_image
import rebin_image

def convert_segm_to_boolean(mask):
    return np.ma.make_mask(mask, shrink=False)



def do_each_polygon(input_image, mask_image, hdu_inp=0, output_image=None, PA=0., polygon_mask = None, sky_image=None):
      if polygon_mask is not None:
          merge_masks.main([mask_image,polygon_mask], 'tmp_mask.fits')
      else:
          shutil.copy(mask_image, 'tmp_mask.fits')
      
      
      if PA!=0.:
          rotate_image.main(input_image, PA, xc=None, yc=None, output_image='rot_tmp.fits', hdu_inp=0, cval=float('nan'), cropping=False)
          rotate_image.main('tmp_mask.fits', PA, xc=None, yc=None, output_image='rot_mask_tmp.fits', hdu_inp=0, cval=float('nan'), cropping=False)
          input_image1 = 'rot_tmp.fits'
          mask_image1 = 'rot_mask_tmp.fits'
      else:
          input_image1 = input_image
          mask_image1 = 'tmp_mask.fits'
      #exit()
      hdulist = pyfits.open(input_image1)
      data = hdulist[hdu_inp].data      
      header = hdulist[hdu_inp].header
      ny,nx = np.shape(data)
      
      new_data = np.copy(data)
      sky_data = np.zeros((ny,nx))

      hdulist_mask = pyfits.open(mask_image1)
      mask = hdulist_mask[0].data     

      mask_astropy = convert_segm_to_boolean(mask)
      #N = 0 
      for k in range(ny):
            #try:
            line = data[k,:]
            line_mask = mask_astropy[k,:]  
            
            for i in range(len(line)):
                if np.isnan(line[i]):
                    line_mask[i]=True
            if not np.all(line_mask):
                mean, median, std = sigma_clipped_stats(line, sigma=3.0, mask=line_mask)
            else:
                #print(line,'here')
                median = 0.
                #N = N + 1
            new_data[k] = line - median
            sky_data[k] = median * np.arange(len(line))
            #print(median)
            #except:
            #  z=1

      if output_image is None:
          output_image = input_image.split('.fits')[0] + '_linesub.fits'
      
      outHDU = pyfits.PrimaryHDU(new_data, header=header)
      outHDU.writeto('tmp_linesub.fits', clobber=True)
      
      if sky_image is not None:
        outHDU = pyfits.PrimaryHDU(sky_data, header=header)
        outHDU.writeto('sky_linesub.fits', clobber=True)          
      
      if PA!=0.:
        rotate_image.main('tmp_linesub.fits', -PA, xc=None, yc=None, output_image='tmp_rot.fits', hdu_inp=0, cval=float('nan'), cropping=False)
        rebin_image.rebin(input_image, 'tmp_rot.fits', output_image, hdu_ref=0, hdu_inp=0, preserve_bad_pixels=True)
        if sky_image is not None:
            rotate_image.main('sky_linesub.fits', -PA, xc=None, yc=None, output_image='tmp_sky_rot.fits', hdu_inp=0, cval=float('nan'), cropping=False)
            rebin_image.rebin(input_image, 'tmp_sky_rot.fits', sky_image, hdu_ref=0, hdu_inp=0, preserve_bad_pixels=True)            
      else:
        shutil.move('tmp_linesub.fits', output_image)
        if sky_image is not None:
            shutil.move('sky_linesub.fits', sky_image)
      
      for file in ['rot_tmp.fits', 'rot_mask_tmp.fits', 'tmp_linesub.fits', 'tmp_mask.fits', 'tmp_rot.fits','sky_linesub.fits']:
        if os.path.exists(file):
            os.remove(file)
      return output_image


    


def main(input_image, mask_image, hdu_inp=0, output_image=None, PA=0., region_file = None, sky_image=None):
    print('Subtracting lines...')
    if region_file is None:
      output_image = do_each_polygon(input_image, mask_image, hdu_inp=hdu_inp, output_image=output_image, PA=PA, polygon_mask = None, sky_image=sky_image)
    else:
      hdulist = pyfits.open(input_image)
      data = hdulist[hdu_inp].data
      ny,nx = np.shape(data)
      polygon_mask = np.ones((ny,nx))
        
      f = open(region_file, "r")
      N = 0
      for line in f: 
        if 'polygon(' in line:
            print('\tPolygon #%i...' % (N+1))
            param = map(float, ((line.split('polygon(')[1]).split(')')[0]).split(','))
            param = list(param)

            param2 = [None]*(int(len(param)/2))
            for i in range(0,len(param2)):  param2[i] = (int(param[2*i]),int(param[2*i+1])) 
            param2.append(param2[0])
            codes = []
            codes.append(Path.MOVETO)
            for i in range(1,len(param2)-1):  codes.append(Path.LINETO)
            codes.append(Path.CLOSEPOLY)
            path = Path(param2, codes)

            coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
            X = []; Y = []
            for k in range(0,len(coords)-1,2):
               X.append(int(float(coords[k])))
               Y.append(int(float(coords[k+1])))

            xmin = min(X)
            ymin = min(Y)
            xmax = max(X)
            ymax = max(Y)
            for i, j in product(range(xmin,xmax), range(ymin,ymax)):
                try:
                    if path.contains_point((i,j)):
                        polygon_mask[j-1][i-1] = 0.
                except:
                    zz=1
            if 'text=' in line:
                PA = float(line.split('text=')[-1].rstrip().split('{')[-1].split('}')[0])
            
            outHDU = pyfits.PrimaryHDU(polygon_mask)
            outHDU.writeto('polygon_mask_tmp.fits', clobber=True)            
            
            output_image = do_each_polygon(input_image, mask_image, hdu_inp=hdu_inp, output_image=output_image, PA=PA, polygon_mask = 'polygon_mask_tmp.fits', sky_image=sky_image)
            shutil.copy(output_image, 'tmp_tmp.fits')
            input_image = 'tmp_tmp.fits'
            N = N + 1
      
      os.remove('polygon_mask_tmp.fits')
      os.remove('tmp_tmp.fits')
    return output_image
      


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky fitting along the rows")
    parser.add_argument("inputImage", help="Input fits image ")
    parser.add_argument("inputRegion", help="Optional: Output region file", type=str, default='general_mask.reg') 
    parser.add_argument("--inputMask", help="Input mask image", type=str, default=None) 
    
    parser.add_argument("--pa", help="Optional: Position angle of the lines, clock-wise, from positive X", type=float, default=0.) 
    parser.add_argument("--output_image", help="Optional: Output image", type=str, default='') 
    parser.add_argument("--output_sky", help="Optional: Output sky", type=str, default=None)
    
    args = parser.parse_args()

    input_image = args.inputImage
    region_file = args.inputRegion
    output_image = args.output_image
    sky_image = args.output_sky
    mask_image = args.inputMask
    PA = args.pa


    main(input_image, mask_image, hdu_inp=0, output_image=output_image, PA=PA, region_file = region_file, sky_image=sky_image)










#output_image = main('NGC2683.phot.1_nonan_someskysub.fits', 'general_mask.fits', hdu_inp=0, output_image=None, PA=-71., sky_image='sky.fits')
#png_image = plot_smoothed_image.main('NGC2683.phot.1_nonan_someskysub_linesub.fits', None, '39492', 20.48, 0.75, SB_bright=24., SB_faint=27., sigma_smooth=5., add_axes=True)
#png_image = plot_smoothed_image.main('NGC2683.phot.1_nonan_someskysub.fits', None, '39492', 20.48, 0.75, SB_bright=24., SB_faint=27., sigma_smooth=5., add_axes=True)