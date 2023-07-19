#! /usr/bin/env python
import os
import sys

from pylab import *
from astropy.io import fits as pyfits

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import glob
from matplotlib.pyplot import cm
import re
from scipy.odr.odrpack import *
import subprocess
import warnings
import argparse

warnings.filterwarnings("ignore")
FNULL = open(os.devnull, 'w')


#*** Colour fonts ***
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

def build_model_image(subcomp_image,model_image_inner,model_image_outer,Rbreak):
  # Read cubcomp_image to get the image of the galaxy:
  hdulist_sub = pyfits.open(model_image_inner)
  image = hdulist_sub[1].data
  
  model_image = np.copy(image)
  residual_image = np.copy(image)

  hdulist_inner = pyfits.open(model_image_inner)
  image_inner = hdulist_inner[2].data  
  center = np.where(image_inner == np.max(image_inner))
  xc,zc = int(center[1]),int(center[0])

  hdulist_outer = pyfits.open(model_image_outer)
  image_outer = hdulist_outer[2].data
  
  ySize, xSize = image_outer.shape

  for k in range(ySize):
    for i in range(xSize):
      if i<=xc-int(Rbreak) or i>=xc+int(Rbreak):
        model_image[k,i] = image_outer[k,i]
        residual_image[k,i] = image[k,i]-image_outer[k,i]
      if i>xc-int(Rbreak) and i<xc+int(Rbreak):
        model_image[k,i] = image_inner[k,i]
        residual_image[k,i] = image[k,i]-image_inner[k,i]
  total_lum = np.sum(image)
  disk_lum = np.sum(model_image)
  
  #hdu = pyfits.PrimaryHDU([image,image,model_image,residual_image])

  hdu = pyfits.HDUList()
  hdu.append(pyfits.ImageHDU(image))
  hdu.append(pyfits.ImageHDU(image))
  hdu.append(pyfits.ImageHDU(model_image))
  hdu.append(pyfits.ImageHDU(residual_image))

  hdu.writeto('model.fits',overwrite=True)  
  
  
  return 'model.fits',disk_lum/total_lum,xc,zc


def get_fits_data(input_image,level):
  hdulist = pyfits.open(input_image)
  data = hdulist[level].data
  header = hdulist[level].header
  return data, header

def make_resid(ref_data,model_data):
        resid_data = np.copy(model_data)
        (dimy,dimx) = ref_data.shape

        for k in range(dimy):
          for i in range(dimx):
            if np.isnan(ref_data[k,i])==False:
                #resid_data[k,i] = fabs(ref_data[k,i] - model_data[k,i]) / fabs(ref_data[k,i])
                resid_data[k,i] = (ref_data[k,i] - model_data[k,i]) / ref_data[k,i]
            else:
                resid_data[k,i] = 0.

        return resid_data

def get_subcomp_ref(input_image):
  hdulist = pyfits.open(input_image)
  data = hdulist[0].data
  if np.sum(data)==0.:
    return 1
  else:
    return 0

def add_keyw_to_header(input_image,level,keyword,value):
  hdulist = pyfits.open(input_image, do_not_scale_image_data=True, mode='update')
  prihdr = hdulist[level].header

  prihdr.append((str(keyword),value),end=True)
  hdulist.flush()  

def main(input_image, model_file, composed_model_file = 'composed_model.fits', comp_names=None, subtract_sky=True, galfitPath='', full_image=False):
  #print bcolors.OKBLUE+'\n\n************ Compiling final GALFIT model image ************' + bcolors.ENDC

  hdulist_input = pyfits.open(input_image)
  data = hdulist_input[0].data
  header = hdulist_input[0].header
  ny,nx = np.shape(data)
  
  # In case of strange Error: Segmentation Fault while providing the sigma image, we copy the input model file and change there the sigma image to none

  f = open(model_file, 'r')
  ff = open('galf_tmp.inp', 'w')
  for line in f:
    if '# Sigma image name (made from data if blank or \"none\")' in line:
      ff.write('C) none          # Sigma image name (made from data if blank or \"none\")\n')
    elif '# Bad pixel mask (FITS image or ASCII coord list)' in line:
       mask_image = line.split()[1]
       ff.write(line)
    elif '# Image region to fit (xmin xmax ymin ymax)' in line:
      if full_image:
        ff.write('H) 1  %i 1  %i # Image region to fit (xmin xmax ymin ymax)' % (nx,ny))
      else:
        ff.write(line)
        xmin = int(line.split()[1])
        xmax = int(line.split()[2])
        ymin = int(line.split()[3])
        ymax = int(line.split()[4])
        if mask_image.lower()!="none":
            hdulist_mask = pyfits.open(mask_image)
            mask = hdulist_mask[0].data
            header_mask = hdulist_mask[0].header
            outHDU = pyfits.PrimaryHDU(mask[ymin-1:ymax,xmin-1:xmax])
            outHDU.writeto('mask_cropped.fits', overwrite=True)                 
            
        
    else:
      ff.write(line)
    
  f.close()
  ff.close()

  model_file = 'galf_tmp.inp'
  subprocess.call("%sgalfit -o1 %s" % (galfitPath, model_file), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
  subprocess.call("%sgalfit -o3 %s" % (galfitPath, model_file), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

  with open(model_file,'r') as fff:
      for line in fff:
        if "B" in line and "# Output data image block\n" in line:
            galf_image = line.split()[1]
  fff.close()

  
  # Read reference image:
  no_ref = get_subcomp_ref('subcomps.fits')

  ref_data, head = get_fits_data('subcomps.fits',no_ref)

  # Subtract sky:
  sky_level = 'None'
  hdulist = pyfits.open('subcomps.fits')
  if subtract_sky==True:
    for k in range(len(hdulist)-(1+no_ref)):
      prihdr = hdulist[k+1+no_ref].header
      if 'sky' in prihdr['OBJECT']:
        sky_level = k+1+no_ref
        break
    try:
      #print sky_level
      ref_data = ref_data - get_fits_data('subcomps.fits',sky_level)[0]
    except:
      z=1
      
  # Read total model image:
  model_data,header_mod = get_fits_data(galf_image,0) 
  
  if subtract_sky==True and sky_level != 'None':
    #print 'here'
    model_data = model_data - get_fits_data('subcomps.fits',sky_level)[0]
  
  # Find residual image:
  resid_data = ref_data - model_data
  
  # Create relative residual image (in percentage):
  resid_data_percent = make_resid(ref_data,model_data)

  # Append these data to the final output cube image:
  hdu = pyfits.HDUList()
  hdu.append(pyfits.ImageHDU(ref_data,header))
  hdu.append(pyfits.ImageHDU(model_data,header_mod))
  hdu.append(pyfits.ImageHDU(resid_data,header))
  hdu.append(pyfits.ImageHDU(resid_data_percent,header))



  
  # Find components of the model:
  name_of_model = []
  levels = []
  #hdulist = pyfits.open('subcomps.fits')
  for k in range(len(hdulist)-(1+no_ref)):
     prihdr = hdulist[k+1+no_ref].header
     if 'sky' not in prihdr['OBJECT']:
        name_of_model.append(prihdr['OBJECT'])
        levels.append(k+1+no_ref)
  #print name_of_model
  for k in range(len(name_of_model)):
    comp_model_data = get_fits_data('subcomps.fits',levels[k])[0]
    hdu.append(pyfits.ImageHDU(comp_model_data,header))

  hdu.writeto(composed_model_file,overwrite=True) 

  # Add names of the layers to the headers
  add_keyw_to_header(composed_model_file,0,'NAME_OF_LAYER','data')
  add_keyw_to_header(composed_model_file,1,'NAME_OF_LAYER','model')
  add_keyw_to_header(composed_model_file,2,'NAME_OF_LAYER','resid')
  add_keyw_to_header(composed_model_file,3,'NAME_OF_LAYER','resid')

  if comp_names is None:
    for k in range(len(name_of_model)):
        add_keyw_to_header(composed_model_file,k+4,'NAME_OF_LAYER',name_of_model[k])
  else:
    for k in range(len(comp_names)):
        add_keyw_to_header(composed_model_file,k+4,'NAME_OF_LAYER',comp_names[k])
        
  # Remove tmp files:
  os.remove(galf_image)
  os.remove('subcomps.fits')
  os.remove(model_file)
  #print 'Done!'
  return composed_model_file


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert region file with a mask to a fits segmentation image")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("model_file", help="Input model file")
    parser.add_argument("--output_image", nargs='?', const=1, help="Optional: Output image", type=str, default='composed_model.fits') 
    parser.add_argument("--full", nargs='?', const=1, help="Optional: Full image (default: false)?", type=bool, default=False) 
    
    args = parser.parse_args()

    input_image = args.input_image
    model_file = args.model_file
    output_image = args.output_image
    full_image = args.full


    
    main(input_image, model_file, composed_model_file = output_image, subtract_sky=True, galfitPath='', full_image=full_image)



