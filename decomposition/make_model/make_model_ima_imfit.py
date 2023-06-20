#! /usr/bin/env python
# EXAMPLE: python ~/MEGA/MyPrograms/IMAN/Decomposition/make_model_ima_imfit.py galaxy.fits imf_output.dat psf_rebin.fits thin_disc,thick_disc,bulge
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
import argparse
import warnings
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


def change_dim(input_image, input_file, full_image=False):
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data    
    ny,nx = np.shape(data)
    hdulist.close()
    f = open(input_file, 'r')

    lines = f.readlines()
    xmin = 1
    xmax = nx
    ymin = 1
    ymax = ny
    mask_image = None
    for line in lines:
            if 'imfit' in line and '-c' in line and '[' in line and ']' in line:
                coords = line.split('[')[1].split(']')[0]
                xmin = int(coords.split(':')[0])
                xmax = int(coords.split(':')[1].split(',')[0])
                ymin = int(coords.split(',')[1].split(':')[0])
                ymax = int(coords.split(':')[2])
            if 'imfit' in line and '-c' in line and '--mask=' in line:
                mask_image = line.split('--mask=')[1].split()[0].split('[')[0]
                
    f.close() 
    if mask_image is not None:
        hdulist_mask = pyfits.open(mask_image)
        mask = hdulist_mask[0].data
        #print(xmin,xmax,ymin,ymax,mask[ymin-1:ymax,xmin-1:xmax])
        outHDU = pyfits.PrimaryHDU(mask[ymin-1:ymax,xmin-1:xmax])
        outHDU.writeto('mask_cropped.fits', clobber=True)         
    
    return xmin,xmax,ymin,ymax



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

  hdu.writeto('model.fits',clobber=True)  
  
  
  return 'model.fits',disk_lum/total_lum,xc,zc


def get_fits_data(input_image,level):
  hdulist = pyfits.open(input_image)
  data = hdulist[level].data
  return data

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

def add_keyw_to_header(input_image,level,keyword,value):
  hdulist = pyfits.open(input_image, do_not_scale_image_data=True, mode='update')
  prihdr = hdulist[level].header

  prihdr.append((str(keyword),value),end=True)
  hdulist.flush()  

def read_imfit_objects(model_file):
    objects = []
    with open(model_file) as f:
        lines = f.readlines()
    for line in lines:
      if 'FUNCTION' in line:
        objects.append(line.split()[1])
    return objects
  

def main(input_image, model_file, psf_image=None, composed_model_file = 'composed_model.fits', comp_names=[], imfitPath='', oned=False):
  print(bcolors.OKBLUE+'\n\n************ Compiling final IMFIT model image ************' + bcolors.ENDC)
  if not oned:
    if psf_image is not None:
        subprocess.call("%smakeimage %s -o total_model.fits --output-functions Comp_ --refimage %s --psf %s" % (imfitPath, model_file,input_image, psf_image), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    else:
        subprocess.call("%smakeimage %s -o total_model.fits --output-functions Comp_ --refimage %s" % (imfitPath, model_file,input_image), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
  else:
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    ny,nx = np.shape(data)
    if nx==1:
        nx = 2
    if ny==1:
        ny = 2
        
    if psf_image is not None:
        subprocess.call("%smakeimage %s -o total_model.fits --output-functions Comp_ --ncols %i --nrows %i --psf %s" % (imfitPath, model_file, nx, ny, psf_image), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    else:
        subprocess.call("%smakeimage %s -o total_model.fits --output-functions Comp_ --ncols %i --nrows %i" % (imfitPath, model_file, nx, ny), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)    

      
  #model_images = glob.glob('Comp_*.fits')
  
  xmin,xmax,ymin,ymax = change_dim(input_image, model_file)
  if oned:
      if xmax==1:
          xmax = 2
      if ymax==1:
          ymax = 2
      
  
  
  #Model_images = [None] * len(model_images)
  
  imfit_objects = read_imfit_objects(model_file)
  
  model_images = []
  for k in range(len(imfit_objects)):
      model_images.append('Comp_%i_%s.fits' % (k+1,imfit_objects[k]))

  print('Model images in order are:',model_images)

  if comp_names!=[]:
    print('They have been replaced by:',comp_names)

  Model_images = []
  Names = []
  for k in range(len(imfit_objects)):
    number_of_model = k
    if 'FlatSky' not in imfit_objects[number_of_model]:
      if comp_names!=[]:
        Model_images.append(model_images[k])
        if '_' in comp_names[k]:
            l = comp_names[k].split('_')
            Names.append(' '.join(l))
        else:
            Names.append(comp_names[k])
      else:
        Model_images.append(model_images[k])
        Names.append(imfit_objects[number_of_model])

  # Read reference image:
  ref_data = get_fits_data(input_image,0)

  # Read total model image:
  model_data = get_fits_data('total_model.fits',0) 
  
  # Find residual image:
  resid_data = ref_data - model_data
  
  # Create relative residual image (in percentage):
  resid_data_percent = make_resid(ref_data,model_data)

  # Append these data to the final output cube image:
  hdu = pyfits.HDUList()
  hdu.append(pyfits.ImageHDU(ref_data[ymin-1:ymax,xmin-1:xmax]))
  hdu.append(pyfits.ImageHDU(model_data[ymin-1:ymax,xmin-1:xmax]))
  hdu.append(pyfits.ImageHDU(resid_data[ymin-1:ymax,xmin-1:xmax]))
  hdu.append(pyfits.ImageHDU(resid_data_percent[ymin-1:ymax,xmin-1:xmax]))
  
  for k in range(len(Model_images)):
    comp_model_data = get_fits_data(Model_images[k],0)
    hdu.append(pyfits.ImageHDU(comp_model_data[ymin-1:ymax,xmin-1:xmax]))

  hdu.writeto(composed_model_file, overwrite=True) 

  # Add names of the layers to the headers
  add_keyw_to_header(composed_model_file,0,'NAME_OF_LAYER','data')
  add_keyw_to_header(composed_model_file,1,'NAME_OF_LAYER','model')
  add_keyw_to_header(composed_model_file,2,'NAME_OF_LAYER','resid')
  add_keyw_to_header(composed_model_file,3,'NAME_OF_LAYER','resid')

  for k in range(len(Model_images)):
    add_keyw_to_header(composed_model_file,k+4,'NAME_OF_LAYER',Names[k])#.lower())
  
  for k in range(len(model_images)):
    os.remove(model_images[k])
    
  # Remove tmp files:
  os.remove('total_model.fits')
  
  print('Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert region file with a mask to a fits segmentation image")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("model_file", help="Input model file")
    parser.add_argument("--output_image", nargs='?', const=1, help="Optional: Output image", type=str, default='composed_model.fits') 
    parser.add_argument("--psf_image", nargs='?', const=1, help="Optional: Psf image", type=str, default=None) 
    parser.add_argument("--comp_names", nargs='?', const=1, help="Optional: Names for the components", type=str, default=None) 
    
    args = parser.parse_args()

    input_image = args.input_image
    model_file = args.model_file
    output_image = args.output_image
    psf_image = args.psf_image
    comp_names = args.comp_names
    if comp_names is not None:
        comp_names = comp_names.split(',')
    else:
        comp_names = []

    main(input_image, model_file, psf_image=psf_image, composed_model_file = output_image, comp_names=comp_names, imfitPath='')
