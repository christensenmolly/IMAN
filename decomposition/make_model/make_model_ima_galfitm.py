#! /usr/bin/env python
# python ~/MyPrograms/IMAN/Decomposition/make_model_ima_galfitm.py model.galfit.01.band
import os
import sys

from pylab import *
import pyfits

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import glob
from matplotlib.pyplot import cm
import re
from scipy.odr.odrpack import *
import subprocess
import shutil

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

  hdu.writeto('model.fits',clobber=True)  
  
  
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
	    if ref_data[k,i]>0. and np.isnan(ref_data[k,i])==False:
	      #resid_data[k,i] = fabs(ref_data[k,i] - model_data[k,i]) / fabs(ref_data[k,i])
	      resid_data[k,i] = (ref_data[k,i] - model_data[k,i]) / fabs(ref_data[k,i])
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

def add_keyw_to_header(hdulist,level,keyword,value):
   hdulist[level].header.append((str(keyword),value),end=True)
   return hdulist

def main(model_file, composed_model_file = 'composed_model.fits', subtract_sky=True, galfitPath=''):
  print bcolors.OKBLUE+'\n\n************ Compiling final GALFITM model image ************' + bcolors.ENDC
 
  # In case of strange Error: Segmentation Fault while providing the sigma image, we copy the input model file and change there the sigma image to none
  f = open(model_file, 'r')
  ff = open('galf_tmp.inp', 'w')
  for line in f:
    if '# Sigma image and min. sigma factor (made from data if blank or \"none\")' in line:
      ff.write('C) none          # Sigma image and min. sigma factor (made from data if blank or \"none\")\n')
    else:
      ff.write(line)
  f.close()
  ff.close()

  model_file = 'galf_tmp.inp'
  sky_comp = None

  with open(model_file,'r') as fff:
      for line in fff:
	  if "B" in line and "# Output data image block\n" in line:
	    galf_image = line.split()[1]
          if "A1" in line and '# Band labels' in line:
            bands = line.split('A1) ')[-1].split('     # Band labels')[0].split(',')
          if "# Component number:" in line:
              comp_number = int(line.split("# Component number: ")[-1])
          if "0) sky" in line  and '#  Component type' in line:
            sky_comp = comp_number - 1

  fff.close()

  
  subprocess.call("%sgalfitm -o3 %s" % (galfitPath, model_file), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
  shutil.move(galf_image, 'subcomps.fits')
  subprocess.call("%sgalfitm -o2 %s" % (galfitPath, model_file), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
  
  n_sky = -1
  NNN = 0
  hdu = pyfits.HDUList()  
  for kk in range(len(bands)):
        # Read in reference image:
        ref_data, header = get_fits_data(galf_image, kk)

        # Subtract sky:
        hdulist = pyfits.open(galf_image)
        if subtract_sky==True and sky_comp!=None:
            prihdr = hdulist[len(bands)].header
            n_sky = n_sky + sky_comp + 1
            ref_data = ref_data - get_fits_data('subcomps.fits', n_sky)[0]

            
        # Read total model image:
        model_data,header_mod = get_fits_data(galf_image,kk+len(bands)) 

        if subtract_sky==True and sky_comp!=None:
            model_data = model_data - get_fits_data('subcomps.fits',n_sky)[0]
        
        # Find residual image:
        resid_data = ref_data - model_data
        
        # Create relative residual image (in percentage):
        resid_data_percent = make_resid(ref_data,model_data)

        # Append these data to the final output cube image:
        hdu.append(pyfits.ImageHDU(ref_data,header))
        hdu.append(pyfits.ImageHDU(model_data,header_mod))
        hdu.append(pyfits.ImageHDU(resid_data))
        hdu.append(pyfits.ImageHDU(resid_data_percent))

  
        # Find components of the model:
        name_of_model = []

        if sky_comp!=None:
            for k in range(kk, kk+comp_number):
                if k!=kk+sky_comp:
                        name_of_model.append(header_mod['COMP_%i' % (k-kk+1)])
                        data_comp,header_comp = get_fits_data('subcomps.fits',k)
                        hdu.append(pyfits.ImageHDU(data_comp,header_comp))
        else:
            for k in range(kk, kk+comp_number):
                        name_of_model.append(header_mod['COMP_%i' % (k-kk+1)])
                        data_comp,header_comp = get_fits_data('subcomps.fits',k)
                        hdu.append(pyfits.ImageHDU(data_comp,header_comp))           
            



        # Add names of the layers to the headers
        hdu = add_keyw_to_header(hdu,NNN+0,'NAME_OF_LAYER','data')
        if NNN==0:
            hdu = add_keyw_to_header(hdu,NNN+0,'number_of_bands',len(bands))
            if sky_comp!=None:
                hdu = add_keyw_to_header(hdu,NNN+0,'number_of_comps',comp_number-1)
            else:
                hdu = add_keyw_to_header(hdu,NNN+0,'number_of_comps',comp_number)

        if bands[kk].strip()=='P100':
           bands[kk]='PACS100' 
        if bands[kk].strip()=='P160':
           bands[kk]='PACS160'         
        if bands[kk].strip()=='S250':
           bands[kk]='SPIRE250' 
        if bands[kk].strip()=='S350':
           bands[kk]='SPIRE350' 
        if bands[kk].strip()=='S500':
           bands[kk]='SPIRE500' 
 
        hdu = add_keyw_to_header(hdu,NNN+0,'BAND',bands[kk])
        hdu = add_keyw_to_header(hdu,NNN+1,'NAME_OF_LAYER','model')
        hdu = add_keyw_to_header(hdu,NNN+2,'NAME_OF_LAYER','resid')
        hdu = add_keyw_to_header(hdu,NNN+3,'NAME_OF_LAYER','resid')

        for k in range(len(name_of_model)):
            hdu = add_keyw_to_header(hdu,NNN+k+4,'NAME_OF_LAYER',name_of_model[k])

        hdu.writeto(composed_model_file,clobber=True)
        
        NNN = NNN+k+4+1
  # Remove tmp files:
  os.remove(galf_image)
  os.remove('subcomps.fits')
  os.remove(model_file)
  print 'Done!'
  return composed_model_file



if __name__ == '__main__':
    model_file = str(sys.argv[1])
    main(model_file)
