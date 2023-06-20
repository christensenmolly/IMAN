#! /usr/bin/env python

import sys
from os.path import exists
from os import remove
from scipy.ndimage.filters import convolve
from numpy import max, sum
import numpy as np
import astropy.io.fits as pyfits



def _add_hdu(hdu_list, data):
        if len(hdu_list) == 0:
            hdu = pyfits.PrimaryHDU(data)
        else:
            hdu = pyfits.ImageHDU(data)
        hdu_list.append(hdu)
        return hdu

def main(input_images, output_image='cube.fits', layer_desc='standard', comps=None):
  # Usually layer=0, layer_desc='standard'

  # Append these data to the final output cube image:
  if layer_desc=='standard':
    hdu = pyfits.HDUList()
        
    for k in range(len(input_images)):    
            hdulist = pyfits.open(input_images[k])
            data = hdulist[0].data
            header = hdulist[0].header
            if comps is not None:
                header['NAME_OF_LAYER'] = comps[k]
            
            hdu.append(pyfits.ImageHDU(data, header))
            hdulist.close()
    hdu.writeto(output_image, overwrite=True) 
    
  else:
        hdu_list = pyfits.HDUList()
        Input_images = []
        for k in range(len(input_images)):    
                hdulist = pyfits.open(input_images[k])
                data = hdulist[0].data
                header = hdulist[0].header
                Input_images.append(data)
        Input_images = np.array(Input_images)
        nimages = len(Input_images)
        im = Input_images[0]
        dtype = float#im.array.dtype
        nx = np.shape(data)[1]
        ny = np.shape(data)[0]
        # Use the first image's wcs and bounds
        #wcs = im.wcs
        #bounds = im.bounds
        # Note: numpy shape is y,x
        array_shape = (nimages, ny, nx)
        cube = np.zeros(array_shape, dtype=dtype)
        for k in range(nimages):
            im = Input_images[k]
            #nx_k = im.xmax-im.xmin+1
            #ny_k = im.ymax-im.ymin+1
            #if nx_k != nx or ny_k != ny:
            #    print('ERROR. exiting...')
            #    exit()
            cube[k,:,:] = Input_images[k]
      
        
        hdu = _add_hdu(hdu_list, cube)
        #if wcs:
        #    wcs.writeToFitsHeader(hdu.header, bounds)

        hdu.writeto(output_image, overwrite=True)
        #header.tofile(output_image, overwrite=True)
      
  


if __name__ == "__main__":
    input_images = sys.argv[1]
    output_image = sys.argv[2]
    layer_desc = str(sys.argv[3])
    comps = str(sys.argv[4])
    
    input_images = input_images.split(',')
    comps = comps.split(',')
    

    main(input_images, output_image, layer_desc, comps)

    
