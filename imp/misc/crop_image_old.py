#!/usr/bin/env python

# Import the necessary modules
import numpy as np
import math
from numpy import *
import sys
from astropy.stats import sigma_clipped_stats

try:
  import astropy.io.fits as pyfits
  from astropy import wcs
except:
  warnings.warn("Astropy is not installed! No WCS has been added to the header!")





def main(input_image, output_image, x_l, y_l, x_r, y_r, hdu=0):
    '''
    Function to crop the input_image based on the given coordinates
    of the bottom left and top right corners. The crop will be
    saved in the output_image.
    This function keeps the WCS!
    
    x_l, y_l, x_r, y_r - pixel coordinates in the ds9 format
    (start from 0.5).
    hdu - Operate on the specified FITS header extension (HDU), int.
    '''

    hdulist = pyfits.open(input_image)
    referenceheader = pyfits.getheader(input_image, hdu)
    inframe = hdulist[hdu].data
    ny,nx = np.shape(inframe)
    
    try:
        if referenceheader['CTYPE2'] == 'DEC---TAN':
            referenceheader['CTYPE2'] = 'DEC--TAN'
        '''
        if referenceheader['CTYPE1'] == 'RA---TAN':
            referenceheader['CTYPE1'] = 'RA--TAN'
        '''

        if 'COMMENT' in referenceheader:
            del referenceheader['COMMENT']
        if 'HISTORY' in referenceheader:
            del referenceheader['HISTORY']
        if '' in referenceheader:
            del referenceheader['']  

        w = wcs.WCS(referenceheader)
        xcen = (x_l+x_r)/2.
        ycen = (y_l+y_r)/2.
        
        pixcrd = np.array([[xcen, ycen]], np.float_)
        world = w.wcs_pix2world(pixcrd, 1)

        xcen_w,ycen_w = world[0,0],world[0,1]

        
        
        ref_pixel = inframe[int(math.floor(ycen-0.5)),int(math.floor(xcen-0.5))]
        inframe[int(math.floor(ycen-0.5)),int(math.floor(xcen-0.5))] = 999999
        w = True
    except:
        w = False
        
    if int(math.floor(y_l-0.5))<0:
        ymin = 0
        print("ymin is set to 0")
    else:
        ymin = int(math.floor(y_l-0.5))
    if int(math.floor(x_l-0.5))<0:
        xmin = 0
        print("xmin is set to 0")
    else:
        xmin = int(math.floor(x_l-0.5))
    if int(math.floor(y_r-0.5))>=ny:
        ymax = ny-1    
        print("ymax is set to dim_y")
    else:
        ymax = int(math.floor(y_r-0.5))
    if int(math.floor(x_r-0.5))>=nx:
        xmax = nx-1
        print("xmax is set to dim_x")
    else:
        xmax = int(math.floor(x_r-0.5))

    try:
      outframe = inframe[ymin:ymax+1, xmin:xmax+1]
    except:
      print('Wrong range to crop the image! Exiting ...')
      exit()


    if w==True:
        ref_coords = np.where(outframe == 999999)

        ycen_new = int(ref_coords[0][0])
        xcen_new = int(ref_coords[1][0])

        outframe[ycen_new, xcen_new] = ref_pixel
        referenceheader['CRPIX1'] = xcen_new+1.0
        referenceheader['CRPIX2'] = ycen_new+1.0
        referenceheader['CRVAL1'] = xcen_w
        referenceheader['CRVAL2'] = ycen_w
    else:
        ny,nx = np.shape(outframe)
        xcen_new = nx/2.
        ycen_new = ny/2.
    hdu = pyfits.PrimaryHDU(outframe, referenceheader)
    hdu.writeto(output_image, clobber=True)
    
    return xcen_new, ycen_new

if __name__ == '__main__':
    # FOR HELP:
    if str(sys.argv[1]) == '-h':
      print('\n\n\nUSAGE:')
      print('python crop_ima.py [input_image] [output_image] [xl] [yl] [xr] [yr]')
      print('EXAMPLE: python crop_image.py galaxy_inp.fits galaxy_out.fits 10 10 100 100')

    input_image = str(sys.argv[1])
    output_image = str(sys.argv[2])
    x_l = float(sys.argv[3])
    y_l = float(sys.argv[4])
    x_r = float(sys.argv[5])
    y_r = float(sys.argv[6])
    
    main(input_image, output_image, x_l, y_l, x_r, y_r)
