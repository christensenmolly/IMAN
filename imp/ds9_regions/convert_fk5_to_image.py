#!/usr/bin/python
# Script to convert a region file written in the image format to fk5 
# EXAMPLE: convert_image_to_fk5('mask.reg', 'galaxy.fits', 'mask_wcs.reg')

# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import argparse
import math

from astropy import wcs

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def resolution(input_image):
    from astropy.io import fits
    hdulist = fits.open(input_image)#, ignore_missing_end=True)
    header = hdulist[0].header

    pixelscale = 1.
    note = ' '

    if 'PIXSCALE' in header:

        pixelscale = header['PIXSCALE']

    elif 'SECPIX' in header:

        pixelscale = header['SECPIX']

    elif 'PFOV' in header:

        pixelscale = header['PFOV']

    elif 'CD1_1' in header and 'CD1_2' in header:

        pixelscale = math.sqrt(header['CD1_1']**2 + header['CD1_2']**2 ) * 3600.0

    elif 'CD1_1' in header:

        pixelscale = abs(header['CD1_1']) * 3600.0

    elif 'CDELT1' in header:

        pixelscale = abs(header['CDELT1']) * 3600.0

    else:

        print("Could not determine the pixel scale from the image header. Set to 1 pix/arcsec.")
        note = '*'

    # Return the pixel scale (in arcseconds)
    hdulist.close()
    return str(pixelscale),note

def main(input_region_file, input_image, output_region_file):
        #print 'Converting the regions to image...'
        hdulist = pyfits.open(input_image)
        data = hdulist[0].data
        header = hdulist[0].header
        
        pix2sec, note = resolution(input_image)
        pix2sec = float(pix2sec)
        
        if 'COMMENT' in header:
            del header['COMMENT']
        if 'HISTORY' in header:
            del header['HISTORY']
        if '' in header:
            del header['']  
        w = wcs.WCS(header)

        f_tmp = open(input_region_file,'r')
        f_reg = open(output_region_file,'w')
        f_reg.write('%s\n' % ('image') )

        for Line in f_tmp:
            if 'polygon' in Line:
              coords = Line.replace('(',',').replace(')',',').split(',')[1:-1]
              polygon = []; world = []
              for k in range(0,len(coords)-1,2):
                world.append([float(coords[k]),float(coords[k+1])])
              
              world = np.array(world, np.float_)

              # Convert pixel coordinates to wcs
              pixcrd = w.wcs_world2pix(world, 1)
      
              for k in range(len(pixcrd)):
                  polygon = str(pixcrd[k][0]) + ',' + str(pixcrd[k][1])
                  if k==0:
                    f_reg.write('polygon(%s,' % (polygon) )
                  elif k>0 and k<len(pixcrd)-1:
                    f_reg.write('%s,' % (polygon) )
                  else:
                    f_reg.write('%s)\n' % (polygon) )    
            if 'ellipse' in Line:
                      params = Line.split(",")
                      cen = Point(float(params[0].split('(')[1]),
                          float(params[1]))
                      ellA = float(params[2].split('\"')[0])
                      ellB = float(params[3].split('\"')[0])
                      ellPA = float(params[4].split(')')[0])
                      if ellA < ellB:
                        ellA, ellB = ellB, ellA
                        ellPA += 90
                      world = np.array([[cen.x,cen.y]], np.float_)
                      pixcrd = w.wcs_world2pix(world, 1)    
                      f_reg.write("ellipse(%f,%f,%f,%f,%f)\n" % (pixcrd[0][0], pixcrd[0][1],ellA/pix2sec, ellB/pix2sec, ellPA))
            if 'circle(' in Line:
                      params = Line.split("(")[-1].split(",")[0]
                      cen = Point(float(Line.split("(")[-1].split(",")[0]),float(Line.split("(")[-1].split(",")[1]))
                      ellA = float(Line.split("(")[-1].split(",")[2].split(')')[0].split('\"')[0])
                      ellB = ellA
                      ellPA = 0.
                      world = np.array([[cen.x,cen.y]], np.float_)
                      pixcrd = w.wcs_world2pix(world, 1)    
                      f_reg.write("circle(%f,%f,%f)\n" % (pixcrd[0][0], pixcrd[0][1],ellA/pix2sec))

        f_tmp.close()
        f_reg.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convertion of the region file from the image format to fk5 format")
    parser.add_argument("input_region_file", help="Input region file in the image format")
    parser.add_argument("input_image", help="Input image with wcs")
    parser.add_argument("output_region_file", help="Output region file")
 
    args = parser.parse_args()

    main(args.input_region_file, args.input_image, args.output_region_file)
