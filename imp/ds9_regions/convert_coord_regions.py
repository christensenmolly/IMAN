# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import argparse
import math
import os
import sys
from astropy import wcs


LOCAL_DIR = "/imp/ds9_regions"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))
import crop_galaxy_image


def resolution(header):


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


    return str(pixelscale),note

def find_crota(input_image):
        hdulist = pyfits.open(input_image)
        data = hdulist[0].data
        header = hdulist[0].header
        if 'CROTA1' not in header and 'CROTA2' not in header:
            CD11 = header["CD1_1"]
            CD12 = header["CD1_2"]
            CD21 = header["CD2_1"]
            CD22 = header["CD2_2"]
            
            alpha1 = np.degrees(math.atan2(CD21, CD11))
            alpha2 = -np.degrees(math.atan2(CD12, CD22))
            
            
            #print(alpha1, alpha2) 
            return alpha1, alpha2, header



def convert_ellipse_image_to_wcs(input_image, input_region_file, output_region_file=None):
    if output_region_file is None:
        output_region_file = input_region_file.split('.reg')[0] + '_wcs.reg'
        
    alpha1, alpha2, header = find_crota(input_image)
    xc,yc,sma,smb,PA = crop_galaxy_image.read_region(input_region_file)
    
    w = wcs.WCS(header)
    pixcrd = np.array([[xc,yc]], np.float_)
    world = w.wcs_pix2world(pixcrd, 1)
    RA = world[0][0]
    DEC = world[0][1]

    
    pix2sec,note = resolution(header)
    pix2sec = float(pix2sec)
    
    sma = sma*pix2sec
    smb = smb*pix2sec
    
    PA = PA - alpha1

    f_reg = open(output_region_file, 'w')
    f_reg.write('%s\n' % ('fk5') )
    f_reg.write("ellipse(%f,%f,%f\",%f\",%f)\n" % (RA, DEC,sma, smb, PA))
    f_reg.close()
        
        

def convert_ellipse_wcs_to_image(input_image, input_region_file, output_region_file=None):        
    if output_region_file is None:
        output_region_file = input_region_file.split('.reg')[0] + '_ima.reg'
        
    alpha1, alpha2, header = find_crota(input_image)
    RA,DEC,sma,smb,PA = crop_galaxy_image.read_region(input_region_file, coord_format='fk5')
    
    w = wcs.WCS(header)
    world = np.array([[RA,DEC]], np.float_)
    
    pixcrd = w.wcs_world2pix(world, 1)

    xc = pixcrd[0][0]
    yc = pixcrd[0][1]
    
    pix2sec,note = resolution(header)
    pix2sec = float(pix2sec)
    
    sma = sma/pix2sec
    smb = smb/pix2sec
    
    PA = PA + alpha1

    f_reg = open(output_region_file, 'w')
    f_reg.write('%s\n' % ('image') )
    f_reg.write("ellipse(%f,%f,%f,%f,%f)\n" % (xc, yc,sma, smb, PA))
    f_reg.close()

    
    
#convert_ellipse_image_to_wcs('new-image.fits', 'test.reg')
#convert_ellipse_wcs_to_image('sky_subtr_galaxy_rot.fits', 'test_wcs.reg')
