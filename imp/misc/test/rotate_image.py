#!/usr/bin/python
# DESCRIPTION:
# Script to rebin image to a reference image.
# MINIMAL USAGE: python rebin_image.py [input_image] [angle]

import os
import numpy as np
from scipy import ndimage
from scipy import misc
import argparse

# Import astronomical modules
import astropy.io.fits as pyfits

# This script rotates a FITS image around the center of the frame over an arbitrary angle.
# It updates the header accordingly.
# Written by Sebastien Viaene and modified by Aleksandr Mosenkov, see PTS/pts/magic/tools/rotation.py

def main(input_image, angle, output_image=None, hdu_inp=0, cval='nan'):
        # Load in image
        hdu = pyfits.open(input_image)
        data = hdu[hdu_inp].data
        header = hdu[hdu_inp].header

        # Get the rotation angle
        #angle = header["CROTA2"]

        # Rotate the header
        new_header = rotate_header(header, angle)

        # Rotate the image
        new_image = rotate_frame(data, angle, cval)

        # Write out rotated image.
        rot_hdu = pyfits.PrimaryHDU(new_image, new_header)
        if output_image is None:
            output_image = input_image.split('.fits')[0] + '_rot.fits'
        rot_hdu.writeto(output_image, clobber=True)

def rotate_header(header, angle):
    
    """
    This function rotates the header
    """
    
    new_header = header
    
    # Check if a rotation matrix element exists
    matrix = True
    try:
        cd1_1 = np.float(header["CD1_1"])
    except:
        matrix = False

    theta = angle * np.pi / 180.
    rot_matrix = np.array( [ [ np.cos(theta), np.sin(theta)],
                            [-1.*np.sin(theta), np.cos(theta)] ] )

    center = np.array([(header['NAXIS1'] - 1)/2., (header['NAXIS2'] - 1)/2. ])
    
    try:
        crpix = np.array([header['CRPIX1'], header['CRPIX2']])
    except:
        crpix = center
        header.append( pyfits.Card('CRPIX1', crpix[0], 'Reference pixel on this axis'), end=True)
        header.append( pyfits.Card('CRPIX2', crpix[1], 'Reference pixel on this axis'), end=True)

    ncrpix = (crpix-1-center).dot(rot_matrix.T) + 1
    ncrpix += center

    new_header["CRPIX1"] = ncrpix[0]
    new_header["CRPIX2"] = ncrpix[1]

    if matrix:

        try:
            cd1_2 = np.float(header["CD1_2"])
        except:
            cd1_2 = 0.
            header.append(pyfits.Card('CD1_2', cd1_2, 'Rotation matrix element 1_2'), end=True)

        try:
            cd2_1 = np.float(header["CD2_1"])
        except:
            cd2_1 = 0.
            header.append(pyfits.Card('CD2_1', cd2_1, 'Rotation matrix element 2_1'), end=True)

        try:
            cd2_2 = np.float(header["CD2_2"])
        except:
            cd2_2 = 0.
            header.append(pyfits.Card('CD2_2', cd2_2, 'Rotation matrix element 2_2'), end=True)

        cd = np.array([[cd1_1, cd1_2], [cd2_1, cd2_2]])

        newcd = rot_matrix.dot(cd)
        new_header["CD1_1"] = newcd[0,0]
        new_header["CD1_2"] = newcd[0,1]
        new_header["CD2_1"] = newcd[1,0]
        new_header["CD2_2"] = newcd[1,1]
    
    else:
        
        #try:
        #    new_header["CROTA1"] = -1.*angle
        #except:
        #    new_header.append(pyfits.Card('CROTA1', -1.*angle, 'Rotation parameter'), end=True)

        #try:
        #    new_header["CROTA2"] = -1.*angle
        #except:
        #   new_header.append( pyfits.Card('CROTA2', -1.*angle, 'Rotation parameter'), end=True)

        new_header["CROTA2"] = -angle #0.0

    return new_header

def rotate_frame(frame, angle, cval='nan'):
    
    # Perform the image rotation and update the fits header
    #frame[np.isnan(frame)] = 0.0
    new_frame = ndimage.interpolation.rotate(frame, angle, reshape=False, order=3, mode='constant', cval=float(cval))
    
    #new_frame = misc.imrotate(frame, angle, interp="bilinear")
    
    # Return the rotated frame
    return new_frame





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Rotation")
    parser.add_argument("input_image", help="Input image")    
    parser.add_argument("angle", help="Angle [deg]", type=float)
    parser.add_argument("--output_image", help="Optional: Output image", type=str, default=None) 
    parser.add_argument("--hdu", help="Optional: HDU layer", type=int, default=0)
    parser.add_argument("--cval", help="Optional: Value used for points outside the boundaries of the input, default is nan", type=str, default='nan')
    args = parser.parse_args()

    input_image = args.input_image
    angle = args.angle
    output_image = args.output_image
    hdu_inp = args.hdu
    cval = args.cval
    
    main(input_image, angle, output_image, hdu_inp=hdu_inp, cval=cval)
