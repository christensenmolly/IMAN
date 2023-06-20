#!/usr/bin/python
# DESCRIPTION:
# Script to rotate an image around the center of the frame over an arbitrary angle. It updates the header accordingly.
# Written by Sebastien Viaene and modified by Aleksandr Mosenkov, see PTS/pts/magic/tools/rotation.py
# MINIMAL USAGE: python rebin_image.py [input_image] [angle]

import os
import numpy as np
from scipy import ndimage
from scipy import misc
import argparse
import astropy.io.fits as pyfits
import math
import numpy as np
try:
  import astropy.io.fits as pyfits
  from astropy import wcs
except:
  warnings.warn("Astropy is not installed! No WCS has been added to the header!")
  
# This script rotates a FITS image around the center of the frame over an arbitrary angle.



def ds9_to_np(x):
    '''
    Function to convert ds9 pixel coordinates to numpy format, i.e
    0.5 pix is 0, 1.5 pix is 1
    '''
    return int(round(x)) - 1

def get_xy_rotation_and_scale(header):
    """
    CREDIT: See IDL code at
    http://www.astro.washington.edu/docs/idl/cgi-bin/getpro/library32.html?GETROT
    """

    def calc_from_cd(cd1_1, cd1_2, cd2_1, cd2_2):

        # TODO: Check if first coordinate in CTYPE is latitude
        # if (ctype EQ 'DEC-') or (strmid(ctype, 1) EQ 'LAT')  then $
        #    cd = reverse(cd,1)

        det = cd1_1*cd2_2 - cd1_2*cd2_1
        if det < 0:
            sgn = -1
        else:
            sgn = 1
        ## if det > 0:
        ##     raise ValueError("Astrometry is for a right-handed coordinate system")

        if (cd2_1 == 0.0) or (cd1_2 == 0.0):
            # Unrotated coordinates?
            xrot = 0.0
            yrot = 0.0
            cdelt1 = cd1_1
            cdelt2 = cd2_2
        else:
            xrot = math.atan2(sgn * cd1_2, sgn * cd1_1)
            yrot = math.atan2(-cd2_1, cd2_2)

            cdelt1 = sgn * math.sqrt(cd1_1**2 + cd1_2**2)
            cdelt2 = math.sqrt(cd1_1**2 + cd2_1**2)

        return xrot, yrot, cdelt1, cdelt2

    def calc_from_crota():
        try:
            crota1 = float(header['CROTA1'])
            xrot = crota1
        except KeyError:
            xrot = None

        try:
            crota2 = float(header['CROTA2'])
            yrot = crota2
        except KeyError:
            yrot = 0.0

        if xrot is None:
            xrot = yrot

        cdelt1 = float(header.get('CDELT1', 1.0))
        cdelt2 = float(header.get('CDELT2', 1.0))

        return xrot, yrot, cdelt1, cdelt2

    # 1st, check for presence of PC matrix
    try:
        pc1_1 = header['PC1_1']
        try:
            pc1_2 = header['PC1_2']
        except:
            pc1_2 = 0.
        try:
            pc2_1 = header['PC2_1']
        except:
            pc2_1 = 0.
        pc2_2 = header['PC2_2']

        cdelt1 = float(header['CDELT1'])
        cdelt2 = float(header['CDELT2'])

        cd1_1, cd1_2 = pc1_1 * cdelt1, pc1_2 * cdelt1
        cd2_1, cd2_2 = pc2_1 * cdelt2, pc2_2 * cdelt2

        xrot, yrot, cdelt1p, cdelt2p = calc_from_cd(pc1_1, pc1_2,
                                                    pc2_1, pc2_2)
        cdelt1 = cdelt1p
        cdelt2 = cdelt2p

    except KeyError:
        # 2nd, check for presence of CD matrix
        try:
            cd1_1 = header['CD1_1']
            cd1_2 = header['CD1_2']
            cd2_1 = header['CD2_1']
            cd2_2 = header['CD2_2']
            xrot, yrot, cdelt1, cdelt2 = calc_from_cd(cd1_1, cd1_2,
                                                      cd2_1, cd2_2)

        except KeyError:
            # 3rd, check for presence of CROTA keyword
            #  (or default is north=up)
            xrot, yrot, cdelt1, cdelt2 = calc_from_crota()

    xrot, yrot = np.rad2deg(xrot), np.rad2deg(yrot)

    return ((xrot, yrot), (cdelt1, cdelt2))


def add_wcs_func(input_image, output_image, xc, yc, angle):  
    referenceheader = pyfits.getheader(input_image)
    w = wcs.WCS(referenceheader)
    
    pixcrd = np.array([[xc, yc]], np.float_)
    world_pixcrd = w.wcs_pix2world(pixcrd, 1)[0]  

    
    CRPIX1 = referenceheader['CRPIX1']
    CRPIX2 = referenceheader['CRPIX2']

    data = pyfits.getdata(output_image)
    header = pyfits.getheader(output_image)

    header['CTYPE1']  = referenceheader['CTYPE1']                                                         
    header['CTYPE2']  = referenceheader['CTYPE2'] 

    header['CRPIX1'] = xc
    header['CRPIX2'] = yc
    header['CRVAL1'] = world_pixcrd[0]
    header['CRVAL2'] = world_pixcrd[1]

    ((xrot, yrot), (cdelt1, cdelt2)) = get_xy_rotation_and_scale(referenceheader)
    header['CDELT1'] = cdelt1
    header['CDELT2'] = cdelt2
    header['CROTA2'] = xrot - angle    
    header['EQUINOX'] = float(2000.)

    hdu = pyfits.PrimaryHDU(data=data, header=header)
    hdu.writeto(output_image, clobber=True)     



def main(input_image, angle, xc=None, yc=None, output_image=None, hdu_inp=0, cval='nan'):
        # Load in image
        hdu = pyfits.open(input_image)
        data = hdu[hdu_inp].data
        header = hdu[hdu_inp].header

        ny,nx = np.shape(data)

        if xc is None:
            xc = nx/2.
            #xc = int(math.floor(xc))

        if yc is None:
            yc = ny/2.
            #yc = int(math.floor(yc))
            
        # Get the rotation angle
        #angle = header["CROTA2"]

        # Rotate the header
        #new_header = rotate_header(header, xc, yc, angle)

        # Rotate the image
        new_image = rotate_frame(data, xc, yc, angle, cval)

        rot_hdu = pyfits.PrimaryHDU(new_image)
        if output_image is None:
                output_image = input_image.split('.fits')[0] + '_rot.fits'
        rot_hdu.writeto(output_image, clobber=True)        
        
        add_wcs_func(input_image, output_image, xc, yc, angle)
     




def rotate_frame(img, xc, yc, angle, cval='nan'):    
    # Perform the image rotation and update the fits header
    #frame[np.isnan(frame)] = 0.0
    xc = ds9_to_np(xc)
    yc = ds9_to_np(yc)
    
    ny,nx = np.shape(img)

    padX = [nx-xc, xc]
    padY = [ny-yc, yc]
    imgP = np.pad(img, [padY, padX], 'constant')
    imgR = ndimage.interpolation.rotate(imgP, angle, reshape=False, order=1, mode='constant', cval=float(cval)) # order=1 ???
    
    # Return the rotated frame
    
    return imgR[padY[0] : -padY[1], padX[0] : -padX[1]]





if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Rebinning")
    parser.add_argument("input_image", help="Input image")    
    parser.add_argument("angle", help="Angle [deg]", type=float)
    parser.add_argument("--xc", help="Optional: X-coordinate, default the center of the X-axis", type=int, default=None)
    parser.add_argument("--yc", help="Optional: Y-coordinate, default the center of the Y-axis", type=int, default=None) 
    parser.add_argument("--output_image", help="Optional: Output image", type=str, default=None) 
    parser.add_argument("--hdu", help="Optional: HDU layer", type=int, default=0)
    parser.add_argument("--cval", help="Optional: Value used for points outside the boundaries of the input, default is nan", type=str, default='nan')
    args = parser.parse_args()

    input_image = args.input_image
    angle = args.angle
    xc = args.xc
    yc = args.yc
    output_image = args.output_image
    hdu_inp = args.hdu
    cval = args.cval
    
    main(input_image, angle, xc = xc, yc = yc, output_image=output_image, hdu_inp=hdu_inp, cval=cval)
