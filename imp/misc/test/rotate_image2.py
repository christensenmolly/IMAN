#!/usr/bin/python
# DESCRIPTION:
# Script to rotate an image around the center of the frame over an arbitrary angle. It updates the header accordingly.
# Written by Sebastien Viaene and modified by Aleksandr Mosenkov, see PTS/pts/magic/tools/rotation.py
# MINIMAL USAGE: python rotate_image.py [input_image] [angle] [coords]

import os
import numpy as np
from scipy import ndimage
from scipy import misc
import argparse
import astropy.io.fits as pyfits

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

    xrot, yrot = degrees(xrot), degrees(yrot)

    return ((xrot, yrot), (cdelt1, cdelt2))

def add_wcs(fileToRotate,outName,xOrig,yOrig,Xnew,Ynew,angle):
    referenceheader = pyfits.getheader(fileToRotate)
    w = wcs.WCS(referenceheader)
    
    pixcrd = np.array([[xOrig, yOrig]], np.float_)
    world = w.wcs_pix2world(pixcrd, 1)
    CRPIX1 = referenceheader['CRPIX1']
    CRPIX2 = referenceheader['CRPIX2']

    data = pyfits.getdata(outName)
    header = pyfits.getheader(outName)
    xOrig_w,yOrig_w = world[0,0],world[0,1]

    header['CTYPE1']  = referenceheader['CTYPE1']                                                         
    header['CTYPE2']  = referenceheader['CTYPE2'] 

    header['CRPIX1'] = Xnew
    header['CRPIX2'] = Ynew
    header['CRVAL1'] = xOrig_w
    header['CRVAL2'] = yOrig_w
    '''
    if 'CD1_1' in referenceheader:
	header['CD1_1'] = referenceheader['CD1_1']
	header['CD1_2'] = referenceheader['CD1_2']
	header['CD2_1'] = referenceheader['CD2_1']
	header['CD2_2'] = referenceheader['CD2_2']
	
    elif 'CDELT1' in referenceheader:
        header['CDELT1'] = referenceheader['CDELT1']
        header['CDELT2'] = referenceheader['CDELT2']
	header['CROTA2'] = float(referenceheader['CROTA2']) - angle
    '''
    ((xrot, yrot), (cdelt1, cdelt2)) = get_xy_rotation_and_scale(referenceheader)
    header['CDELT1'] = cdelt1
    header['CDELT2'] = cdelt2
    header['CROTA2'] = xrot - angle    
    header['EQUINOX'] = float(2000.)

    hdu = pyfits.PrimaryHDU(data=data, header=header)
    hdu.writeto(outName,clobber=True) 

'''
def rotate_frame(img, xc, yc, angle, cval='nan'):    
    # Perform the image rotation and update the fits header
    #frame[np.isnan(frame)] = 0.0
    ny,nx = np.shape(img)

    padX = [nx-xc, xc]
    padY = [ny-yc, yc]
    imgP = np.pad(img, [padY, padX], 'constant')
    imgR = ndimage.interpolation.rotate(imgP, angle, reshape=False, order=3, mode='constant', cval=float(cval))
    
    # Return the rotated frame
    return imgR[padY[0] : -padY[1], padX[0] : -padX[1]]
'''

def rotate_frame(image, xy, angle):
    im_rot = ndimage.interpolation.rotate(image, angle, reshape=False, order=3, mode='constant') 
    org_center = (np.array(image.shape[:2][::-1])-1)/2.
    rot_center = (np.array(im_rot.shape[:2][::-1])-1)/2.
    org = xy-org_center
    a = np.deg2rad(angle)
    new = np.array([org[0]*np.cos(a) + org[1]*np.sin(a),
            -org[0]*np.sin(a) + org[1]*np.cos(a) ])
    return im_rot, new+rot_center

def main(input_image, angle, xc=None, yc=None, output_image=None, hdu_inp=0, cval='nan'):
        # Load in image
        hdu = pyfits.open(input_image)
        data = hdu[hdu_inp].data
        header = hdu[hdu_inp].header

        ny,nx = np.shape(data)

        if xc is None:
            xc = nx/2.
            xc = int(math.floor(xc))

        if yc is None:
            yc = nx/2.
            yc = int(math.floor(yc))
            
    
        # Rotate the image
        #new_image = rotate_frame(data, xc, yc, angle, cval)
        new_image, (x1,y1) = rotate_frame(data, np.array([xc,yc]), angle)
        print(x1,y1)
        hdu = pyfits.PrimaryHDU(data=new_image)
        hdu.writeto(output_image, clobber=True) 

        if True:#set_wcs==True:
            try:
                add_wcs(input_image, output_image,xc,yc,x1,y1,angle)
            except:
                print('Failed! WCS was not added!')

main('galaxy.fits', 30., xc=30., yc=30., output_image='galaxy_rot.fits', hdu_inp=0, cval='nan')