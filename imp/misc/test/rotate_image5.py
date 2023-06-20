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

    return data, header  


def pad_image(input_image, xc=None, yc=None, output_image=None, hdu_inp=0):
        # Load in image
        hdu = pyfits.open(input_image)
        data = hdu[hdu_inp].data
        header = hdu[hdu_inp].header

        ny,nx = np.shape(data)

        if xc is None:
            xc = nx/2.

        if yc is None:
            yc = ny/2.
            
        # Pad image:
        data_pad,delta_x,delta_y = pad_frame(data, xc, yc)
        #ny_pad,nx_pad = np.shape(data_pad)
        header_pad = pad_header(header, delta_x, delta_y)


        rot_hdu = pyfits.PrimaryHDU(data_pad, header_pad)
        if output_image is None:
                output_image = input_image.split('.fits')[0] + '_pad.fits'
        rot_hdu.writeto(output_image, clobber=True)
        return output_image

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
        header.append( fits.Card('CRPIX1', crpix[0], 'Reference pixel on this axis'), end=True)
        header.append( fits.Card('CRPIX2', crpix[1], 'Reference pixel on this axis'), end=True)

    ncrpix = (crpix-1-center).dot(rot_matrix.T) + 1
    ncrpix += center

    new_header["CRPIX1"] = ncrpix[0]
    new_header["CRPIX2"] = ncrpix[1]

    if matrix:

        try:
            cd1_2 = np.float(header["CD1_2"])
        except:
            cd1_2 = 0.
            header.append(fits.Card('CD1_2', cd1_2, 'Rotation matrix element 1_2'), end=True)

        try:
            cd2_1 = np.float(header["CD2_1"])
        except:
            cd2_1 = 0.
            header.append(fits.Card('CD2_1', cd2_1, 'Rotation matrix element 2_1'), end=True)

        try:
            cd2_2 = np.float(header["CD2_2"])
        except:
            cd2_2 = 0.
            header.append(fits.Card('CD2_2', cd2_2, 'Rotation matrix element 2_2'), end=True)

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
        #    new_header.append(fits.Card('CROTA1', -1.*angle, 'Rotation parameter'), end=True)

        #try:
        #    new_header["CROTA2"] = -1.*angle
        #except:
        #   new_header.append( fits.Card('CROTA2', -1.*angle, 'Rotation parameter'), end=True)

        new_header["CROTA2"] = 0.0

    return new_header



def rotate_image(data, header, angle, cval='nan'):
        ny,nx = np.shape(data)

        xc = nx/2.
        yc = ny/2.
        
        '''
        # Provides worse WCS alignment!
        rot_data = ndimage.interpolation.rotate(data, angle, reshape=False, order=1, mode='constant', cval=float(cval))

        rot_hdu = pyfits.PrimaryHDU(rot_data)
        rot_hdu.writeto('rot_tmp.fits', clobber=True)          
        
        rot_data, header_rot = add_wcs_func('rot_tmp.fits', output_image, xc, yc, angle)
        os.remove('rot_tmp.fits')
        '''

        rot_data = ndimage.interpolation.rotate(data, angle, reshape=False, order=1, mode='constant', cval=float(cval))
        header_rot = rotate_header(header, angle)

        return rot_data, header_rot
        
    

def main(input_image, angle, xc=None, yc=None, output_image=None, hdu_inp=0, cval='nan'):       
        hdu = pyfits.open(input_image)
        data = hdu[hdu_inp].data
        header = hdu[hdu_inp].header
        new_header = header.copy()

        if xc is None:
            xc = nx/2.

        if yc is None:
            yc = ny/2.
            
        xc = ds9_to_np(xc)
        yc = ds9_to_np(yc)
        
        ny,nx = np.shape(data)

        padX = [nx-xc, xc]
        padY = [ny-yc, yc]

        # Pad image:
        data_pad,delta_x,delta_y = pad_frame(data, xc, yc)
        #ny_pad,nx_pad = np.shape(data_pad)
        header_pad = pad_header(new_header, delta_x, delta_y)
        
        #hdu = pyfits.PrimaryHDU(data_pad, header_pad)       
        #hdu.writeto('pad.fits', clobber=True)  
        #exit()

        # Rotate image
        #hdu_pad = pyfits.open('pad.fits')
        #data_pad = hdu_pad[hdu_inp].data
        #header_pad = hdu_pad[hdu_inp].header         
        rot_data, header_rot = rotate_image(data_pad, header_pad, angle, cval=cval) 

        hdu1 = pyfits.PrimaryHDU(rot_data, header_rot)       
        hdu1.writeto('rot.fits', clobber=True)  
        exit()
        
        # Crop back
        crop_data,crop_header = crop_back(rot_data, header_rot, padX, padY)
        
        # Save
        hdu = pyfits.PrimaryHDU(crop_data, crop_header)
        if output_image is None:
                output_image = input_image.split('.fits')[0] + '_rot.fits'
        
        hdu.writeto(output_image, clobber=True)             



def crop_back(data, header, padX, padY):       
    data_crop = data[padY[0] : -padY[1], padX[0] : -padX[1]]

    #header['CRPIX1'] = header['CRPIX1'] - max(padX)
    #header['CRPIX2'] = header['CRPIX2'] - max(padY)
    header.update(CRPIX1=header['CRPIX1'] - max(padX), CRPIX2=header['CRPIX2'] - max(padY))

    return data_crop, header         
        


def pad_header(hdr, delta_x, delta_y):
    #hdr['CRPIX1'] = hdr['CRPIX1'] + delta_x
    #hdr['CRPIX2'] = hdr['CRPIX2'] + delta_y
    hdr.update(CRPIX1=hdr['CRPIX1'] + delta_x, CRPIX2=hdr['CRPIX2'] + delta_y)

    return hdr



def pad_frame(img, xc, yc):
    xc = ds9_to_np(xc)
    yc = ds9_to_np(yc)
    
    ny,nx = np.shape(img)

    padX = [nx-xc, xc]
    padY = [ny-yc, yc]
    imgP = np.pad(img, [padY, padX], 'constant')
    return imgP,max(padX),max(padY)






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
