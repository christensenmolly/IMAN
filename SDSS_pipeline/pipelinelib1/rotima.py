#! /usr/bin/env python

from math import *
import sys
import os

from ctypes import *
import warnings
import numpy as np
import math
try:
  import astropy.io.fits as pyfits
  from astropy import wcs
except:
  warnings.warn("Astropy is not installed! No WCS has been added to the header!")

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

def rotPoint(x, y, x0, y0, angle):
    """ performs rotation of a point (x,y) around point (x0,y0)
        by angle radians"""
    x1 = cos(angle)*(x-x0) - sin(angle)*(y-y0) + x0
    y1 = sin(angle)*(x-x0) + cos(angle)*(y-y0) + y0
    return x1, y1


def crotima(name, outName, xOrig, yOrig, angle, layer=1, set_wcs=True):
    #name = b'%s' % (name)
    #outName = b'%s' % (outName)
    dllpath = "%s/crotima.so" % (os.path.dirname(__file__))
    crotlib = CDLL(dllpath)
    crotlib.crotima_crop.restype = c_int
    newCoordsOfOrigin = (c_double * 2)(0.0, 0.0)
    try:
        r = crotlib.crotima(bytes(name,encoding='utf8'), bytes(outName,encoding='utf8'), c_double(xOrig), c_double(yOrig),
                        c_double(-angle), c_int(layer), newCoordsOfOrigin)
    except:
        r = crotlib.crotima(name, outName, c_double(xOrig), c_double(yOrig),
                        c_double(-angle), c_int(layer), newCoordsOfOrigin)
    xNew = newCoordsOfOrigin[0]
    yNew = newCoordsOfOrigin[1]
    if set_wcs==True:
        try:
            add_wcs(name,outName,xOrig,yOrig,xNew,yNew,angle)
        except:
            print('Failed! WCS was not added!')
    return xNew, yNew

def cstretchima(fitsName, outName, scale, workHDU=1):
    dllpath = "%s/crotima.so" % (os.path.dirname(__file__))
    crotlib = CDLL(dllpath)
    crotlib.stretchima.restype = c_int
    r = crotlib.stretchima(fitsName, outName, c_double(scale), c_int(workHDU))
    return r

def cdeproject(infile, outfile, xCen, yCen, incl, posang, workHDU):
    dllpath = "%s/crotima.so" % (os.path.dirname(__file__))
    crotlib = CDLL(dllpath)
    crotlib.stretchima.deproject = c_int
    r = crotlib.deproject(infile, outfile, c_double(xCen), c_double(yCen), c_double(incl), c_double(posang), c_int(workHDU))
    return r



if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage:")
        print("./rotima.py in out xCen yCen angle")
        exit()

    fileToRotate = sys.argv[1]
    outName = sys.argv[2]
    xOrig = float(sys.argv[3])
    yOrig = float(sys.argv[4])
    angle = float(sys.argv[5])
    coords = crotima(name=fileToRotate,
                     outName=outName,
                     xOrig=xOrig,
                     yOrig=yOrig,
                     angle=angle)

    print("New coordinates:", coords)
    
