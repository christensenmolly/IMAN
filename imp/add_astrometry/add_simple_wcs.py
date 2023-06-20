#!/usr/bin/python
# DESCRIPTION:
# Script to add rough wcs to an image if we know pixel coordinates of the reference point (x_ref, y_ref)
# along with its WCS coordinates (RA_ref, DEC_ref), pixelscale (CDELT1, CDELT2 -in deg) and rotation angle CROTA2 (otherwise 0).
# NOTE: For adding robust wcs you should then feed the output image into nova.astrometry.net.
# MINIMAL USAGE: python add_simple_wcs.py [input_image] [x_ref] [y_ref] [RA_ref] [DEC_ref] [cdelt1] [cdelt2]
# EXAMPLE: python add_simple_wcs.py HCG041_11Dec2018_ud.fits 1917 2013 149.456624 45.246035 -0.000233333 0.000233333

# Import the necessary modules
from astropy.io import fits as pyfits
import argparse
from astropy import coordinates
import astropy.units as u
from astropy import wcs



def main(input_image, x_ref, y_ref, RA_ref, DEC_ref, cdelt1, cdelt2, angle=None, output_image=None):
    if output_image is None:
        output_image = input_image.split('.fits')[0] + '_wcs.fits'
    
    hdulist = pyfits.open(input_image)
    referenceheader = pyfits.getheader(input_image, 0)
    outframe = pyfits.getdata(input_image, 0)


    referenceheader['EQUINOX'] = 2.000000000000E+03
    referenceheader['RADECSYS'] = 'FK5'
    
    referenceheader['CTYPE1'] = 'RA---TAN'
    referenceheader['CUNIT1'] = 'deg'    
    referenceheader['CRVAL1'] = RA_ref
    referenceheader['CRPIX1'] = x_ref

    
    referenceheader['CTYPE2'] = 'DEC--TAN'	#### Was:  'DEC---TAN'
    referenceheader['CUNIT2'] = 'deg'  
    referenceheader['CRVAL2'] = DEC_ref
    referenceheader['CRPIX2'] = y_ref

    referenceheader['CDELT1'] = cdelt1
    referenceheader['CDELT2'] = cdelt2
    
    if angle is not None:
        referenceheader['CROTA2'] = angle        
    
    hdu = pyfits.PrimaryHDU(outframe, referenceheader)
    hdu.writeto(output_image, clobber=True)
    
    



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add world coordinate system to image")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("x_ref", help="Input x-coordinate (pix)", type=float)  
    parser.add_argument("y_ref", help="Input y-coordinate (pix)", type=float)      

    parser.add_argument("RA_ref", help="Input RA-coordinate (deg) of x_ref", type=float)  
    parser.add_argument("DEC_ref", help="Input DEC-coordinate (deg) of y_ref", type=float)    

    parser.add_argument("cdelt1", help="Input CDELT1 (in deg)", type=float)  
    parser.add_argument("cdelt2", help="Input CDELT2 (in deg)", type=float)  

    parser.add_argument("--angle", help="Input CROTA2 (in deg)", type=float, default=None)  
    
    parser.add_argument("--output_image", help="Output fits image", type=str, default=None)  

    args = parser.parse_args()

    input_image = args.input_image
    x_ref = args.x_ref
    y_ref = args.y_ref
    RA_ref = args.RA_ref
    DEC_ref = args.DEC_ref
    cdelt1 = args.cdelt1
    cdelt2 = args.cdelt2
    angle = args.angle
    output_image = args.output_image


    main(input_image, x_ref, y_ref, RA_ref, DEC_ref, cdelt1, cdelt2, angle=angle, output_image=output_image)