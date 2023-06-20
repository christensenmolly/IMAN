#!/usr/bin/python
# DESCRIPTION:
# Script to flip image in the left/right or up/down direction.
# MINIMAL USAGE: python flip_image.py [input_image]
# EXAMPLE: python flip_image.py HCG041_11Dec2018.fits

# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import argparse
import warnings

warnings.filterwarnings("ignore")



def main(input_image, output_image=None, flip='ud'):
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    ny,nx = np.shape(data)
    header = hdulist[0].header
    
    
    if flip=='lr':
        new_data = np.fliplr(data)
    elif flip=='ud':
        new_data = np.flipud(data)
        
    hdu = pyfits.PrimaryHDU(new_data, header)
    if output_image is None:
        output_image = input_image.split('.fits')[0]+'_%s.fits' % (flip)
    hdu.writeto(output_image, clobber=True)    



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Flip image in the left/right or up/down direction")
    parser.add_argument("input_image", help="Input fits image")
    
    parser.add_argument("--output_image", help="Output fits image", type=str, default=None)  
    parser.add_argument("--flip", help="Direction of flipping (ud for up/down and lr for left/right). Dafault: ud", type=str, default='ud')  

    args = parser.parse_args()

    input_image = args.input_image
    output_image = args.output_image
    flip = args.flip



    main(input_image, output_image=output_image, flip=flip)