#!/usr/bin/python
# DESCRIPTION:
# Script to flip and median image 
# MINIMAL USAGE: python flip_mean_image.py [input_image]
# EXAMPLE: python flip_mean_image.py HCG041_11Dec2018.fits
# NOTE: The image of the object should be centred! The major axis should be horizontal.

# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import argparse
import warnings
import os
import sys

warnings.filterwarnings("ignore")


LOCAL_DIR = "/imp"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/rotate'))
import rotate_image

def main(input_image, output_image=None, PA=90., layer=0):

    hdulist = pyfits.open(input_image)
    data = hdulist[layer].data
    ny,nx = np.shape(data)
    header = hdulist[layer].header
    #if layer!=0:
    #    data = data[layer]
    
    if PA!=90.:
        angle = 90-PA
        rotate_image.main(input_image, -angle, xc=None, yc=None, output_image='tmp_rot.fits', hdu_inp=layer, cval=float('nan'), cropping=False, verbosity=True)
        hdulist = pyfits.open('tmp_rot.fits')
        data = hdulist[0].data    
    
    data1 = np.copy(data)
    
    data2 = np.fliplr(data)
    
    data3 = np.flipud(data)
    
    data4 = np.flipud(data2)
    
    #data_out = np.zeros((ny,nx))
    #for k in range(ny):
    #    for i in range(nx):
    #        data_out[k,i] = np.median([data1[k,i],data2[k,i],data3[k,i],data4[k,i]])
    
    data_out = np.median([data1,data2,data3,data4], axis=0)
    
    hdu = pyfits.PrimaryHDU(data_out, header)
    if output_image is None:
        output_image = input_image.split('.fits')[0]+'_%s.fits' % ('flipmean')
    hdu.writeto(output_image, clobber=True) 
    if PA!=90.:
        os.remove('tmp_rot.fits')




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Flip-mean image")
    parser.add_argument("input_image", help="Input fits image")
    
    parser.add_argument("--o", help="Output fits image", type=str, default=None)  
    parser.add_argument("--pa", help="Position angle of the galaxy: CCW,y: 0 is up, 90 is left", type=float, default=90.)
    parser.add_argument("--layer", help="Layer in the cube image", type=int, default=0)

    args = parser.parse_args()

    input_image = args.input_image
    output_image = args.o
    PA = args.pa
    layer = args.layer

    main(input_image, output_image=output_image, PA=PA, layer=layer)
