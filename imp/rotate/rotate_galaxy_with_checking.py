#! /usr/bin/env python

import sys
import subprocess
import os
import shutil
import argparse
from os.path import exists
from math import hypot, pi, log10
import numpy  as np
from numpy import ma, mean, std, zeros_like, where, arange, exp, zeros
from numpy import copy as npcopy
from numpy import sum as npsum
from pylab import plot, show, savefig, xlabel, ylabel, vlines, clf
import astropy.io.fits as pyfits
from scipy.odr.odrpack import *
import warnings

import rotate_image


warnings.filterwarnings("ignore")

FNULL = open(os.devnull, 'w')


def main(input_images, xc, yc, PA=0., output_images = None, ask_user=False,cropping=True):

    optPosAng = PA
    print("Input position angle = %1.1f" % (optPosAng))

    
    if output_images is None:
        output_images = []
        for k in range(len(input_images)):
            output_images.append(input_images[k].split('.fits')[0]+'_rot.fits')


    xCenRot, yCenRot = rotate_image.main(input_images[0], optPosAng, xc=xc, yc=yc, output_image=output_images[0], hdu_inp=0, cval=float('nan'),cropping=cropping)
    
    deltaPosAng = 0.
    if ask_user:
        # Now lets show one of rotated images, so user can see if rotation was ok
        # and make some corrections if it is nessesary
        deltaPosAng = 0.
        answer='no'
        while answer=='no':
                subprocess.call("ds9 %s -scale histequ %s -scale histequ" % (input_images[0], output_images[0]), shell=True)  
                answer = 'yes'
                answer = str(raw_input('Are you happy with this (YES/no)?') or 'yes')
                if answer=='no':
                    inputString = raw_input("Correction to position angle (deg): ").strip()
                    if inputString:
                        deltaPosAng = float(inputString)
                    else:
                        deltaPosAng = 0.0

                    xCenRot, yCenRot = rotate_image.main(input_images[0], optPosAng-deltaPosAng, xc=xc, yc=yc, output_image=output_images[0], hdu_inp=0, cval=float('nan'),cropping=cropping)
  

    # Rotate image using the mean position angle
    for k in range(len(input_images)):
        xCenRot, yCenRot = rotate_image.main(input_images[k], optPosAng-deltaPosAng, xc=xc, yc=yc, output_image=output_images[k], hdu_inp=0, cval=float('nan'),cropping=cropping)




    print('New coordinates of the center: %f,%f' % (xCenRot, yCenRot))
    print("Output position angle = %1.1f" % (optPosAng-deltaPosAng))
    return xCenRot, yCenRot, optPosAng-deltaPosAng


 
#main(['NGC3628.phot.1_nonan.fits','NGC3628.1.finmask_nonan.fits','NGC3628_sigma2014.fits','new_mask.fits'], 1120.95, 1441.58, I_DN_min=2.45, I_DN_max=None, output_images = None)  

#main(['NGC4302.phot.1_nonan.fits','NGC4302.1.finmask_nonan.fits','NGC4302_sigma2014.fits'], 879.0, 1227.0, I_DN_min=0.1, I_DN_max=None, output_images = None)  

#main(['NGC891_coadd.fits','NGC891_coadd_mask.fits','NGC891_coadd_sigma.fits'], 1285.0, 1585.0, I_DN_min=0.2, I_DN_max=None, output_images = None) 

#main(['model.fits','mask.fits','sigma.fits'], 879.0, 1227.0, I_DN_min=1.6, I_DN_max=None, output_images = None) 

'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Rebinning")
    parser.add_argument("input_images", help="Input images, separated by comma: first galaxy image, then mask image, then others.")    
    parser.add_argument("xc", help="X-center of the galaxy", type=float)
    parser.add_argument("yc", help="Y-center of the galaxy", type=float)
    parser.add_argument("--output_images", help="Optional: Output images, separated by comma: first galaxy image, then mask image, then others.", type=str, default=None) 
    parser.add_argument("--min_DN", help="Optional: Outermost isophote. By default 5 sigma.", type=float, default=None)
    parser.add_argument("--max_DN", help="Optional: Innermost isophote. By default None", type=float, default=None)
    args = parser.parse_args()

    input_images = args.input_images
    xc = args.xc
    yc = args.yc
    output_images = args.output_images
    I_DN_min = args.min_DN
    I_DN_max = args.max_DN
    
    input_images = input_images.split(',')
    output_images = output_images.split(',')
    main(input_images, xc, yc, I_DN_min=I_DN_min, I_DN_max=I_DN_max, output_images = output_images)
'''