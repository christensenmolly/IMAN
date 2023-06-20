#!/usr/bin/python
# -*- coding:  cp1251 -*-
# EXAMPLE: python3 ~/MEGA/MyPrograms/IMAN/misc_funcs/download_fits.py NGC891 --s 'DSS1 Blue' --r 5 --o NGC891_dss.fits
# python3 ~/MEGA/MyPrograms/IMAN/misc_funcs/download_fits.py M81 --s 'SDSSr' --r 10 --o M81.fits --scale 2.0
import sys
import math
import numpy as np
from numpy import *
from pylab import *
import os
import urllib
from time import sleep
import subprocess
import shutil
import glob
import argparse
from astroquery.skyview import SkyView
import astropy.units as u
FNULL = open(os.devnull, 'w')



def download_fits(position, coordinates=None, radius=None, output_files=['g.fits'], survey=['SDSSg'], scale=1.):
    if output_files is None:
        output_files = []
        for surv in survey:
            output_files.append('%s_%s.fits' % (position,surv))

    # radius - in arcmin
    Radius = radius
    if radius is None:
        radius = 2.*1./60. * u.deg
    else:
        radius = 2.*radius/60. * u.deg
    # Now radius - in degrees

    pixels_x = int(math.ceil(2.*Radius*60./scale))
    pixels_y = pixels_x
    
    paths = SkyView.get_images(position=position, coordinates=coordinates, survey=survey, radius = radius, pixels='%i,%i' % (pixels_x,pixels_y))
    for k in range(len(paths)):
        paths[k].writeto(output_files[k], clobber=True)  

    print('Done!')
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download fits image of an object with (RA,DEC)")
    parser.add_argument("position", help="Name of the object (SIMBAD) or coordinates", type=str)
    parser.add_argument("--c", help="Coordinate system", type=str, default=None)
    parser.add_argument("--r", nargs='?', const=1, help="Optional: Radius of the image (in arcmin)", type=float, default=1.)
    parser.add_argument("--o", help="Optional: Output fits files, separated by comma (should be placed in the same directory!)", type=str, default=None)
    parser.add_argument("--scale", nargs='?', const=1, help="Optional: Scale of the output image (arcsec/pix), default 1.0", type=float, default=1.)
    parser.add_argument("--s", help="Optional: Surveys from the lust, default SDSSg", type=str, default='SDSSg')    

    
    args = parser.parse_args()

    position = args.position
    coordinates = args.c
    radius = args.r
    output_files = args.o
    survey = args.s
    scale = args.scale
    
    if output_files is not None:
        output_files = output_files.split(',')
    
    survey = survey.split(',')

    download_fits(position, coordinates=coordinates, radius=radius, output_files=output_files, survey=survey, scale=scale)
