#! /usr/bin/env python

import subprocess
import math
from time import time, sleep
import urllib.request
import sys
import os
import glob
import shutil
import argparse
import numpy as np

#https://www.legacysurvey.org/viewer/jpeg-cutout?ra=190.1086&dec=1.2005&width=100&layer=galex&pixscale=0.3&bands=nf
#FUV: mAB = -2.5 x log10(CPS) + 18.82 # For 1.5 arcsec/pix
#NUV: mAB = -2.5 x log10(CPS) + 20.08 # For 1.5 arcsec/pix

def main(name, RA, DEC, R, bands='nf', pixscale=1.5):
    # R in arcmin
    print('Downloading...')
    RA = np.array(RA, float)
    DEC = np.array(DEC, float)
    R = np.array(R, float)
    
    
    
    # convert R to pixels
    
    for k in range(len(RA)):
        if True:
           RR = int(math.ceil(R[k]*60./pixscale))
           output_file = '%s.fits' % (name[k])
           url="http://legacysurvey.org/viewer/fits-cutout?ra=%f&dec=%f&width=%i&height=%i&layer=galex&pixscale=%.3f&bands=%s" % (RA[k], DEC[k], 2*RR, 2*RR, pixscale,bands) 
           print('%i %s' % (k+1, name[k]))
           urllib.request.urlretrieve(url, output_file)
        else:
            zz=1
    print('Done!')
    



#main(['PGC8299'], [2.17155*15.], [-10.32115], [2.*6*10**1.13/60.], [0.1], pixscale=0.262) 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Download from Legacysurvey") 
    parser.add_argument("ra", help="Right ascension [deg]", type=float)
    parser.add_argument("dec", help="Declination [deg]", type=float)
    parser.add_argument("width", help="Width of the image [arcmin]", type=float)
    parser.add_argument("--name", help="Optional: Name of the object", type=str, default=None)
    parser.add_argument("--bands", help="Optional: Bands to be downloaded, e.g. grz", type=str, default='grz') 
     

    args = parser.parse_args()

    ra = args.ra
    dec = args.dec
    width = args.width
    
    name = args.name
    bands = args.bands
    
    R = width/2.
    
    main([name], [ra], [dec], [R], bands, pixscale=1.5)
