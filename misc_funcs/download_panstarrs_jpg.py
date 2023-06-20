#!/usr/bin/python
# -*- coding:  cp1251 -*-

import sys
import math
import numpy as np
from numpy import *
from pylab import *
import os
import urllib
from time import sleep
from panstamps.downloader import downloader
from panstamps.image import image
import subprocess
import shutil
import glob
import argparse
FNULL = open(os.devnull, 'w')



def download_panstarrs(RA, DEC, width=None, output_file=None):
    if output_file is not None:
        file = output_file.split('/')[-1]
        output_folder = output_file.split(file)[0]
        if output_folder=='':
            output_folder = './'
    else:
        output_folder = './'
    #print(output_folder)
    #exit()
    # width - in arcmin
    if width is None:
        width = 2
    else:
        width = int(math.ceil(width))
    callString = 'panstamps -FJc --width=%i --filters=gri --downloadFolder=%s stack %.6f %.6f' % (width, output_folder, RA, DEC)
    subprocess.call(callString, shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)

    try:
        file = glob.glob(output_folder+'color__ra%.6f_dec%.6f_*.jpeg' % (RA,DEC))[0]
    except:
        print('Image not found!')
        return 1
    if output_file is not None:
        shutil.move(file, output_file)
    print('Done!')
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download jpeg image of an object with (RA,DEC)")
    parser.add_argument("RA", help="Right ascension (in deg, decimal)", type=float)
    parser.add_argument("DEC", help="Declination (in deg, decimal)", type=float)
    parser.add_argument("--w", nargs='?', const=1, help="Optional: Width of the image (in arcmin)", type=float, default=2.)
    parser.add_argument("--o", help="Optional: Output jpeg file", type=str, default=None) 

    
    args = parser.parse_args()

    RA = args.RA
    DEC = args.DEC
    width = args.w
    output_file = args.o

    download_panstarrs(RA, DEC, width=width, output_file=output_file)
