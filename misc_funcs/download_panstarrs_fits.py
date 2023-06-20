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



def download_panstarrs(RA, DEC, width=None, output_files=None, filters='gri'):
    if output_files is not None:
        file = output_files[0].split('/')[-1]
        output_folder = output_files[0].split(file)[0]
        if output_folder=='':
            output_folder = './'
    else:
        output_folder = './'

    # width - in arcmin
    if width is None:
        width = 2
    else:
        width = int(math.ceil(width))
    callString = 'panstamps --width=%i --filters=%s --downloadFolder=%s stack %s %s' % (width, filters, output_folder, str(RA), str(DEC))
    subprocess.call(callString, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)


    files = []
    for band in filters:
            files.append(glob.glob(output_folder+'stack_%s_ra%s*_dec%s*_*.fits' % (band,str(RA),str(DEC)))[0])
    
    if files==[]:
        print('Image lis is empty!')
        return 1
    
    if output_files is not None:
        for k in range(len(files)):
            shutil.move(files[k], output_files[k])

    print('Done!')
    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download fits image of an object with (RA,DEC)")
    parser.add_argument("RA", help="Right ascension (in deg, decimal)", type=str)
    parser.add_argument("DEC", help="Declination (in deg, decimal)", type=str)
    parser.add_argument("--w", nargs='?', const=1, help="Optional: Width of the image (in arcmin)", type=float, default=2.)
    parser.add_argument("--o", help="Optional: Output fits files, separated by comma (should be placed in the same directory!)", type=str, default=None)
    parser.add_argument("--f", help="Optional: Filters to download, default gri", type=str, default='gri')    

    
    args = parser.parse_args()

    RA = args.RA
    DEC = args.DEC
    width = args.w
    output_files = args.o
    filters = args.f
    
    if output_files is not None:
        output_files = output_files.split(',')

    download_panstarrs(RA, DEC, width=width, output_files=output_files, filters=filters)
