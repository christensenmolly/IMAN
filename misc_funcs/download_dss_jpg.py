#!/usr/bin/python
# -*- coding:  cp1251 -*-
# Read: http://aladin.u-strasbg.fr/java/FAQ.htx
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
import argparse
FNULL = open(os.devnull, 'w')


def download_dss(object_names, width, folder=None):
    # width - in arcmin
    p = subprocess.Popen(['java -jar /home/amosenko/Downloads/Aladin/Aladin.jar -nogui'],
                        shell=True,
                        stdin=subprocess.PIPE)#, stdout=FNULL, stderr=subprocess.STDOUT)
    #p.stdin.write('grid on\n')
    for k in range(len(object_names)):
        obj = object_names[k]
        #p.stdin.write('reset; get aladin '+obj+'\n')#';get STScI(dss1); get ESO(dss2)\n')
        #p.stdin.write('zoom '+str(int(math.ceil(width[k]/2.)))+'arcmin;rgb;save '+obj+'.jpg\n')
        p.stdin.write("get hips(P/DSS2/color) " + obj + "\n")
        p.stdin.write("zoom "+str(int(math.ceil(width[k])))+"arcmin;save " + obj + ".jpg\n")
        

    p.stdin.write('quit\n')
    p.wait()
    for obj in (object_names):
        if folder is not None:
            shutil.move('/home/amosenko/Downloads/Aladin/%s.jpg' % (obj), folder+'/%s.jpg' % (obj))

download_dss(['NGC891'], [12.], folder='/home/amosenko/MyCurrentWork/Edge_on_S4G/test')
            
'''            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download jpeg image of an object with (RA,DEC)")
    parser.add_argument("RA", help="Right ascension (in deg, decimal)", type=float)
    parser.add_argument("DEC", help="Declination (in deg, decimal)", type=float)
    parser.add_argument("--w", nargs='?', const=1, help="Optional: Width of the image (in arcmin)", type=float, default=2.)
    parser.add_argument("--o", help="Optional: Output jpeg file", type=str, default=None) 

    
    args = parser.parse_args()

    RA = args.RA
    DEC = args.dec
    width = args.w
    output_file = args.o

    download_panstarrs(RA, DEC, width=width, output_file=output_file)
'''
