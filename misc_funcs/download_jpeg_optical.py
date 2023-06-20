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
FNULL = open(os.devnull, 'w')


def download_sdss(RA, DEC, radius=None, output_file=None):
        # RA,DEC - decimals
        # radius - in arcmin
        if radius is None:
            scale=1.0
        else:
            scale=1.1*(2.*radius/10.)
            
        try:
           # scale=1 is for a box of 10'x10'
           pic="http://skyservice.pha.jhu.edu/DR14/ImgCutout/getjpeg.aspx?ra=%s&dec=%s&scale=%.1f&width=600&height=600&opt=G" % (str(RA), str(DEC), scale) # Was DR9
           f = urllib.urlopen(pic)
           image = f.read()
           fout = open(output_file, 'w')
           fout.write(image)
           fout.close()
        except IOError:
           sleep(1)
           print('Try again!')
           good=0
        finally:
           if int(os.path.getsize(output_file))>10000:
              good = 1
           else:
              good = 0
              sleep(1)
              try:
                 os.remove(output_file)
                 print('Removed')
                 exit()
              except:
                 print('Not found')
                 
                 
def download_dss(object_names, radius, folder=None):
    # radius - in arcmin
    p = subprocess.Popen(['java -jar /home/amosenko/Downloads/Aladin/Aladin.jar -nogui'],
                        shell=True,
                        stdin=subprocess.PIPE, stdout=FNULL, stderr=subprocess.STDOUT
                        )
    #p.stdin.write('grid on\n')
    for k in range(len(object_names)):
        obj = object_names[k]
        p.stdin.write('reset; get aladin '+obj+'\n')#';get STScI(dss1); get ESO(dss2)\n')
        p.stdin.write('zoom '+str(int(math.ceil(1.1*radius[k])))+'arcmin;rgb;save '+obj+'.jpg\n')


    p.stdin.write('quit\n')
    p.wait()
    for obj in (object_names):
        if folder is not None:
            shutil.move('/home/amosenko/Downloads/Aladin/%s.jpg' % (obj), folder+'/%s.jpg' % (obj))

def download_panstarrs(RA, DEC, radius=None, output_file=None):
    if output_file is not None:
        file = output_file.split('/')[-1]
        output_folder = output_file.split(file)[0]
        if output_folder=='':
            output_folder = './'
    else:
        output_folder = './'
    #print(output_folder)
    #exit()
    # radius - in arcmin
    if radius is None:
        width = 2.
    else:
        width = 2*int(math.ceil(radius))
    callString = 'panstamps -FJc --width=%i --filters=gri --downloadFolder=%s stack %.6f %.6f' % (width, output_folder, RA, DEC)
    subprocess.call(callString, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    if output_file is not None:
        ff = glob.glob(output_folder+'color__ra%.6f_dec%.6f_*.jpeg' % (RA,DEC))
        shutil.move(ff[0], output_file)
    print('Done!')

'''
def download_panstarrs(RA, DEC, radius=None, output_file=None):
    fitsPaths, jpegPaths, colorPath = downloader(
        log=log,
        settings=False,
        downloadDirectory=False,
        fits=False,
        jpeg=True,
        arcsecSize=600,
        filterSet='gri',
        color=True,
        singleFilters=True,
        ra="70.60271",
        dec="-21.72433",
        imageType="stack",  # warp | stack
        mjdStart=False,
        mjdEnd=False,
        window=False
    ).get()

    for j in jpegPaths:

        myimage = image(
            log=log,
            settings=False,
            imagePath=j,
            arcsecSize=120,
            crosshairs=True,
            transient=False,
            scale=True,
            invert=False,
            greyscale=False
        ).get()    

'''




#!/usr/bin/env python

RA = 6.011821
DEC = +16.486389
output_file = '/Users/mosenkov/Downloads/NGC891.jpg'
radius = 5.
download_panstarrs(RA, DEC, radius, output_file)

'''
download_dss(['NGC891'], [10.], folder='/home/amosenko/MyCurrentWork/Edge_on_S4G/test')

exit()
RA = 6.011821
DEC = +16.486389
output_file = 'NGC891.jpg'
download_panstarrs(RA, DEC, radius=None, output_file=None)
'''
#download_sdss(RA, DEC, radius=4., output_file=output_file)
