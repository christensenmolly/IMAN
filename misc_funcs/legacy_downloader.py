#! /usr/bin/env python

#
# Author: Sergey Savchenko (savchenko.s.s@gmail.com)
# Example: python ~/programs/Sergey/sdss_downloader-master/sdss_downloader.py gri -i sample.dat -a -s -c --ps -t -f --scatter

from threading import Thread
import subprocess
from math import hypot, log10, sqrt
from time import time, sleep
import urllib.request
import sys
import os
import glob
import shutil
import argparse
import bz2

import numpy as np

import parse_cube_fits



def define_swarpName():
        # Check what name has SWarp package on this system
        rCode = subprocess.call("which swarp >/dev/null", shell=True)
        if rCode == 0:
            swarpName = "swarp"
        else:
            rCode = subprocess.call("which SWarp >/dev/null", shell=True)
            if rCode == 0:
                swarpName = "SWarp"
            else:
                print("\033[31m Error: SWarp was not found on your system.\033[0m")
                print("\033[31m The command has to be either 'swarp' or 'SWarp'\033[0m")
                print("\033[31m Intall SWarp package or try to run this script without -s option.\033[0m")
                exit(1)
        return swarpName

def swarping(name, files, band, swarpName):
    print("Running SWarp for %s band..." % (band))
    callSt = "%s -verbose_type quiet -BACK_TYPE MANUAL " % (swarpName)
    callSt += " ".join(["%s[0]" % (s) for s in files])
    subprocess.call(callSt, shell="True")   
    move("coadd.fits", "%s_%s.fits" % (name, band))
    os.remove("coadd.weight.fits")
    os.remove("swarp.xml")


def main(name, RA, DEC, R):
    print('Downloading...')
    R = R / 60. #  now in deg
    
    raa = np.arange(RA-1.5*R, RA+1.5*R, 450.*0.262/3600.)
    decc = np.arange(DEC-1.5*R, DEC+1.5*R, 450.*0.262/3600.)

    N = 0
    g_files = []
    r_files = []
    z_files = []
    for x in raa:
        for y in decc:
           N = N + 1
           output_file = '%i.fits' % (N)
           url="http://legacysurvey.org/viewer/fits-cutout?ra=%f&dec=%f&size=512&layer=dr8&pixscale=0.262&bands=grz" % (x,y) 
           print('\t'+url)
           '''
           urllib.request.urlretrieve(url, output_file)
           exit()
           parse_cube_fits.main(output_file, '%s_g.fits' % (N), 0, layer_desc='standard')
           parse_cube_fits.main(output_file, '%s_r.fits' % (N), 1, layer_desc='standard')
           parse_cube_fits.main(output_file, '%s_z.fits' % (N), 2, layer_desc='standard')
           
           os.remove(output_file)
           g_files.append('%s_g.fits' % (N))
           r_files.append('%s_r.fits' % (N))
           z_files.append('%s_z.fits' % (N))
           '''
    exit()

    
    print('SWarping...')    
    swarpName = define_swarpName()      
    swarping(name, g_files, 'g', swarpName)
    swarping(name, r_files, 'r', swarpName)
    swarping(name, z_files, 'z', swarpName)
    print('Done!')
    
    
    
main('NGC4452', 187.180455, 11.755032, 3.)    
    
    