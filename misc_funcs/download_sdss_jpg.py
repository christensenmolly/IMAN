#!/usr/bin/python
# -*- coding:  cp1251 -*-

import sys
import math
import numpy as np
from numpy import *
from pylab import *
import os
#import urllib
import urllib.request
from time import sleep
#from panstamps.downloader import downloader
#from panstamps.image import image
import subprocess
import shutil
import glob
import argparse
import matplotlib.image as mpimg

FNULL = open(os.devnull, 'w')



def add_text(input_image, text, output_image):
    hor_pos=0.03; vert_pos=0.9
    fig =  plt.figure(0)
    ax = plt.subplot()
    fsize = 13
    
    img = mpimg.imread(input_image)
    mg2 = ax.imshow(img)

    

    ax.text(hor_pos, vert_pos, text, fontsize=fsize, color='white',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline')

    ax.set_frame_on(False)
    ax.set_xticks([]); ax.set_yticks([])
    plt.axis('off')

    plt.savefig(output_image, bbox_inches='tight', pad_inches=0)
    plt.clf()
    plt.close()
    

def download_sdss(RA, DEC, width=None, output_file=None, add_bar=True, text=None):
        # RA,DEC - decimals
        # width - in arcmin
        print('Downloading ra=%.6f, dec=%.6f' % (RA,DEC))
        if output_file is None:
            output_file = 'ra%.6f_dec%.6f.jpg' % (RA,DEC)
            
        if width is None:
            scale=1.0
        else:
            scale=width/10.
            
        if True:
           # scale=1 is for a box of 10'x10'
           if add_bar:
                pic="http://skyservice.pha.jhu.edu/DR14/ImgCutout/getjpeg.aspx?ra=%s&dec=%s&scale=%.1f&width=600&height=600&opt=G" % (str(RA), str(DEC), scale) # Was DR9
           else:
                pic="http://skyservice.pha.jhu.edu/DR14/ImgCutout/getjpeg.aspx?ra=%s&dec=%s&scale=%.1f&width=600&height=600" % (str(RA), str(DEC), scale) # Was DR9               
           '''
           f = urllib.urlopen(pic)
           image = f.read()
           fout = open(output_file, 'w')
           fout.write(image)
           fout.close()
           '''
           urllib.request.urlretrieve(pic, output_file)
        '''
        except IOError:
           print('Try again!')
        finally:
           if int(os.path.getsize(output_file))>15000:
              print('Done!')
              return 0
           else:
              try:
                 os.remove(output_file)
                 print('Empty! Removed.')
                 return 1
              except:
                 print('Not found')
                 return 1
        '''
        
        if text is not None:
            add_text(output_file, text, 'tmp.jpg')
            os.rename('tmp.jpg', output_file)
        
        
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

    download_sdss(RA, DEC, width=width, output_file=output_file)
