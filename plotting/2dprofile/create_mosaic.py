#! /usr/bin/env python
#!/usr/bin/python
# -*- coding:  cp1251 -*-

# EXAMPLE: python ~/CurrentWork/ImaPrep/IMAN/SKIRT/plot_fitskirt_images.py imf_output.fski /
# /home/amosenko/CurrentWork/HEROES/models/IC2531/reference /home/amosenko/CurrentWork/HEROES/models/IC2531/RESULTS/fit3

import sys
import math
import numpy as np
from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import NullFormatter
from mpl_toolkits.axes_grid1 import AxesGrid
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import random
import re
import glob
from matplotlib.colors import LogNorm
import astropy.io.fits as pyfits
import argparse
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.image as mpimg
import argparse

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True


#warnings.filterwarnings("ignore")

fsize = 12

def plot_indiv_frame(grid, number, image_file, borders='0,0,0,0', min_int=0., max_int=0., name = ''):
        img = mpimg.imread(image_file)
        h, w, c = img.shape
        '''
        xaxis = data.shape[1]
        yaxis = data.shape[0]
        if borders!='0,0,0,0':    
            borders = borders.split(',')
            x0 = int(borders[0])
            y0 = int(borders[1])
            x1 = int(borders[2])
            y1 = int(borders[3])
        else:       
            x0 = 0
            y0 = 0
            x1 = xaxis
            y1 = yaxis
        '''
        #vmax = np.max(data)
        #vmin = min_int#np.min(data)
        #if name=='Tilted envelope': 
        #    vmax = 50.

        im = grid[number].imshow(img)

        #grid[number].text(0.03, 0.95, name, color='black', transform=grid[number].transAxes, fontsize=fsize, fontweight='bold', va='top')
        #return xaxis,yaxis


def create_black_image(dimx, dimy, extension="jpg", color=0):
    from PIL import Image

    # create a black image
    img = np.zeros((dimy, dimx, 3), dtype = np.uint8)
    img.fill(color) # numpy array!
    im = Image.fromarray(img) #convert numpy array to image
    im.save('whh.%s' % (extension))




def main(folder, ratio="4:3", extension='png', min_int=0., borders='0,0,0,0', output_file='mosaic.png'):

    files = glob.glob('%s/*.%s' % (folder, extension))
    #files = files.sort(key=lambda f: int(''.join(filter(str.isdigit, f))))

    Number_of_files = len(files)
    print('Number of files:', Number_of_files)
    #exit()

    Ratio = float(ratio.split(':')[0]) / float(ratio.split(':')[1])

    dim_y = int(math.ceil(math.sqrt(Number_of_files/Ratio)))
    dim_x = int(math.ceil(Number_of_files/dim_y))

    print('Number of rows:', dim_y)
    print('Number of columns:', dim_x)
    #exit()
    inch_per_image = 3

    dim_y_inch = dim_y * inch_per_image
    dim_x_inch = dim_x * inch_per_image

    fig = plt.figure(0, (dim_x_inch, dim_y_inch))
    fig.subplots_adjust(left=0.05, right=0.95)

    grid = AxesGrid(fig, 111,
                        nrows_ncols=(dim_y, dim_x),
                        axes_pad=0.02,
                        label_mode="L",
                        share_all=True
                        ) # cbar_mode="single"

    for cax in grid.cbar_axes:
        cax.toggle_label(False)

    img = mpimg.imread(files[0])
    h, w, c = img.shape

    #print(w,h)
    #exit()

    number = 0
    for i in range(len(files)):
        plot_indiv_frame(grid, number, files[i], borders='0,0,0,0', min_int=0., max_int=0., name = '')
        number +=1

    if number<dim_x*dim_y:
        # Fill in the rest with blank images:
        create_black_image(w, h, extension=extension, color=0) # Color can be changed: 0 is black, 255 is white
        for i in range(dim_x*dim_y-number):
            plot_indiv_frame(grid, number, 'whh.%s' % (extension), borders='0,0,0,0', min_int=0., max_int=0., name = '')
            number +=1
        os.remove('whh.%s' % (extension))



    
    #grid.axes_llc.set_xlim(0,xmax)
    #grid.axes_llc.set_ylim(0,ymax)
    grid.axes_llc.set_xticklabels([])
    grid.axes_llc.set_yticklabels([])
    grid.axes_llc.get_xaxis().set_ticks([])	# To remove ticks
    grid.axes_llc.get_yaxis().set_ticks([])	# To remove ticks


    plt.draw() 
    #plt.show()
    plt.savefig(output_file, bbox_inches='tight', pad_inches=0.01, dpi = 300)
    plt.clf()
    plt.close()   
    #plt.show()
    return output_file


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A script to create a mosaic of images. All input images must have the same dimensions")
    parser.add_argument("folder", help="Directory with all images to be placed on the mosaic.")
    parser.add_argument("--e", help="Extension for the input images: jpg", type=str, default="jpg")
    parser.add_argument("--o", help="Output image: mosaic.jpg", type=str, default="mosaic.jpg")
    parser.add_argument("--r", help="Image ratio for the output image: 4:3", type=str, default="4:3")

    args = parser.parse_args()

    folder = args.folder
    extension = args.e
    output_file = args.o
    ratio = args.r

    main(folder, ratio=ratio, extension=extension, output_file=output_file)
