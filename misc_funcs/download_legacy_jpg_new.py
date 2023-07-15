#! /usr/bin/env python
''' 
Example of input file:

UGC10043 237.171708 21.869389 3 0.148
SDSSJ151223.37+013823.9 228.0973 1.6401 1 0.572
PRCA-1 24.730255 -7.765600 1.5 0.373
PRCA-6 228.984784 43.166636 2 0.365
PRCB-17 222.810048 35.542269 1.5 0.088
SPRC-7 118.143038 29.347166 0.8 1.160
SPRC-69 312.023607 0.068844 1 0.498
'''

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
from PIL import Image, ImageEnhance, ImageFilter, ImageOps # pip3 install pillow
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

LOCAL_DIR = "/misc_funcs"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'plotting/2dprofile'))
import plot_smoothed_image

def create_composite(Img, img):
    thresh = 150
    fn = lambda x : 255 if x > thresh else 0
    r = img.convert('L').point(fn, mode='1')


    width, height = img.size 
    
    for y in range(height):
        for x in range(width):
            coordinate = x, y
            rgb = Img.getpixel(coordinate)

            if r.getpixel(coordinate)== 0:
                img.putpixel( (x, y), (rgb[0], rgb[1], rgb[2], 255) ) 


    return img

def resize_deep_images(deep_image, output_image, R=0., PA=0., resolution=600, brightness_factor=1., contrast_factor=1., sharpness_factor=1.,sigma_smooth=5., invert=True, composite=True, text=None, pix2sec=1., kpc_per_arc=1., hor_pos=0.03, vert_pos=0.8, L_bar=60., plot_scale_bar=True):
    basewidth = resolution
    img = Image.open(deep_image)
    NX = img.size[0]
    NY = img.size[1]

    if PA!=0.:
        img = img.rotate(-PA, resample=Image.CUBIC)

        width, height = img.size
        new_r = R/math.sqrt(2.)

        xc = width/2.
        yc = height/2.
        border = (int(xc-new_r),int(yc-new_r),int(xc+new_r),int(yc+new_r))
        img = img.crop(border)
        Img = img
    else:
        Img = img

    if brightness_factor!=1.:
        enhancer = ImageEnhance.Brightness(img)
        img = enhancer.enhance(brightness_factor)

    if contrast_factor!=1.:
        enhancer = ImageEnhance.Contrast(img)
        img = enhancer.enhance(contrast_factor)


    
    if sigma_smooth is not None:
        img = img.filter(ImageFilter.GaussianBlur(radius = sigma_smooth)) 

    if sharpness_factor!=1.:
        enhancer = ImageEnhance.Sharpness(img)
        img = enhancer.enhance(sharpness_factor)
    
    if invert:
        img = ImageOps.invert(img)
    
    if composite:
        img = create_composite(Img, img)

    wpercent = (basewidth/float(img.size[0]))
    hsize = int((float(img.size[1])*float(wpercent)))
    img = img.resize((basewidth,hsize), Image.ANTIALIAS)
    img = img.convert("RGB") # Remove this if you want png!!!!
    #img.save(output_image) 
    

    # Plot final image
    my_dpi = 100
    fig =  plt.figure(0, figsize=(6,6))
    #fig.set_size_inches(6, 6)
    ax = plt.subplot()
    fsize = 13
    nx = img.size[0]
    ny = img.size[1]

    
    ax.imshow(img, origin="upper", aspect='equal')
    ax.set_frame_on(False)
    ax.set_xticks([]); ax.set_yticks([])
    plt.axis('off')    
    ax.text(hor_pos, vert_pos, text, fontsize=fsize, color='black',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline',backgroundcolor='whitesmoke')
    
    if plot_scale_bar:
        plot_smoothed_image.add_scale_bar(ax, nx, ny, (float(NX)/float(nx))*pix2sec, L_bar=L_bar, kpc_per_arc=kpc_per_arc)
        plt.draw() 
    plt.savefig(output_image, bbox_inches='tight', pad_inches=0.01, format='png', dpi=my_dpi)
    plt.clf()
    plt.close()

    
    


def main(k, name, RA, DEC, R, PA, kpc_per_arc, pixscale=0.262, resolution=600, brightness_factor=4.0, contrast_factor=15.,sharpness_factor=0.01, invert=True, composite=True, output_dir='./legacy',L_bar=30., output_file=None):
    # R in arcmin
    
    if output_file is None:
        output_file = '%s.png' % (name)

    try:
           RR = int(math.ceil(R*60./pixscale))
           url="http://legacysurvey.org/viewer/jpeg-cutout?ra=%f&dec=%f&width=%i&height=%i&layer=dr8&pixscale=%.3f&bands=grz" % (RA, DEC, 2*RR, 2*RR, pixscale) 
           urllib.request.urlretrieve(url, 'tmp_%s.jpg' % (name))

           resize_deep_images('tmp_%s.jpg' % (name), output_file,RR, PA, resolution=resolution, brightness_factor=brightness_factor, contrast_factor=contrast_factor, sharpness_factor=sharpness_factor,sigma_smooth=5., invert=invert, composite=composite, text=name, pix2sec=pixscale, kpc_per_arc=kpc_per_arc, hor_pos=0.03, vert_pos=0.93, L_bar=L_bar) # TODO: brightness_factor, contrast_factor, sharpness_factor can be tuned!
           # brightness_factor=3, contrast_factor=5., sharpness_factor=0.01
           os.remove('tmp_%s.jpg' % (name))
           shutil.move(output_file, '%s/%s' % (output_dir,output_file))
           print('%i %s Done!' % (k+1, name))
    except:
           print('%i %s Failed' % (k+1, name))

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download jpeg images of objects with (RA,DEC)")
    parser.add_argument("--i", help="Input file in the following format: Name RA(deg) DEC(deg) Radius(arcmin)  Scale(kpc/arcsec)", type=str)
    parser.add_argument("--s", nargs='?', const=1, help="Optional: pixscale (in arcsec)", type=float, default=0.262)

    parser.add_argument("--r", nargs='?', const=1, help="Optional: image resolution (in pixels)", type=int, default=600)
    parser.add_argument("--b", nargs='?', const=1, help="Optional: brightness factor", type=float, default=4.0)
    parser.add_argument("--c", nargs='?', const=1, help="Optional: contrast factor", type=float, default=15.0)
    parser.add_argument("--h", nargs='?', const=1, help="Optional: sharpness factor", type=float, default=0.01)
    parser.add_argument("--o", nargs='?', const=1, help="Optional: output directory", type=str, default='./legacy')
    parser.add_argument("--l", nargs='?', const=1, help="Optional: length of the scale bar in arcsec", type=float, default=30.)    
    parser.add_argument("--proc", nargs='?', const=1, help="Optional: number of processes", type=int, default=1)

    
    args = parser.parse_args()
    
    input_file = args.i
    pixscale = args.s
    resolution = args.r
    brightness_factor = args.b
    contrast_factor = args.c
    sharpness_factor = args.h
    output_dir = args.o
    L_bar = args.l
    n_jobs = args.proc
    
    RA,DEC,R,Scale = np.loadtxt(input_file, usecols=[1,2,3,4], dtype=float, unpack=True)
    name = np.loadtxt(input_file, usecols=[0], dtype=str, unpack=True)

    # R in arcmin
    print('Downloading...')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    kpc_per_arc = [float('nan')] * len(RA)
    PA = [0.01] * len(RA)
    
    Parallel(n_jobs=n_jobs)(delayed(main)(k, name[k], RA[k], DEC[k], R[k], PA[k], kpc_per_arc[k], pixscale=0.262, resolution=600, brightness_factor=brightness_factor, contrast_factor=contrast_factor,sharpness_factor=sharpness_factor, invert=True, output_dir=output_dir,L_bar=L_bar) for k in range(len(name)))
    
    #main(k, name, RA, DEC, R, PA, kpc_per_arc, pixscale=0.262, resolution=600, brightness_factor=4.0, contrast_factor=15.,sharpness_factor=0.01, output_dir='./legacy',L_bar=30.)
    

    print('Done!')


