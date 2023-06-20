import os
import sys

from pylab import *
import astropy.io.fits as pyfits

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import glob
from matplotlib.pyplot import cm
import re
#from scipy.odr.odrpack import *
import argparse
from scipy.ndimage import rotate
from astropy.stats import sigma_clipped_stats
#matplotlib.use('TkAgg')
import warnings
import itertools as it
import shutil
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
from joblib import Parallel, delayed
from astropy.stats import sigma_clip
from scipy import special
import scipy.ndimage as ndi
import pickle
import collections
import subprocess
import tarfile
from astropy.convolution import Gaussian2DKernel
from scipy import interpolate


LOCAL_DIR = "/imred"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'iraf_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'decomposition/simple_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))


#sys.path.append('/home/byu.local/mosav/MyGit/IMAN/iraf_fitting')
#sys.path.append('/home/byu.local/mosav/MyGit/IMAN/decomposition/simple_fitting')
#sys.path.append('/home/byu.local/mosav/MyGit/IMAN/imp/cropping')

#import iraf_ellipse
#import sersic_fitting
#import get_galaxy_polygon
#import add_keyw_galfit

from astropy.convolution import convolve
from astropy.modeling import models, fitting
import time
from scipy import stats
# Try: https://photometrypipeline.readthedocs.io/en/latest/index.html

import statistics
import warnings
warnings.filterwarnings("ignore")


def get_master_bias(frames, output_image='superbias.fits'):
    print('\n\nBias frames reduced:')
    for frame in frames:
        print('\t\t%s' % (frame))
            
    imgs = []
    for frame in frames:
        hdulist = pyfits.open(frame)
        img = hdulist[0].data
        header = hdulist[0].header
        
        imgs.append(img)
    imgs = np.stack(imgs)

    median_img = np.median(imgs, axis=0)

    hdu_out = pyfits.PrimaryHDU(data=median_img, header=header)
    hdu_out.writeto(output_image, overwrite=True, output_verify='ignore')
    return output_image


def get_master_dark(frames, superbias = None, output_image='superdark_onesec.fits', exptime=None):
    print('\n\nDark frames reduced:')
    for frame in frames:
        print('\t\t%s' % (frame))
        
    imgs = []
    for frame in frames:
        hdulist = pyfits.open(frame)
        img = hdulist[0].data
        header = hdulist[0].header
        
        imgs.append(img)
    imgs = np.stack(imgs)

    median_img = np.median(imgs, axis=0)
    
    if exptime is None:
        exptime = float(header['EXPTIME'])
    
    if superbias is None:
        hdu_out = pyfits.PrimaryHDU(data=median_img/exptime, header=header)
    else:
        hdulist1 = pyfits.open(superbias)
        superbias_img = hdulist1[0].data
        header['EXPTIME'] = 1.0
        hdu_out = pyfits.PrimaryHDU(data=(median_img- superbias_img)/exptime, header=header)
    
    
    hdu_out.writeto(output_image, overwrite=True, output_verify='ignore')
    return output_image

def find_mode(array):
    vals,counts = np.unique(array, return_counts=True)
    index = np.argmax(counts)
    return(vals[index])

def get_master_flatfield(flatfield_frames, master_bias_frame, master_dark_frame, exptimes=None, output_image='superflatfield.fits'):
    print('\n\nFlat field frames reduced:')
    for flatfield_frame in flatfield_frames:
        print('\t\t%s' % (flatfield_frame))

    hdulist = pyfits.open(master_bias_frame)
    superbias_img = hdulist[0].data    
    
    hdulist1 = pyfits.open(master_dark_frame)
    superdark_onesec_img = hdulist1[0].data        
    
    flatfield_imgs = []
    for k in range(len(flatfield_frames)):
        hdulist2 = pyfits.open(flatfield_frames[k])
        flatfield_img = hdulist2[0].data
        header = hdulist2[0].header    

        if exptimes is None:
            exptime = float(header['EXPTIME'])
        else:
            exptime = exptimes[k]
        
        flatfield_imgs.append(flatfield_img - superbias_img - superdark_onesec_img*exptime)
        #ny,nx = np.shape(flatfield_img)
        #print(ny,nx)
    
    combined_img = np.array(flatfield_imgs).ravel()
    
    #mode_img = find_mode(combined_img)
    #print(mode_img)
    #exit() 
    
    #plt.hist(combined_img, bins=200)
    #plt.show() 
    
    
    #exit()
    
    
    #print(stats.mode(combined_img))
    #exit()
    

    mode_all = np.median(combined_img)   # Must be mode!!!

    

    for k in range(len(flatfield_imgs)):
        mode_img = np.median(flatfield_imgs[k].ravel())  # Must be mode!!!
        flatfield_imgs[k] = flatfield_imgs[k] * mode_all/mode_img

    imgs = np.stack(flatfield_imgs)

    img_stacked = np.median(imgs, axis=0) 
    
    median_img_stacked = np.median(img_stacked.ravel())
    
    median_img_norm = img_stacked / median_img_stacked
    
    hdu_out = pyfits.PrimaryHDU(data=median_img_norm)
    hdu_out.writeto(output_image, overwrite=True)    
    
    return output_image
    

def get_reduced_science_frames(science_frames, masterbias, masterdark, masterflatfield, exptime=None):
    print('\n\nScience frames reduced:')
    for science_frame in science_frames:
        print('\t\t%s' % (science_frame))
    
    
    hdulist1 = pyfits.open(masterbias)
    img_bias = hdulist1[0].data

    hdulist2 = pyfits.open(masterdark)
    img_dark_onesec = hdulist2[0].data

    hdulist3 = pyfits.open(masterflatfield)
    img_flat = hdulist3[0].data
    
    output_images = []
    for science_frame in science_frames:
        hdulist = pyfits.open(science_frame)
        img = hdulist[0].data
        header = hdulist[0].header
        if exptime is None:
            exptime = float(header['EXPTIME'])
      
        output_image = science_frame.split('/')[-1].split('.fits')[0] + '_reduced.fits'
        img = (img - img_bias - img_dark_onesec*exptime)/img_flat
        hdu_out = pyfits.PrimaryHDU(data=img, header=header)
        hdu_out.writeto(output_image, overwrite=True)
        output_images.append(output_image)
    
    return output_images
            
        


def main(cals_path, science_path, science_prefix, bias_prefix='*bias*.fits', dark_prefix='*dark*.fits', flat_prefix='*flat.g*'):
    print('***DATA REDUCTION***')
    
    # Bias
    bias_frames = glob.glob(cals_path + '/%s' % (bias_prefix))
    masterbias = get_master_bias(bias_frames)
    
    # Dark
    dark_frames = glob.glob(cals_path + '/%s' % (dark_prefix))
    masterdark = get_master_dark(dark_frames, superbias = masterbias, output_image='superdark_onesec.fits', exptime=None)
    
    # Flat field
    flatfield_frames = glob.glob(cals_path + '/%s' % (flat_prefix))

    masterflatfield = get_master_flatfield(flatfield_frames, masterbias, masterdark, exptimes=None, output_image='superflatfield.fits')
    
    # Final reduction
    science_frames = glob.glob(science_path + '/%s' % (science_prefix))
    output_images = get_reduced_science_frames(science_frames, masterbias, masterdark, masterflatfield, exptime=None)
    
    print('Done!')    
    return output_images
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Initial data reduction")
    parser.add_argument("cals_path", help="Path to calibrations")
    parser.add_argument("science_path", help="Path to science images")
    parser.add_argument("science_prefix", help="science prefix") 
    parser.add_argument("--bias_prefix", nargs='?', const=1, help="Optional: default *bias*.fits", type=str, default='*bias*.fits')     
    parser.add_argument("--dark_prefix", nargs='?', const=1, help="Optional: default *dark*.fits", type=str, default='*dark*.fits')      
    parser.add_argument("--flat_prefix", nargs='?', const=1, help="Optional: default *flat.g*", type=str, default='*flat.g*')        
    

    
    args = parser.parse_args()

    cals_path = args.cals_path
    science_path = args.science_path
    science_prefix = args.science_prefix
    bias_prefix = args.bias_prefix
    dark_prefix = args.dark_prefix
    flat_prefix = args.flat_prefix

    main(cals_path, science_path, science_prefix, bias_prefix='*bias*.fits', dark_prefix='*dark*.fits', flat_prefix=flat_prefix)

#main('/media/mosav/MY_DATA_DRIVE/APO_observations_red/UT211009/UT211009/cals/','/media/mosav/MY_DATA_DRIVE/APO_observations_red/UT211009/UT211009','UGC10043.g.*.fits')
