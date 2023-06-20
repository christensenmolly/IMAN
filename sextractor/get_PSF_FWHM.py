#!/usr/bin/python
# DESCRIPTION:
# Script to run psfex.
# EXAMPLE: python3 ~/MyGit/IMAN/sextractor/get_PSF_FWHM.py sky_subtr_someskysub.fits --pix2sec 0.24 --xc 1590 --yc 1590

import sys
import subprocess
from pylab import *
import itertools
import os
from os.path import exists
from os import remove
from scipy.spatial import cKDTree
from scipy.optimize import fmin_tnc, fmin
#import arithm_operations
import argparse
from astropy.modeling import models, fitting
from datetime import datetime
import shutil
# Disable astropy logging except for warnings and errors
from astropy import log
from astropy.io import fits as fits
log.setLevel("WARNING")
import tempfile
import warnings
warnings.filterwarnings("ignore")

LOCAL_DIR = "/sextractor"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]


sys.path.append(os.path.join(IMAN_DIR, 'imp/misc'))

import run_SExtractor
import read_SEcat
import norm_image



def main(input_image, pix2sec, xc, yc, size=25, output_psf_file='psf.fits',output_moffatpsf_file='moffat_psf.fits'):
    tmp_dir = tempfile.mkdtemp()
    current_dir = os.path.abspath(os.getcwd())
    shutil.copy(os.path.abspath(input_image), tmp_dir + '/' + os.path.abspath(input_image).split('/')[-1])
    
    os.chdir(tmp_dir)
    
    # Run
    run_SExtractor.call_SE(input_image, snr=None, min_pix=None, sextr_dir=None, sextr_setup='cold.sex', sextr_param='default.param', output_cat='field.cat', checkimage_type='SEGMENTATION',checkimage_name='segm.fits', sextr_add_string=None,determine_sky=False, m0=28.0, GAIN=10., pix2sec=1., fwhm=4., verbosity=True)

    columns = read_SEcat.find_sex_column('field.cat', ['CLASS_STAR','FWHM_IMAGE'], dtypes=str)
    CLASS_STAR = np.array(columns[0], dtype=float)
    FWHM = np.array(columns[1], dtype=float)
    FWHM_STARS = FWHM[CLASS_STAR>0.9] * pix2sec # in arcsec!
    
    print('FWHM of stars (in arcsec) = %.2f +/- %.2f' % (np.mean(FWHM_STARS),np.std(FWHM_STARS)))    
    
    if output_psf_file is not None:
        # Run Sextractor
        run_SExtractor.call_SE(input_image, snr=None, min_pix=None, sextr_dir=None, sextr_setup='default_psfex.sex', sextr_param='prepsfex.param', output_cat='field.cat', checkimage_type='SEGMENTATION',checkimage_name='segm.fits', sextr_add_string=None,determine_sky=False, m0=28.0, GAIN=10., pix2sec=pix2sec, fwhm=np.mean(FWHM_STARS), verbosity=True)
    
        # Run psfEX
        subprocess.call('psfex field.cat -CHECKPLOT_DEV NULL -PSF_SIZE %i,%i -VERBOSE_TYPE QUIET' % (size,size), shell=True)
    
        # Create PSF for the specific position
        subprocess.call('psfex-rec field.psf psf_tmp.fits %i %i' % (xc,yc), shell=True)
    
        # Norm PSF:
        norm_image.main('psf_tmp.fits', output_psf_file)
        os.remove('psf_tmp.fits')
        shutil.copy(output_psf_file, current_dir + '/' + output_psf_file)
        
        
        # Fit with a Moffat function:
        data = fits.getdata(output_psf_file, ext=0)
        fit_w = fitting.LevMarLSQFitter()
        y0, x0 = np.unravel_index(np.argmax(data), data.shape)
        sigma = np.mean(FWHM_STARS)/2.35
        bbeta = 2.5
        ggama = np.mean(FWHM_STARS)/(2.*math.sqrt(2**(1./bbeta)-1.))
        
        
        amp = np.max(data) # Maximal intensity
        w = models.Moffat2D(amp, x0, y0, ggama, bbeta)

        yi, xi = np.indices(data.shape)
        g = fit_w(w, xi, yi, data)
        print('FWHM=%.2f pix' % (g.fwhm))
        print('FWHM=%.2f arcsec' % (float(g.fwhm)*pix2sec))

        print('gamma=%.2f' % (g.gamma[0]))
        print('alpha=%.2f' % (g.alpha[0]))
        
        
        g.x_0 = g.x_0 + (51-27)/2
        g.y_0 = g.y_0 + (51-27)/2
        yi, xi = np.indices((51,51))
        model_data = g(xi, yi)
        hdu = fits.PrimaryHDU(model_data)
        hdu.writeto('model.fits', overwrite=True)
        shutil.copy('model.fits', current_dir + '/' + output_moffatpsf_file)
        
        
    shutil.rmtree(tmp_dir)
    return np.mean(FWHM_STARS),np.std(FWHM_STARS)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PSF estimation")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("--pix2sec", nargs='?', const=1, help="Optional: Pixel scale",type=float, default=1.)
    parser.add_argument("--xc", nargs='?', const=1, help="Optional: x-Center to create PSF",type=int, default=1)
    parser.add_argument("--yc", nargs='?', const=1, help="Optional: y-Center to create PSF",type=int, default=1)
    parser.add_argument("--size", nargs='?', const=1, help="Optional: Size of the output psf file",type=int, default=25)
    parser.add_argument("--out_psf", nargs='?', const=1, help="Optional: Output psf file",type=str, default='psf.fits')
    parser.add_argument("--out_moffatpsf", nargs='?', const=1, help="Optional: Output moffat psf file",type=str, default='moffat_psf.fits')
    
    #parser.add_argument("--snr", nargs='?', const=1, help="Optional: Signal-to-noise ratio of detected pixels",type=float,default=2.) 
    #parser.add_argument("--min_pix", nargs='?', const=1, help="Optional: Number of joint pixels",type=int,default=5)
    #parser.add_argument("--sextr_dir", nargs='?', const=1, help="Optional:Sextractor directory",type=str,default=None)
    #parser.add_argument("--sextr_setup", nargs='?', const=1, help="Optional: Sextractor setup file",type=str,default='cold.sex')
    #parser.add_argument("--sextr_param", nargs='?', const=1, help="Optional: Sextractor parameter file",type=str,default='default.param')
    #parser.add_argument("--sky", help="Determine sky level", action="store_true", default=False)
    args = parser.parse_args()

    input_image = args.inputImage
    pix2sec = args.pix2sec
    xc = args.xc
    yc = args.yc
    size = args.size
    output_psf_file = args.out_psf
    output_moffatpsf_file = args.out_moffatpsf
    #snr = args.snr
    #min_pix = args.min_pix
    #sextr_dir = args.sextr_dir
    #sextr_setup = args.sextr_setup
    #sextr_param = args.sextr_param
    #determine_sky = args.sky
    
    main(input_image, pix2sec, xc, yc, size, output_psf_file, output_moffatpsf_file)
