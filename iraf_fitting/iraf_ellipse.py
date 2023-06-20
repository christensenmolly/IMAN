#!/usr/bin/python
# DESCRIPTION:
# STSDAS IRAF/ELLIPSE fitting. IRAF is required. Do not create your own login.cl where you launch this script. 
# MINIMAL USAGE: python  iraf_ellipse.py [input_image]
# EXAMPLE:  python iraf_ellipse.py galaxy_clean_galf.fits --xc 1122 --yc 1021 --maxsma 850 --ZP 29.5962 --pix2sec 0.833 --mask mask.fits

import sys
import math
import numpy as np
#from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import astropy.io.fits as pyfits
from astropy.stats import sigma_clipped_stats
import argparse
import glob
import warnings
warnings.filterwarnings("ignore")

import plot_iraf_results

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

tmp_out = sys.stdout


fsize=16
FNULL = open(os.devnull, 'w')


def create_reg_with_ellipses(x, y, sma, ell, PA):
    # Function to create region file with IRAF/ELLIPSE ellipses
    f = open('iraf_ellipses.reg', 'w')
    f.write('%s\n' % ('image') )
    for k in range(len(x)):
        f.write('ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # color=red\n' % (x[k], y[k], sma[k], sma[k]*(1.-ell[k]), PA[k]+90.))
    f.close()    


def crea_ell(input_image, xc, yc, ellip=0.2, pa=20., sma=10., step=0.03, minsma=0.0, maxsma=-1., ell_file='ellipse.txt', fits_mask="", fflag=1.0, olthresh=0.0, linear='no', hcenter='no', model_file=None, hellip='no', hpa='no'):
        if os.path.exists('gal.tab'):
            os.remove('gal.tab')

        
        
        # Function to create cl script to run IRAF/ELLIPSE
        f = open("ell.cl", "w") 
        sys.stdout = f

        print("# Script: ellipse")
        print("!rm -f %s" % (ell_file))
        print("stsdas")
        print("analysis")
        print("isophote")
        
        print("geompar.linear=%s" % (linear))
        print("geompar.step=%.2f" % (step))
        print("geompar.minsma=%.2f" % (minsma))
        if maxsma==-1.:
            print("geompar.maxsma=INDEF")
        else:
            print("geompar.maxsma=%.2f" % (maxsma))
        
        print("geompar.x0=%.3f" % (xc))
        print("geompar.y0=%.3f" % (yc))
        print("geompar.ellip0=%.3f" % (ellip))
        print("geompar.pa0=%.3f" % (pa))
        print("geompar.sma0=%.3f" % (sma))
        
        print("samplepar.fflag=%.1f" % (fflag)) # Acceptable fraction of flagged data points in intensity sample.
        
        print("controlpar.hcenter=%s" % (hcenter))
        print("controlpar.hellip=%s" % (hellip)) 
        print("controlpar.hpa=%s" % (hpa)) 
        print("controlpar.olthresh=%.1f" % (olthresh)) # Threshold for the object locator algorithm. By lowering this value the object locator becomes less strict, in the sense that it will accept lower signal-to-noise data. If set to zero, the x0, y0 values found in the geompar pset are used without questioning.
        #if fits_mask!="no":
        #    print("geompar.dqf=%s" % (fits_mask))
        print("controlpar.minit=10")
        print("controlpar.maxit=100")
        
        print("ellipse %s gal.tab dqf=%s" % (input_image, fits_mask))
        print("tprint gal.tab pwidth=600 plength=5000 > %s" % (ell_file))
        if model_file is not None:
            if os.path.exists(model_file):
                os.remove(model_file)
            print("bmodel table=gal.tab output=%s parent=%s fulltable=yes" % (model_file, input_image))

        #print("!rm -f gal.tab")
        print("logout")

        sys.stdout = tmp_out
        f.close()






def main_ell(input_image, xc, yc, ellip=0.2, pa=20., sma0=10., m0=28., pix2sec=1.0, step=0.03, minsma=0., maxsma=None, outp_format='png', ell_file='ellipse.txt', fits_mask="", fflag=1.0, olthresh=0.0, linear='no', hcenter='no', layers=['all'], verbosity=True, model_file=None, hellip='no', hpa='no'):
        
        # IRAF/ELLIPSE fitting
        if verbosity:
            print('IRAF/ELLIPSE fitting:')

        #if pa>=180.:
        #    pa = pa - 180.
        
        if pa>=0.:
            pa = pa % 360
        else:
            pa = pa % -360
        
        
        if pa>=90.:
            pa = pa-180.
        elif pa<=-90.:
            pa = pa+180.
        
        
        if round(pa)==90.:
            pa = 89.9
                      
        if round(pa)==-90.:
            pa = -89.9


        if np.isnan(pa) or np.isinf(pa):
            pa = 0.
            

        
        if ellip<=0.05:
            ellip=0.06
        elif ellip>=1.:
            ellip=0.99

        sma0 = math.ceil(sma0)
        
        if fits_mask=="":
            fits_mask = "no"

        hdulist = pyfits.open(input_image)
        
        if layers == ['all']:
            layers = range(len(hdulist))
        
        layer_names = []
        ellipse_files = []
        for layer in layers:
            layer = int(layer)
            if verbosity:
                print('\t Layer %i ...' % (layer))
            
            data = hdulist[layer].data
            header = hdulist[layer].header  #### WARNING!
            
            layer_name = None
            if 'NAME_OF_LAYER' in header:
                layer_name = header['NAME_OF_LAYER']
                
            tmp_image = input_image.split('.fits')[0] + '_%i.fits' % (layer)
            
            if layer!=0:
                if layer_name is not None:
                    layer_names.append(layer_name)
                    if layer_names.count(layer_name)>1:
                        layer_name = '%s%i' % (layer_name, layer_names.count(layer_name))
                    ellipse_file = ell_file.split('.txt')[0]+'_%s.txt' % (layer_name)
                else:
                    ellipse_file = ell_file.split('.txt')[0]+'_%i.txt' % (layer)
            else:
                ellipse_file = ell_file
            ellipse_files.append(ellipse_file)
            ny,nx = np.shape(data)
            
            if xc==None and yc==None:
                xc = nx/2.
                yc = ny/2.
            else:
                xc = float(xc)
                yc = float(yc)
                
            if maxsma==None:
                maxsma = math.sqrt( (nx/2.)**2 + (ny/2.)**2 )
            else:
                maxsma = float(maxsma)            
        
            outHDU = pyfits.PrimaryHDU(data)
            outHDU.writeto(tmp_image, overwrite=True)          
            
            if layer!=0:
                mask_image = "no"
            else:
                mask_image = fits_mask

            crea_ell(tmp_image, xc, yc, ellip=ellip, pa=pa, sma=sma0, step=step, minsma=minsma, maxsma=float(maxsma), ell_file=ellipse_file, fits_mask=mask_image, fflag=fflag, olthresh=olthresh, linear=linear, hcenter=hcenter, model_file=model_file, hellip=hellip, hpa=hpa)
            os.chmod(r"ell.cl",0o777)
            if not verbosity:
                subprocess.call("cl < ell.cl -o", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
            else:
                subprocess.call("cl < ell.cl -o", shell=True)
            
            os.remove(tmp_image)
            
            if layer==0:
                sma, inten, inten_err, ell, errell, PA, errPA, x0, y0, B4, errB4 = plot_iraf_results.read_ell(ellipse_file)
                create_reg_with_ellipses(x0, y0, sma, ell, PA)
        
            files = glob.glob('upar*.par')
            for file in files:
                os.remove(file)
            os.remove('ell.cl')

        
        
        plot_iraf_results.main(ellipse_files, m0, pix2sec, show_error_bars=True)       
                
        
        print('\nDone!')        
        
        
        
        
        


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="IRAF/ELLIPSE fitting")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("--outp_format", nargs='?', const=1, help="Optional: Input the name of the format for the output picture ", type=str, default='eps') 
    parser.add_argument("--mask", nargs='?', const=1, help="Optional: Input the name of the mask file",type=str,default='')

    parser.add_argument("--xc", nargs='?', const=1, help="Optional: Input the x-coordinate of the centre",type=str,default=None) 
    parser.add_argument("--yc", nargs='?', const=1, help="Optional: Input the y-coordinate of the centre",type=str,default=None) 

    parser.add_argument("--ellip", nargs='?', const=1, help="Optional: Input initial ellipticity",type=float, default=0.2)
    parser.add_argument("--pa", nargs='?', const=1, help="Optional: Input initial PA",type=float, default=20.)
    parser.add_argument("--sma", nargs='?', const=1, help="Optional: Input initial sma",type=float, default=10.)
    
    parser.add_argument("--maxsma", nargs='?', const=1, help="Optional: Input the maximum radius",type=float, default=None)
    parser.add_argument("--minsma", nargs='?', const=1, help="Optional: Input the minimum radius",type=float, default=1.)
    parser.add_argument("--step", nargs='?', const=1, help="Optional: Input the step",type=float, default=0.03)    
    parser.add_argument("--ZP", nargs='?', const=1, help="Optional: Input the y-coordinate of the centre",type=float,default=28.) 
    parser.add_argument("--pix2sec", nargs='?', const=1, help="Optional: Input the y-coordinate of the centre",type=float,default=1.) 
    parser.add_argument("--fflag", nargs='?', const=1, help="Optional: Acceptable fraction of flagged data points in intensity sample.", type=float, default=1.0) 
    parser.add_argument("--olthresh", nargs='?', const=1, help="Optional: Threshold for the object locator algorithm. See http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?controlpar.hlp", type=float, default=0.0)
    parser.add_argument("--layers", nargs='?', const=1, help="Optional: Layers of the fits file to be done (ALL)", type=str, default='all')     
    
  
    parser.add_argument("--linear", action="store_true", default=False,
                        help="Linear geometric step? Default False.")  

    parser.add_argument("--fix_cen", action="store_true", default=False,
                        help="Fix the centre? Default False.")  

    parser.add_argument("--fix_ell", action="store_true", default=False,
                        help="Fix the ellipticity? Default False.")  

    parser.add_argument("--fix_pa", action="store_true", default=False,
                        help="Fix the PA? Default False.")

    parser.add_argument("--model_file", nargs='?', const=1, help="Optional: Model file",type=str,default=None)    
    
    args = parser.parse_args()

    input_image = args.inputImage
    outp_format = args.outp_format
    xc = args.xc
    yc = args.yc
    ellip = args.ellip
    pa = args.pa
    sma = args.sma
    maxsma = args.maxsma
    minsma = args.minsma
    ZP = args.ZP
    pix2sec = args.pix2sec
    mask = args.mask
    step = args.step
    fflag = args.fflag
    olthresh = args.olthresh
    linear = args.linear
    if linear:
        linear='yes'
    else:
        linear='no'
    layers = args.layers.split(',')
    
    fix_cen = args.fix_cen
    fix_ell = args.fix_ell
    fix_pa = args.fix_pa
    if fix_cen:
        hcenter = 'yes'
    else:
        hcenter = 'no'

    if fix_ell:
        hellip = 'yes'
    else:
        hellip = 'no'    

    if fix_pa:
        hpa = 'yes'
    else:
        hpa = 'no'
    
    model_file = args.model_file
    
    main_ell(input_image, xc, yc, ellip=ellip, pa=pa, sma0=sma, m0=ZP, pix2sec=pix2sec, step=step, minsma=minsma, maxsma=maxsma, outp_format=outp_format, ell_file='ellipse.txt', fits_mask=mask, fflag=fflag, olthresh=olthresh, linear=linear, layers=layers, hcenter=hcenter, hellip=hellip, hpa=hpa, model_file=model_file)


