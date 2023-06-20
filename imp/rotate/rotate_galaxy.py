#! /usr/bin/env python

import sys
import subprocess
import os
import shutil
import argparse
from os.path import exists
from math import hypot, pi, log10
import math
import numpy  as np
from numpy import ma, mean, std, zeros_like, where, arange, exp, zeros
from numpy import copy as npcopy
from numpy import sum as npsum
from pylab import plot, show, savefig, xlabel, ylabel, vlines, clf
import astropy.io.fits as pyfits
from scipy.odr.odrpack import *
import warnings

import rotate_image
import run_SExtractor

warnings.filterwarnings("ignore")

FNULL = open(os.devnull, 'w')


def gauss(B, x):
    a = B[0]
    x0 = B[1]
    sigma = B[2]
    return a * exp(-(x-x0)**2 / (2*sigma**2))

def check_ds9(output_image, outer_level, inner_level):
    subprocess.call("ds9 %s -scale log -cmap b -regions system image -contour yes -contour scale log -contour limits %.1f %.1f -contour smooth %i -contour nlevels %i -contour" % (output_image,outer_level,inner_level,5,10), shell=True)
    

def fit_by_gauss(xList, yList):
    aGuess = max(yList)
    x0Guess = (xList[0] + xList[-1]) / 2
    sigmaGuess = x0Guess - xList[0]  # FIXME
    fitting = ODR(RealData(xList, yList), 
                  Model(gauss),
                  [aGuess, x0Guess, sigmaGuess])
    fitting.set_job()
    result = fitting.run()
    B1 = (result.beta[0], result.beta[1], result.beta[2])
    # second interation
    dyLust = abs(yList-gauss(B1, xList)) ** 0.5
    fitting = ODR(RealData(xList, yList, sy=dyLust), 
                  Model(gauss),
                  B1)
    fitting.set_job()
    result = fitting.run()
    B2 = (result.beta[0], result.beta[1], result.beta[2])
    return B2

'''
def call_SE(fitsFile):
    # detect the name of SExtractor
    if subprocess.call("which sex >/dev/null", shell=True) == 0:
        callString = "sex %s " % fitsFile
    elif subprocess.call("which sextractor >/dev/null", shell=True) == 0:
        callString = "sextractor %s " % fitsFile
    else:
        print("SExtractor was not found. Exiting...")
        exit(1)
    callString += "-c /home/amosenko/MEGA/MyPrograms/IMAN/SDSS_pipeline/pipelinelib1/default.sex"
    callString += " -VERBOSE_TYPE=QUIET"
    subprocess.call(callString, shell=True)
'''

def get_galaxy_params(fitsFile, xCenField, yCenField):
    """ Function finds object nearest to the center
    of the field."""
    hdu = pyfits.open(fitsFile)
    ySize, xSize = hdu[0].data.shape
    hdu.close()

    minArea = 100  # Minimal area of interested objects [pix^2]
    minCenterDist = 1e10
    for line in open("field.cat"):
        if line.startswith("#"):
            continue
        params = line.split()
        N = int(params[0])
        xCen = float(params[1])
        yCen = float(params[2])
        kron = float(params[8])
        ellA = kron * float(params[4])
        ellB = kron * float(params[5])
        PA = float(params[6])
        ellArea = pi * ellA * ellB
        if ellArea > minArea:
            centerDist = hypot(xCen-xCenField, yCen-yCenField)
            if centerDist < minCenterDist:
                minCenterDist = centerDist
                galN = N
                galXCen = xCen
                galYCen = yCen
                galEllA = ellA
                galEllB = ellB
                galPA = PA
    return galN, galXCen, galYCen, galEllA, galEllB, galPA


def get_backgr_params(fitsFile):
    """Function finds the mean value and the sigma 
    of the backgound. Only pixels outside of all objects
    are being used (based on SE segmentation map)"""
    imageHDU = pyfits.open(fitsFile)
    imageData = imageHDU[0].data
    segmHDU = pyfits.open("segm.fits")
    segmData = segmHDU[0].data
    imageDataMasked = ma.masked_array(imageData, segmData)
    backMean = mean(imageDataMasked)
    backSTD = std(imageDataMasked)
    imageHDU.close()
    segmHDU.close()
    return backMean, backSTD


def get_pa(galaxy_image, xCen, yCen, sePA, ellA, ellB):
    # Finds PA as an angle where biggest number of galaxy pixels are located at
    # a narrow horisontal line
    sliceWidth = min(ellA, ellB) / 2.0
    paList = []
    nPixList = []
    for curPA in arange(sePA-10.0, sePA+10.0, 0.25):
        xCenRot, yCenRot = rotate_image.main(galaxy_image, curPA, xc=xCen, yc=yCen, output_image="rot_tmp.fits", hdu_inp=0, cval=float('nan'))
        rotHDU = pyfits.open("rot_tmp.fits")
        rotData = rotHDU[0].data
        xSizeRot, ySizeRot = rotData.shape
        # Find number of galaxy pixels located inside of horisontal slice of given width
        galPixelsInSlice = 0
        upperEdge = int(yCenRot + sliceWidth)
        lowerEdge = int(yCenRot - sliceWidth)
        for sliceLoc in range(lowerEdge, upperEdge+1):
            row = rotData[sliceLoc]
            galPixelsInSlice += len(row[where(row != 0)])
        paList.append(curPA)
        nPixList.append(galPixelsInSlice)
        os.remove("rot_tmp.fits")
    params = fit_by_gauss(paList, nPixList)
    optPosAng = params[1]
    plot(paList, nPixList)
    plot(paList, gauss(params, paList))
    vlines(optPosAng, min(nPixList), max(nPixList)*1.1, linestyles=":")
    xlabel("posang")
    ylabel("pixels in slice")
    savefig("pa_gauss.png")
    clf()
    return optPosAng



def mask_all_except_galaxy(xc, yc, input_image, segm_image, mask_image, I_DN_min, I_DN_max=None, output_image='galaxy_mask_tmp.fits'):
    """ Function sets values of all background pixels
    and pixels of all objects except of galaxy equal to zero"""
    imageHDU = pyfits.open(input_image)
    imageData = imageHDU[0].data
    
    maskHDU = pyfits.open(mask_image)
    maskData = maskHDU[0].data

    segmHDU = pyfits.open(segm_image)
    segmData = segmHDU[0].data
    
    nObj = segmData[int(yc),int(xc)]
    
    maskedData = npcopy(imageData)
    
    if I_DN_max is None:
        I_DN_max = np.max(imageData)
    
    maskedData[where((maskedData <= I_DN_min) | (maskedData >= I_DN_max))] = 0.0
    
    maskedData[where((maskData !=0) & (segmData != nObj))] = 0.0
    
    maskedHDU = pyfits.PrimaryHDU(data=maskedData)
    maskedHDU.writeto(output_image, overwrite=True)


def rotate_manual(input_image, method = 'zoom'):
    if method == 'line':
        print('Please use a line region in DS9 to show the orientation of the major axis of the galaxy')
    if method=='points':
        print('Please use two point regions in DS9 to show the orientation of the major axis of the galaxy')
    
    
    if method == 'line' or method=='points':
        if not os.path.exists('rot.reg'):  
          open("rot.reg", "w").close()

        #p = subprocess.Popen(["ds9",input_image,"-scale","histequ","-cmap","Cool","-regions","load","rot.reg","-contour","limits"])
        #p.wait()
        subprocess.call("ds9 %s -scale histequ -cmap Cool -regions %s -contour limits" % (input_image, "rot.reg"), shell=True)  
        
        f = open("rot.reg", "r") 
        lines = f.readlines()
        if method=='line':
            for line in lines:
              if 'line(' in line:
                coords = np.array(line.split('(')[1].split(')')[0].split(','), float)
                x_l,y_l,x_r,y_r = coords
                break
        if method=='points':
            side='left'
            for line in lines:
              if 'point(' in line:
                coords = np.array(line.split('(')[1].split(')')[0].split(','), float)
                if side=='left':
                  x_l,y_l = coords
                  side='right'
                else:
                  x_r,y_r = coords
                
        if x_l>x_r:
          a = x_r
          b = y_r
          x_r = x_l
          y_r = y_l
          x_l = a
          y_l = b
        
        

        left_x = int(math.ceil(x_l+0.5))
        left_y = int(math.ceil(y_l+0.5))
        right_x = int(math.ceil(x_r+0.5))
        right_y = int(math.ceil(y_r+0.5))


        xc = (left_x+right_x) / 2.
        yc = (left_y+right_y) / 2.
        PosAng = math.degrees(math.atan(float(left_y - right_y)/float(left_x - right_x)))
        f.close()


    if method=='zoom':
        #ds9Proc = subprocess.Popen(["ds9",input_image,"-scale","histequ","-cmap","Cool","-contour","limits"])
        subprocess.call("ds9 %s -scale histequ -cmap Cool -contour limits" % (input_image), shell=True) 
      
        inputString = str(input("Input xc,yc,PA (In DS9: Frame 1, deg value): "))
        if inputString=='':
            print('No data was given. Exiting...')
            exit()
        inputString = inputString.split(',')
        if inputString:
            xc = float(inputString[0])
            yc = float(inputString[1])
            PosAng = float(inputString[2])
            print('here')
        else:
            PosAng = 0.0
        '''
        if PosAng>=180.:
            PosAng = PosAng - 180.
        
        if PosAng>270.:
          PosAng = 360. - PosAng
        if PosAng>90. and PosAng<180.:
          PosAng = 180. - PosAng
        '''
        if PosAng>90.:
            PosAng = 90. - PosAng
            
        PosAng = PosAng-90.
    return xc,yc,PosAng

def main(input_images, xc, yc, I_DN_min=None, I_DN_max=None, output_images = None, regime = 'manual'):
    # 1. Launching of the SE for background subtracted images
    # to obtain some geometric parameters of the galaxy
    
    print('Rotation will be done counterclockwise from the x-axis')
    galaxy_image = input_images[0]
    if len(input_images)>1:
        mask_image = input_images[1]
    else:
        mask_image = None


    if regime == 'auto':
        print('Computing the position angle to rotate the galaxy...')
        #call_SE(galaxy_image)
        run_SExtractor.call_SE(galaxy_image, snr=None, min_pix=None, sextr_dir=None, sextr_setup='default.sex', sextr_param='default.param', output_cat='field.cat', checkimage_type='SEGMENTATION',checkimage_name='segm.fits', sextr_add_string=None,determine_sky=False, m0=28.0, GAIN=10., pix2sec=1., fwhm=4.)

        galN, xCen, yCen, ellA, ellB, ellPA = get_galaxy_params(galaxy_image, xc, yc)

        if I_DN_min is None:
            backg_mean,backg_std = get_backgr_params(galaxy_image)
            I_DN_min = backg_mean + 5.*backg_std
    
        
        # Now we want to estimate position angle according to the SE results
        mask_all_except_galaxy(xc, yc, galaxy_image, 'segm.fits', mask_image, I_DN_min, I_DN_max, output_image='galaxy_mask_tmp.fits')
    
        
        # Use max pixels method
        optPosAng = ellPA#get_pa('galaxy_mask_tmp.fits', xc, yc, ellPA, ellA, ellB)
        print("SE Position angle = %1.1f" % (optPosAng))
        #optPosAng = ellPA
        #exit()
    else:
        xc,yc,optPosAng = rotate_manual(galaxy_image, method = 'line')
    
    
    
    
    
    
    if output_images is None:
        output_images = []
        for k in range(len(input_images)):
            output_images.append(input_images[k].split('.fits')[0]+'_rot.fits')



    xCenRot, yCenRot = rotate_image.main(galaxy_image, optPosAng, xc=xc, yc=yc, output_image=output_images[0], hdu_inp=0, cval=float('nan'))


    # Now lets show one of rotated images, so user can see if rotation was ok
    # and make some corrections if it is nessesary
    deltaPosAng = 0.
    answer='no'
    while answer=='no':
            subprocess.call("ds9 %s -scale histequ %s -scale histequ" % (galaxy_image, output_images[0]), shell=True)  
            answer = 'yes'
            answer = str(input('Are you happy with this (YES/no)?') or 'yes')
            if answer=='no':
                inputString = input("Correction to position angle (deg): ").strip()
                if inputString:
                    deltaPosAng = float(inputString)
                else:
                    deltaPosAng = 0.0

                xCenRot, yCenRot = rotate_image.main(galaxy_image, optPosAng-deltaPosAng, xc=xc, yc=yc, output_image=output_images[0], hdu_inp=0, cval=float('nan'))
  

    # Rotate image using the mean position angle
    for k in range(len(input_images)):
        xCenRot, yCenRot = rotate_image.main(input_images[k], optPosAng-deltaPosAng, xc=xc, yc=yc, output_image=output_images[k], hdu_inp=0, cval=float('nan'))




    print('New coordinates of the center: %f,%f' % (xCenRot, yCenRot))
    return xCenRot, yCenRot


 
#main(['NGC3628.phot.1_nonan.fits','NGC3628.1.finmask_nonan.fits','NGC3628_sigma2014.fits','new_mask.fits'], 1120.95, 1441.58, I_DN_min=2.45, I_DN_max=None, output_images = None)  

#main(['NGC4302.phot.1_nonan.fits','NGC4302.1.finmask_nonan.fits','NGC4302_sigma2014.fits'], 879.0, 1227.0, I_DN_min=0.1, I_DN_max=None, output_images = None)  

#main(['NGC891_coadd.fits','NGC891_coadd_mask.fits','NGC891_coadd_sigma.fits'], 1285.0, 1585.0, I_DN_min=0.2, I_DN_max=None, output_images = None) 

#main(['model.fits','mask.fits','sigma.fits'], 879.0, 1227.0, I_DN_min=1.6, I_DN_max=None, output_images = None) 

'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Rebinning")
    parser.add_argument("input_images", help="Input images, separated by comma: first galaxy image, then mask image, then others.")    
    parser.add_argument("xc", help="X-center of the galaxy", type=float)
    parser.add_argument("yc", help="Y-center of the galaxy", type=float)
    parser.add_argument("--output_images", help="Optional: Output images, separated by comma: first galaxy image, then mask image, then others.", type=str, default=None) 
    parser.add_argument("--min_DN", help="Optional: Outermost isophote. By default 5 sigma.", type=float, default=None)
    parser.add_argument("--max_DN", help="Optional: Innermost isophote. By default None", type=float, default=None)
    args = parser.parse_args()

    input_images = args.input_images
    xc = args.xc
    yc = args.yc
    output_images = args.output_images
    I_DN_min = args.min_DN
    I_DN_max = args.max_DN
    
    input_images = input_images.split(',')
    output_images = output_images.split(',')
    main(input_images, xc, yc, I_DN_min=I_DN_min, I_DN_max=I_DN_max, output_images = output_images)
'''
