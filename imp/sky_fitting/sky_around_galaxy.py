#!/usr/bin/python
# DESCRIPTION:
# A script to fit sky background in an annulus around galaxy. An outer galaxy ellipse --galaxy_ellipse (where outermost isophotes end) is required, either as a DS9 region file (ellipse region in images coordinates), or in the format xc,yc,sma,smb,PA. The width of the annulus is controlled by the key --annulus_width (in pix). To check the annulus in DS9 before starting the fitting the key --manual should be given. For now, only constant background level (--degree 0) within the annulus is computed.
# MINIMAL USAGE: python sky_around_galaxy.py [input_image] --galaxy_ellipse [region_file OR xc,yc,sma,smb,PA]

from astropy.modeling import models, fitting
import warnings
warnings.filterwarnings("ignore")
from astropy.stats import sigma_clipped_stats
import sys
import subprocess
from astropy.io import fits as pyfits
from pylab import *
import itertools
import os
from os.path import exists
from os import remove
from scipy.spatial import cKDTree
from scipy.optimize import fmin_tnc, fmin
import argparse
from astropy.modeling import models, fitting
from datetime import datetime
import shutil
from astropy import log
import warnings
import numpy as np
from photutils import EllipticalAnnulus
import matplotlib.pyplot as plt
from types import SimpleNamespace
import random
import scipy as sp

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


class ObjParams:
    def __init__(self, xCen, yCen, ellA, ellB, PA):
        self.cen = Point(xCen, yCen)
        self.ellA = ellA
        self.ellB = ellB
        self.PA = PA


def read_region(reg_file):
    f = open(reg_file, "r")    
    for line in f:
        if "ellipse" in line:
            params = line.split(",")
            cen = [float(params[0].split('(')[1]),float(params[1])]
            ellA = float(params[2])
            ellB = float(params[3])
            ellPA = float(params[4].split(')')[0])
            if ellA < ellB:
                ellA, ellB = ellB, ellA
                ellPA += 90    
            break
    f.close()
    return cen[0],cen[1],ellA,ellB,ellPA

def rot_point(p, orig, angle):
    """ performs rotation of a point (x,y) around point (x0,y0)
        by angle radians"""
    x1 = cos(angle)*(p.x-orig.x) - sin(angle)*(p.y-orig.y) + orig.x
    y1 = sin(angle)*(p.x-orig.x) + cos(angle)*(p.y-orig.y) + orig.y
    return Point(x1, y1)

def crea_annulus(input_image,inner_ellipse,outer_ellipse):
    [cen_in,sma_in,smb_in,PA_in] = inner_ellipse
    [cen_out,sma_out,smb_out,PA_out] = outer_ellipse

    fout = open("annulus.reg", "w")
    fout.truncate(0)
    fout.write("# Region file format: DS9 version 4.1\n")
    fout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 ')
    fout.write('normal roman" select=1 highlite=1 dash=0 fixed=0 ')
    fout.write('edit=1 move=1 delete=1 include=1 source=1\n')
    fout.write('image\n')
    fout.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # ellipse annulus\n" % (cen_in.x, cen_in.y,
                                                                 sma_in, smb_in, sma_out, smb_out,
                                                                 PA_in))    
    
    
    
    fout.close()      
    ds9Proc = subprocess.Popen(["ds9", input_image,
                                "-regions", "annulus.reg",
                                "-scale", "histequ"])
    ds9Proc.wait()      

def read_annulus(annulus_file="annulus.reg"):
   ellipses = []
   for line in open(annulus_file):
        if "ellipse" in line:
          params = line.split(",")
          cen = Point(float(params[0].split('(')[1]),
		      float(params[1]))
          sma_in = float(params[2])
          smb_in = float(params[3])
          sma_out = float(params[4])
          smb_out = float(params[5])
          PA = float(params[6].split(')')[0])
   inner_ellipse = [cen,sma_in,smb_in,PA]
   outer_ellipse = [cen,sma_out,smb_out,PA]
   ann_width = abs(sma_out-sma_in)
   return inner_ellipse,outer_ellipse,ann_width



def angles_in_ellipse(
        num,
        a,
        b):
    assert(num > 0)
    #assert(a > b)
    angles = 2 * np.pi * np.arange(num) / num
    if a != b:
        e = (1.0 - b ** 2.0 / b ** 2.0) ** 0.5
        tot_size = sp.special.ellipeinc(2.0 * np.pi, e)
        arc_size = tot_size / num
        arcs = np.arange(num) * arc_size
        res = sp.optimize.root(
            lambda x: (sp.special.ellipeinc(x, e) - arcs), angles)
        angles = res.x 
    return angles




def estimate_sky_uncertainty(ellipse, data, mask=None, width=30.):

    # Function to estimate difference of the median sky within a row of equal boxes around the galaxy
    # TODO: re-write function
    '''
    ini_exc_anomaly = np.linspace(0, 2*pi, int(10*ellipse.sma))
    exc_anomaly = random.sample(list(ini_exc_anomaly), 100)

    sinpa = sin(ellipse.pa)
    cospa = cos(ellipse.pa)
    cose = np.cos(exc_anomaly)
    sine = np.sin(exc_anomaly)
    q_value = 1 - ellipse.eps
    centers_x = []
    centers_y = []
    if mask is None:
      for i in range(len(exc_anomaly)):
        fx, ix = modf(ellipse.x0 + ellipse.sma * cose[i] * cospa - ellipse.sma * q_value * sine[i] * sinpa)
        fy, iy = modf(ellipse.y0 + ellipse.sma * q_value * sine[i] * cospa + ellipse.sma * cose[i] * sinpa)
        iix = int(ix)
        iiy = int(iy)
        centers_x.append(iix)
        centers_y.append(iiy)

    else:
      for i in range(len(exc_anomaly)):
        fx, ix = modf(ellipse.x0 + ellipse.sma * cose[i] * cospa - ellipse.sma * q_value * sine[i] * sinpa)
        fy, iy = modf(ellipse.y0 + ellipse.sma * q_value * sine[i] * cospa + ellipse.sma * cose[i] * sinpa)
        iix = int(ix)
        iiy = int(iy)
        try:
          if mask[iiy][iix]==0.: #### Take into account only those pixels which are not masked
            centers_x.append(iix)
            centers_y.append(iiy)
        except:
            zz=1
    '''
    smb = ellipse.sma * (1.-ellipse.eps)
    '''
    n=100
    t = np.arange(0, 2.*math.pi, 2.*math.pi/n)
    e = 1. - smb**2/((ellipse.sma)**2)
    
    phi = t - (e**2/8. + e**4/16. + 71.*e**6/2048.)*sin(2.*t) + (5.*e**4/256. + 5.*e**6/256.)*sin(4.*t) + 29.*e**6/6144.*sin(6.*t)
    xx = ellipse.sma * np.cos(phi)
    yy = smb * np.sin(phi)
    centers_x = np.array(ellipse.x0 + xx * cos(ellipse.pa) - yy * sin(ellipse.pa), int)
    centers_y = np.array(ellipse.y0 + yy * cos(ellipse.pa) + xx * sin(ellipse.pa), int)    

    
    fig = plt.figure()
    ax = fig.gca()
    ax.axes.set_aspect('equal')
    #ax.plot(centers_x, centers_y, 'o', color='red', markersize=1.)
    '''
    n=100
    smb = ellipse.sma * (1.-ellipse.eps)
    phi = angles_in_ellipse(n, ellipse.sma, smb)
    xx = ellipse.sma * np.cos(phi)
    yy = smb * np.sin(phi)
    centers_x = np.array(ellipse.x0 + xx * cos(ellipse.pa) - yy * sin(ellipse.pa), int)
    centers_y = np.array(ellipse.y0 + yy * cos(ellipse.pa) + xx * sin(ellipse.pa), int)
    #print(centers_x)
    #print(centers_y)
    #exit()
    #ax.plot(centers_x, centers_y, 'o', alpha=0.5, markersize=1.)
    #plt.show()
    #exit()
    
    medians = []
    mask_astropy = (mask!=0.)
    #plt.plot(centers_x, centers_y, 'o')
    #plt.show()
    #exit()
    for k in range(len(centers_x)):
      try:
        mean_sigclip, median_sigclip, std_sigclip = sigma_clipped_stats(data[centers_y[k]-int(width/2.):centers_y[k]+int(width/2.), centers_x[k]-int(width/2.):centers_x[k]+int(width/2.)], mask=mask_astropy[centers_y[k]-int(width/2.):centers_y[k]+int(width/2.), centers_x[k]-int(width/2.):centers_x[k]+int(width/2.)], sigma=3)
        medians.append(median_sigclip)
      except:
          zz=1
    
    unc = std(medians)
    if np.isnan(unc):
        unc = 0.
    
    return unc





def ellipse_annulus(output_image, inner_ellipse, outer_ellipse, inframe, mask, xSize, ySize, header, sky_subtr, degree, show_running=True, ann_width=30.):
        
        cen_in, ellA_in, ellB_in, ellPA_in = inner_ellipse
        cen_out, ellA_out, ellB_out, ellPA_out = outer_ellipse
        
        annulus_aperture = EllipticalAnnulus([(cen_in.x,cen_in.y)], ellA_in, ellA_out, ellB_out, theta=np.radians(ellPA_in))
        annulus_mask = annulus_aperture.to_mask(method='center')

        
        annulus_data = annulus_mask[0].multiply(inframe)
        annulus_data_1d = annulus_data[annulus_data != 0]


        
        mask_new = np.ones((ySize, xSize))
        mask_new[mask > 0.]=0.
        
        annulus_mask_mask = annulus_mask[0].multiply(mask_new)
        annulus_mask_mask_1d = annulus_mask_mask[annulus_data != 0]
        mask_astropy = (annulus_mask_mask_1d==0.)

        
        
        mean_sigclip, median_sigclip, std_sigclip = sigma_clipped_stats(annulus_data_1d, mask=mask_astropy, sigma=3)  
        bkg_median = median_sigclip
        print('Sky inside of the annulus (mean,med,std): %.8f, %.8f, %.8f' % (mean_sigclip, median_sigclip, std_sigclip))


        if sky_subtr==True:
            sframe = inframe - bkg_median
        else:
            sframe = inframe
        
        
        header['Sky_medg'] = bkg_median
        header['Sky_stdg'] = std_sigclip

        header['Ann_x0'] = cen_in.x
        header['Ann_y0'] = cen_in.y
        header['Ann_sma'] = ellA_in
        header['Ann_smb'] = ellB_out*(ellA_in/ellA_out)
        header['Ann_PA'] = ellPA_in
        
        
        ellipse = SimpleNamespace(x0=cen_in.x, y0=cen_in.y, sma=(ellA_in+ellA_out)/2., intens=0,
                                  eps=1.-ellB_in/ellA_in, pa=np.radians(ellPA_in), grad=0, a3=0, b3=0, a4=0, b4=0)
        

        bkg_median_std = estimate_sky_uncertainty(ellipse, sframe, mask=mask, width=ann_width)

        
        header['Sky_medg_std'] = bkg_median_std
        print('Std of the median background around galaxy: %.8f' % (bkg_median_std))
        

        if output_image is None:
            output_image = 'sky_subtr_galaxy.fits'
        outHDU = pyfits.PrimaryHDU(sframe, header=header)
        outHDU.writeto(output_image, overwrite=True)  

        outHDU = pyfits.PrimaryHDU(bkg_median*np.ones((ySize, xSize)), header=header)
        outHDU.writeto('sky_galaxy.fits', overwrite=True)
        print('Done!')
        return bkg_median,std_sigclip,bkg_median_std


def ellipse_annulus1(output_image, inner_ellipse, outer_ellipse, inframe, mask, xSize, ySize, header, sky_subtr, degree, show_running=True):
        # This approach is more time consuming than above!
        def define_ellipse(cen, ellA, ellB, ellPA):
            cospa = cos(radians(ellPA))
            sinpa = sin(radians(ellPA))
            # Check if whole ellipse is inside of the image
            # and obtain size of the ellipse in xy plane
            xMax = 0.0
            yMax = 0.0
            xMin = 1e10
            yMin = 1e10
            for e in linspace(0, 4*pi, 1000):
                cose = cos(e)
                sine = sin(e)
                x = cen.x + ellA * cose * cospa - ellB * sine * sinpa
                y = cen.y + ellB * sine * cospa + ellA * cose * sinpa
                if x > xMax:
                    xMax = x
                if y > yMax:
                    yMax = y
                if x < xMin:
                    xMin = x
                if y < yMin:
                    yMin = y
            xMin = max(0, int(round(xMin)))
            xMax = min(xSize, int(round(xMax)))
            yMin = max(0, int(round(yMin)))
            yMax = min(ySize, int(round(yMax)))
            
            
            focusR = (ellA ** 2.0 - ellB ** 2.0) ** 0.5
            focus10 = Point(cen.x + focusR, cen.y)  # Unrotated
            focus20 = Point(cen.x - focusR, cen.y)  #
            focus1 = rot_point(focus10, cen, radians(ellPA))
            focus2 = rot_point(focus20, cen, radians(ellPA))
            return xMin,xMax,yMin,yMax,focus1,focus2

        cen_in, ellA_in, ellB_in, ellPA_in = inner_ellipse
        cen_out, ellA_out, ellB_out, ellPA_out = outer_ellipse

        ## Inner ellipse
        xMin_in,xMax_in,yMin_in,yMax_in,focus1_in,focus2_in = define_ellipse(cen_in, ellA_in, ellB_in, ellPA_in)

        ## Outer ellipse
        xMin_out,xMax_out,yMin_out,yMax_out,focus1_out,focus2_out = define_ellipse(cen_out, ellA_out, ellB_out, ellPA_out)

        xMin = min([xMin_in,xMin_out])
        yMin = min([yMin_in,yMin_out])

        xMax = max([xMax_in,xMax_out])
        yMax = max([yMax_in,yMax_out])

        # Find pixels inside of the ellipse annulus
        dEll_in = 2 * ellA_in
        dEll_out = 2 * ellA_out
        
        I_an = []
        I_new_mask = (mask>0.)
        I_dan = np.copy(inframe)
        X = []
        Y = []
        Z = []
        for x in range(xMin, xMax):
            for y in range(yMin, yMax):
                dFocus1_in = hypot(x-focus1_in.x, y-focus1_in.y)
                dFocus2_in = hypot(x-focus2_in.x, y-focus2_in.y)
                dPoint_in = dFocus1_in + dFocus2_in
                if dPoint_in <= dEll_in:
                    continue

                dFocus1_out = hypot(x-focus1_out.x, y-focus1_out.y)
                dFocus2_out = hypot(x-focus2_out.x, y-focus2_out.y)
                dPoint_out = dFocus1_out + dFocus2_out                
                if dPoint_out < dEll_out and mask[y,x]==0. and np.isnan(inframe[y,x])==False:
                  try:
                    I_an.append(inframe[y,x])
                    I_dan[y,x] = 1000.
                    X = np.append(X, float(x))
                    Y = np.append(Y, float(y))
                    Z = np.append(Z, float(inframe[y,x]))
                  except:
                    zz=1
        if show_running==True:
            print('Number of pixels in the annulus is', len(I_an))
        outHDU = pyfits.PrimaryHDU(I_dan)
        outHDU.writeto('annulus.fits',overwrite=True)

        mean, median, std = sigma_clipped_stats(I_an, sigma=3.0)
        sky_median = median
        print(median, std)
        '''
        #### TODO: DON'T USE THIS OPTION BELOW
        if degree!=0:
          degree = degree + 1
          if degree>0:
            # Fit a 3rd degree, 2d polynomial
            m = polyfit2d(X,Y,Z,degree,True)
            
            # Evaluate it on a grid...
            nx, ny = 20, 20
            xx, yy = np.meshgrid(np.linspace(0, xSize, nx),
				np.linspace(0, ySize, ny))
            zz = polyval2d(xx, yy, m)

            # Make a sky array and subtract it from the object frame
            sky = np.copy(inframe)
            strip = np.arange(1.0,np.float(xSize+1.0))
            for i in range(0,ySize):
                sky[i,:] = polyval2d(strip,i,m)
        else:
            sky = sky_median
        '''
        sky = sky_median
        
        
        if sky_subtr==True:
            sframe = inframe - sky
        else:
            sframe = inframe

        header['Sky_medg'] = sky_median
        header['Sky_stdg'] = std

        header['Ann_x0'] = cen_in.x
        header['Ann_y0'] = cen_in.y
        header['Ann_sma'] = ellA_in
        header['Ann_smb'] = ellB_in
        header['Ann_PA'] = ellPA_in

        if output_image is None:
            output_image = 'sky_subtr_galaxy.fits'
        outHDU = pyfits.PrimaryHDU(sframe, header=header)
        outHDU.writeto(output_image, overwrite=True)  

        outHDU = pyfits.PrimaryHDU(sky, header=header)
        outHDU.writeto('sky_galaxy.fits', overwrite=True)
        print('Done!')
        return sky_median,std



def main(input_image, output_image, mask_image, manual=False, degree=0, annulus_width=None, galaxy_ellipse=None, Npix=100000, offset_size=1., show_running=True, galaxy_annulus=None):
    print('Estimating sky within an elliptical annulus...')
    hdulist_ima = pyfits.open(input_image)
    data_image = hdulist_ima[0].data
    header = hdulist_ima[0].header
    ySize,xSize = np.shape(data_image)
    
    if mask_image is not None:
        hdulist_mask = pyfits.open(mask_image)
        data_mask = hdulist_mask[0].data
    else:
        data_mask = np.zeros(shape=(ySize, xSize))
 
    
    '''
    mask_astropy = np.zeros_like(mask_infr, dtype=bool)
    nx, ny =mask_infr.shape[1], mask_infr.shape[0]
    for k in range(ny):
      for i in range(nx):
	if mask_infr[k,i]>0:
	  mask_astropy[k,i] = True
	else:
	  mask_astropy[k,i] = False
    '''
    if galaxy_annulus is None:
        if galaxy_ellipse is not None:
            if os.path.exists(galaxy_ellipse):
                xc_in,yc_in,sma_in,smb_in,PA_in = read_region(galaxy_ellipse)
            else:
                [xc_in, yc_in, sma_in, smb_in, PA_in] = galaxy_ellipse
            
            sma_in = sma_in*offset_size
            smb_in = smb_in*offset_size
            inner_ellipse = [Point(xc_in,yc_in),sma_in,smb_in,PA_in]
        else:
            print('ERROR! No galaxy ellipse is given! Exiting...')
            exit()
    

    

        if annulus_width is None:  
            # Find width of the annulus: the annulus should include Npix pixels
            ann_width = ceil(0.5*(-(sma_in+smb_in)+sqrt((sma_in+smb_in)**2+4.*Npix*2./math.pi)))
        else:
            ann_width = float(annulus_width)
            NN = int(ceil(math.pi/8. * (4.*ann_width**2 + 4.*ann_width*(sma_in+smb_in))))
            if NN<1000:
                NN = 1000
                ann_width = ceil(0.5*(-(sma_in+smb_in)+sqrt((sma_in+smb_in)**2+4.*NN*2./math.pi)))
                if show_running==True:
                    print('Enetered width of the annulus is too low to contain 1000 pixels, therefore it was changed!')
      

        if show_running==True: 
            print('Width of the annulus is', ann_width)
    
        # Outer ellipse (calculated)
        xc_out = xc_in
        yc_out = yc_in
        sma_out = sma_in + ann_width
        smb_out = smb_in + ann_width #*1.5 for NGC4013
        PA_out = PA_in
        outer_ellipse = [Point(xc_out,yc_out),sma_out,smb_out,PA_out]
        
        if manual==True:
                crea_annulus(input_image,inner_ellipse,outer_ellipse)
                inner_ellipse,outer_ellipse,ann_width = read_annulus()
    else:
        if manual==True:    
            ds9Proc = subprocess.Popen(["ds9", input_image,
                                        "-regions", galaxy_annulus,
                                        "-scale", "histequ"])
            ds9Proc.wait()               
        inner_ellipse,outer_ellipse,ann_width = read_annulus(galaxy_annulus)
    
    
    sky,std,bkg_median_std = ellipse_annulus(output_image, inner_ellipse, outer_ellipse, data_image, data_mask, xSize, ySize, header, True, degree, show_running=show_running, ann_width=ann_width)

    return sky,std,bkg_median_std

    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky background estimation")
    parser.add_argument("inputImage", help="Input fits image with the object")
    parser.add_argument("--maskImage", nargs='?', const=1, help="Optional: Input fits image with masked objects", type=str, default=None) 
    parser.add_argument("--outputImage", nargs='?', const=1, help="Optional: Output image with subtracted sky", type=str, default=None) 
    parser.add_argument("--degree", nargs='?', const=1, help="Optional: Input the polinomial order to fit the sky within the annulus. NOT WORKING PROPERLY!", type=int, default=0)
    parser.add_argument("--manual", nargs='?', const=1, help="Optional: Do you want to manually change the annulus in DS9", type=bool, default=False)
    parser.add_argument("--annulus_width", nargs='?', const=1, help="The width of the annulus in pix.", type=float, default=30)
    parser.add_argument("--galaxy_ellipse", nargs='?', const=1, help="Galaxy ellipse, either as a DS9 region file with an ellipse region, or in the format xc,yc,sma,smb,PA", type=str, default=None)
    args = parser.parse_args()

    input_image = args.inputImage
    mask_image = args.maskImage
    output_image = args.outputImage
    degree = args.degree
    manual = args.manual
    annulus_width = args.annulus_width
    galaxy_ellipse = args.galaxy_ellipse
    main(input_image, output_image, mask_image, manual=manual, degree=degree, annulus_width=annulus_width, galaxy_ellipse=galaxy_ellipse, Npix=100000, show_running=True)
