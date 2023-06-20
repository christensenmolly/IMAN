#!/usr/bin/python
# DESCRIPTION:
# Script to do photometric calibration. This script finds good PSF stars
# and searches them in a catalogue.
# MINIMAL USAGE: python calibration.py [input_image]

# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import math
import sys

from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import os
#import subprocess
from photutils import CircularAnnulus, EllipticalAnnulus
from photutils import CircularAperture, EllipticalAperture
from photutils import aperture_photometry
from astroquery.simbad import Simbad
from astropy import coordinates
import astropy.units as u
from astropy import wcs
import collections
from astroquery.vizier import Vizier
from astropy.stats import sigma_clipped_stats
from photutils import find_peaks
from astropy.stats import sigma_clip

# Import additional modules
LOCAL_DIR = "/imp/phot_calibration"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'sextractor'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))
sys.path.append(os.path.join(IMAN_DIR, 'decomposition/simple_fitting'))

import run_SExtractor
import rebin_image
import add_keyw_to_header
import photometric_conversions

from astropy.stats import sigma_clipped_stats

#PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
#PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
#sys.path.append(PATH_TO_PACKAGE+'/Sextractor')


FNULL = open(os.devnull, 'w')

tmp_out = sys.stdout

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx],idx

def find_in_cat(x_image,y_image,xc_s,yc_s,FWHM):
    # Find nearest star in the SExtractor catalogue
    r = np.sqrt( (x_image-xc_s)**2 + (y_image-yc_s)**2 )
    R_near,idx_near = find_nearest(r, 0.)
    if R_near>FWHM:
      print('The star (%f,%f) was not found!' % (xc_s,yc_s))
      return float('nan')
    else:
      print('The star (%f,%f) corresponds to %i in the SExtractor catalogue.' % (xc_s,yc_s,idx_near+1))
      return idx_near
    
    

def pointInEllipse(x,y,xp,yp,d,D,angle):
    #tests if a point[xp,yp] is within
    #boundaries defined by the ellipse
    #of center[x,y], diameter d D, and tilted at angle
    angle = np.radians(angle)
    cosa=math.cos(angle)
    sina=math.sin(angle)
    dd=d/2*d/2
    DD=D/2*D/2

    a =math.pow(cosa*(xp-x)+sina*(yp-y),2)
    b =math.pow(sina*(xp-x)-cosa*(yp-y),2)
    try:
        ellipse=(a/dd)+(b/DD)

        if ellipse <= 1:
            return True
        else:
            return False
    except:
        return True

def get_true_ra_dec(RA,DEC,Radius,catalog,band,Table=None):
        if Table is None:
            viz = Vizier(keywords=["stars", "optical"])
            viz.ROW_LIMIT = -1
            try:
                c = coordinates.SkyCoord(float(RA),float(DEC),unit=('deg','deg'),frame='fk5')
                result = viz.query_region(c, radius=Radius*u.deg, catalog=catalog)
                table = result[0]
            except:
                print('ERROR: The reference coordinate is not found in the catalogue! Exiting ...')
                exit()
            Table = []; RAA = []; DECC = []; MAG = []
            for i in range(len(table)):
                RAA.append(float(table["RAJ2000"][i])) #### WARNING: was _RAJ2000
                DECC.append(float(table["DEJ2000"][i])) #### DEJ2000
                if catalog=='NOMAD':
                  MAG.append(float(table[band + "mag"][i]))
                elif catalog=='V/139':
                  u_mag = float(table["umag"][i])
                  g_mag = float(table["gmag"][i])
                  r_mag = float(table["rmag"][i])
                  i_mag = float(table["imag"][i])
                  z_mag = float(table["zmag"][i])
                  B,V,R,I = photometric_conversions.transformation_from_sdss_to_BVRI(u_mag,g_mag,r_mag,i_mag,z_mag)
                  MAG.append(float(photometric_conversions.get_mag(u_mag,g_mag,r_mag,i_mag,z_mag,B,V,R,I,band)))
            RAA = np.array(RAA, float)
            DECC = np.array(DECC, float)
            MAG = np.array(MAG, float)
            Table = np.array([RAA,DECC,MAG], float)

        dist = np.sqrt( (RA-Table[0])**2 + (DEC-Table[1])**2) * 3600.
        dist = np.array(dist, float)

        min_dist = dist.min()
        index = np.where(dist == min_dist)[0][0]
        return Table,Table[0][index],Table[1][index],Table[2][index],index

def rot_point(p, orig, angle):
    """ performs rotation of a point (x,y) around point (x0,y0)
        by angle radians"""
    x1 = cos(angle)*(p.x-orig.x) - sin(angle)*(p.y-orig.y) + orig.x
    y1 = sin(angle)*(p.x-orig.x) + cos(angle)*(p.y-orig.y) + orig.y
    return Point(x1, y1)

# Create mask image based on the Sextractor catalogue
def ellipse_mask(cen, ellA, ellB, ellPA, num_ellipse, inframe, xSize, ySize):
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
        # Find pixels inside of the ellipse
        dEll = 2 * ellA
        

        for x in range(xMin, xMax+1):
            for y in range(yMin, yMax+1):
                  dFocus1 = hypot(x-focus1.x, y-focus1.y)
                  dFocus2 = hypot(x-focus2.x, y-focus2.y)
                  dPoint = dFocus1 + dFocus2
                  if dPoint < dEll:
                    try:
                      inframe[y,x] = inframe[y,x] + num_ellipse
                    except:
                      zz=1


# Find mean background in annulus:
def back_in_annulus(data, mask, cen, ellA, ellB, ellPA, num_ellipse, xSize, ySize):
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
        # Find pixels inside of the ellipse
        dEll = 2 * ellA
        
        I_back = []
        for x in range(xMin, xMax+1):
           for y in range(yMin, yMax+1):
                  dFocus1 = hypot(x-focus1.x, y-focus1.y)
                  dFocus2 = hypot(x-focus2.x, y-focus2.y)
                  dPoint = dFocus1 + dFocus2
                  if dPoint < dEll and mask[y,x]==0.:
                    try:
                      I_back.append(data[y,x])
                    except:
                      zz=1
        mean, median, std = sigma_clipped_stats(I_back, sigma=3.0)
        return median, std


# Select non-overlaping ellipses
def select_good_ellipses(cen, ellA, ellB, ellPA, num_ellipse, mask, xSize, ySize):
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
        # Find pixels inside of the ellipse
        dEll = 2 * ellA
        
        I_ell_area = 0.
        I_best = 0.
        for x in range(xMin, xMax+1):
           for y in range(yMin, yMax+1):
                  dFocus1 = hypot(x-focus1.x, y-focus1.y)
                  dFocus2 = hypot(x-focus2.x, y-focus2.y)
                  dPoint = dFocus1 + dFocus2
                  if dPoint < dEll:
                    try:
                      I_ell_area = I_ell_area + mask[y,x]
                      I_best = I_best + num_ellipse
                    except:
                      zz=1
        if I_best==I_ell_area:
          return True
        else:
          return False


def ellipse_borders(xy,width,height,angle):
  width = width 
  height = height 
  angle = np.radians(angle)
  
  X = math.sqrt( (width*math.cos(angle))**2 + (height*math.sin(angle))**2 )
  Y = math.sqrt( (width*math.sin(angle))**2 + (height*math.cos(angle))**2 )
  x_max = xy[0] + X
  x_min = xy[0] - X
  y_max = xy[1] + Y
  y_min = xy[1] - Y
  return [[x_min,y_min],[x_max,y_max]]



def main(input_image, output_image, band='R',catalogue='NOMAD', m0=None, sex_cat=None, region_file=None, peaks_file=None, mag_min=15.6, mag_max=18.3, verbosity=True):
  if verbosity: print('Using catalogue %s' % (catalogue))
  
  # Open the input image:
  hdulist = pyfits.open(input_image)
  data = hdulist[0].data
  header = hdulist[0].header
  ny,nx = np.shape(data)

  data_copy = np.copy(data)
  
  if peaks_file==None or os.path.exists(peaks_file)==False:
    x_reg = open('x.reg','w')
    mean_data, median_data, std_data = sigma_clipped_stats(data_copy, sigma=3.0)
    threshold = median_data + (3.0 * std_data)
    tbl = find_peaks(data_copy, threshold, box_size=5, subpixel=False)
    x_peaks = []; y_peaks = []
    f_peaks = open('peaks.txt', 'w')
    for k in range(len(tbl)):
          x_peaks.append(float(tbl['x_peak'][k])+1.)
          y_peaks.append(float(tbl['y_peak'][k])+1.)
          f_peaks.write('%f\t%f\n' % (x_peaks[k],y_peaks[k]))
          x_reg.write('point(%f,%f) # point=x\n' % (x_peaks[k],y_peaks[k]))
  else:
    x_peaks,y_peaks = np.loadtxt(peaks_file, usecols=[0,1], dtype=float, unpack=True, skiprows = 0, delimiter='\t')

  hdulist.close()


  ####NOTE #######
  delta_R = 1. #Add pixels to mask
  dx = nx/10. # Stay away from the border
  dy = ny/10. # Stay away from the border
  class_star_lim = 0.95
  GAIN = 10
  fwhm = 4.
  #mag_min = 15.6
  #mag_max = 18.3
  ###########  

  # Create zero-mask array
  mask = np.zeros(shape=(ny,nx))




  # Determine image resolution (pix/arcsec)
  pix2sec,note = rebin_image.resolution(input_image)
  pix2sec = float(pix2sec)
  if note=='*':
    if verbosity: print('WCS is not found. Exiting ...')
    exit()
  else:
    if verbosity: print('Resolution is: %.3f arcsec/pix' % (pix2sec))
  
  Radius = pix2sec*math.sqrt(nx**2+ny**2)/2. / 3600. # in deg


  if sex_cat==None or os.path.exists(sex_cat)==False:
    # Run Sextractor
    run_SExtractor.call_SE(input_image, snr=None, min_pix=None, sextr_dir=None, sextr_setup='default_fwhm.sex', sextr_param='default_fwhm.param', output_cat='field.cat', checkimage_type='SEGMENTATION',checkimage_name='segm.fits', sextr_add_string=None,determine_sky=False, m0=20.0, GAIN=GAIN, pix2sec=pix2sec, fwhm=fwhm, verbosity=verbosity)
    sex_cat = 'field.cat'
    #exit()
  # Read Sextractor catalogue
  number,x_image,y_image,x_world,y_world,XPEAK_IMAGE,YPEAK_IMAGE,XPEAK_WORLD,YPEAK_WORLD,flux_radius25,flux_radius50,flux_radius75,flux_radius99,mag_auto,a_image,b_image,theta_image,ellipticity,kron_radius,backgr,class_star,fwhm = np.loadtxt(sex_cat, usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21], unpack=True, skiprows = 19, dtype=float)
  

  

  # Fill in the mask for ALL objects
  for k in range(len(number)):
    num_ellipse = float(number[k])
    xCen = float(x_image[k])
    yCen = float(y_image[k])
    cen = Point(xCen,yCen)
    kron = float(kron_radius[k])
    ellA = kron * float(a_image[k]+delta_R)
    ellB = kron * float(b_image[k]+delta_R)
    ellPA = float(theta_image[k])  
    ellipse_mask(cen, ellA, ellB, ellPA, num_ellipse, mask, nx, ny)
  
  hdu = pyfits.PrimaryHDU(mask, header)
  hdu.writeto('mask.fits',clobber=True)  

  if region_file==None or os.path.exists(region_file)==False:
        # Select stars:
        cl_stars = []
        '''
        # By statistics:
        class_star_filter = sigma_clip(class_star, sigma=3, iters=5)

        cl_med = np.median(class_star_filter)
        cl_std = np.std(class_star_filter)

        rad_med = np.median(flux_radius99)
        rad_std = np.std(flux_radius99) 
        
        for k in range(len(class_star)):
	  if class_star[k]<=cl_med+cl_std and class_star[k]>=cl_med-cl_std and flux_radius99[k]<=rad_med+rad_std and flux_radius99[k]>=rad_med-rad_std:
	    cl_stars.append(k)
        '''
        # Or by the limit value of class_star:  
        for k in range(len(class_star)):
          if class_star[k]>=class_star_lim and float(kron_radius[k])>0.:
            cl_stars.append(k)  

        if verbosity: print('Number of good stars after class_star filter: %i' % (len(cl_stars)))


        try:
          fwhm_med = np.median(fwhm[cl_stars])
        except:
          fwhm_med = np.mean(fwhm[cl_stars])
        fwhm_std = np.std(fwhm[cl_stars]) 

        try:
          ell_med = np.median(ellipticity[cl_stars])
        except:
          ell_med = np.mean(ellipticity[cl_stars])
        ell_std = np.std(ellipticity[cl_stars])   

        try:
          PA_med = np.median(theta_image[cl_stars])
        except:
          PA_med = np.mean(theta_image[cl_stars])
        PA_std = np.std(theta_image[cl_stars]) 

        if verbosity: print('FWHM [pix] = %.3f+/-%.3f' % (fwhm_med,fwhm_std)) 
        if verbosity: print('Ellipticity = %.3f+/-%.3f' % (ell_med,ell_std)) 
        if verbosity: print('PA [deg] = %.3f+/-%.3f' % (PA_med,PA_std)) 
        
        # Select stars by FWHM and ellipticity
        fwhm_stars = []
        for k in range(len(cl_stars)):
          K = cl_stars[k]
          if fwhm[K]<=fwhm_med+fwhm_std and fwhm[K]>=fwhm_med-fwhm_std and ellipticity[K]<=ell_med+ell_std and ellipticity[K]>=ell_med-ell_std: # and ellipticity[K]<0.2: #and ellipticity[K]<=ell_med+ell_std and ellipticity[K]>=ell_med-ell_std:
            fwhm_stars.append(K)

        if verbosity: print('Number of good stars after FWHM filter: %i' % (len(fwhm_stars)))
        mag_stars = fwhm_stars



        BEST_STARS = []
        fout = open("stars.reg", "w")
        fout.truncate(0)
        fout.write("# Region file format: DS9 version 4.1\n")
        fout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 ')
        fout.write('normal roman" select=1 highlite=1 dash=0 fixed=0 ')
        fout.write('edit=1 move=1 delete=1 include=1 source=1\n')
        fout.write('image\n')

        for k in range(len(mag_stars)):
          K = mag_stars[k]
          num_ellipse = float(number[K])
          xCen = float(x_image[K])
          yCen = float(y_image[K])
          cen = Point(xCen,yCen)
          kron = float(kron_radius[K])
          ellA = kron * float(a_image[K]+delta_R)
          ellB = kron * float(b_image[K]+delta_R)
          ellPA = float(theta_image[K]) 
          res = select_good_ellipses(cen, ellA, ellB, ellPA, num_ellipse, mask, nx, ny)

          if res==True:
            borders = ellipse_borders([xCen,yCen],ellA,ellB,ellPA)
            if borders[0][0]>dx and borders[0][1]>dy and borders[1][0]<nx-dx and borders[1][1]<ny-dy:
              BEST_STARS.append(K)
              fout.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # text={%i}\n" % (xCen, yCen,
									ellA, ellB,
									ellPA, K))      
        fout.close()

        if verbosity: print('Number of good non-overlaping stars away from the frame borders: %i' % (len(BEST_STARS)))

  else:
        mag_max = 1000.
        mag_min = -1000.
        f_reg = open(region_file, 'r')
        BEST_STARS = []
        for line in f_reg:
          if "ellipse" in line or "circle" in line:
                if 'text={' in line:
                  BEST_STARS.append(int(line.split('text={')[-1].split('}')[0]))
                else:
                  #print 'All regions should be subscribed!'
                  if "circle" in line:
                    xc_s = float(line.split('circle(')[-1].split(',')[0])
                    yc_s = float(line.split('circle(')[-1].split(',')[1])
                  if "ellipse" in line:
                    xc_s = float(line.split('ellipse(')[-1].split(',')[0])
                    yc_s = float(line.split('ellipse(')[-1].split(',')[1])
                  ind = find_in_cat(x_image,y_image,xc_s,yc_s,np.median(fwhm))
                  if np.isnan(ind)==False:
                    BEST_STARS.append(int(ind)) 

        f_reg.close()
        if verbosity: print('Number of found stars in the region file is: %i' % (len(BEST_STARS)))

  try:
          fwhm_med = np.median(fwhm[BEST_STARS])
  except:
          fwhm_med = np.mean(fwhm[BEST_STARS])
  fwhm_std = np.std(fwhm[BEST_STARS]) 

  try:
          rad_med = np.median(flux_radius99[BEST_STARS])
  except:
          rad_med = np.mean(flux_radius99[BEST_STARS])
  rad_std = np.std(flux_radius99[BEST_STARS])


  try:
          theta_med = np.median(theta_image[BEST_STARS])
  except:
          theta_med = np.mean(theta_image[BEST_STARS])
  theta_std = np.std(theta_image[BEST_STARS]) 

  try:
          ell_med = np.median(ellipticity[BEST_STARS])
  except:
          ell_med = np.mean(ellipticity[BEST_STARS])
  ell_std = np.std(ellipticity[BEST_STARS]) 





  # CHECK IF TWO OR MORE PEAKS LIE WITHIN THE BEST_STARS ELLIPSES AND REMOVE SUCH STARS
  FINAL_BEST_STARS = []

  for i in range(len(BEST_STARS)):
        I = BEST_STARS[i]
        xCen = float(x_image[I])
        yCen = float(y_image[I])

        kron = float(kron_radius[I])
        ellA = kron * float(a_image[I]+delta_R)
        ellB = kron * float(b_image[I]+delta_R)
        ellPA = float(theta_image[I])
        count = 0
        for k in range(len(x_peaks)):
          res_pix = pointInEllipse(xCen,yCen,x_peaks[k],y_peaks[k],ellA*2.,ellB*2.,ellPA)
  
          if res_pix==True:
            count = count + 1
          if count>=2:
            break
        if count<2:
          FINAL_BEST_STARS.append(I)

  if len(FINAL_BEST_STARS)>1 and region_file==None:
    BEST_STARS = FINAL_BEST_STARS
    comment = ''
  else:
    # Too little FINAL_BEST_STARS
    comment = ''	#'*'





  N_try = 0
  repeat = False
  while N_try<=2:
    N_try = N_try + 1  
    fout = open("best_stars.reg", "w")
    fout.truncate(0)
    fout.write("# Region file format: DS9 version 4.1\n")
    fout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 ')
    fout.write('normal roman" select=1 highlite=1 dash=0 fixed=0 ')
    fout.write('edit=1 move=1 delete=1 include=1 source=1\n')
    fout.write('image\n')
    
    fout_wcs = open("best_stars_wcs.reg", "w")
    fout_wcs.truncate(0)
    fout_wcs.write("# Region file format: DS9 version 4.1\n")
    fout_wcs.write('global color=green dashlist=8 3 width=1 font="helvetica 10 ')
    fout_wcs.write('normal roman" select=1 highlite=1 dash=0 fixed=0 ')
    fout_wcs.write('edit=1 move=1 delete=1 include=1 source=1\n')
    fout_wcs.write('fk5\n')  

    # Find background for them in the annulus with mask:
    Fluxes = []; Mags = []; M0 = []; FWHMS = []; ELLS = []; PAS = []
    if verbosity: print("Zero-point for each star:")
    f_ind = open("ZP_results.dat", "w")
    VERY_FINAL_STARS = []
    for k in range(len(BEST_STARS)):
        #try:
        K = BEST_STARS[k]

        num_ellipse = float(number[K])
        xCen = float(x_image[K])
        yCen = float(y_image[K])
        cen = Point(xCen,yCen)
        kron = float(kron_radius[K])
        ellA = kron * float(a_image[K]+delta_R)
        ellB = kron * float(b_image[K]+delta_R)
        ellPA = float(theta_image[K]) 
        
        if verbosity: print(a_image[K], b_image[K], kron, number[K])

        #aperture = CircularAperture([(x[k],y[k])], r=R_in[k])
        aperture = EllipticalAperture([(x_image[K],y_image[K])], a=kron * float(a_image[K]), b=kron * float(b_image[K]), theta=float(theta_image[K]))
        #annulus_aperture = CircularAnnulus([(x[k],y[k])], r_in=R_in[k], r_out=R_out[k])
        annulus_aperture = EllipticalAnnulus([(x_image[K],y_image[K])], a_in=kron * float(a_image[K]), a_out=1.5*kron * float(a_image[K]), b_out=1.5*kron * float(b_image[K]), theta=float(theta_image[K]))
        
        '''
        #### REDO?
        phot_table = aperture_photometry(data, annulus_aperture)#, segm_mask)
        bkg_mean = phot_table['aperture_sum'] / annulus_aperture.area()
        bkg_sum = bkg_mean * aperture.area()
        #print bkg_mean
        #### REDO?    
        '''

     
        bkg_mean,bkg_std = back_in_annulus(data, mask, cen, ellA*1.5, ellB*1.5, ellPA, num_ellipse, nx, ny)
        try:
            bkg_sum = bkg_mean * aperture.area()
        except:
            bkg_sum = bkg_mean * aperture.area
        
        phot_table = aperture_photometry(data, aperture)
        final_sum = phot_table['aperture_sum'] - bkg_sum
        
        flux = float(final_sum[0])
        
        Fluxes.append(flux)
    
        if m0 is None:
            if k==0:
                Table,RA_prec,DEC_prec,Mag_prec,ind_prec = get_true_ra_dec(x_world[K],y_world[K],Radius,catalogue,band,Table=None)
            else:
                Table,RA_prec,DEC_prec,Mag_prec,ind_prec = get_true_ra_dec(x_world[K],y_world[K],Radius,catalogue,band,Table=Table)
            if Mag_prec>=mag_min and Mag_prec<=mag_max or N_try==2:
                zp = Mag_prec + 2.5*math.log10(flux)
                Mags.append(Mag_prec)
                M0.append(zp)      

      
        else:
            if Mag_prec>=mag_min and Mag_prec<=mag_max or N_try==2:
                Mags.append(m0-2.5*math.log10(flux))
                M0.append(m0)

        if (Mag_prec>=mag_min and Mag_prec<=mag_max) or N_try==2:
            fout.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # text={%i}\n" % (xCen, yCen,
                                                                    ellA, ellB, ellA*1.5, ellB*1.5,
                                                                    ellPA, K))   

            fout_wcs.write("ellipse(%.7f,%.7f,%.7f,%.7f,%.7f,%.7f,%1.1f) # text={%i}\n" % (x_world[K],y_world[K],
                                                                    ellA*pix2sec/3600., ellB*pix2sec/3600., ellA*1.5*pix2sec/3600., ellB*1.5*pix2sec/3600.,
                                                                    ellPA, K)) 

            FWHMS.append(fwhm[K])
            ELLS.append(1.-ellB/ellA)
            PAS.append(ellPA)
            VERY_FINAL_STARS.append(K)
        #except:
        #  z=1
    fout.close()
    fout_wcs.close()
    if N_try==2:
        repeat=True
        if verbosity: print('WARNING: The list of best PSF stars is empty! Using good PSF stars.')
    if VERY_FINAL_STARS!=[]:
        N_try=3
    #print(N_try, repeat)  
  
  
  
  #print(VERY_FINAL_STARS)
  #exit()
  
  #exit()
  #print(repeat)
  if not repeat:
    filtered_data = sigma_clip(M0, sigma=2)
    good_only = filtered_data.data[~filtered_data.mask]
    good_M0 = np.in1d(M0,good_only)
  else:
    good_M0 = len(M0) * [True] 
  
  #print(M0, good_M0)
  #exit()
  for k in range(len(VERY_FINAL_STARS)):
    K = VERY_FINAL_STARS[k]
    if good_M0[k]==True:
      comment_star = "GOOD"
    else:
      comment_star = "BAD"
    f_ind.write('STAR %i: %.5f %.3f %.2f %.1f %s\n' % (K+1,M0[k],FWHMS[k],ELLS[k],PAS[k],comment_star))
    if verbosity: print("Star #%i: %.5f %.3f %.2f %.1f %s" % (K+1,M0[k],FWHMS[k],ELLS[k],PAS[k],comment_star))

  
  f_ind.close()


  
  if m0 is None:
    m0_mean, m0_median, m0_std = sigma_clipped_stats(M0, sigma=2)
    fwhm_mean, fwhm_med, fwhm_std = sigma_clipped_stats(FWHMS, sigma=2)
    ell_mean, ell_med, ell_std = sigma_clipped_stats(ELLS, sigma=2)
    theta_mean, theta_med, theta_std = sigma_clipped_stats(PAS, sigma=2)
    
    if verbosity: print('FINAL NUMBER OF PSF STARS IS: %i' % (len(M0)))
    if verbosity: print("PSF:")
    if verbosity: print('\tFWHM [pix] = %.3f%s +/- %.3f' % (fwhm_med,comment,fwhm_std)) 
    if verbosity: print('\tEllipticity = %.3f%s +/- %.3f' % (ell_med,comment,ell_std)) 
    if verbosity: print('\tP.A. = %.3f%s +/- %.3f' % (theta_med,comment,theta_std)) 


    if verbosity: print("Photometric calibration for NOMAD (DO NOT USE IT!!!):")
    if verbosity: print("\t Zero-point [mag] %.5f%s +/- %.5f" % (m0_median,comment,m0_std))
    
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data

    header = hdulist[0].header

    add_keyw_to_header.write_keyw(header,'ZP', m0_median)
    add_keyw_to_header.write_keyw(header,'PSF_FWHM', fwhm_med)
    add_keyw_to_header.write_keyw(header,'PSF_ELL', ell_med)
    add_keyw_to_header.write_keyw(header,'PSF_PA', theta_med)

    hdulist.close()
    hdu = pyfits.PrimaryHDU(data, header)
    if output_image is None:
        output_image = input_image
    hdu.writeto(output_image, clobber=True)

    #return m0_mean, m0_median, m0_std, fwhm_med, ell_med, theta_med, comment
  
  if os.path.exists('segm.fits'):
    shutil.move('segm.fits', 'segm_ini.fits')
    
  for file in ['x.reg', 'stars.reg', 'peaks.txt', 'models.fits', 'best_stars_wcs.reg', 'ZP_results.dat']:
      if os.path.exists(file):
        os.remove(file)  
  if verbosity: print('Done!')      
  if not np.isnan(fwhm_med):       
        return fwhm_med, ell_med, theta_med
  else:
        return fwhm_mean, ell_mean, theta_mean      



'''
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Do photometric calibration")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("--output_image", help="Output fits image", type=str, default=None)  
    parser.add_argument("--band", help="Photometric band in the catalogue (e.g. B, V, R, J, H, K)", type=str, default='R') 
    parser.add_argument("--catalogue", help="Photometric catalogue of stars (NOMAD, UCAC4, II/246 (2MASS), V/139 (SDSS))", type=str, default='NOMAD')
    parser.add_argument("--sex_cat", help="Sextractor catalogue", type=str, default=None)     
    parser.add_argument("--mag_min", help="Minimum apparent magnitude of stars which should be taken into account", type=float, default=15.6)   
    parser.add_argument("--mag_max", help="Maximum apparent magnitude of stars which should be taken into account", type=float, default=18.3)       
    parser.add_argument("--region_file", help="Region file with selected PSF stars", type=str, default=None) 
    args = parser.parse_args()

    input_image = args.input_image
    output_image = args.output_image
    sex_cat = args.sex_cat
    catalogue = args.catalogue
    band = args.band
    mag_min = args.mag_min
    mag_max = args.mag_max
    region_file = args.region_file

    main(input_image, output_image, band=band, catalogue=catalogue, m0=None, sex_cat=sex_cat, region_file=region_file, peaks_file=None, mag_min=mag_min, mag_max=mag_max)
'''

'''
input_image = 'new-image.fits'
band = 'R'
catalogue = "NOMAD"
main(input_image,None,band,catalogue, sex_cat='field.cat')
'''
