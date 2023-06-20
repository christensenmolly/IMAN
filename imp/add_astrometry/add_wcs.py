#!/usr/bin/python
# DESCRIPTION:
# Script to add rough wcs to an image if a region file with reference stars (more than 3) is given.
# For each star a circle must be given which encompasses the whole star, with the following text label:
# ra,dec , where ra and dec are give in degrees (decimals).
# NOTE: Sometimes it produces wrong resuts. For adding robust wcs you should then feed the output image into nova.astrometry.net.
# MINIMAL USAGE: python add_wcs.py [input_image] [input_region_with_ref_stars]


# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import math
import sys
import shutil
import argparse
import os
import astroquery
from astroquery.simbad import Simbad
from astropy import coordinates
import astropy.units as u
from astropy import wcs
import collections
from astroquery.vizier import Vizier
from photutils import CircularAperture
from astropy.stats import sigma_clipped_stats
from photutils import centroid_2dg
import warnings
from photutils import CircularAnnulus, EllipticalAnnulus
from photutils import CircularAperture, EllipticalAperture
from photutils import aperture_photometry
warnings.filterwarnings("ignore")



import add_simple_wcs


def get_true_ra_dec(RA,DEC,extent,catalog,band,table=None):
        if table is None:
            viz = Vizier(keywords=["stars", "optical"])
            viz.ROW_LIMIT = -1
            try:
                c = coordinates.SkyCoord(float(RA),float(DEC),unit=('deg','deg'),frame='icrs')
                result = viz.query_region(c, radius=(extent/2./3600.)*u.deg, catalog=catalog)
                table = result[0]
            except:
                print('ERROR: The reference coordinate is not found in the catalogue! Exiting ...')
                exit()

        # Find the star nearest to the given RA,DEC
        min_dist = 10.*3600.
        for i in range(len(table)):
            star_ra = float(table["RAJ2000"][i]) ####WARNING: was _RAJ2000
            star_dec = float(table["DEJ2000"][i]) ####WARNING: was DEJ2000
            star_mag = float(table[band + "mag"][i])
            dist = math.sqrt( (RA-star_ra)**2 + (DEC-star_dec)**2) * 3600.
            if dist<min_dist:
              RA_prec = star_ra
              DEC_prec = star_dec
              Mag_prec = star_mag
              ind_prec = i
              min_dist = dist
        print("True coordinates of the reference star %f,%f are: %f,%f (%f arcsec)" % (RA, DEC, RA_prec, DEC_prec, min_dist))
        return table,RA_prec,DEC_prec,Mag_prec,ind_prec

  





def precise_center(data,nx,ny,x,y,R):
  ymin = max([0,int(math.floor(y-R))])
  ymax = min([ny,int(math.floor(y+R))])
  xmin = max([0,int(math.floor(x-R))])
  xmax = min([nx,int(math.floor(x+R))])
  x_cen, y_cen = centroid_2dg(data[ymin:ymax,xmin:xmax])
  x = x_cen+xmin; y = y_cen+ymin
  return x+1.,y+1.



def read_reg_file(region_file,data,nx,ny):
    f = open(region_file, "r")
    f_new = open('tmp.reg', "w")
    lines = f.readlines()
    x = []; y = []; R_in = []; R_out = []; RA = []; DEC = []; mag = []

    for line in lines:
      if 'circle(' in line:
        xc = float(line.split('(')[-1].split(',')[0])
        yc = float(line.split(',')[1])
        rad_in = float(line.split(',')[2].split(')')[0])

        xc,yc = precise_center(data,nx,ny,xc,yc,rad_in)
        rad_out = 1.5*rad_in
        x.append(xc)
        y.append(yc)
        R_in.append(rad_in)
        R_out.append(rad_out)
        if 'text={' in line:
          ra = float(line.split('text={')[-1].split(',')[0])
          try:
            #xx = float(line.split('text={')[-1].split(',')[2].split('}')[0])
            dec = float(line.split('text={')[-1].split(',')[1])
            m = float(line.split('text={')[-1].split(',')[-1].split('}')[0])
          except:
            dec = float(line.split('text={')[-1].split(',')[-1].split('}')[0])
            m = float('nan')
          RA.append(ra)
          DEC.append(dec)
          mag.append(m)
          f_new.write('annulus(%f,%f,%f,%f) # text={%f,%f,%f}\n' % (xc,yc,rad_in,rad_out,ra,dec,m))
        else:
          RA.append(float('nan'))
          DEC.append(float('nan'))
          mag.append(float('nan'))  
          f_new.write('annulus(%f,%f,%f,%f)\n' % (xc,yc,rad_in,rad_out))
      if 'ellipse(' in line or 'annulus(' in line:
        xc = float(line.split('(')[-1].split(',')[0])
        yc = float(line.split(',')[1])
        rad_in = float(line.split(',')[2])
        xc,yc = precise_center(data,nx,ny,xc,yc,rad_in)
        rad_out = 1.5*rad_in
        x.append(xc)
        y.append(yc)
        R_in.append(rad_in)
        R_out.append(rad_out)
        if 'text={' in line:  
          ra = float(line.split('text={')[-1].split(',')[0])

          try:
            #xx = float(line.split('text={')[-1].split(',')[2].split('}')[0])
            dec = float(line.split('text={')[-1].split(',')[1])
            m = float(line.split('text={')[-1].split(',')[-1].split('}')[0])
          except:
            dec = float(line.split('text={')[-1].split(',')[-1].split('}')[0])
            m = float('nan')
          RA.append(ra)
          DEC.append(dec)
          mag.append(m)
          f_new.write('annulus(%f,%f,%f,%f) # text={%f,%f,%f}\n' % (xc,yc,rad_in,rad_out,ra,dec,m))
        else:
          RA.append(float('nan'))
          DEC.append(float('nan'))
          mag.append(float('nan'))  
          f_new.write('annulus(%f,%f,%f,%f)\n' % (xc,yc,rad_in,rad_out))
    f.close()
    f_new.close()
    
    if len(x)<3:
      print('The number of reference stars in your region file is insufficient to create WCS! Exiting ...') 
    else:
      print('Number of stars in the DS9 region files is %i' % (len(x)))
      shutil.move('tmp.reg',region_file)
    return x,y,R_in,R_out,RA,DEC,mag
      



def define_matrix(pts_pix_ref, pts_wcs_ref, pts_pix, pts_wcs):
  # pts_pix_ref - reference pixel coordinates
  # pts_wcs_red - reference wcs coordinates
  I0 = pts_pix_ref[0]
  J0 = pts_pix_ref[1]

  X0 = pts_wcs_ref[0]
  Y0 = pts_wcs_ref[1]

  X = []; Y = []; I = []; J = []
  for k in range(len(pts_pix)):
    X.append(pts_wcs[k][0])
    Y.append(pts_wcs[k][1])  
    I.append(pts_pix[k][0])
    J.append(pts_pix[k][1])
  X = np.array(X, float)
  Y = np.array(Y, float)
  I = np.array(I, float)
  J = np.array(J, float)

  # Let's re-designate the variables:  
  # Solving z = a*x + b*y using lsq
  
  # First equation: X = X0 + c11*(I-I0) + c12*(J-J0)
  z = X - X0
  x = I - I0
  y = J - J0
  
  zy = sum(z*y)
  xz = sum(x*z)
  xy = sum(x*y)
  x2 = sum(x**2)
  y2 = sum(y**2)
  
  
  c12 = ( zy - xz*xy/x2 ) / ( y2 - xy*xy/x2 )	# b
  c11 = ( xz - c12*xy ) / x2	# a
  
  #c11 = a; c12 = b
  

  # Second equation: Y = Y0 + c21*(I-I0) + c22*(J-J0)
  z = Y - Y0
  x = I - I0
  y = J - J0
  
  zy = sum(z*y)
  xz = sum(x*z)
  xy = sum(x*y)
  x2 = sum(x**2)
  y2 = sum(y**2)
  
  
  c22 = ( zy - xz*xy/x2 ) / ( y2 - xy*xy/x2 )	#b
  c21 = ( xz - c22*xy ) / x2	#a
  
  #c21 = a; c22 = b
  return c11,c12,c21,c22




def add_wcs_to_header(x_ref,y_ref,RA_ref,DEC_ref,c11,c12,c21,c22,inputname,outputname):
    hdulist = pyfits.open(inputname)
    referenceheader = pyfits.getheader(inputname, 0)
    outframe = pyfits.getdata(inputname, 0)


    referenceheader['EQUINOX'] = 2.000000000000E+03
    referenceheader['RADECSYS'] = 'FK5'
    
    referenceheader['CTYPE1'] = 'RA---TAN'
    referenceheader['CUNIT1'] = 'deg'    
    referenceheader['CRVAL1'] = RA_ref
    referenceheader['CRPIX1'] = x_ref

    
    referenceheader['CTYPE2'] = 'DEC--TAN'	#### Was:  'DEC---TAN'
    referenceheader['CUNIT2'] = 'deg'  
    referenceheader['CRVAL2'] = DEC_ref
    referenceheader['CRPIX2'] = y_ref

    referenceheader['CD1_1'] = c11
    referenceheader['CD1_2'] = c12
    referenceheader['CD2_1'] = c21
    referenceheader['CD2_2'] = c22  
    
    hdu = pyfits.PrimaryHDU(outframe,referenceheader)
    hdu.writeto(outputname,clobber=True)




def calculate_angle(pts_pix, pts_wcs):
    x1 = pts_pix[0][0]
    y1 = pts_pix[0][1]


    x2 = pts_pix[1][0]
    y2 = pts_pix[1][1]

    RA1 = pts_wcs[0][0]*3600.
    DEC1 = pts_wcs[0][1]*3600.
    print(x1,y1,RA1/3600.,DEC1/3600.)
    RA2 = pts_wcs[1][0]*3600.
    DEC2 = pts_wcs[1][1]*3600.
    print(x2,y2,RA2/3600.,DEC2/3600.)
    #exit()

    
    alpha = np.degrees(math.atan2(DEC2-DEC1, RA2-RA1))
    beta = np.degrees(math.atan2(y2-y1, x2-x1))
    #print(alpha, beta)
    cdelt = math.sqrt( (RA2-RA1)**2 + (DEC2-DEC1)**2)/(math.sqrt( (x2-x1)**2 + (y2-y1)**2))
    PA = alpha-beta
    return PA,cdelt




def main(input_image, scale, band, region_file, output_image=None):
  if output_image is None:
    output_image = input_image.split('.fits')[0]+'_wcs.fits'
    
  # Read region file, change circles on circular annuli if needed:
  hdulist = pyfits.open(input_image)
  data = hdulist[0].data
  ny,nx = np.shape(data)
  x,y,R_in,R_out,RA,DEC,mags = read_reg_file(region_file,data,nx,ny)


  # Download stars from the catalogue:
  # Find true coordinates for the reference star:
  extent = max([ny,nx]) * scale # in arcsec
  close_stars = []
  final_stars = []
  ref_star_index = []
  table = None
  for k in range(len(RA)):
    close_stars.append([])
    final_stars.append([])
    if np.isnan(RA[k])==False and np.isnan(DEC[k])==False:
      table,RA_prec,DEC_prec,Mag_prec,ind_prec = get_true_ra_dec(RA[k],DEC[k],extent,"NOMAD",band,table=table)
      RA[k] = RA_prec
      DEC[k] = DEC_prec
      if np.isnan(mags[k])==True:
        mags[k] = Mag_prec
      ref_star_index.append(k)
      
  # Check if we have more than three stars with defined WCS coordinates:
  count = 0
  pts_pix = []; pts_wcs = []; circular_annuli = []; mag_prec = []
  for k in range(len(x)):
    if np.isnan(RA[k])==False:
      if count == 0:
        x_ref = x[k]; y_ref = y[k]; RA_ref = RA[k]; DEC_ref = DEC[k]
        pts_pix_ref = [x[k],y[k]]
        pts_wcs_ref = [RA[k],DEC[k]]
      else:
        pts_pix.append([x[k],y[k]])
        pts_wcs.append([RA[k],DEC[k]])
      count = count + 1
      
      if np.isnan(mags[k])==False:
          circular_annuli.append([x[k],y[k],R_in[k]])
          mag_prec.append(mags[k])
  
  #print(pts_pix, pts_wcs)
  #PA,cdelt = calculate_angle(pts_pix, pts_wcs)
  #print(PA,cdelt)
  
  #add_simple_wcs.main(input_image, x_ref, y_ref, RA_ref, DEC_ref, scale/3600., scale/3600., angle=PA, output_image=None)
  

  if count >=3:
    c11,c12,c21,c22 = define_matrix(pts_pix_ref, pts_wcs_ref, pts_pix, pts_wcs)
    add_wcs_to_header(x_ref,y_ref,RA_ref,DEC_ref,c11,c12,c21,c22,input_image,output_image)
    print('Done!')
  else:
      print('Error. In the region file less than 3 stars can be used to create wcs. Exiting!')

  
  
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add world coordinate system to image")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("region_file", help="Input region file with more than 3 reference stars")  
    parser.add_argument("--output_image", help="Output fits image", type=str, default=None)  
    parser.add_argument("--scale", help="Rough estimate on the pixel scale in arcsec/pix", type=float, default=1.0)  
    parser.add_argument("--band", help="Photometric band in the NOMAD catalogue (B, V, R, J, H, K)", type=str, default='R') 

    args = parser.parse_args()

    input_image = args.input_image
    region_file = args.region_file
    output_image = args.output_image
    scale = args.scale
    band = args.band


    main(input_image, scale, band, region_file, output_image = output_image)


