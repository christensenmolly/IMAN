#!/usr/bin/python
# -*- coding:  cp1251 -*-
# Import the necessary modules
import pyfits
import numpy as np
import math
import itertools
import matplotlib.pyplot as plt
from scipy import ndimage
import sys
from matplotlib.path import Path
from itertools import product
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import imp_setup
import os
import subprocess

from astroquery.simbad import Simbad
from astropy import coordinates
import astropy.units as u
from astropy import wcs
import collections
from astroquery.vizier import Vizier
from photutils import aperture_photometry
from photutils import CircularAperture
from astropy.stats import sigma_clipped_stats
from photutils import centroid_2dg
import warnings
warnings.filterwarnings("ignore")
PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]








def backgr_around_star(scidata,xc,yc,R):
  nx, ny =scidata.shape[1], scidata.shape[0]
  Rin = int(ceil(1.1*R))
  Rout = int(ceil(2.*R))
  backgr = []
  X = []; Y = []; I = []
  for y in range(int(yc)-Rout,int(yc)+Rout,1):
    for x in range(int(xc)-Rout,int(xc)+Rout,1):
      try:
	I.append(scidata[y,x])
	X.append(x)
	Y.append(y)
	if (x-int(xc))**2+(y-int(yc))**2>Rin**2 and (x-int(xc))**2+(y-int(yc))**2<=Rout**2 and x>=0 and y>=0 and x<nx and y<ny:
	  backgr.append(scidata[y,x])
      except:
	zz = 1
  Imax = max(I)-np.median(backgr)

  g_init = models.Gaussian2D(amplitude=Imax, x_mean=xc, y_mean=yc, x_stddev=R/2., y_stddev=R/2., theta=0.)
  f = fitting.LevMarLSQFitter()
  mean, median, std = sigma_clipped_stats(backgr, sigma=3.0, iters=5)
  p = f(g_init, X, Y, I-median)
  
  
  return median,std,float(p.x_mean[0]),float(p.y_mean[0]),float(p.x_stddev[0]),float(p.y_stddev[0])


























def get_true_ra_dec(RA,DEC,extent,catalog,band,table=None):
	#print type(table)
	if table is None:
	    viz = Vizier(keywords=["stars", "optical"])
	    viz.ROW_LIMIT = -1
	    try:
		c = coordinates.SkyCoord(float(RA),float(DEC),unit=('deg','deg'),frame='icrs')
		
		result = viz.query_region(c, radius=(extent/2./3600.)*u.deg, catalog=catalog)
		table = result[0]
	    except:
		print 'ERROR: The reference coordinate is not found in the catalogue! Exiting ...'
		exit()

	# Find the star nearest to the given RA,DEC
	min_dist = 10.*3600.
	#print(table)
	#exit()
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
        print "True coordinates of the reference star are:", RA_prec,DEC_prec
        return table,RA_prec,DEC_prec,Mag_prec,ind_prec
        













'''
        number_of_stars = 0
        number_of_stars_in_frame = 0
        number_of_new_stars = 0

	MAGS = []
	MAGS_ERR = []
	dist_to_star = []
	
        # Loop over all entries in the table
        for i in range(len(table)):
	    magnitudes = {}
	    magnitude_errors = {}

	    # Get the magnitude in different bands
	    for name in table.colnames:

		# If this column name does not end with "mag", skip it
		if not name.endswith("mag"): continue

		# If the column name contains more than one character before "mag", skip it
		if len(name.split("mag")[0]) > 1: continue

		# Get the name of the band
		band = name.split("mag")[0]

		# Create empty lists for the magnitudes and errors
		magnitudes[band] = []
		magnitude_errors[band] = []

            # -- General information --

            # Get the ID of this star in the catalog
            if catalog == "UCAC4": star_id = table["UCAC4"][i]
            elif catalog == "NOMAD": star_id = table["NOMAD1"][i]
            elif catalog == "II/246": star_id = table["_2MASS"][i]
            elif catalog == "V/139": star_id = table["SDSS9"][i]
            else: raise ValueError("Catalogs other than 'UCAC4', 'NOMAD' or 'II/246' are currently not supported")

            # -- Positional information --

            # Get the position of the star as a SkyCoord object and as pixel coordinate
            position = coordinates.SkyCoord(ra=table["RAJ2000"][i], dec=table["DEJ2000"][i], unit=(u.deg, u.deg), frame='fk5')
            #pixel_position_x, pixel_position_y = position.to_pixel(frame.wcs, origin=0, mode='wcs')
            #pixel_position = Position(pixel_position_x, pixel_position_y)

            # Get the right ascension and declination for the current star
            star_ra = table["RAJ2000"][i]
            star_dec = table["DEJ2000"][i]
	    #dist_to_star.append(sqrt((float(star_ra)-float(coord[0]))**2+(float(star_dec)-float(coord[1]))**2) / 3600.) # find distance to the star from the given one
	    Dist = sqrt((float(star_ra)-float(coord[0]))**2+(float(star_dec)-float(coord[1]))**2) / 3600.
	    

            number_of_stars += 1

            # If this star does not lie within the frame, skip it
            #if not frame.contains(position): continue

            number_of_stars_in_frame += 1

            # Get the mean error on the right ascension and declination
            if catalog == "UCAC4" or catalog == "NOMAD":

                ra_error = table["eRAJ2000"][i] * u.mas
                dec_error = table["eDEJ2000"][i] * u.mas

            elif catalog == "II/246":

                error_maj = table["errMaj"][i] * u.arcsec
                error_min = table["errMin"][i] * u.arcsec
                error_theta = Angle(table["errPA"][i], u.deg)

                # Temporary: use only the major axis error (convert the error ellipse into a circle)
                ra_error = error_maj.to("mas")
                dec_error = error_maj.to("mas")

            elif catalog == "V/139":

                error_maj = 1#table["errMaj"][i] * u.arcsec
                error_min = 1#table["errMin"][i] * u.arcsec
                error_theta = 1#Angle(table["errPA"][i], u.deg)

                # Temporary: use only the major axis error (convert the error ellipse into a circle)
                ra_error = 1#error_maj.to("mas")
                dec_error = 1#error_maj.to("mas")

            else: raise ValueError("Catalogs other than 'UCAC4', 'NOMAD' or 'II/246' are currently not supported")

            # -- Magnitudes --

            # Loop over the different bands for which a magnitude is defined
            for band in magnitudes:

                # Determine the column name
                column_name = band + "mag"

                value = table[column_name][i]
		#print band,value

                if isinstance(value, np.ma.core.MaskedConstant):

                    magnitudes[band].append(None)
                    magnitude_errors[band].append(None)

                else:

                    # Add the magnitude value
                    magnitudes[band].append(value)

                    # Check for presence of error on magnitude
                    error_column_name = "e_" + column_name
                    if error_column_name in table.colnames:
                        error = table[error_column_name][i]
                        if isinstance(error, np.ma.core.MaskedConstant): magnitude_errors[band].append(None)
                        else: magnitude_errors[band].append(error)
                    else: magnitude_errors[band].append(None)
	    MAGS.append(magnitudes)
	    MAGS_ERR.append(magnitude_errors)
	    
	dist_to_star = np.array(dist_to_star)
	index_closest_star = int(np.where( dist_to_star == np.min(dist_to_star))[0])
	mags.append(MAGS[index_closest_star])
	mags_err.append(MAGS_ERR[index_closest_star])

    return mags,mags_err
'''  
  
  
  
  
  
  
  
  
  
  

'''
def get_stars(coords,Radius,catalog):
    mags = []
    mags_err = []
    for coord in coords:
	#catalog="NOMAD"
	#catalog="V/139"
	
	viz = Vizier(keywords=["stars", "optical"])
	viz.ROW_LIMIT = -1
	c = coordinates.SkyCoord(float(coord[0]),float(coord[1]),unit=('deg','deg'),frame='icrs')


	try:
	    result = viz.query_region(c, radius=(Radius/3600.)*u.deg, catalog=catalog)			#### 3!!!
	    table = result[0]
	except:
	    mags.append(float('nan'))
	    mags_err.append(float('nan'))
	    continue



        number_of_stars = 0
        number_of_stars_in_frame = 0
        number_of_new_stars = 0










	MAGS = []
	MAGS_ERR = []
	dist_to_star = []
	
        # Loop over all entries in the table
        for i in range(len(table)):
	    magnitudes = {}
	    magnitude_errors = {}

	    # Get the magnitude in different bands
	    for name in table.colnames:

		# If this column name does not end with "mag", skip it
		if not name.endswith("mag"): continue

		# If the column name contains more than one character before "mag", skip it
		if len(name.split("mag")[0]) > 1: continue

		# Get the name of the band
		band = name.split("mag")[0]

		# Create empty lists for the magnitudes and errors
		magnitudes[band] = []
		magnitude_errors[band] = []

            # -- General information --

            # Get the ID of this star in the catalog
            if catalog == "UCAC4": star_id = table["UCAC4"][i]
            elif catalog == "NOMAD": star_id = table["NOMAD1"][i]
            elif catalog == "II/246": star_id = table["_2MASS"][i]
            elif catalog == "V/139": star_id = table["SDSS9"][i]
            else: raise ValueError("Catalogs other than 'UCAC4', 'NOMAD' or 'II/246' are currently not supported")

            # -- Positional information --

            # Get the position of the star as a SkyCoord object and as pixel coordinate
            position = coordinates.SkyCoord(ra=table["RAJ2000"][i], dec=table["DEJ2000"][i], unit=(u.deg, u.deg), frame='fk5')
            #pixel_position_x, pixel_position_y = position.to_pixel(frame.wcs, origin=0, mode='wcs')
            #pixel_position = Position(pixel_position_x, pixel_position_y)

            # Get the right ascension and declination for the current star
            star_ra = table["RAJ2000"][i]
            star_dec = table["DEJ2000"][i]
	    dist_to_star.append(sqrt((float(star_ra)-float(coord[0]))**2+(float(star_dec)-float(coord[1]))**2) / 3600.) # find distance to the star from the given one

            number_of_stars += 1

            # If this star does not lie within the frame, skip it
            #if not frame.contains(position): continue

            number_of_stars_in_frame += 1

            # Get the mean error on the right ascension and declination
            if catalog == "UCAC4" or catalog == "NOMAD":

                ra_error = table["eRAJ2000"][i] * u.mas
                dec_error = table["eDEJ2000"][i] * u.mas

            elif catalog == "II/246":

                error_maj = table["errMaj"][i] * u.arcsec
                error_min = table["errMin"][i] * u.arcsec
                error_theta = Angle(table["errPA"][i], u.deg)

                # Temporary: use only the major axis error (convert the error ellipse into a circle)
                ra_error = error_maj.to("mas")
                dec_error = error_maj.to("mas")

            elif catalog == "V/139":

                error_maj = 1#table["errMaj"][i] * u.arcsec
                error_min = 1#table["errMin"][i] * u.arcsec
                error_theta = 1#Angle(table["errPA"][i], u.deg)

                # Temporary: use only the major axis error (convert the error ellipse into a circle)
                ra_error = 1#error_maj.to("mas")
                dec_error = 1#error_maj.to("mas")

            else: raise ValueError("Catalogs other than 'UCAC4', 'NOMAD' or 'II/246' are currently not supported")

            # -- Magnitudes --

            # Loop over the different bands for which a magnitude is defined
            for band in magnitudes:

                # Determine the column name
                column_name = band + "mag"

                value = table[column_name][i]
		#print band,value

                if isinstance(value, np.ma.core.MaskedConstant):

                    magnitudes[band].append(None)
                    magnitude_errors[band].append(None)

                else:

                    # Add the magnitude value
                    magnitudes[band].append(value)

                    # Check for presence of error on magnitude
                    error_column_name = "e_" + column_name
                    if error_column_name in table.colnames:
                        error = table[error_column_name][i]
                        if isinstance(error, np.ma.core.MaskedConstant): magnitude_errors[band].append(None)
                        else: magnitude_errors[band].append(error)
                    else: magnitude_errors[band].append(None)
	    MAGS.append(magnitudes)
	    MAGS_ERR.append(magnitude_errors)
	    
	dist_to_star = np.array(dist_to_star)
	index_closest_star = int(np.where( dist_to_star == np.min(dist_to_star))[0])
	mags.append(MAGS[index_closest_star])
	mags_err.append(MAGS_ERR[index_closest_star])

    return mags,mags_err
'''

def calc_from_cd(cd1_1, cd1_2, cd2_1, cd2_2):

        # TODO: Check if first coordinate in CTYPE is latitude
        # if (ctype EQ 'DEC-') or (strmid(ctype, 1) EQ 'LAT')  then $
        #    cd = reverse(cd,1)

        det = cd1_1*cd2_2 - cd1_2*cd2_1
        if det < 0:
            sgn = -1
        else:
            sgn = 1
        ## if det > 0:
        ##     raise ValueError("Astrometry is for a right-handed coordinate system")

        if (cd2_1 == 0.0) or (cd1_2 == 0.0):
            # Unrotated coordinates?
            xrot = 0.0
            yrot = 0.0
            cdelt1 = cd1_1
            cdelt2 = cd2_2
        else:
            xrot = math.atan2(sgn * cd1_2, sgn * cd1_1)
            yrot = math.atan2(-cd2_1, cd2_2)

            cdelt1 = sgn * math.sqrt(cd1_1**2 + cd1_2**2)
            cdelt2 = math.sqrt(cd1_1**2 + cd2_1**2)

        return xrot, yrot, cdelt1, cdelt2

def convert_pix_to_wcs(I,J,I0,J0,X0,Y0,cdelt1,cdelt2,crota2=0.0):
  if crota2!=0.0:
    crota2 = crota2 * math.pi / 180	# in radians

  sinrot = math.sin(crota2)
  cosrot = math.cos(crota2)
  pcmatrix = [cosrot, -sinrot, sinrot, cosrot]
  cd1_1 = cdelt1*pcmatrix[0]
  cd1_2 = cdelt2*pcmatrix[1]
  cd2_1 = cdelt1*pcmatrix[2]
  cd2_2 = cdelt2*pcmatrix[3]
  
  X = X0 + cd1_1* (I - I0) + cd1_2* (J - J0)
  Y = Y0 + cd2_1* (I - I0) + cd2_2* (J - J0)
  return X,Y

def define_wcs(input_image, pix_coords,cel_coords):
  from astropy.io import fits as pyfits
  shutil.copy(input_image, 'galaxy.fits')
  hdulist = pyfits.open('galaxy.fits')
  data = hdulist[0].data
  
  #http://www.atnf.csiro.au/people/mcalabre/WCS/Intro/WCS06.html
  #print cel_coords
  I0 = pix_coords[0][0]
  J0 = pix_coords[0][1]

  I1 = pix_coords[1][0]
  J1 = pix_coords[1][1]

  I2 = pix_coords[2][0]
  J2 = pix_coords[2][1]

  X0 = cel_coords[0][0]
  Y0 = cel_coords[0][1]
  
  X1 = cel_coords[1][0]
  Y1 = cel_coords[1][1]

  X2 = cel_coords[2][0]
  Y2 = cel_coords[2][1]
  
  #(I1-I0) * C11 + (J1-J0) * C12 = X1-X0
  #(I2-I0) * C11 + (J2-J0) * C12 = X2-X0

  a = np.array([[(I1-I0),(J1-J0)], [(I2-I0),(J2-J0)]])
  b = np.array([X1-X0,X2-X0])
  
  [C11,C12] = np.linalg.solve(a, b)

  #(I1-I0) * C21 + (J1-J0) * C22 = Y1-Y0
  #(I2-I0) * C21 + (J2-J0) * C22 = Y2-Y0
  
  #a = np.array([[(I1-I0),(J1-J0)], [(I2-I0),(J2-J0)]])
  b = np.array([Y1-Y0,Y2-Y0])
  
  [C21,C22] = np.linalg.solve(a, b)


  xrot, yrot, cdelt1, cdelt2 = calc_from_cd(C11, C12, C21, C22)
  
  from astropy.io import fits as pyfits
  from astropy import wcs

  w = wcs.WCS(naxis=2)

  # what is the center pixel of the XY grid.
  w.wcs.crpix = [I0,J0]

  # what is the galactic coordinate of that pixel.
  w.wcs.crval = [X0, Y0]

  # what is the pixel scale in lon, lat.
  w.wcs.cdelt = np.array([cdelt1, cdelt2])

  # you would have to determine if this is in fact a tangential projection. 
  w.wcs.ctype = ["RA-TAN", "DEC-TAN"]

  pc1_1 = C11 / cdelt1
  pc1_2 = C12 / cdelt1
  pc2_1 = C21 / cdelt2
  pc2_2 = C22 / cdelt2



  # write the HDU object WITH THE HEADER
  header = w.to_header()

  header['PC1_1'] = pc1_1
  header['PC1_2'] = pc1_2
  header['PC2_1'] = pc2_1
  header['PC2_2'] = pc2_2
  
  #print header
  #exit()
  #header['CDELT1'] 
  #header['CDELT2']

  hdu = pyfits.PrimaryHDU(data, header=header)
  hdu.writeto('galaxy.fits', overwrite=True)



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
	#print nx,ny,xc,yc,rad_in
	#exit()
	xc,yc = precise_center(data,nx,ny,xc,yc,rad_in)
	rad_out = 1.5*rad_in
	x.append(xc)
	y.append(yc)
	R_in.append(rad_in)
	R_out.append(rad_out)
	if 'text={' in line:
	  ra = float(line.split('text={')[-1].split(',')[0])
	  try:
	    xx = float(line.split('text={')[-1].split(',')[2].split('}')[0])
	    dec = float(line.split('text={')[-1].split(',')[1])
	    m = float(line.split('text={')[-1].split(',')[-1].split('}')[0])
	  except:
	    dec = float(line.split('text={')[-1].split(',')[-1].split('}')[0])
	    m = float('nan')
	  RA.append(ra)
	  DEC.append(dec)
	  mag.append(m)
	  print >>f_new, 'annulus(%f,%f,%f,%f) # text={%f,%f,%f}' % (xc,yc,rad_in,rad_out,ra,dec,m)
	else:
	  RA.append(float('nan'))
	  DEC.append(float('nan'))
	  mag.append(float('nan'))  
	  print >>f_new, 'annulus(%f,%f,%f,%f)' % (xc,yc,rad_in,rad_out)
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
	    xx = float(line.split('text={')[-1].split(',')[2].split('}')[0])
	    dec = float(line.split('text={')[-1].split(',')[1])
	    m = float(line.split('text={')[-1].split(',')[-1].split('}')[0])
	  except:
	    dec = float(line.split('text={')[-1].split(',')[-1].split('}')[0])
	    m = float('nan')
	  RA.append(ra)
	  DEC.append(dec)
	  mag.append(m)
	  print >>f_new, 'annulus(%f,%f,%f,%f) # text={%f,%f,%f}' % (xc,yc,rad_in,rad_out,ra,dec,m)
	else:
	  RA.append(float('nan'))
	  DEC.append(float('nan'))
	  mag.append(float('nan'))  
	  print >>f_new, 'annulus(%f,%f,%f,%f)' % (xc,yc,rad_in,rad_out)
    f.close()
    f_new.close()
    
    if len(x)<3:
      print 'The number of reference stars in your region file is insufficient to create WCS! Exiting ...' 
      #exit()
    else:
      print 'Number of stars in the DS9 region files is %i' % (len(x))
      shutil.move('tmp.reg',region_file)
    return x,y,R_in,R_out,RA,DEC,mag

def find_mags(fluxes, mags):
  Mags = []; m0 = []
  for k in range(len(fluxes)):
    if np.isnan(mags[k])==False and np.isnan(fluxes[k])==False:
      m0.append(mags[k] + 2.5*math.log10(fluxes[k]))
  M0 = np.mean(m0)
  M0_std = np.std(m0)

  for k in range(len(fluxes)):
    if np.isnan(fluxes[k])==False:
      Mags.append(M0 - 2.5*math.log10(fluxes[k]))
    else:
      Mags.append(float('nan'))
  return M0,M0_std,Mags       



'''
def find_stars():
  delta_R = 5. # arcsec 
  # Consider each star individually
  for k in range(len(x)):
   if k!= ref_star_index:
    R = sqrt((x[k]-x_ref)**2 + (y[k]-y_ref)**2) * scale # distance between the main reference star and the given star, in arcsec    
    
    # Select stars from the catalogue within the distances of R +/- deltaR:    
    for i in range(len(RA_cat)):
      Dist = sqrt((RA_cat[k]-RA_ref)**2 + (DEC_cat[k]-DEC_ref)**2)
      if Dist-delta_R>=R and Dist+delta_R<=R:
	print RA_cat[k],DEC_cat[k]
'''




def define_matrix(pts_pix_ref, pts_wcs_ref, pts_pix, pts_wcs):
  # pts_pix_ref - reference pixel coordinates
  # pts_wcs_red - reference wcs coordinates
  I0 = pts_pix_ref[0]
  J0 = pts_pix_ref[1]
  #print I0, J0
  X0 = pts_wcs_ref[0]
  Y0 = pts_wcs_ref[1]
  #print X0,Y0
  #exit()
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
  #print X,Y
  #print I,J
  #exit()
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

    #w = wcs.WCS(referenceheader)
    #print referenceheader
    #exit()



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
    hdu.writeto(outputname,overwrite=True)




def get_angle(p0, p1=np.array([0,0]), p2=None):
    # http://stackoverflow.com/questions/13226038/calculating-angle-between-two-lines-in-python
    ''' compute angle (in degrees) for p0p1p2 corner
    Inputs:
        p0,p1,p2 - points in the form of [x,y]
    '''
    if p2 is None:
        p2 = p1 + np.array([1, 0])
    v0 = np.array(p0) - np.array(p1)
    v1 = np.array(p2) - np.array(p1)

    angle = np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))
    return np.degrees(angle)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_wcs_angle(pts_wcs_ref, wcs_stars, pts_wcs_wrong):
  # pts_wcs_wrong - arbitrary converted wcs coordinates of the selected points
  # wcs_stars - selected stars from the catalogue which correspond to the given points
  #print pts_wcs_ref
  #print wcs_stars
  #print pts_wcs_wrong
  #exit()
  Angles = []
  #for k in range(len(wcs_stars)):
  #for i in range(len(pts_wcs_wrong[k])):
  for k in range(len(pts_wcs_wrong)):
    for i in range(len(wcs_stars[k])):
      Angles.append(get_angle(np.array(wcs_stars[k][i]), p1=np.array(pts_wcs_ref), p2=np.array(pts_wcs_wrong[k])))
  #print Angles
  #exit()
  
  
  n, bins,patches = plt.hist(Angles, bins=360., normed=True, fc='k', alpha=0.3)
  #plt.show()
  elem = np.argmax(n)
  PA =  bins[elem]
  print Angles
  PA_prec = find_nearest(np.array(Angles,float),bins[elem])
  print PA_prec
  #exit()
  return PA_prec

def add_wcs_rot(inputname,outputname,x_ref,y_ref,RA_ref,DEC_ref,cdelt1,cdelt2,angle):
    hdulist = pyfits.open(inputname)
    header = pyfits.getheader(inputname, 0)
    outframe = pyfits.getdata(inputname, 0)


    header['EQUINOX'] = 2.000000000000E+03
    header['RADECSYS'] = 'FK5'
    
    header['CTYPE1'] = 'RA---TAN'
    header['CUNIT1'] = 'deg'    
    header['CRVAL1'] = RA_ref
    header['CRPIX1'] = x_ref

    
    header['CTYPE2'] = 'DEC---TAN'
    header['CUNIT2'] = 'deg'  
    header['CRVAL2'] = DEC_ref
    header['CRPIX2'] = y_ref

    header['CDELT1'] = cdelt1
    header['CDELT2'] = cdelt2
    header['CROTA2'] = angle    


    hdu = pyfits.PrimaryHDU(data=outframe, header=header)
    hdu.writeto(outputname,overwrite=True)   


def main(input_image, scale, band, region_file=None, output_image=None, check_reg=True):
  # Create region file and edit it (if needed):

  if region_file==None:
      region_file = 'calibr_stars.reg'
      open(region_file, "a").close()
  
  if check_reg==True:
    print 'Please select at least three good isolated PSF (unsaturated, with high S/N) stars using the circle regions in DS9.\
    For one star, specify RA, DEC and mag (optional) in the text region field separated with comma, e.g. 1.43242,-65.2313' 
    p = subprocess.Popen(["ds9",input_image,"-scale","histequ","-invert","-regions","load",region_file])
    p.wait()  
  
  # Read region file, change circles on circular annuli if needed:
  hdulist = pyfits.open(input_image)
  data = hdulist[0].data
  ny,nx = np.shape(data)
  x,y,R_in,R_out,RA,DEC,mags = read_reg_file(region_file,data,nx,ny)
 
  '''
  # Check if we have more than three stars with defined WCS coordinates:
  count = 0
  pts_pix = []; pts_wcs = []
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

  if count >=3:
    #pts_pix_ref = [x_ref,y_ref]
    #pts_wcs_ref = [RA_ref,DEC_ref]
    c11,c12,c21,c22 = define_matrix(pts_pix_ref, pts_wcs_ref, pts_pix, pts_wcs)
    print c11,c12,c21,c22
    add_wcs_to_header(x_ref,y_ref,RA_ref,DEC_ref,c11,c12,c21,c22,input_image,input_image.split('.fits')[0]+'_wcs.fits')
    #exit()
  '''
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #print RA
  #exit()

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
  pts_pix = []; pts_wcs = []
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

  if count >=3:
    #pts_pix_ref = [x_ref,y_ref]
    #pts_wcs_ref = [RA_ref,DEC_ref]
    print pts_pix_ref, pts_wcs_ref, pts_pix, pts_wcs
    #exit()
    c11,c12,c21,c22 = define_matrix(pts_pix_ref, pts_wcs_ref, pts_pix, pts_wcs)
    print c11,c12,c21,c22
    add_wcs_to_header(x_ref,y_ref,RA_ref,DEC_ref,c11,c12,c21,c22,input_image,input_image.split('.fits')[0]+'_wcs.fits')
    print 'Done!'
    exit()
  
  # Find (precise) scale
  
  if count >=2:
    Scales = []
    for k in range(len(ref_star_index)-1):
      for i in range(k+1,len(ref_star_index)):
	K = ref_star_index[k]
	I = ref_star_index[i]
	Scales.append( math.sqrt((RA[K]-RA[I])**2+(DEC[K]-DEC[I])**2)/math.sqrt((x[K]-x[I])**2+(y[K]-y[I])**2) )
    scale = np.mean(np.array(Scales))*3600.
  #exit()
  #print RA_prec,DEC_prec,Mag_prec
  #print mags
  #exit()
  
  
  from photutils import CircularAnnulus
  from photutils import CircularAperture
  Fluxes = []
  for k in range(len(x)):
    aperture = CircularAperture([(x[k],y[k])], r=R_in[k])
    annulus_aperture = CircularAnnulus([(x[k],y[k])], r_in=R_in[k], r_out=R_out[k])
    apers = [aperture, annulus_aperture]
    phot_table = aperture_photometry(data, apers)
    bkg_mean = phot_table['aperture_sum_1'] / annulus_aperture.area()
    bkg_sum = bkg_mean * aperture.area()
    final_sum = phot_table['aperture_sum_0'] - bkg_sum
    #phot_table['residual_aperture_sum'] = final_sum
    #print final_sum[0]
    Fluxes.append(float(final_sum[0]))
    #print(phot_table['residual_aperture_sum'])      

  M0,M0_std,Mags = find_mags(Fluxes, mags)
  #print M0,M0_std
  #exit()
  #print Mags[ref_star_index]
  #exit()
  RA_cat = []; DEC_cat = []; MAG_cat = []
  for i in range(len(table)):
    RA_cat.append(float(table["RAJ2000"][i]))
    DEC_cat.append(float(table["DEJ2000"][i]))
    MAG_cat.append(float(table[band + "mag"][i]))

  delta_R = 5.#np.min(R_in)*scale#5. # search withn +/- arcsec
  #print delta_R
  #exit()
  delta_mag = 0.1
  
  # Consider each star individually
  x_ref = x[ref_star_index[0]]; y_ref = y[ref_star_index[0]]
  RA_ref = RA[ref_star_index[0]]; DEC_ref = DEC[ref_star_index[0]]
  
  RA_unique = []; DEC_unique = []; MAGS_unique = []
  for k in range(len(x)):
   if k not in [ref_star_index]:
    R = math.sqrt((x[k]-x_ref)**2 + (y[k]-y_ref)**2) * scale # distance between the main reference star and the given star, in arcsec    
    # Select stars from the catalogue within the distances of R +/- deltaR:    

    RA_found = []; DEC_found = []; MAGS_found = []    
    if np.isnan(Mags[k]) == True:
      for i in range(len(RA_cat)):
	Dist = math.sqrt((RA_cat[i]-RA_ref)**2 + (DEC_cat[i]-DEC_ref)**2) * 3600.

	if Dist>=R-delta_R and Dist<=R+delta_R:
	  RA_found.append(RA_cat[i])
	  DEC_found.append(DEC_cat[i])
	  MAGS_found.append(MAG_cat[i])
    else:
      #print k,'here'
      for i in range(len(RA_cat)):
	Dist = math.sqrt((RA_cat[i]-RA_ref)**2 + (DEC_cat[i]-DEC_ref)**2) * 3600.

	if Dist>=R-delta_R and Dist<=R+delta_R:# and MAG_cat[i]>=Mags[k]-delta_mag and MAG_cat[i]<=Mags[k]+delta_mag:
	  RA_found.append(RA_cat[i])
	  DEC_found.append(DEC_cat[i])
	  MAGS_found.append(MAG_cat[i])
    #print k,RA_found
    #print k,len(RA_found)
    if len(RA_found) != 0:
      RA_unique.append(RA_found)
      DEC_unique.append(DEC_found)
      MAGS_unique.append(MAGS_found)
    else:
      RA_unique.append([float('nan')])
      DEC_unique.append([float('nan')])
      MAGS_unique.append([float('nan')])
   else:
      RA_unique.append([k])
      DEC_unique.append([k])
      MAGS_unique.append([Mags[k]])
      
  pts_wcs_wrong = []
  wcs_stars = []
  for k in range(len(x)):
    if k not in ref_star_index and np.isnan(RA_unique[k][0])==False:
      RA_wrong,DEC_wrong = convert_pix_to_wcs(x[k],y[k],x_ref,y_ref,RA_ref,DEC_ref,scale/3600.,scale/3600.,crota2=0.0)
      #print RA_wrong,DEC_wrong
      pts_wcs_wrong.append([RA_wrong,DEC_wrong])
      points = []
      for i in range(len(RA_unique[k])):
	points.append([RA_unique[k][i],DEC_unique[k][i]])
      wcs_stars.append(points)
  #exit()
  #PA_prec = find_wcs_angle([RA_ref,DEC_ref], wcs_stars, pts_wcs_wrong)	### Commented!!!
  PA_prec = 180.
  #print 'here'
  add_wcs_rot(input_image,input_image.split('.fits')[0]+'_wcs.fits',x_ref,y_ref,RA_ref,DEC_ref,-scale/3600.,scale/3600.,PA_prec)
  '''
  print x[3]
  for k in range(len(RA_unique[3])):
      #if RA_unique[3][k]>192.98 and RA_unique[3][k]<192.99:
      print RA_unique[3][k],DEC_unique[3][k],MAGS_unique[3][k]#,Mags[3]
      print Mags
  '''






  exit()
  pts_pix = []; pts_wcs = []
  for k in range(len(x)):
    if len(RA_unique[k])==1 and k!=ref_star_index:
      if np.isnan(RA_unique[k])==False:
	 pts_pix.append([x[k],y[k]])
	 pts_wcs.append([RA_unique[k][0],DEC_unique[k][0]])
	 print x[k],y[k], RA_unique[k][0],DEC_unique[k][0]
  
  exit()
  
  #print pts_pix
  #print pts_wcs
  #exit()
  pts_pix_ref = [x_ref,y_ref]
  pts_wcs_ref = [RA_ref,DEC_ref]
  c11,c12,c21,c22 = define_matrix(pts_pix_ref, pts_wcs_ref, pts_pix, pts_wcs)
  print c11,c12,c21,c22
  add_wcs_to_header(x_ref,y_ref,RA_ref,DEC_ref,c11,c12,c21,c22,input_image,input_image.split('.fits')[0]+'_wcs.fits')

  '''
  print x
  print y
  from photutils import centroid_2dg
  for k in range(len(x)):
    ymin = max([0,int(math.floor(y[k]-R[k]))])
    ymax = min([ny,int(math.floor(y[k]+R[k]))])
    xmin = max([0,int(math.floor(x[k]-R[k]))])
    xmax = min([nx,int(math.floor(x[k]+R[k]))])
    x_cen, y_cen = centroid_2dg(data[ymin:ymax,xmin:xmax])
    x[k] = x_cen+xmin; y[k] = y_cen+ymin
  print x
  print y
  '''  
  # and determine their fluxes and magnitudes:
  
  
  # xx

  exit()
  # Estimate the size of the field in arcsec:

  extent = max([ny,nx]) * scale # in arcsec

  # Download stars from the catalogue:
  # Find true coordinates for the reference star:
  close_stars = []
  final_stars = []
  for k in range(len(RA)):
    close_stars.append([])
    final_stars.append([])
    if np.isnan(RA[k])==False and np.isnan(DEC[k])==False:
      table,RA_prec,DEC_prec,ind_prec = get_true_ra_dec(RA[k],DEC[k],extent,"NOMAD")
      ref_star_index = k
  
  
  
  for k in range(len(RA)):
    if k!=ref_star_index:
      dist_between_ref_and_other = math.sqrt( (x[k]-x[ref_star_index])**2 + (y[k]-y[ref_star_index])**2 ) * scale # in arcsec

      for i in range(len(table)):
            star_ra = float(table["RAJ2000"][i])
            star_dec = float(table["DEJ2000"][i])
            dist = math.sqrt( (RA_prec-star_ra)**2 + (DEC_prec-star_dec)**2) * 3600.
            if abs(dist-dist_between_ref_and_other)<3.:
	      close_stars[k].append(i)


  '''
  final_stars = []
  for k in range(len(RA)):
    if k!=ref_star_index:
      new = []
      for i in range(k+1,len(RA)):
	if i!=ref_star_index:
	  new.append('no')
      final_stars.append(new)
  #print final_stars
  #exit()
  '''
  #final_stars = [[[]]*len(RA) for i in range(len(RA))]
  
  Stars = collections.OrderedDict()
  for k in range(len(RA)):
    if k!=ref_star_index:
      for i in range(k+1,len(RA)):
	if i!=ref_star_index:
	  true_dist = math.sqrt( (x[k]-x[i])**2 + (y[k]-y[i])**2) * scale # in arcsec
	  #print k,i,x[k],x[i],true_dist
	  stars_k = close_stars[k]
	  stars_i = close_stars[i]
	  #print stars_k, stars_i
	  #exit()
	  s = []
	  for kk in range(len(stars_k)):
	    for ii in range(len(stars_i)):
	      star_ra_kk = float(table["RAJ2000"][stars_k[kk]])
	      star_dec_kk = float(table["DEJ2000"][stars_k[kk]])
	      star_ra_ii = float(table["RAJ2000"][stars_i[ii]])
	      star_dec_ii = float(table["DEJ2000"][stars_i[ii]])
	      #print true_dist, math.sqrt( (star_ra_ii-star_ra_kk)**2 + (star_dec_ii-star_dec_kk)**2) * 3600.
	      if abs(math.sqrt( (star_ra_ii-star_ra_kk)**2 + (star_dec_ii-star_dec_kk)**2) * 3600. - true_dist) < 3.:
		#print x[k],x[i],math.sqrt( (star_ra_ii-star_ra_kk)**2 + (star_dec_ii-star_dec_kk)**2) * 3600.,true_dist,close_stars[k][kk],close_stars[i][ii]
		#final_stars[k].append(str(k) + '_' + str(i) + '_' + str(close_stars[k][kk])+'_'+str(close_stars[i][ii]))
		s.append(str(close_stars[k][kk])+'_'+str(close_stars[i][ii]))
	  Stars[str(k) + '_' + str(i)] = s
	  #exit()
  #print Stars
  #exit()
  FinalStars = collections.OrderedDict()
  add_str1 = ''; add_str2 = ''

  for k, v in Stars.items():
    #print k, v
    star1 = k.split('_')[0] 
    star2 = k.split('_')[1] 
    if star1 not in FinalStars:
      FinalStars[star1] = []
      #add_str1 = ''
    else:
      add_str1 = add_str1 #+ '*'
    if star2 not in FinalStars:
      FinalStars[star2] = []
      #add_str2 = ''
    else:
      add_str2 = add_str2 #+ '*'   
      
    values1 = []; values2 = []
    for value in v:
      st1 = value.split('_')[0]
      st2 = value.split('_')[1]
      values1.append(st1)
      values2.append(st2)    
    values1 = list(set(values1))
    values2 = list(set(values2))
    #print values1
    #exit()
    
    for value in values1:
     if value in FinalStars[star1]:
      FinalStars[star1].append( value + '*')
     else:
      FinalStars[star1].append( value )       

    for value in values2:
     if value in FinalStars[star2]:
      FinalStars[star2].append( value + '*')
     else:
      FinalStars[star2].append( value )   
  #print FinalStars,'\n'

  for k, v in FinalStars.items():
    #print k,v
    for value in v:
      print value
      if '*' in value:
	numb = value.split('*')[0]
	print 'Corresponding star for: %i,%i is %f,%f' % (x[int(k)],y[int(k)],float(table["RAJ2000"][int(numb)]), float(table["DEJ2000"][int(numb)]))
    

#main('3351subsky1.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('3379subsky7.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('4214sky2.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('4565subsky12.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('4638subsky3.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('5719subsky11.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('6946subsky6.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('7331subsky4.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('891subsky3.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('M101subsky6.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('M104subsky12.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('M51subsky3.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('M63subsky3.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
##main('MasterNGC7463.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('MasterNGC7479.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#main('SkyFitMasterNGC3628.fits', 1.7, 'R', region_file='stars.reg', output_image=None, check_reg=False)
#import glob
#arr = glob.glob('SkyfitMaster*.fits')
#arr = glob.glob('SkyfitMasterNGC247_wcs_crop.fits')
#print arr[0]
#exit()
#main(arr[0], 1.7, 'R', region_file='ds9.reg', output_image=None, check_reg=False)


#main('FinalNGC4302.fits', 0.83, 'R', region_file='astrom_stars.reg', output_image=None, check_reg=False)
#main('FinalNGC4096.fits', 0.83, 'R', region_file='astrom_stars.reg', output_image=None, check_reg=False)
#main('NGC1400&1407 copy.fits', 0.83, 'R', region_file='astrom_stars.reg', output_image=None, check_reg=False)
#main('NGC2683.fits', 0.83, 'R', region_file='astrom_stars.reg', output_image=None, check_reg=False)
#main('NGC4038.fits', 0.83, 'R', region_file='astrom_stars.reg', output_image=None, check_reg=False)
#main('SkyfitMasterNGC596 copy.fits', 0.83, 'R', region_file='astrom_stars.reg', output_image=None, check_reg=False)
#main('NGC4395.fits', 0.83, 'R', region_file='astrom_stars.reg', output_image=None, check_reg=False)
#main('NGC3079.fits', 0.83, 'R', region_file='astrom_stars.reg', output_image=None, check_reg=False)
#main('NGC9257oct16x2y1 copy.fits', 0.83, 'R', region_file='astrom_stars.reg', output_image=None, check_reg=False)
main('NGC3115-final.fits', 0.83, 'R', region_file='astrom_stars.reg', output_image=None, check_reg=False)