#!/usr/bin/python
# DESCRIPTION:
# Script to do photometric calibration. This script finds good PSF stars
# and searches them in a catalogue.
# NOTE: Run frist calibration.py to create a region file with PSF stars!!!
# In the region file only the following region is alowed (annulus):
# ellipse(xc,yc,sma_in,smb_in,sma_out,smb_out,PA) - sky will be estimated within the annulus.
# NOTE: COLOUR is not applicable at the moment. Use just band!
# CAUTION: ugriz are in AB-system! UBVRI are Landolt (VEGA)!
# MINIMAL USAGE: python calibration.py [input_image]
# python ~/MEGA/MyPrograms/IMAN/imp/phot_calibration/calibration_with_colour.py new-image.fits --mask_image ini_segm.fits  --catalogue NOMAD --region_file best_stars.reg --band R
# python ~/MEGA/MyPrograms/IMAN/imp/phot_calibration/calibration_with_colour.py new-image.fits --mask_image ini_segm.fits  --catalogue NOMAD,V/139,II/349/ps1,I/345/gaia2 --region_file best_stars.reg --band R,r,r,r
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
import astropy.units as un 
from astroquery import ned
from astropy.coordinates import SkyCoord
import extinction
from dustmaps.bayestar import BayestarQuery
from scipy.interpolate import interp1d

# Import additional modules
LOCAL_DIR = "/imp/phot_calibration"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'sextractor'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))
sys.path.append(os.path.join(IMAN_DIR, 'decomposition/simple_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'extinction_maps'))


#import run_SExtractor
import rebin_image
import add_keyw_to_header
import photometric_conversions
import Arenou_model


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def retrieve_l_b(RA, DEC):
    c = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='fk5').galactic
    return c.l.degree,c.b.degree




def show_complete(i,N):
  percent = 100. * float(i) / float(N)
  sys.stdout.write("\r%2d%%" % percent)
  sys.stdout.flush()



def retrieve_extinction(RA, DEC, D):
        l,b = retrieve_l_b(RA, DEC)
        Aexts = schlafly_extinction.main([[l,b]], wavelengths=[wavelength], frame='galactic', units='deg', output_file=None)
        return Aexts[0][0]



def get_cross_star(RA, DEC, Radius=30./3600., catalog='NOMAD'):
            viz = Vizier(keywords=["stars", "optical"])
            viz.ROW_LIMIT = -1

            # Get table with all stars found within Radius
            try:
                c = coordinates.SkyCoord(float(RA),float(DEC),unit=('deg','deg'),frame='fk5')
                if catalog=='I/345/gaia2':
                    viz = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS', 
                             'RPmag', 'BPmag', 'Gmag', 'Plx']) 
                    c = c.transform_to('icrs')
                    result_dist = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS', 
                             'rest']).query_region(c, radius=Radius*un.deg, catalog='I/347/gaia2dis') # Take distance from Bailer-Jones
                    table_dist = result_dist[0]
                    #print(table_dist)
                
                result = viz.query_region(c, radius=Radius*un.deg, catalog=catalog)
                table = result[0]
                #print(table)
                
            except:
                print('ERROR: The reference coordinate is not found in the catalogue!')
                return None
            
            #exit()
            #print(len(table))    
            
            # Cross retrieved stars with the given star (by minimal distance):
            RAA = []; DECC = [] #; MAG = []; COLOUR = []
            for ii in range(len(table)):               
                if catalog=='NOMAD': # in VEGA                
                  raa = float(table["RAJ2000"][ii])
                  decc = float(table["DEJ2000"][ii])
                  RAA.append(raa) #### WARNING: was _RAJ2000
                  DECC.append(decc) #### DEJ2000  
               
                elif catalog=='V/139': # SDSS in AB-system
                  gc = coordinates.SkyCoord(ra=float(table["RA_ICRS"][ii]), dec=float(table["DE_ICRS"][ii]), unit='deg', frame='icrs')
                  gc = gc.transform_to('fk5') 
                  raa = float(gc.ra.degree) #### WARNING: was _RAJ2000
                  decc = float(gc.dec.degree) #### DEJ2000
                  RAA.append(raa) #### WARNING: was _RAJ2000
                  DECC.append(decc) #### DEJ2000

                elif catalog=='II/349/ps1': # PANSTARRS in AB-system
                  raa = float(table["RAJ2000"][ii])  #### WARNING: was _RAJ2000
                  decc = float(table["DEJ2000"][ii]) #### DEJ2000
                  RAA.append(raa) #### WARNING: was _RAJ2000
                  DECC.append(decc) #### DEJ2000

                elif catalog=='I/345/gaia2': # GAIA DR2
                  gc = coordinates.SkyCoord(ra=float(table["RA_ICRS"][ii]), dec=float(table["DE_ICRS"][ii]), unit='deg', frame='icrs')
                  gc = gc.transform_to('fk5') 
                  raa = float(gc.ra.degree) #### WARNING: was _RAJ2000
                  decc = float(gc.dec.degree) #### DEJ2000
                  RAA.append(raa) #### WARNING: was _RAJ2000
                  DECC.append(decc) #### DEJ2000


                    
            RAA = np.array(RAA, float)
            DECC = np.array(DECC, float)

            dist = np.sqrt( (RA-RAA)**2 + (DEC-DECC)**2) * 3600.
            dist = np.array(dist, float)
            

            
            try:
                min_dist = dist.min()
                index = np.where(dist == min_dist)[0][0]  # number of star in Table which we need!
            except:
                return None


            if min_dist>2.:
                # NOTE: In this case the minimum distance is too large and porbably this star does not have a good match with the catalogue. OR problem with WCS? 
                return None            
            

            if catalog=='NOMAD': # in VEGA
                    U =float('nan')
                    B =float(table["Bmag"][index])
                    V =float(table["Vmag"][index])                  
                    R =float(table["Rmag"][index])                  
                    I =float('nan')
                    return [U,B,V,R,I]

            elif catalog=='V/139': # SDSS in AB-system
                    u = float(table["umag"][index])
                    g = float(table["gmag"][index])
                    r = float(table["rmag"][index])
                    i = float(table["imag"][index])
                    z = float(table["zmag"][index])
                    y = float('nan')
                    return [u,g,r,i,z,y]

            elif catalog=='II/349/ps1': # PANSTARRS in AB-system
                    g = float(table["gmag"][index])
                    r = float(table["rmag"][index])
                    i = float(table["imag"][index])
                    z = float(table["zmag"][index])
                    y = float(table["ymag"][index])
                    return [g,r,i,z,y]
                
            elif catalog=='I/345/gaia2': # GAIA DR2
                    GBP = float(table["BPmag"][index])
                    G = float(table["Gmag"][index])
                    GRP = float(table["RPmag"][index])
                    parallax = float(table["Plx"][index])
                    try:
                        ind = list(table_dist["Source"]).index(str(table["Source"][index]))
                        D = float(table_dist["rest"][ind])
                    except:
                        D = math.fabs(1000./(parallax+0.029)) # See Lindegren et al. 2018
                    return [GBP, G, GRP, D]










def back_in_annulus(inner_ellipse, outer_ellipse, inframe, mask, xSize, ySize):

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

        
        
        mean_sigclip, median_sigclip, std_sigclip = sigma_clipped_stats(annulus_data_1d, mask=mask_astropy, sigma=3, iters=5)  
        return median_sigclip, std_sigclip






def ext_func(Av, Rv, wavelength, ext_map):
    # Wavelength in um

    wavelength = float(wavelength)*10000. # in A
    if ext_map=='Arenou':
        return extinction.ccm89(np.array([wavelength*10000.]), Av, Rv)[0]
    if ext_map=='Green':
        # Here Rv_* is Ab/E(B-V)SFD
        lambda_eff, Rv_21, Rv_31, Rv_41, Rv_51 = np.loadtxt(
            os.path.join(IMAN_DIR, 'extinction_maps/Table_6_SF2011.dat'), usecols=[1, 2, 3, 4, 5],
            unpack=True, skiprows=1, dtype=float, delimiter='\t')
        if Rv==2.1:
            Ab_EBV = Rv_21

        if Rv==3.1:
            Ab_EBV = Rv_31

        if Rv==4.1:
            Ab_EBV = Rv_41        

        if Rv==5.1:
            Ab_EBV = Rv_51        

        inds = lambda_eff.argsort()
        Ab_EBV = Ab_EBV[inds]
        lambda_eff = lambda_eff[inds]

        f = interp1d(lambda_eff, Ab_EBV)  
        egr = Av/2.682
        return egr*f(wavelength)
 
 





def main(input_image, region_file, mask_file, bands=['R'], catalogues=['NOMAD'], star_number=None, obs_wavelength=0.655):
        # WARNING: At the moment, colour is not taken into account, so this parameters should be omitted.
        
        print('Photometric calibration by the region file with stars...')
 
        if 'I/345/gaia2' not in catalogues:
               catalogues = ['I/345/gaia2'] + catalogues 
 

 
        # Read in the input image
        hdulist = pyfits.open(input_image)
        header1 = hdulist[0].header
        scidata = hdulist[0].data
        ySize,xSize = np.shape(scidata)
        
        # Read in the mask image
        hdulist_mask = pyfits.open(mask_file)
        mask = hdulist_mask[0].data
        if mask is None:
            mask = hdulist_mask[1].data
        
        # Read in important keywords EXPTIME and NCOMBINE:
        if 'EXPTIME' in header1:
            exptime = header1['EXPTIME']
        else:
            exptime = float('nan')

        if 'NCOMBINE' in header1:
            ncombine = header1['NCOMBINE']        
        else:
            ncombine = float('nan')

        # Read in WCS:
        w = wcs.WCS(header1)

        # Determine image resolution (pix/arcsec):
        pix2sec,note = rebin_image.resolution(input_image)
        pix2sec = float(pix2sec)


        # Read in region file:
        ff = open(region_file, "r")
        lines = ff.readlines()
        
        # Read region file line by line:
        results = collections.OrderedDict()
        for kk in range(len(catalogues)):
            results['m0:%s:%s:%s' % (catalogues[kk],bands[kk],str(None))] = []

            
        # Define extinction map from Green et al. (2019)
        bayestar = BayestarQuery(max_samples=1)
        
        N = 0
        Number_of_lines = len(lines)
        
        
        ZP = []
        for k in range(Number_of_lines):
           # Do each star:
           show_complete(k+1, Number_of_lines) # Progress status
           
           if 'ellipse(' in lines[k] or 'annulus(' in lines[k]:
                N = N + 1
                star = lines[k].split(',')
                xc_star = float(star[0].split('(')[1])
                yc_star = float(star[1])
                sma_star_in = float(star[2])
                smb_star_in = float(star[3])
                sma_star_out = float(star[4])
                smb_star_out = float(star[5])
                PA_star_out = float(star[6].split(')')[0])
                PA_star_in = PA_star_out
     
                if 'text={' in lines[k]:
                    N = int(lines[k].split('text={')[1].split('}')[0])
                
                if star_number is not None:
                    if N!=star_number:
                        continue
           
           else: continue  
                    
              
           # Find total flux of the star taking into account the background:     
           cen = Point(xc_star, yc_star)
           aperture = EllipticalAperture([(xc_star,yc_star)], a=sma_star_in, b=smb_star_in, theta=PA_star_in)
                
           bkg_mean,bkg_std = back_in_annulus([cen, sma_star_in, smb_star_in, PA_star_in], [cen, sma_star_out, smb_star_out, PA_star_out], scidata, mask, xSize, ySize)
           try:
                bkg_sum = bkg_mean * aperture.area()
           except:
                bkg_sum = bkg_mean * aperture.area
           
           phot_table = aperture_photometry(scidata, aperture)
           final_sum = phot_table['aperture_sum'] - bkg_sum
   
           pixcrd = np.array([[xc_star, yc_star]], np.float_)
           world = w.wcs_pix2world(pixcrd, 1)
           RA,DEC = world[0,0],world[0,1]
           l,b = retrieve_l_b(RA, DEC)
           #print(l,b)
           
           flux = float(final_sum[0])

           MAG = []; COLOUR = []

           egr = bayestar(coords = SkyCoord(l*un.deg, b*un.deg, distance=1000.*un.pc, frame='galactic'), mode='best')

           if np.isnan(egr):
                Av = Arenou_model.Arenou([[l,b,1000.]])[0][0]
                ext_map = 'Arenou'
           else:
                Av = egr * 2.682
                ext_map = 'Green'

           Rv = 3.1
           flux_corr = flux * (10**(0.4*ext_func(Av,Rv,obs_wavelength,ext_map)))
           
           zp = []

           for kk in range(len(catalogues)):
                cross_res = get_cross_star(RA, DEC, Radius=30./3600., catalog=catalogues[kk])

                if cross_res is not None:
                    if catalogues[kk] == 'I/345/gaia2': # Should go first as we need to know distances for other catalogues
                        [mag_GBP, mag_G, mag_GRP, D] = cross_res

                        egr = bayestar(coords = SkyCoord(l*un.deg, b*un.deg, distance=D*un.pc, frame='galactic'), mode='best')

                        if np.isnan(egr):
                            Av = Arenou_model.Arenou([[l,b,D]])[0][0]
                            ext_map = 'Arenou'
                        else:
                            Av = egr * 2.682
                            ext_map = 'Green'
                        #print(l,b,D,"%.5f" % (egr), 'HEREEEEEEEEE')                            

                        
                        mag_GBP = mag_GBP - ext_func(Av,Rv,0.532,ext_map)
                        mag_G = mag_G - ext_func(Av,Rv,0.673,ext_map)
                        mag_GRP = mag_GRP - ext_func(Av,Rv,0.797,ext_map)

                        mag_g,mag_r,mag_i = photometric_conversions.transformation_from_gaia_to_sdss(mag_G, mag_GRP, mag_GBP)
                        mag_y = float('nan')
                        mag_u = float('nan')
                        mag_z = float('nan')
                            
                        mag_U = float('nan')
                        mag_I = float('nan')
                        mag_V = mag_G
                        mag_B = mag_GBP
                        mag_R = mag_GRP

                        flux_corr = flux * (10**(0.4*ext_func(Av,Rv,obs_wavelength,ext_map)))
                        #print(ext_func(Av,Rv,obs_wavelength))
                        mag = photometric_conversions.get_mag(mag_u, mag_g, mag_r, mag_i, mag_z, mag_y, mag_U, mag_B, mag_V, mag_R, mag_I, bands[kk])
                        zp.append(mag + 2.5*math.log10(flux_corr))         
                    
                    if catalogues[kk] == 'NOMAD':
                        [mag_U,mag_B,mag_V,mag_R,mag_I] = cross_res
                        mag_U = mag_U - ext_func(Av,Rv,0.365,ext_map)
                        mag_B = mag_B - ext_func(Av,Rv,0.445,ext_map)
                        mag_V = mag_V - ext_func(Av,Rv,0.551,ext_map)
                        mag_R = mag_R - ext_func(Av,Rv,0.658,ext_map)
                        mag_I = mag_I - ext_func(Av,Rv,0.806,ext_map)
                        
                        mag_u,mag_g,mag_r,mag_i,mag_z = photometric_conversions.transformation_from_UBVRI_to_sdss(mag_U, mag_B, mag_V, mag_R, mag_I) # The output values are in the AB-system
                        mag_y = float('nan')
                        
                        
                        mag = photometric_conversions.get_mag(mag_u, mag_g, mag_r, mag_i, mag_z, mag_y, mag_U, mag_B, mag_V, mag_R, mag_I, bands[kk])
                        zp.append(mag + 2.5*math.log10(flux_corr))
                        
                    if catalogues[kk] == 'V/139':
                        [mag_u,mag_g,mag_r,mag_i,mag_z,mag_y] = cross_res
                        mag_u = mag_u - ext_func(Av,Rv,0.354,ext_map)
                        mag_g = mag_g - ext_func(Av,Rv,0.475,ext_map)
                        mag_r = mag_r - ext_func(Av,Rv,0.622,ext_map)
                        mag_i = mag_i - ext_func(Av,Rv,0.763,ext_map)
                        mag_z = mag_z - ext_func(Av,Rv,0.905,ext_map)
                        mag_y = mag_y - ext_func(Av,Rv,0.962,ext_map)

                        mag_U,mag_B,mag_V,mag_R,mag_I = photometric_conversions.transformation_from_sdss_to_UBVRI(mag_u, mag_g, mag_r, mag_i, mag_z)
                        
                        mag = photometric_conversions.get_mag(mag_u, mag_g, mag_r, mag_i, mag_z, mag_y, mag_U, mag_B, mag_V, mag_R, mag_I, bands[kk])
                        zp.append(mag + 2.5*math.log10(flux_corr))

                    
                    if catalogues[kk] == 'II/349/ps1':
                        [mag_g,mag_r,mag_i,mag_z,mag_y] = cross_res
                        mag_g = mag_g - ext_func(Av,Rv,0.481,ext_map)
                        mag_r = mag_r - ext_func(Av,Rv,0.617,ext_map)
                        mag_i = mag_i - ext_func(Av,Rv,0.752,ext_map)
                        mag_z = mag_z - ext_func(Av,Rv,0.866,ext_map)
                        mag_y = mag_y - ext_func(Av,Rv,0.962,ext_map)
                        
                        mag_u,mag_g,mag_r,mag_i,mag_z,mag_y = photometric_conversions.transformation_from_panstarrs_to_sdss(float('nan'), mag_g, mag_r, mag_i, mag_z, mag_y)
                        mag_U,mag_B,mag_V,mag_R,mag_I = photometric_conversions.transformation_from_sdss_to_UBVRI(mag_u, mag_g, mag_r, mag_i, mag_z)
                        
                        mag = photometric_conversions.get_mag(mag_u, mag_g, mag_r, mag_i, mag_z, mag_y, mag_U, mag_B, mag_V, mag_R, mag_I, bands[kk])
                        zp.append(mag + 2.5*math.log10(flux_corr))                                      


                else:
                    zp.append(float('nan'))


           ZP.append(zp)
        
        if ZP==[]:
            print('No match of the selected stars with the catalogue. Check WCS! Exiting...')
            #return None
            exit()
            
           
        for k in range(len(ZP)):
            for kk in range(len(catalogues)):
                if not np.isnan(ZP[k][kk]):
                    results['m0:%s:%s:%s' % (catalogues[kk],bands[kk],str(None))].append(ZP[k][kk])

           
        RESULTS = collections.OrderedDict()

        for kk in range(len(catalogues)):
               m0_mean, m0_median, m0_std = sigma_clipped_stats(results['m0:%s:%s:%s' % (catalogues[kk],bands[kk],str(None))], sigma=2, iters=5)
               if m0_median==0.0:
                  m0_median = float('nan')
                  m0_std = float('nan')
               print('%s:%s:%s  %f %f' % (catalogues[kk],bands[kk],str(None),m0_median, m0_std))
               RESULTS['M0:%s:%s:%s' % (catalogues[kk],bands[kk],str(None))] = m0_median
               RESULTS['M0_STD:%s:%s:%s' % (catalogues[kk],bands[kk],str(None))] = m0_std
        RESULTS['N_stars_calibr'] = len(results['m0:%s:%s:%s' % ('I/345/gaia2','r',str(None))])
        print('Number of stars used: %i' % (RESULTS['N_stars_calibr']))
        RESULTS['EXPTIME'] = exptime
        RESULTS['NCOMBINE'] = ncombine
        print('Done!')
        return RESULTS

                


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Do photometric calibration based on the file with stars given by annuli")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("--mask_image", help="Mask fits image", type=str, default=None)  
    parser.add_argument("--bands", help="Photometric band in the catalogue (e.g. B, V, R, J, H, K)", type=str, default='R') 
    parser.add_argument("--catalogues", help="Photometric catalogue of stars (NOMAD, UCAC4, II/246 (2MASS), V/139 (SDSS))", type=str, default='NOMAD')
    parser.add_argument("--colours", help="Colour (NOT WORKING CORRECTLY)", type=str, default=None)     
    parser.add_argument("--star_number", help="Star number from the region file", type=int, default=None)       
    parser.add_argument("--region_file", help="Region file with selected PSF stars", type=str, default=None) 
    args = parser.parse_args()

    input_image = args.input_image
    mask_file = args.mask_image
    catalogues = args.catalogues.split(',')
    bands = args.bands.split(',')
    try:
        colours = args.colours.split(',')
    except:
        colours = len(bands)* [None]
    region_file = args.region_file
    star_number = args.star_number

    main(input_image, region_file, mask_file, bands=bands, colours=colours, catalogues=catalogues, star_number=star_number)       
           
#main('new-image.fits', 'best_stars.reg', 'ini_segm.fits', band='r', colour=None, catalogue='V/139', star_number=None)
#main('new-image.fits', 'best_stars.reg', 'ini_segm.fits', band='R', colour=None, catalogue='NOMAD', star_number=None)  
