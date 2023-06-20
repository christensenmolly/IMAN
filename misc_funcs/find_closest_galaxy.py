#! /usr/bin/env python

import sys
from os.path import exists
from os import remove
from scipy.ndimage.filters import convolve
from numpy import max, sum
import numpy as np
import astropy.io.fits as pyfits
import os
import math
import warnings
warnings.filterwarnings("ignore")


LOCAL_DIR = "/misc_funcs"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'detect_objects'))
sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))
import get_galaxy_center_ned
import read_data





sample_file = '/home/amosenko/MyWork/Edge_on_HERON/HyperLeda_meandata_1574360197.txt'




def projected_density(neighbour_dists, n=5):
    # up to n=10!
    r_n = neighbour_dists[n-1]
    return float(n) / ( math.pi * r_n**2 )    

def volume_density(neighbour_dists, n=5):
    # up to n=10!
    r_n = neighbour_dists[n-1]
    return float(n) / ( (4./3.) * math.pi * r_n**3 )




def calc_scale(Dist):
    return Dist*1000./206265. # kpc per arcsec

def convert_modbest_to_D(modbest):    
    return 10**((modbest-25.)/5.) 

def find_distance_between_galaxies(D_A, D_B, alpha):
    return np.sqrt(D_A**2 + D_B**2 - 2.*D_A*D_B*np.cos(np.radians(alpha)))



def calculate_alpha(ra1, dec1, ra2, dec2):
    return math.sqrt((ra1-ra2)**2 + (dec1-dec2)**2)







def main_volume(galaxy_name, D_galaxy=None, diam_ratio=10., min_dist=2.0):
    ra,dec = get_galaxy_center_ned.get_coords_by_name(galaxy_name)

    data, units = read_data.main(sample_file, delimiter='|', header_line=0, units_line=None, skip_lines = [], comment='#')

    
    pgc_HL = np.array(data['pgc'], int)
    ra_HL = np.array(data['al2000'], float)*15.
    dec_HL = np.array(data['de2000'], float)
    logd25_HL = np.array(data['logd25'], float)
    modbest_HL = np.array(data['modbest'], float)
    cz_HL = np.array(data['v'], float)

    # get distance to galaxy:
    for k in range(len(ra_HL)):
        if math.sqrt((ra_HL[k]-ra)**2 + (dec_HL[k]-dec)**2)*3600.<10.:
            D_galaxy = convert_modbest_to_D(modbest_HL[k])
            D25_galaxy = 6.*10**(logd25_HL[k]) * calc_scale(D_galaxy)
            pgc_galaxy = pgc_HL[k]
            print(6.*10**(logd25_HL[k]))
            break
    
    if D_galaxy is None:
        print('ERROR: No distance found for galaxy %s. Exiting ...' % (galaxy_name))
        exit()
    

    closest_galaxies = []
    distances_between = []
    
    print('\nAll galaxies within R=%.2f Mpc from %s and with a ratio of their diameters larger than %.2f:' % (min_dist,galaxy_name,diam_ratio))
    for k in range(len(ra_HL)):
        if ra_HL[k]>=ra-15. and ra_HL[k]<=ra+15. and dec_HL[k]>=dec-15. and dec_HL[k]<=dec+15.:
            alpha = math.sqrt((ra_HL[k]-ra)**2 + (dec_HL[k]-dec)**2)
            
            Dist = find_distance_between_galaxies(D_galaxy, convert_modbest_to_D(modbest_HL[k]), alpha)
            D25_gal = calc_scale(convert_modbest_to_D(modbest_HL[k])) * 6.*10**(logd25_HL[k])
            
            ratio = max([D25_gal,D25_galaxy])/min([D25_gal,D25_galaxy])

            if Dist<min_dist and pgc_galaxy!=pgc_HL[k] and ratio<=diam_ratio and not np.isnan(modbest_HL[k]):
                distances_between.append(Dist)
                closest_galaxies.append(k)
                print('\t PGC%s: Dist=%.2f Mpc, D25=%.2f kpc' % (pgc_HL[k], Dist, D25_gal))
                #print(alpha*60.)
                print(6.*10**(logd25_HL[k]), ra_HL[k], dec_HL[k],alpha)
    
    
    if closest_galaxies == []:
        print('None!')
        return float('nan'),'---'
    else:
        closest_galaxy = distances_between.index(min(distances_between))
        print('\nThe closest galaxy to %s is PGC%s (%f,%f):' % (galaxy_name, pgc_HL[closest_galaxies[closest_galaxy]],ra_HL[closest_galaxies[closest_galaxy]],dec_HL[closest_galaxies[closest_galaxy]]))
        print('\tDistance between them %.2f Mpc:' % (distances_between[closest_galaxy]))
        print('\tTheir distances to us %.2f Mpc and %.2f Mpc:' % (D_galaxy, convert_modbest_to_D(modbest_HL[closest_galaxies[closest_galaxy]])))
        print('\tTheir D25 sizes: %.2f and %.2f kpc\n' % (D25_galaxy,calc_scale(convert_modbest_to_D(modbest_HL[closest_galaxies[closest_galaxy]]))* 6.*10**(logd25_HL[closest_galaxies[closest_galaxy]])))

        return distances_between[closest_galaxy], pgc_HL[closest_galaxies[closest_galaxy]]











def main_nearest_neighbors(galaxy_name, D_galaxy=None, delta_cz=None, method='projected', n_neighbors=5):
    # Get coodinates by galaxy name:
    ra,dec = get_galaxy_center_ned.get_coords_by_name(galaxy_name)

    # Read in the HyperLeda database
    data, units = read_data.main(sample_file, delimiter='|', header_line=0, units_line=None, skip_lines = [], comment='#')
    
    pgc_HL = np.array(data['pgc'], int)
    ra_HL = np.array(data['al2000'], float)*15.
    dec_HL = np.array(data['de2000'], float)
    logd25_HL = np.array(data['logd25'], float)
    modbest_HL = np.array(data['modbest'], float)
    cz_HL = np.array(data['v'], float)


    # Get galaxy distance, diamtere, pgc name, cz and Scale (kpc/arcsec):
    for k in range(len(ra_HL)):
        if math.sqrt((ra_HL[k]-ra)**2 + (dec_HL[k]-dec)**2)*3600.<10.:
            D_galaxy = convert_modbest_to_D(modbest_HL[k])
            D25_galaxy = 6.*10**(logd25_HL[k]) * calc_scale(D_galaxy)
            pgc_galaxy = pgc_HL[k]
            cz_galaxy = cz_HL[k]
            Scale_galaxy = calc_scale(D_galaxy)
            break


    if method=='projected':
        
        neighbour_dists = []
       
        if delta_cz is None:
            for k in range(len(ra_HL)):
                if ra_HL[k]>=ra-15. and ra_HL[k]<=ra+15. and dec_HL[k]>=dec-15. and dec_HL[k]<=dec+15.:
                    # Calculate angular distance between two galaxies, in arcsec
                    alpha = calculate_alpha(ra_HL[k], dec_HL[k], ra, dec)*3600.
                    
                    # Convert angular to projected distance, taking that the distances to the other galaxy is the same as to the target galaxy:
                    projected_distance = Scale_galaxy * alpha / 1000. # Now in Mpc
                    
                    neighbour_dists.append(projected_distance)
        else:
            for k in range(len(ra_HL)):
                if ra_HL[k]>=ra-15. and ra_HL[k]<=ra+15. and dec_HL[k]>=dec-15. and dec_HL[k]<=dec+15. and abs(cz_galaxy-cz_HL[k])<delta_cz:
                    # Calculate angular distance between two galaxies, in arcsec
                    alpha = calculate_alpha(ra_HL[k], dec_HL[k], ra, dec)*3600.
                    
                    # Convert angular to projected distance, taking that the distances to the other galaxy is the same as to the target galaxy:
                    projected_distance = Scale_galaxy * alpha / 1000. # Now in Mpc
                    
                    neighbour_dists.append(projected_distance)
                    
                    

        # Sort the list of projected distances:        
        neighbour_dists_ordered = sorted(neighbour_dists)
        
        # Find the projected density:
        projected_sigma = projected_density(neighbour_dists_ordered, n=n_neighbors)
        
        print('Projected density for galaxy %s: %.1f Mpc-2' % (galaxy_name, projected_sigma))
        
        return projected_sigma

    elif method=='volume':
        
        neighbour_dists = []
       
        for k in range(len(ra_HL)):
                if ra_HL[k]>=ra-15. and ra_HL[k]<=ra+15. and dec_HL[k]>=dec-15. and dec_HL[k]<=dec+15. and abs(cz_galaxy-cz_HL[k])<delta_cz:
                    # Calculate angular distance between two galaxies, in deg
                    alpha = calculate_alpha(ra_HL[k], dec_HL[k], ra, dec)
                    
                    # Calculate physical distance between the galaxies
                    Dist = find_distance_between_galaxies(D_galaxy, convert_modbest_to_D(modbest_HL[k]), alpha)
                                        
                    neighbour_dists.append(Dist)
                    


        # Sort the list of projected distances:        
        neighbour_dists_ordered = sorted(neighbour_dists)
        
        
        # Find the 3D density:
        volume_sigma = volume_density(neighbour_dists_ordered, n=n_neighbors)
        
        print('Volume density for galaxy %s: %.3f Mpc-3' % (galaxy_name, volume_sigma))
        
        return volume_sigma




    

if __name__ == "__main__":
    galaxy_name = sys.argv[1]
    #main(galaxy_name, D_galaxy=None)
    
    main_nearest_neighbors(galaxy_name, D_galaxy=None, delta_cz=500., method='volume', n_neighbors=5)

    
