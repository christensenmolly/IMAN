import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage as ndi
from astropy.visualization import simple_norm
from astropy.modeling import models
import photutils
import time
import statmorph
import sys
import argparse
from astropy.io import fits as pyfits
import os
from statmorph.utils.image_diagnostics import make_figure

from photutils import make_source_mask
from astropy.stats import sigma_clipped_stats

LOCAL_DIR = "/non-param_indicators"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))

import get_galaxy_polygon

tmp_out = sys.stdout

def save_res_statmorph(morph, file):
    f = open(file, "w") 
    sys.stdout = f
    print('xc_centroid =', morph.xc_centroid)
    print('yc_centroid =', morph.yc_centroid)
    print('ellipticity_centroid =', morph.ellipticity_centroid)
    print('elongation_centroid =', morph.elongation_centroid)
    print('orientation_centroid =', morph.orientation_centroid)
    print('xc_asymmetry =', morph.xc_asymmetry)
    print('yc_asymmetry =', morph.yc_asymmetry)
    print('ellipticity_asymmetry =', morph.ellipticity_asymmetry)
    print('elongation_asymmetry =', morph.elongation_asymmetry)
    print('orientation_asymmetry =', morph.orientation_asymmetry)
    print('rpetro_circ =', morph.rpetro_circ)
    print('rpetro_ellip =', morph.rpetro_ellip)
    print('rhalf_circ =', morph.rhalf_circ)
    print('rhalf_ellip =', morph.rhalf_ellip)
    print('r20 =', morph.r20)
    print('r80 =', morph.r80)
    print('Gini =', morph.gini)
    print('M20 =', morph.m20)
    print('F(G, M20) =', morph.gini_m20_bulge)
    print('S(G, M20) =', morph.gini_m20_merger)
    print('sn_per_pixel =', morph.sn_per_pixel)
    print('C =', morph.concentration)
    print('A =', morph.asymmetry)
    print('S =', morph.smoothness)
    print('sersic_amplitude =', morph.sersic_amplitude)
    print('sersic_rhalf =', morph.sersic_rhalf)
    print('sersic_n =', morph.sersic_n)
    print('sersic_xc =', morph.sersic_xc)
    print('sersic_yc =', morph.sersic_yc)
    print('sersic_ellip =', morph.sersic_ellip)
    print('sersic_theta =', morph.sersic_theta)
    print('sky_mean =', morph.sky_mean)
    print('sky_median =', morph.sky_median)
    print('sky_sigma =', morph.sky_sigma)
    print('flag =', morph.flag)
    print('flag_sersic =', morph.flag_sersic)
    sys.stdout = tmp_out
    f.close()


def main(input_image, mask_image, xc=None, yc=None, I_min=None, sigma_image=None, psf_image=None, output_file='statmorph_indicators.dat', verbosity=True):
   if verbosity: print('StatMorph analysis...') 
    
   hdulist = pyfits.open(input_image)
   data = hdulist[0].data
   ny,nx = np.shape(data)
   if xc is None or yc is None:
        xc = nx/2.
        yc = ny/2.
   
   
   if mask_image is not None:
        hdulist_mask = pyfits.open(mask_image)
        mask = hdulist_mask[0].data
   else:
       mask = None
   
   if sigma_image is not None:
        hdulist_sigma = pyfits.open(sigma_image)
        sigma = hdulist_sigma[0].data
   else:
        sigma = None
   
   if psf_image is not None:
        hdulist_psf = pyfits.open(psf_image)
        psf = hdulist_psf[0].data
   else:
       psf = None
     
   
   if sigma_image is None:
        gain = 1000.0
   
   
   x0,y0,sma,smb,PA,Flux,segmap = get_galaxy_polygon.main(input_image, mask_image, output_region='galaxy_polygon.reg', xc=xc, yc=yc, I_min=I_min, min_radius=10., verbosity=verbosity)
   os.remove('galaxy_polygon.reg')
   os.remove('galaxy_polygon_ell.reg')
   if sigma is None:
        source_morphs = statmorph.source_morphology(data, segmap, gain=gain, psf=psf, mask=mask)
   else:
        source_morphs = statmorph.source_morphology(data, segmap, psf=psf, weightmap=sigma, mask=mask)


   morph = source_morphs[0]
   save_res_statmorph(morph, output_file)
   
   fig = make_figure(morph)
   fig.savefig('statmorph_results.png', dpi=300)
   
   if verbosity: print('Done!')
   return 
  
  
  
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="StatMorph")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("--sigma_image", nargs='?', const=1, help="Optional: Input sigma image",type=str, default=None)     
    parser.add_argument("--psf_image", nargs='?', const=1, help="Optional: Input psf image",type=str, default=None)   
    parser.add_argument("--mask_image", nargs='?', const=1, help="Optional: Input mask image",type=str, default=None)       
    parser.add_argument("--o", nargs='?', const=1, help="Optional: Output file with results",type=str, default='statmorph_indicators.dat')       
    
    parser.add_argument("--xcen", nargs='?', const=1, help="Optional: Galaxy center x",type=float,default=None)  
    parser.add_argument("--ycen", nargs='?', const=1, help="Optional: Galaxy center y",type=float,default=None)     

    parser.add_argument("--I_min", nargs='?', const=1, help="Optional: SB of the outer isophote (DN)", type=float, default=None)  
    
    args = parser.parse_args()

    input_image = args.input_image
    sigma_image = args.sigma_image
    psf_image = args.psf_image
    mask_image = args.mask_image
    output_file = args.o
    xcen = args.xcen
    ycen = args.ycen
    I_min = args.I_min

    main(input_image, mask_image, xc=xcen, yc=ycen, I_min=I_min, sigma_image=sigma_image, psf_image=psf_image, output_file=output_file, verbosity=True)  
  
  
    
