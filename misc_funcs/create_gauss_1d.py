#!/usr/bin/env python
# EXAMPLE:

from math import sin
from math import cos
from math import pi
from math import modf
from math import sqrt
from math import ceil
import numpy as np
import argparse
from types import SimpleNamespace

from astropy.io import fits
from photutils.isophote import IsophoteList
from photutils.isophote import build_ellipse_model
from astropy.stats import sigma_clipped_stats
from astropy.modeling.models import Gaussian1D


def main(FWHM, radius_fwhm=5., output_file='gauss_profile.txt'):
    s1 = Gaussian1D(amplitude=1.,mean=0.,stddev=FWHM/2.355)
    r = np.arange(0., radius_fwhm*FWHM, 1.)
    I = s1(r)
    I_sum = np.sum(I)

    f_res = open(output_file, 'w') #### Output model is saved in a text-file as well
    f_res.truncate(0)
    f_res.write("# sma[pix]\tflux[DN]\tflux_err[DN]\n")


    for k in range(len(r)):
        f_res.write("%7.2f\t%15.8e\t%15.8e\n" % (r[k], I[k]/I_sum, 0.))
 
    f_res.close()

    return r,I/I_sum



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("FWHM", help="FWHM (in pixels)",
                        type=float)
    parser.add_argument("--radius", help="Radius (in FWHMs)",
                        type=float, default=5.)
    parser.add_argument("--output_file", help="Name of output file",
                        type=str, default='gauss_profile.txt')



    args = parser.parse_args()

    FWHM = args.FWHM
    radius_fwhm = args.radius
    output_file = args.output_file
    
    
    
    main(FWHM, radius_fwhm, output_file)
