#!/usr/bin/env python

from math import sin
from math import cos
from math import pi
from math import modf
from math import sqrt
import numpy as np
import argparse
from types import SimpleNamespace

from astropy.io import fits
from photutils.isophote import IsophoteList
from photutils.isophote import build_ellipse_model
from astropy.stats import sigma_clipped_stats



def main(azim_model_file, output_model, ell=0., posang=0.):
        data = fits.getdata('sky_subtr.fits')
        sma,inten,inten_err = np.loadtxt(azim_model_file, usecols=[0,1,2], dtype=float, unpack=True, skiprows = 1, delimiter='\t')
        ellipses = IsophoteList([])
        
        nx = int(2*np.max(sma))
        ny = nx
        
        xc = float(nx)/2.
        yc = float(ny)/2.
        
        
        for k in range(len(sma)):
            ellipse = SimpleNamespace(x0=int(xc), y0=int(yc), sma=sma[k], intens=0,
                                    eps=ell, pa=np.radians(posang), grad=0, a3=0, b3=0, a4=0, b4=0, fit_ok=False)       
            
            ellipse.intens = inten[k]
            
            ellipses.append(ellipse)


        model = build_ellipse_model(data.shape, ellipses)
        fits.PrimaryHDU(data=model).writeto(output_model, overwrite=True)    