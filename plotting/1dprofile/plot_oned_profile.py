import os
import sys

from pylab import *
import astropy.io.fits as pyfits

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import glob
from matplotlib.pyplot import cm
import re
from scipy.odr.odrpack import *
import argparse
from scipy.ndimage import rotate
from astropy.stats import sigma_clipped_stats
#matplotlib.use('TkAgg')
import warnings
#PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
#PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
#sys.path.append(PATH_TO_PACKAGE+'/FindFluxes')

#import iraf_fit_ellipse
import radial_profile
#import ds9_contour

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

#matplotlib.use('agg',warn=False, force=True)   ####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

warnings.filterwarnings("ignore")

COLORS = ['b','g','m','y','c','lime','orange']
line_styles = ['-.','--',':','-.','-.',':','--']


def main(input_image, m0, pix2sec, output_file='azim_profile.png', geom_units='arcsec', SB_units='mag/arcsec2', legend_size=10, x_scale='linear'):

      hdu = pyfits.open(input_image)
      number_of_layers = len(hdu) 

      fig = plt.figure(0,figsize=(5, 5))
      ax = fig.add_subplot(111)
          
      if x_scale=='linear':
            plt.xlabel(r' $r$ (%s) ' % (geom_units), fontsize=15)
      elif x_scale=='1/4':
            plt.xlabel(r' $r^{1/4}$ (%s$^{1/4}$) ' % (geom_units), fontsize=15)
      elif x_scale=='log':
            plt.xlabel(r' $\log\,r$ (%s) ' % (geom_units), fontsize=15)
          
      if SB_units=='mag/arcsec2':
            plt.ylabel(r' $\mu$ (mag arcsec$^{-2}$) ', fontsize=15)
      elif SB_units=='ADU/pix2':
            plt.ylabel(r' Intensity (ADU) ', fontsize=15)
      plt.gca().invert_yaxis()
      
      for layer in range(number_of_layers):

            if layer>3:
              color = COLORS[layer-4]
              lyne_style = line_styles[layer-4]


            data = hdu[layer].data
            prihdr = hdu[layer].header
            ny,nx = np.shape(data)
            

            try:
              Label = prihdr['NAME_OF_LAYER']
            except:
              Label = 'galaxy'

            if ny==2:
                    I = data[0,:]
                    r = np.arange(nx) + 1
            if nx==2:
                    I = data[:,0]
                    r = np.arange(ny) + 1


            if ny==1:
                    I = data[0,:]
                    r = np.arange(nx) + 1
            if nx==1:
                    I = data[:,0]
                    r = np.arange(ny) + 1
            
            if geom_units=='arcsec':
                r = r * pix2sec
            
            if x_scale=='1/4':
                r = r**(0.25)
            
            if x_scale=='log':
                r = np.log10(r)
                
            if SB_units=='mag/arcsec2':
                I = m0 - 2.5*np.log10(I/(pix2sec**2))
            
            if layer==0:
                ax.plot(r, I, 'o', color='gray', markeredgecolor='gray', markersize=3., label = Label)
            elif layer==1:
                ax.plot(r, I, ls='-', color='red', lw=2, label = Label)
            else:
                if layer>5:
                    ax.plot(r, I, lyne_style, color=color, lw=2, label = Label)
      ax.legend(loc=1, borderaxespad=0., fontsize=legend_size, numpoints=1)
      plt.savefig(output_file, transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)    
      plt.clf()
      plt.close()                
                
                  
              
          