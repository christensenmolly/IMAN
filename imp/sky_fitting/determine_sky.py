#!/usr/bin/python
# DESCRIPTION:
# Script to fit sky background with a polynomial of n degree. A mask is required.
# MINIMAL USAGE: python determine_sky.py [input_image] [input_mask]

import numpy as np
from astropy.io import fits as pyfits
from astropy.modeling import models, fitting
import warnings
warnings.filterwarnings("ignore")
from astropy.stats import sigma_clipped_stats
import argparse
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import subprocess
import sys
import os
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from scipy import stats


LOCAL_DIR = "/imp/sky_fitting"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))
sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))
import rebin_image
import arithm_operations

def remove_files(files):
    for file in files:
        if os.path.exists(file):
            os.remove(file)

def convert_segm_to_boolean(mask):
    return np.ma.make_mask(mask, shrink=False)

def split_xyz(data, mask=None, arrays=False):
    # Initialize lists to contain the x, y and z values
    x_values = []
    y_values = []
    z_values = []

    # Loop over all x and y values
    for x in range(data.shape[1]):
        for y in range(data.shape[0]):

            # If no mask is specified or the pixel is not masked, add the coordinates and value to the appropriate lists
            if mask is None or not mask[y,x]:
                x_values.append(x)
                y_values.append(y)
                z_values.append(data[y,x])

    if arrays: return np.array(x_values), np.array(y_values), np.array(z_values)
    else: return x_values, y_values, z_values

def evaluate_model(model, x_min, x_max, y_min, y_max, x_delta=1, y_delta=1):
    # Create x and y meshgrid for evaluating the model
    y_plotvalues, x_plotvalues = np.mgrid[y_min:y_max:y_delta, x_min:x_max:x_delta]

    # Evaluate the model
    evaluated_model = model(x_plotvalues, y_plotvalues)

    # Return the evaluated data
    return evaluated_model



def fit_polynomial(data, degree, mask=None):

    ySize,xSize = np.shape(data)

    # Fit the data using astropy.modeling
    poly_init = models.Polynomial2D(degree=degree)

    fit_model = fitting.LevMarLSQFitter()

    # Split x, y and z values that are not masked
    x_values, y_values, z_values = split_xyz(data, mask=mask, arrays=True)

    # Ignore model linearity warning from the fitter
    with warnings.catch_warnings():

        warnings.simplefilter('ignore')
        poly = fit_model(poly_init, x_values, y_values, z_values)  # What comes out is the model with the parameters set

    # Return the polynomial model
    sky = evaluate_model(poly, 0, xSize, 0, ySize)    
    
    return sky



def do_sky(input_image, mask_image, stats, output_image='sky_sub_interp.fits', sky_image=None, sigma_smooth=None, sampling_factor=None, polynomial_degree=1): # sigma_smooth=30., sampling_factor=5.
    [mean, median, std] = stats

    
    hdulist0 = pyfits.open(input_image)
    data0 = hdulist0[0].data
    header0 = hdulist0[0].header 
    ny0,nx0 = np.shape(data0)

    # Replace nan and inf 
    data0[np.isnan(data0)] = median
    data0[np.isinf(data0)] = median
    
    if mask_image is None:
        mask = np.zeros((ny0,nx0))
        outHDU = pyfits.PrimaryHDU(mask, header0)
        outHDU.writeto('tmp_mask.fits', overwrite=True)           
        mask_image ='tmp_mask.fits'

    hdulist_mask = pyfits.open(mask_image)
    mask0 = hdulist_mask[0].data
    input_image0 = input_image

    
    if sampling_factor is not None:
        rebin_image.downsample(input_image, sampling_factor, output_image='tmp_rebin.fits', set_wcs=True, print_mes=True, norm=False)
        rebin_image.downsample(mask_image, sampling_factor, output_image='tmp_mask_rebin.fits', set_wcs=True, print_mes=True, norm=False, no_interp=True)
        input_image = 'tmp_rebin.fits'
        mask_image = 'tmp_mask_rebin.fits'

    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    header = hdulist[0].header
    ny,nx = np.shape(data)

    # Replace nan and inf 
    data[np.isnan(data)] = median
    data[np.isinf(data)] = median


    hdulist1 = pyfits.open(mask_image)
    mask = hdulist1[0].data  
    mask_astropy = convert_segm_to_boolean(mask)

    
    if sigma_smooth is not None:
        data[mask > 0.] = np.nan
        kernel = Gaussian2DKernel(x_stddev=sigma_smooth)
        astropy_conv = convolve(data, kernel)
    else:
        astropy_conv = data

    sky = fit_polynomial(astropy_conv, int(polynomial_degree), mask=mask_astropy)

    # Add keywords about sky
    header['Sky_med'] = median
    header['Sky_mean'] = mean
    header['Sky_std'] = std      
    
    outHDU = pyfits.PrimaryHDU(sky, header)
    outHDU.writeto('fit_tmp_rebin.fits', overwrite=True)    
    
    if nx!=nx0:
        #rebin_image.downsample('fit_tmp_rebin.fits', float(nx)/float(nx0), output_image='interp_tmp.fits', set_wcs=True, print_mes=True, norm=False)
        rebin_image.rebin(input_image0, 'fit_tmp_rebin.fits', output_image='interp_tmp.fits', hdu_ref=0, hdu_inp=0, preserve_bad_pixels=True)
    else:
        os.rename('fit_tmp_rebin.fits', 'interp_tmp.fits')
    
    arithm_operations.main(input_image0, 'interp_tmp.fits', 'sub', output_image)
    if sky_image is not None:
        os.rename('interp_tmp.fits', sky_image)
    
    
    remove_files(['tmp_rebin.fits','tmp_mask_rebin.fits','fit_tmp_rebin.fits','interp_tmp.fits','tmp_mask.fits'])
       




def sky_subtraction(input_image, mask_image, polynomial_degree=5, output_image='sky_subtr.fits', output_sky=None, hdu_inp=0, sampling_factor=5., sigma_smooth=None, verbosity=True, sky_value='mode'):   
      if verbosity: print('Fitting the background with a polynomial of %i degree...' % (polynomial_degree))
      hdulist = pyfits.open(input_image)
      data = hdulist[hdu_inp].data      
      header = hdulist[hdu_inp].header
      ny,nx = np.shape(data)
      
      
      
      hdulist_mask = pyfits.open(mask_image)
      mask = hdulist_mask[0].data
      mask[np.isnan(data)]=1
      
      # Create astropy mask (boolean)
      mask_astropy = convert_segm_to_boolean(mask)
      

      

      # Find average sky
      mean, median, std = sigma_clipped_stats(data, sigma=3.0, mask=mask_astropy)
      #mode = stats.mstats.mode(data, axis=None) # Does not work properly!
      mode = 3.*median-2.*mean #### WARNING!
      print(mean, median, mode, std)
      print('Mean: %.10f +/- %.10f' % (mean, std))
      print('Median: %.10f' % (median))
      print('Mode: %.10f' % (mode))
      
      # Replace nan and inf 
      data[np.isnan(data)] = mode
      data[np.isinf(data)] = mode
      
      
      # If there are inf values, replace them
      


      MaskedArray = np.ma.array(np.array(data), mask=mask_astropy)
      #mean = np.ma.mean(MaskedArray)
      #median = np.ma.median(MaskedArray)
      
      min_data = MaskedArray.min(axis=None)
      max_data = MaskedArray.max(axis=None)
      
      data_hist = []
      for k in range(ny):
          for i in range(nx):
              if mask[k,i]==0:
                  data_hist.append(data[k,i])
      
      #hist, bins = np.histogram(data, bins=100, range=(min_data,max_data), density=False, weights=1*mask_astropy) # Does not provide correct results!
      hist, bins = np.histogram(np.array(data_hist), bins=50, range=(min_data,max_data), density=False)
      bincentres = [(bins[i]+bins[i+1])/2. for i in range(len(bins)-1)]
      plt.step(bincentres,hist,where='mid',color='b',linestyle='--')
      
      #mode = bincentres[list(hist).index(np.max(hist))] # Rough mode. Better fit with a gaussian (TODO)
    
      plt.axvline(x=mean,ls='-', color='green', label='mean')
      plt.axvline(x=median,ls='-.', color='blue', label='median')
      plt.axvline(x=mode,ls='--', color='red', label='mode')
      plt.legend()
      #plt.show()
      plt.savefig('sky_histogram.png', transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)    
      plt.clf()
      plt.close()
      plt.close('all')
      
      
      # Fit the sky
      if polynomial_degree==0:
        if sky_value=='mode':
            sky = mode
        else:
            sky = median
        
        # Add keywords about sky
        header['Sky_med'] = median
        header['Sky_mean'] = mean
        header['Sky_std'] = std      
        
        # Save the output images
        outHDU = pyfits.PrimaryHDU(data-sky, header=header)
        outHDU.writeto(output_image, overwrite=True)
        
        if output_sky is not None:
            outHDU = pyfits.PrimaryHDU(sky, header=header)
            outHDU.writeto(output_sky, overwrite=True)
        
      else:
        stats1 = [mean, median, std]
        do_sky(input_image, mask_image, stats1, output_image=output_image, sky_image=output_sky, sigma_smooth=sigma_smooth, sampling_factor=sampling_factor, polynomial_degree=polynomial_degree)
      
      if verbosity: print('Done!') 
      return output_image, mean, median, std

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Fitting sky background with a polynomial")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("inputMask", help="Input mask image",type=str) 
    parser.add_argument("--n", help="Optional: Polynomial degree (0 is constant level)", type=int, default=0) 
    parser.add_argument("--output_image", help="Optional: Output image", type=str, default=None) 
    parser.add_argument("--output_sky", help="Optional: Output sky", type=str, default=None)
    parser.add_argument("--hdu", help="Optional: HDU layer", type=int, default=0)
    
    args = parser.parse_args()

    input_image = args.inputImage
    mask_image = args.inputMask
    output_image = args.output_image
    output_sky = args.output_sky
    hdu = args.hdu
    n = args.n
    
    if output_image is None:
        output_image='sky_subtr.fits'


    sky_subtraction(input_image, mask_image, polynomial_degree=n, output_image=output_image, output_sky=output_sky, hdu_inp=hdu)
