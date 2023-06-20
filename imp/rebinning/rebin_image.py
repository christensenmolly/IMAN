#!/usr/bin/python
# DESCRIPTION:
# Script to rebin image to a reference image.
# MINIMAL USAGE: python rebin_image.py [input_image] [reference_image]

# Import standard modules
import os
import os.path
import math
import numpy as np
import shutil
import astropy.io.fits as pyfits
from astropy import wcs
import glob
import argparse
from astropy import log
import warnings
from scipy import ndimage
warnings.filterwarnings("ignore")
import sys

LOCAL_DIR = "/imp/rebinning"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))

from hcongrid import hcongrid


# Disable astropy logging except for warnings and errors
log.setLevel("WARNING")


def zoomdown(in_array, zoom_ratio, nan_lim=0.5, fix_nan_flux=False, pad=True):
    # Determine the shape of the output array
    in_y_size, in_x_size = in_array.shape
    out_x_size = int(math.ceil(in_x_size * zoom_ratio))
    out_y_size = int(math.ceil(in_y_size * zoom_ratio))

    # Pad input array
    pad_width = ((0, math.ceil(1/zoom_ratio)),
                 (0, math.ceil(1/zoom_ratio)))
    in_array = np.pad(in_array, pad_width=pad_width, mode="edge")

    out_array = np.zeros((out_y_size, out_x_size))
    # Iterate over indices of the resulting (i.e. smaller) array
    for x in np.arange(out_x_size, dtype=float):
        print(x)
        for y in np.arange(out_y_size, dtype=float):
            # Compute the coordinates of this out array pixel inside of the
            # input array
            x_in_lb = x / zoom_ratio  # Left bottom corner
            y_in_lb = y / zoom_ratio  #
            x_in_rb = (x+1) / zoom_ratio  # Right bottom corner
            y_in_rb = y_in_lb             #
            x_in_lt = x_in_lb             # Left top corner
            y_in_lt = (y+1) / zoom_ratio  #
            x_in_rt = x_in_rb  # Right top corner
            y_in_rt = y_in_lt  #
            # Calculate the total intensiy in this square
            intensities = []
            weights = []
            # Add left bottom corner
            weights.append((1 - x_in_lb % 1) * (1 - y_in_lb % 1))
            intensities.append(in_array[int(y_in_lb), int(x_in_lb)])
            # Add right bottom corner
            weights.append((x_in_rb % 1) * (1 - y_in_rb % 1))
            intensities.append(in_array[int(y_in_rb), int(x_in_rb)])
            # Add left top corner
            weights.append((1 - x_in_lt % 1) * (y_in_lt % 1))
            intensities.append(in_array[int(y_in_lt), int(x_in_lt)])
            # Add right top corner
            weights.append((x_in_rt % 1) * (y_in_rt % 1))
            intensities.append(in_array[int(y_in_rt), int(x_in_rt)])
            # Add left border (without corners)
            left_border_pixel_area = 1 - x_in_lb % 1
            if left_border_pixel_area != 0:
                weights.extend([left_border_pixel_area] * (int(y_in_lt) - int(y_in_lb+1)))
                intensities.extend(in_array[int(y_in_lb+1): int(y_in_lt), int(x_in_lb)])
            # Add top border (without corners)
            top_border_pixel_area = y_in_lt % 1
            if top_border_pixel_area != 0:
                weights.extend([top_border_pixel_area] * (int(x_in_rt) - int(x_in_lt+1)))
                intensities.extend(in_array[int(y_in_lt), int(x_in_lt+1): int(x_in_rt)])
            # Add right border (without corners)
            right_border_pixel_area = x_in_rt % 1
            if right_border_pixel_area != 0:
                weights.extend([right_border_pixel_area] * (int(y_in_rt) - int(y_in_rb+1)))
                intensities.extend(in_array[int(y_in_rb+1): int(y_in_rt), int(x_in_rb)])
            # Add bottom border (without cirners)
            bottom_border_pixel_area = 1 - y_in_lb % 1
            if bottom_border_pixel_area != 0:
                weights.extend([bottom_border_pixel_area] * (int(x_in_rb) - int(x_in_lb+1)))
                intensities.extend(in_array[int(y_in_lb), int(x_in_lb+1): int(x_in_rb)])
            # Add central region (the whole pixels)
            central = in_array[int(y_in_lb+1): int(y_in_rt), int(x_in_lb+1): int(x_in_rt)].ravel()
            weights.extend(np.ones_like(central))
            intensities.extend(central)
            # Write the intensity into the output array
            weights = np.array(weights)
            intensities = np.array(intensities)
            nan_fraction = np.sum(weights[np.isnan(intensities)]) / np.sum(weights)
            if nan_fraction > nan_lim:
                # The total area of nan pixels is greater than the given limit,
                # so we write nan in the output array
                out_array[int(y), int(x)] = np.nan
            else:
                weights = np.array(weights)
                intensities = np.array(intensities)
                out_array[int(y), int(x)] = np.nansum(intensities * weights)
                if (fix_nan_flux is True) and (nan_fraction > 0.0):
                    out_array[int(y), int(x)] /= (1-nan_fraction)
    if pad is False:
        if in_x_size * zoom_ratio > int(in_x_size*zoom_ratio):
            out_array = out_array[:, :-1]
        if in_y_size * zoom_ratio > int(in_y_size*zoom_ratio):
            out_array = out_array[:-1, :]
    return out_array



def zoom_fits(fitsfile, scalefactor, preserve_bad_pixels=True, set_wcs=True, no_interp=False, **kwargs):
    """
    Zoom in on a FITS image by interpolating using `~scipy.ndimage.interpolation.zoom`

    Parameters
    ----------
    fitsfile : str
        FITS file name
    scalefactor : float
        Zoom factor along all axes
    preserve_bad_pixels : bool
        Try to set NAN pixels to NAN in the zoomed image.  Otherwise, bad
        pixels will be set to zero
    """

    arr = pyfits.getdata(fitsfile)
    h = pyfits.getheader(fitsfile)
    
    if set_wcs==True:
      try:
        h['CRPIX1'] = (h['CRPIX1']-1)*scalefactor + scalefactor/2. + 0.5
        h['CRPIX2'] = (h['CRPIX2']-1)*scalefactor + scalefactor/2. + 0.5
        if 'CD1_1' in h:
            for ii in (1, 2):
                for jj in (1, 2):
                    k = "CD%i_%i" % (ii, jj)
                    if k in h:  # allow for CD1_1 but not CD1_2
                        h[k] = h[k]/scalefactor
        elif 'CDELT1' in h:
            h['CDELT1'] = h['CDELT1']/scalefactor
            h['CDELT2'] = h['CDELT2']/scalefactor
      except:
        print('ERROR: WCS was not added!')
    bad_pixels = np.isnan(arr) + np.isinf(arr)

    arr[bad_pixels] = 0
    
    #'''
    if no_interp==False:
        upscaled = ndimage.zoom(arr, scalefactor, **kwargs)
    else:
        upscaled = ndimage.zoom(arr, scalefactor, mode='constant', order=0)
        
    if preserve_bad_pixels:
        bp_up = ndimage.zoom(bad_pixels, scalefactor,
                                   mode='constant', cval=np.nan, order=0)
        try:
            upscaled[bp_up] = np.nan #######???
        except:
            z=1
    upscaled = upscaled / (scalefactor**2)
    #'''
    #upscaled = zoomdown(arr, scalefactor, nan_lim=0.5, fix_nan_flux=False, pad=True) # Very slow!!!
    
    
    up_hdu = pyfits.PrimaryHDU(data=upscaled, header=h)
    return up_hdu

## This function downsamples the frame by a specified zooming factor
def downsample(input_image, factor, output_image=None, set_wcs=True, print_mes=True, norm=False, no_interp=False):
        if print_mes==True:
            print('Rebining the input image by a factor of: %s' % str(factor)) 
        # Open the input file and read it
        hdulist_in = pyfits.open(input_image)
        data_in = hdulist_in[0].data
        header_in = hdulist_in[0].header

        # Use zoom_fits function from FITS_tools:  http://fits-tools.readthedocs.org/en/latest/_modules/FITS_tools/hcongrid.html
        hdu_out = zoom_fits(input_image, 1.0/factor, set_wcs=set_wcs, no_interp=no_interp)
        
        data = hdu_out.data
        if norm:
            hdu_out.data = data/np.sum(data)
        else:
            hdu_out.data = data
        
        # Save the HDU
        if output_image is None:
            output_image = input_image.split('.fits')[0] + '_sampl.fits'
        hdu_out.writeto(output_image, overwrite=True)
        if print_mes==True:
            print('Done!')
        return output_image


def rebin(reference_image, input_image, output_image=None, hdu_ref=0, hdu_inp=0, preserve_bad_pixels=True, no_interp=False):
        '''
        This function rebins the input image
        based on a certain reference FITS file.
        '''
        
        # Get pixels scale for each image
        pix2sec_ref, note_ref = resolution(reference_image)
        pix2sec_inp, note_inp = resolution(input_image)        
        pix2sec_ref = float(pix2sec_ref)
        pix2sec_inp = float(pix2sec_inp)


        # Open the HDU list for the reference FITS file
        hdulist = pyfits.open(reference_image)

        # Get the primary image
        hdu = hdulist[hdu_ref]

        referenceheader = hdu.header

        referenceheader["NAXIS"] = 2
        referenceheader.pop("NAXIS3", None)

        # Obtain the coordinate system from the reference header
        coordinates = wcs.WCS(referenceheader)

        hdulist_in = pyfits.open(input_image)
        data_in = hdulist_in[hdu_inp].data
        header_in = hdulist_in[hdu_inp].header

        if 'NAXIS3' in header_in:
          # For CUBE images
          del header_in['NAXIS3']
          del header_in['NAXIS']
          hdu_tmp = pyfits.PrimaryHDU(data_in[0], header_in)
          hdu_tmp.writeto('tmp.fits',overwrite=True)
          hdulist_in = pyfits.open('tmp.fits')
          data_in = hdulist_in[0].data
          header_in = hdulist_in[0].header
          os.remove('tmp.fits')

        # Do the rebinning based on the header of the reference image
        data_out = hcongrid(data_in, header_in, referenceheader, preserve_bad_pixels=preserve_bad_pixels, no_interp=no_interp)


        # Use the coordinate system of the reference image to create the new header for this image
        header_out = coordinates.to_header()
        
        # Add keywords
        res_factor = pix2sec_ref/pix2sec_inp
        header_out.append(('PIXSCALE_OLD', pix2sec_inp),end=True)
        header_out.append(('PIXSCALE_NEW', pix2sec_ref),end=True)
        header_out.append(('FACTOR_M0', res_factor),end=True) 
    
        hdu_out = pyfits.PrimaryHDU(data_out, header_out)
        if output_image is None:
            output_image = input_image.split('.fits')[0] + '_rebin.fits'

        hdu_out.writeto(output_image, overwrite=True)



def resolution(input_image):
    hdulist = pyfits.open(input_image)#, ignore_missing_end=True)
    header = hdulist[0].header

    pixelscale = 1.
    note = ' '

    if 'PIXSCALE' in header:
        pixelscale = header['PIXSCALE']

    elif 'SECPIX' in header:
        pixelscale = header['SECPIX']

    elif 'PFOV' in header:
        pixelscale = header['PFOV']

    elif 'CD1_1' in header and 'CD1_2' in header:
        pixelscale = math.sqrt(header['CD1_1']**2 + header['CD1_2']**2 ) * 3600.0
        

    elif 'CD1_1' in header:
        pixelscale = abs(header['CD1_1']) * 3600.0

    elif 'CDELT1' in header:
        pixelscale = abs(header['CDELT1']) * 3600.0
        if 'PC1_1' in header:
            pixelscale = abs(header['PC1_1']) * 3600.0

    else:
        print("Could not determine the pixel scale from the image header. Set to 1 pix/arcsec.")
        note = '*'

    if 'CDELT1' in header:
        pixelscale = abs(header['CDELT1']) * 3600.0
        if 'PC1_1' in header:
            pixelscale = abs(header['PC1_1']) * 3600.0

    #print(pixelscale)
    #exit()
    # Return the pixel scale (in arcseconds)
    hdulist.close()
    return str(pixelscale),note


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Rebinning")
    parser.add_argument("input_image", help="Input image")    
    parser.add_argument("reference_image", help="Input reference image")
    parser.add_argument("--o", help="Optional: Output image",type=str,default=None)
    parser.add_argument("--nointerp", action="store_true", default=False,
                        help="Do not interpolate the rebinned image (important for masking)")    
    args = parser.parse_args()

    reference_image = args.reference_image
    input_image = args.input_image
    output_image = args.o
    no_interp = args.nointerp

    
    rebin(reference_image, input_image, output_image, hdu_ref=0, hdu_inp=0, preserve_bad_pixels=True, no_interp=no_interp)

