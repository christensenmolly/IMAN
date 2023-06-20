# img_interp.py
import os
import sys
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from PIL import Image
from astropy.io import fits as pyfits

from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator
import random
import scipy.ndimage as ndimage
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel


LOCAL_DIR = "/imp/sky_fitting"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))

import rebin_image

def astropy_smoothing(input_image, mask_image, output_image='interpolated.fits', sigma=10):
    hdulist0 = pyfits.open(input_image)
    data0 = hdulist0[0].data 
    ny0,nx0 = np.shape(data0)
    
    rebin_image.downsample(input_image, 5., output_image='tmp_rebin.fits', set_wcs=True, print_mes=True, norm=False)
    rebin_image.downsample(mask_image, 5., output_image='tmp_mask_rebin.fits', set_wcs=True, print_mes=True, norm=False)
    
    
    hdulist = pyfits.open('tmp_rebin.fits')
    data = hdulist[0].data 
    header = hdulist[0].header 
    ny,nx = np.shape(data)

    hdulist1 = pyfits.open('tmp_mask_rebin.fits')
    mask = hdulist1[0].data    
    
    data[mask > 0.] = np.nan
    
    kernel = Gaussian2DKernel(x_stddev=sigma)
    
    astropy_conv = convolve(data, kernel)
    
    for k in range(ny):
        for i in range(nx):
            #if mask[k,i]>0.:
            data[k,i] = astropy_conv[k,i]
    
    
    outHDU = pyfits.PrimaryHDU(data, header=header)
    outHDU.writeto('interp_tmp_rebin.fits', overwrite=True)    
    
    rebin_image.downsample('interp_tmp_rebin.fits', float(nx)/float(nx0), output_image=output_image, set_wcs=True, print_mes=True, norm=False)





def make_interpolated_image2(im, X, Y, nsamples):
    """Make an interpolated image from a random selection of pixels.

    Take nsamples random pixels from im and reconstruct the image using
    scipy.interpolate.griddata.

    """

    ix = np.random.randint(im.shape[1], size=nsamples)
    iy = np.random.randint(im.shape[0], size=nsamples)
    samples = im[iy,ix]
    int_im = griddata((iy, ix), samples, (Y, X), method='cubic')
    return int_im


def make_interpolated_image1(im, mask, X, Y, nsamples):
    """Make an interpolated image from a random selection of pixels.

    Take nsamples random pixels from im and reconstruct the image using
    scipy.interpolate.griddata.

    """
    im[mask>0] = np.nan
    im = np.ma.masked_invalid(im)

    
    ix = []; iy = []
    N = 0
    while N<=nsamples:
        xx = random.randint(0, im.shape[1]-1)
        yy = random.randint(0, im.shape[0]-1)
        if not np.isnan(im[yy,xx]):
            ix.append(xx)
            iy.append(yy)
            N = N + 1
    
    ix = np.array(ix)
    iy = np.array(iy)
    
    samples = im[iy,ix]
    int_im = griddata((iy, ix), samples, (Y, X), method='cubic')
    return int_im

def make_interpolated_image(im, mask, X, Y):
    """Make an interpolated image from a random selection of pixels.

    Take nsamples random pixels from im and reconstruct the image using
    scipy.interpolate.griddata.

    """
    im[mask>0] = np.nan
    im = np.ma.masked_invalid(im)
    
    iy = Y[~im.mask]
    ix = X[~im.mask]
    
    int_im = griddata((iy, ix), im[~im.mask], (Y, X), method='cubic')
    return int_im


def main1():
    # Read in image and convert to greyscale array object
    img_name = sys.argv[1]
    #im = Image.open(img_name)
    #im = np.array(im.convert('L'))

    hdulist = pyfits.open(img_name)
    im = hdulist[0].data   

    # A meshgrid of pixel coordinates
    nx, ny = im.shape[1], im.shape[0]
    X, Y = np.meshgrid(np.arange(0, nx, 1), np.arange(0, ny, 1))

    # Create a figure of nrows x ncols subplots, and orient it appropriately
    # for the aspect ratio of the image.
    nrows, ncols = 2, 2
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(6,4), dpi=100)
    if nx < ny:
        w, h = fig.get_figwidth(), fig.get_figheight()
        fig.set_figwidth(h), fig.set_figheight(w)

    # Convert an integer i to coordinates in the ax array
    get_indices = lambda i: (i // nrows, i % ncols)

    # Sample 100, 1,000, 10,000 and 100,000 points and plot the interpolated
    # images in the figure
    for i in range(4):
        nsamples = 10**(i+2)
        axes = ax[get_indices(i)]
        axes.imshow(make_interpolated_image1(im, X, Y, nsamples),
                            cmap=plt.get_cmap('Greys_r'))
        axes.set_xticks([])
        axes.set_yticks([])
        axes.set_title('nsamples = {0:d}'.format(nsamples))
    filestem = os.path.splitext(os.path.basename(img_name))[0]
    plt.savefig('{0:s}_interp.png'.format(filestem), dpi=100)
    
    
def main(input_image, mask_image=None, nsamples=100000, output_image='interp.fits',sigma_smooth=50.):
    # Read in image and convert to greyscale array object


    hdulist = pyfits.open(input_image)
    data = hdulist[0].data   
    nx, ny = data.shape[1], data.shape[0]
    
    if mask_image is not None:
        hdulist1 = pyfits.open(mask_image)
        mask = hdulist1[0].data
    else:
        mask = np.zeros((ny,nx))

    #data = ndimage.uniform_filter(data, size=sigma_smooth)#ndimage.gaussian_filter(data, sigma=sigma_smooth, order=0)

    # A meshgrid of pixel coordinates
    
    X, Y = np.meshgrid(np.arange(0, nx, 1), np.arange(0, ny, 1))

    new_image = make_interpolated_image(data, mask, X, Y)
    #new_image = make_interpolated_image1(data, mask, X, Y, nsamples)
    
    
    outHDU = pyfits.PrimaryHDU(new_image)#, header=header)
    outHDU.writeto(output_image, clobber=True)       



def main2(input_image, mask_image=None, output_image='interp.fits'):
    # Too many points: WARNING DOES NOT WORK!!!!!!!!!!!!!!
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data   
    
    nx, ny = data.shape[1], data.shape[0]

    if mask_image is not None:
        hdulist1 = pyfits.open(mask_image)
        mask = hdulist1[0].data
    else:
        mask = np.zeros((ny,nx))
    
    z = np.ma.array(data, mask=mask)

    x, y = np.mgrid[0:z.shape[0], 0:z.shape[1]]
    x1 = x[~z.mask]
    y1 = y[~z.mask]
    z1 = z[~z.mask]
    new_data = interpolate.interp2d(x1, y1, z1)(np.arange(z.shape[0]), np.arange(z.shape[1]))    

    outHDU = pyfits.PrimaryHDU(new_data)#, header=header)
    outHDU.writeto(output_image, clobber=True)   
    


def main3(input_image, mask_image=None, output_image='interp.fits'):
    # Too many points: WARNING DOES NOT WORK!!!!!!!!!!!!!!
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data   
    
    nx, ny = data.shape[1], data.shape[0]

    if mask_image is not None:
        hdulist1 = pyfits.open(mask_image)
        mask = hdulist1[0].data
    else:
        mask = np.zeros((ny,nx))
    
    z_masked = np.ma.array(data, mask=mask)
    
    z = z_masked.filled(np.nan)
    
    zinterp = RegularGridInterpolator((np.arange(0, nx, 1), np.arange(0, ny, 1)), z.T)

    X, Y = np.meshgrid(np.arange(0, nx, 1), np.arange(0, ny, 1))
    newpoints = np.array((X, Y)).T
    
    z2 = zinterp(newpoints)
    z2_masked = np.ma.array(z2, mask=np.isnan(z2))
  

    outHDU = pyfits.PrimaryHDU(z2)#, header=header)
    outHDU.writeto(output_image, clobber=True)   



def main4(input_image, mask_image=None, output_image='interp.fits'):
    # Too many points: WARNING DOES NOT WORK!!!!!!!!!!!!!!
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data   
    
    nx, ny = data.shape[1], data.shape[0]

    if mask_image is not None:
        hdulist1 = pyfits.open(mask_image)
        mask = hdulist1[0].data
    else:
        mask = np.zeros((ny,nx))
    
    z_masked = np.ma.array(data, mask=mask)
    
    z = z_masked.filled(np.nan)
    
    zinterp = RegularGridInterpolator((np.arange(0, nx, 1), np.arange(0, ny, 1)), z.T)
    minterp = RegularGridInterpolator((np.arange(0, nx, 1), np.arange(0, ny, 1)), (mask+0.).T)

    X, Y = np.meshgrid(np.arange(0, nx, 1), np.arange(0, ny, 1))
    newpoints = np.array((X, Y)).T

    # actual interpolation
    z2 = zinterp(newpoints)
    mask2 = minterp(newpoints) > 0  # apply threshold, e.g. 0.5 is considered contaminated and will be removed.
    z2[mask2] = np.nan  # fill with nans or whatever missing data flag
  

    outHDU = pyfits.PrimaryHDU(z2)#, header=header)
    outHDU.writeto(output_image, clobber=True)        
    
#main4('35840_r_trim.fits',mask_image='general_mask.fits')
main('combined.fits',mask_image='mask.fits')
#astropy_smoothing('28498_r_trim.fits', 'general_mask_r.fits', output_image='interpolated.fits', sigma=10)