#!/usr/bin/python
# Written by Sviatoslav Borisov, MSU
from mpl_toolkits.mplot3d import *
import matplotlib.pyplot as plt
from astropy.io import fits
import scipy as sci
import numpy as np
import itertools
import warnings
import timeit
import math
import copy
import os
import sys
import argparse

warnings.filterwarnings("ignore")



def read_image(image_file, extension):
    print("READING THE IMAGE")
    hdulist = fits.open(image_file)
    header = hdulist[extension].header
    image = hdulist[extension].data
    size = image.shape
    image_orig = copy.copy(image)
    return image, image_orig, size, header

def mask_image(image, square_nan_size, n_iter, n_sigma, verbosity):
    if verbosity: print("MASKING OUTLIERS IN THE IMAGE")
    for i in range(n_iter):
                if i == 1:
                    image_temp = image
                if i >= 1:
                    image = image_temp
                    image_temp = image
                    med_temp = np.nanmedian(image_temp)
                    sig_temp = np.nanstd(image_temp)
                    ind_nan = np.where(image_temp >= med_temp+n_sigma*sig_temp)
                    num_nan = (ind_nan[1]).shape[0]
                    image_temp[ind_nan] = np.NaN
                    for j in range(num_nan):
                                x = (ind_nan[0])[j]
                                y = (ind_nan[1])[j]
                                image[x-square_nan_size:x+square_nan_size, y-square_nan_size:y+square_nan_size] = np.NaN
    return image

    
def mask_plot(image_orig, image):
    llim = np.median(image_orig)
    ulim = np.mean(image_orig[(image_orig>np.mean(image_orig)+1*np.std(image_orig)) &
        (image_orig<np.mean(image_orig)+3*np.std(image_orig))])
    fig, ax = plt.subplots(1, 2)
    ax[0].imshow(image_orig, vmin=llim, vmax=ulim)
    ax[1].imshow(image, vmin=llim, vmax=ulim)
    plt.show(block=False)

def image_rebin(square_size, size, nan_thr, image, verbosity):
    if verbosity: print("REBINNING THE IMAGE")
    height = int(np.floor(size[0]/square_size))
    width = int(np.floor(size[1]/square_size))
    image_trunc = image[0:height*square_size, 0:width*square_size]
    image_rebin = image_trunc.reshape((height, image_trunc.shape[0]//height, width, -1))
    image_rebin_nanmean = np.zeros(shape=(height*width, 1))
    x_center = np.zeros(shape=(height*width, 1))
    y_center = np.zeros(shape=(height*width, 1))
    num = 0
    for i in range(image_rebin.shape[0]):
        for j in range(image_rebin.shape[2]):
            x_center[num, 0] = j*square_size+np.floor(square_size/2)+1
            y_center[num, 0] = i*square_size+np.floor(square_size/2)+1
            square = image_rebin[i, 0:square_size, j, 0:square_size]
            n_nan = (np.argwhere(np.isnan(square))).shape[0]
            if (float(n_nan)/square_size**2 <= nan_thr):
                image_rebin_nanmean[num, 0] = np.nanmean(square)
            else:
                image_rebin_nanmean[num, 0] = np.NaN
            num += 1
    ind_good = np.where(~np.isnan(image_rebin_nanmean))
    num_good = ind_good[0].shape[0]
    x = x_center[ind_good]
    y = y_center[ind_good]
    z = image_rebin_nanmean[ind_good]
    data = np.zeros(shape=(num_good, 3))
    for i in range(num_good):
        data[i, :] = [x[i], y[i], z[i]]
    return data

def polyfit2d(x, y, z, pol_degree, verbosity):
    if verbosity: print("SURFACE FITTING")
    ncols = (pol_degree + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(pol_degree+1), range(pol_degree+1))
    for k, (i, j) in enumerate(ij):
        G[:, k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

def polyval2d(x, y, m, verbosity):
    if verbosity: print("GETTING A SURFACE")
    pol_degree = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(pol_degree+1), range(pol_degree+1))
    z = np.zeros_like(x)
    for a, (i, j) in zip(m, ij):
        z += a * x**i * y**j
    return z

def get_surf(size, nx, ny, poly_coeff, verbosity):
    X, Y = np.meshgrid(np.linspace(0, size[1], nx),
        np.linspace(0, size[0], ny))
    Z = polyval2d(X, Y, poly_coeff, verbosity)
    return X, Y, Z

def surf_plot(data, size, nx_plot, ny_plot, poly_coeff, verbosity):
    X, Y, Z = get_surf(size, nx_plot, ny_plot, poly_coeff, verbosity)
    llim = np.nanmedian(data[:, 2])-5*np.nanstd(data[:, 2])
    ulim = np.nanmedian(data[:, 2])+5*np.nanstd(data[:, 2])
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.9)
    ax.scatter(data[:, 0], data[:, 1], data[:, 2], c='r', s=20)
    ax.set_zlim3d(llim, ulim)
    plt.xlabel('X')
    plt.ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

def write_fits(filename, image2write, header, verbosity):
    if verbosity: print("WRITTING THE DATA TO A NEW FILE")
    if os.path.exists(filename):
        os.remove(filename)
    fits.writeto(filename, image2write, header)
    if verbosity: print("NEW FILE '"+filename+"'")





def main(input_image, output_file=None, mask_image=None, sky_image=None, nan_thr = 0.3, pol_degree = 3, extension=0, plot_mask = 0, plot_surf = 0, verbosity=True):
    start = timeit.default_timer()
    
    if verbosity:
        print('')
        print("BACKGROUND SUBTRACTION")
        print("----------------------")
    
    
    if output_file is None:
        output_file = input_image.split('.fits')[0] + '_bg.fits'

    nx_plot, ny_plot = 50, 50   #number of bins of surface (used only for display)


    #########PROCEDURE EXECUTION##########
    ######################################
    #number in a comment indicate the number of a block of variables used in this function(s) (see above)
    image, image_orig, size, header = read_image(input_image, extension)  #1
    square_nan_size = int(np.ceil(np.sqrt(size[0]*size[1])/150.))   
    square_size = int(square_nan_size*1.5)
    nx, ny = size[1], size[0]
    
    if mask_image is not None:
        hdulist = fits.open(mask_image)
        mask = hdulist[0].data        
        image[mask!=0]=float('nan')
    
    
    if plot_mask:
        mask_plot(image_orig, image)
    #print(square_size)
    data = image_rebin(square_size, size, nan_thr, image, verbosity)   #3

    #outHDU = fits.PrimaryHDU(data)
    #outHDU.writeto('rebin_out.fits', overwrite=True)    
    #exit()
    
    poly_coeff = polyfit2d(data[:, 0], data[:, 1], data[:, 2], pol_degree, verbosity)

    if plot_surf:
        surf_plot(data, size, nx_plot, ny_plot, poly_coeff, verbosity) #4
    X, Y, Z = get_surf(size, nx, ny, poly_coeff, verbosity)

    image_orig -= Z

    write_fits(output_file, image_orig, header, verbosity)    #5
    
    if sky_image is not None:
        outHDU = fits.PrimaryHDU(Z)
        outHDU.writeto(sky_image, overwrite=True)    
    
    if verbosity:
        print("===========================")
        print("BACKGROUND SUBTRACTION DONE")
        print("TIME ELAPSED "+str(int((timeit.default_timer()-start)*100)/100.)+" SECONDS")
    


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

    
    main(input_image, output_file=output_image, mask_image=mask_image, sky_image=output_sky, nan_thr = 0.3, pol_degree = n, extension=hdu, plot_mask = 0, plot_surf = 0)




#main('28498_r_trim.fits',mask_image='general_mask_r.fits')

#check_background.main('28498_r_trim_bg.fits', 'general_mask_r.fits', sigma_smooth=10)
#check_background.main('sky_subtr.fits', 'general_mask_r.fits', sigma_smooth=10)