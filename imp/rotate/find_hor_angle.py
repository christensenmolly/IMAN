
import os
import numpy as np
from scipy import ndimage
from scipy import misc
import argparse
import astropy.io.fits as pyfits
import math
import numpy as np
import shutil
import warnings
import matplotlib.pyplot as plt
try:
  import astropy.io.fits as pyfits
  from astropy import wcs
except:
  warnings.warn("Astropy is not installed! No WCS has been added to the header!")

warnings.filterwarnings("ignore")  
# This script rotates a FITS image around the center of the frame over an arbitrary angle.

from sklearn.linear_model import LinearRegression
from scipy.ndimage.filters import gaussian_filter1d

import rotate_image


def count_vert_pixels(input_image, mask_image, m0, pix2sec, SBbright, SBfaint):
    SBbright = 10**(0.4*(m0-SBbright)) * (pix2sec**2)
    SBfaint = 10**(0.4*(m0-SBfaint)) * (pix2sec**2)
    print(SBfaint)
    print(SBbright)
    
    # Load in image
    hdu = pyfits.open(input_image)
    data = hdu[0].data
    ny,nx = np.shape(data)
    yc = ny/2.
    
    if mask_image is not None:
        hdu_mask = pyfits.open(mask_image)
        mask = hdu_mask[0].data
    else:
        mask = np.copy((ny,nx))

    y = np.arange(ny) - yc
    N_y = [0.] * ny
    for k in range(ny):
        for i in range(nx):
            if mask[k,i]==0 and data[k,i]>SBfaint and data[k,i]<SBbright:
                N_y[k] = N_y[k] + 1
    
    #plt.plot(y, np.log10(N_y), '-')
    y = gaussian_filter1d(y, sigma=2)
    #plt.plot(y, N_y, '-')
    #plt.plot(N_y,np.array(y)*np.array(N_y))
    
    #plt.show()
    non_zero_inds = [i for i, e in enumerate(N_y) if e != 0]
    non_zero_inds = np.array(non_zero_inds)
    
    
    #median_y_distance = np.median(np.array(N_y)[non_zero_inds]*np.fabs(y)[non_zero_inds])
    return np.max(N_y)


def find_regression(input_image, mask_image, m0, pix2sec, SBbright, SBfaint):
    SBbright = 10**(0.4*(m0-SBbright)) * (pix2sec**2)
    SBfaint = 10**(0.4*(m0-SBfaint)) * (pix2sec**2)
    
    # Load in image
    hdu = pyfits.open(input_image)
    data = hdu[0].data
    ny,nx = np.shape(data)
    yc = ny/2.
    
    if mask_image is not None:
        hdu_mask = pyfits.open(mask_image)
        mask = hdu_mask[0].data
    else:
        mask = np.copy((ny,nx))    
    
    x = []; y = []
    for k in range(ny):
        for i in range(nx):
            if mask[k,i]==0 and data[k,i]>SBfaint and data[k,i]<SBbright:
                x.append(i)
                y.append(k)
    x = np.array(x).reshape(-1, 1)
    y = np.array(y).reshape(-1, 1)

    linear_regressor = LinearRegression()  # create object for the class
    linear_regressor.fit(x, y)  # perform linear regression
    y_pred = linear_regressor.predict(x)  # make predictions
    params = linear_regressor.get_params()
    print(params)
    exit()
    plt.plot(x,y)
    plt.plot(x, y_pred)
    plt.show()
    exit()


def do_each_step(input_image, mask_image, angle, m0, pix2sec, SBbright, SBfaint):
    if mask_image is not None:
        output_mask_image = 'tmp_mask_%s.fits' % (angle)
        rotate_image.rotate_image(mask_image, angle, output_image=output_mask_image, hdu_inp=0, cval='nan')    
    else:
        output_mask_image = None
        
    output_image = 'tmp_%s.fits' % (angle)
    rotate_image.rotate_image(input_image, angle, output_image=output_image, hdu_inp=0, cval='nan')    
    median_y_distance = count_vert_pixels(output_image, output_mask_image, m0, pix2sec, SBbright, SBfaint)
    
    #find_regression(output_image, output_mask_image, m0, pix2sec, SBbright, SBfaint)
    
    
    os.remove(output_image)
    if mask_image is not None:
        os.remove(output_mask_image)
    
    return median_y_distance


def main(input_image, m0, pix2sec, mask_image=None, angle_limits=None, angle_step=1., SBbright=None, SBfaint=26.):
        if SBbright is None:
            SBbright = -99.
    
        #median_y_distance_0 = count_vert_pixels(input_image, mask_image, m0, pix2sec, SBbright, SBfaint)
        #angle_0 = 0.

        
        if angle_limits is None:
            angle_limits = [-90.,90.]
        
        angles = np.arange(angle_limits[0],angle_limits[1],angle_step)
        
        median_y_distances = []
        for angle in angles:
            median_y_distance = do_each_step(input_image, mask_image, angle, m0, pix2sec, SBbright, SBfaint)
            median_y_distances.append(median_y_distance)
            print(angle, median_y_distance)
        
        plt.plot(angles, median_y_distances)
        plt.show()
        
        
        angles = list(angles)
        median_y_distances = median_y_distances
        #plt.plot(angles, median_y_distances, 'o')
        #plt.show()
        best_angle = angles[median_y_distances.index(max(median_y_distances))]
        print('Inclination angle is: %.2f deg' % (best_angle))


main('galaxy_rot_crop.fits', 29.61, 0.833, mask_image='mask_rot_crop.fits', angle_limits=[-10.,11], angle_step=1., SBbright=None, SBfaint=24.)