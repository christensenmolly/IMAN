#!/bin/python3.7

"""

Use spline interpolation to on grid of x,y,z value where z is either xdiff or ydiff for use as imagemagick 2D displacement maps

"""

import numpy as np
from scipy.interpolate import SmoothBivariateSpline
from skimage import io
from astropy.io import fits as pyfits

import determine_sky


def define_arrays(data, mask):
    ny,nx = np.shape(data)
    ax = []; ay = []; az = []
    for k in range(ny):
        for i in range(nx):
            if mask[k,i]==0:
                ax.append(i)
                ay.append(k)
                az.append(data[k,i])
    return np.array(ax),np.array(ay),np.array(az)


def main(input_image, mask_image):
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data 
    ny,nx = np.shape(data)

    header = hdulist[0].header

    hdulist1 = pyfits.open(mask_image)
    mask = hdulist1[0].data 

    ax, ay, az = determine_sky.split_xyz(data, mask=mask, arrays=True)
    #ax, ay, az = define_arrays(data, mask)
    
    
    #print(x_values, y_values, z_values)
    #exit()


    # convert python lists to numpy arrays
    #ax = np.asarray(x_values)
    #ay = np.asarray(y_values)
    #az = np.asarray(z_values)
    ax = np.array(ax)[0:10000]
    ay = np.array(ay)[0:10000]
    az = np.array(az)[0:10000]
    #print(len(ax), len(ay), len(az))
    #exit()
    #print(az)
    #exit()

    # define bbox of interpolated data
    # bbox=[minx, maxx, miny, maxy]
    bbox=[0, 129, 0, 129]

    # convert bbox to numpy array
    abbox = np.asarray(bbox)

    # do interpolations
    xd = SmoothBivariateSpline(ax, ay, az, w=None, kx=3, ky=3)


    # define integer grid onto which to interpolate
    grid_x=np.arange(0,100)#np.linspace(0, 129, num=129)
    grid_y=np.arange(0,100)#np.linspace(0, 129, num=129)


    # evaluate at grid points
    xdisplace = xd.__call__(grid_x, grid_y, grid=True)
    #ydisplace = yd.__call__(grid_x, grid_y, grid=True)

    outHDU = pyfits.PrimaryHDU(xdisplace)
    outHDU.writeto('test_out.fits', clobber=True)
    
    print(xdisplace)
    exit()

    # save output using skimage
    io.imsave("xdimgs.png", xdisplace.astype('uint8'))
    #io.imsave("ydimgs.png", ydisplace.astype('uint8'))

    # view output using skimage
    io.imshow(xdisplace.astype('uint8')) 
    io.show()
    #io.imshow(ydisplace.astype('uint8')) 
    #io.show()
    
main('28498_r_trim.fits', 'general_mask_r.fits')