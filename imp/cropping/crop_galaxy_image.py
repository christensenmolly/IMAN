#!/usr/bin/python

import math
import numpy as np
from astropy.io import fits as pyfits
import argparse
import os
import subprocess


def read_region(reg_file, coord_format='image'):
    f = open(reg_file, "r")
    for line in f:
        if "ellipse" in line:
            params = line.split(",")
            cen = [float(params[0].split('(')[1]), float(params[1])]
            if coord_format == 'image':
                ellA = float(params[2])
                ellB = float(params[3])
            else:
                ellA = float(params[2].split('\"')[0])
                ellB = float(params[3].split('\"')[0])

            ellPA = float(params[4].split(')')[0])
            if ellA < ellB:
                ellA, ellB = ellB, ellA
                ellPA += 90
            break
    f.close()
    return cen[0], cen[1], ellA, ellB, ellPA


def ds9_to_np(x):
    '''
    Function to convert ds9 pixel coordinates to numpy format, i.e
    0.5 pix is 0, 1.5 pix is 1
    '''
    return int(round(x)) - 1


def np_to_ds9(x):
    '''
    Function to convert numpy coordinates to ds9 format, i.e
    0 is between 0.5 and 1.5 (not incl), 1 is between 1.5 (incl) and
    2.5 (not incl)
    '''
    return x + 1.


def ellipse_borders(xy, width, height, angle, nx=None, ny=None, square=False):
    '''
    This function determines coordinates for the cropping box
    (such that the provided ellipse is inscribed in this box).

    Pixel coordinates are in the ds9 format.
    xy = [xc,yc] - coordinates of the ellipse (pix),
    width - semi-major axis of the ellipse (pix),
    height - semi-minor axis of the ellipse (pix),
    angle - position angle of the ellipse (in deg),
    the rotation goes in a counterclockwise direction, starting
    from the positive x-axis (i.e. 0 deg is right, 90 deg is up).
    '''

    # Convert pixels to numpy
    xy = [ds9_to_np(xy[0]), ds9_to_np(xy[1])]

    width = width
    height = height
    angle = np.radians(angle)

    X = math.sqrt((width * math.cos(angle)) ** 2 + (height * math.sin(angle)) ** 2)
    Y = math.sqrt((width * math.sin(angle)) ** 2 + (height * math.cos(angle)) ** 2)
    if square == True:
        # In this case the output image will be have a square frame
        X = max([X, Y])
        Y = X

    x_max = xy[0] + int(math.ceil(X))
    x_min = xy[0] - int(math.ceil(X))
    y_max = xy[1] + int(math.ceil(Y))
    y_min = xy[1] - int(math.ceil(Y))

    # Check that we are inside the image, and if not, change the found coordinates of the rectangle
    if nx is not None and ny is not None:
        deltax = [];
        deltay = []
        if x_min < 0:
            deltax.append(abs(x_min))
        if x_max >= nx:
            deltax.append(abs(x_max - nx))

        if y_min < 0:
            deltay.append(abs(y_min))
        if y_max >= ny:
            deltay.append(abs(y_max - ny))

        if len(deltax) != 0:
            x_min = x_min + max(deltax)
            x_max = x_max - max(deltax)

        if len(deltay) != 0:
            y_min = y_min + max(deltay)
            y_max = y_max - max(deltay)

        if len(deltax) != 0 or len(deltay) != 0:
            print('Warning: The galaxy ellipse is larger than the image. It will be reduced to fit the image!')

    return [[x_min, y_min], [x_max, y_max]]


def crop_manual(input_image: str, output_image: str, offset_size=1., region_file=None, hdu=0, method='Montage'):
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    nx, ny = data.shape[1], data.shape[0]

    if region_file is not None:
        xc, yc, sma, smb, PA = read_region(region_file)
        [[x_l, y_l], [x_r, y_r]] = ellipse_borders([xc, yc], 2. * offset_size * sma, 2. * offset_size * smb, PA, nx=nx,
                                                   ny=ny, square=False)
        ff = open("crop.reg", "w")
        ff.write('image\n')
        ff.write('box(%f,%f,%f,%f,0)\n' % (xc, yc, x_r - xc, y_r - yc))
        ff.close()
    else:
        if not os.path.exists('crop.reg'):
            open("crop.reg", "w").close()

    # p = subprocess.Popen(["ds9",input_image,"-scale","histequ","-cmap","Cool","-regions","load","crop.reg","-contour","limits"])
    # p.wait()
    # subprocess.call("ds9 %s -scale histequ -cmap Cool -regions %s -contour limits" % (input_image, "crop.reg"), shell=True)
    subprocess.call("ds9 %s -scale histequ -cmap Cool -regions %s" % (input_image, "crop.reg"), shell=True)

    f = open("crop.reg", "r")
    lines = f.readlines()
    for line in lines:
        if 'box(' in line:
            box = line.split(',')
            xc_box = float(box[0].split('(')[1])
            yc_box = float(box[1])
            a_box = float(box[2]) / 2.
            b_box = float(box[3].split(')')[0]) / 2.
            x_l = xc_box - a_box
            y_l = yc_box - b_box
            x_r = xc_box + a_box
            y_r = yc_box + b_box
            break
    f.close()

    if output_image is None:
        output_image = input_image.split('.fits')[0] + '_crop.fits'

    xstartpix = np_to_ds9(x_l)
    ystartpix = np_to_ds9(y_l)
    xendpix = np_to_ds9(x_r)
    yendpix = np_to_ds9(y_r)

    xpixsize = x_r - x_l
    ypixsize = y_r - y_l

    # Apply one of the methods: Montage or Astropy (both tested and give the same results).
    if method == 'Montage':
        import montage_wrapper
        montage_wrapper.commands.mSubimage_pix(input_image, output_image, xstartpix, ystartpix, xpixsize, debug=False,
                                               hdu=hdu, status_file=None, ypixsize=ypixsize)
    else:
        import crop_image
        crop_image.main(input_image, xstartpix, ystartpix, xendpix, yendpix, output_image, hdu=hdu)

    return output_image


def main(input_image, sma=None, smb=None, xc=None, yc=None, output_image=None, PA=0., hdu=0, method='Astropy', square=False,
         galaxy_region=None, offset_size=1.0, crop_box=False):
    '''
    Function to cut out the galaxy image using its outer ellipse.
    Pixel coordinates are in the ds9 format.
    xc,yc - coordinates of the ellipse (pix),
    sma - semi-major axis of the ellipse (pix),
    smb - semi-minor axis of the ellipse (pix),
    PA - position angle of the ellipse (in deg),
    the rotation goes in a counterclockwise direction, starting
    from the positive x-axis (i.e. 0 deg is right, 90 deg is up).
    square - the output image will be a square, otherwise - a rectangle.
    method - Monatge or Astropy (both give the same results).
    '''
    if galaxy_region is not None:
        xc, yc, sma, smb, PA = read_region(galaxy_region)

    sma = sma * offset_size
    smb = smb * offset_size

    if crop_box:
        smb = sma

    # Determine the borders of the ellipse
    hdulist = pyfits.open(input_image)
    data = hdulist[0].data
    nx, ny = data.shape[1], data.shape[0]
    if xc is None:
        xc = nx / 2.
    if yc is None:
        yc = ny / 2.

    if output_image is None:
        output_image = input_image.split('.fits')[0] + '_crop.fits'

    [[x_l, y_l], [x_r, y_r]] = ellipse_borders([xc, yc], sma, smb, PA, nx=nx, ny=ny, square=square)
    xstartpix = np_to_ds9(x_l)
    ystartpix = np_to_ds9(y_l)
    xendpix = np_to_ds9(x_r)
    yendpix = np_to_ds9(y_r)

    xpixsize = x_r - x_l
    ypixsize = y_r - y_l

    # Apply one of the methods: Montage or Astropy (both tested and give the same results).
    if method == 'Montage':
        import montage_wrapper
        montage_wrapper.commands.mSubimage_pix(input_image, output_image, xstartpix, ystartpix, xpixsize, debug=False,
                                               hdu=hdu, status_file=None, ypixsize=ypixsize)
    else:
        import crop_image
        crop_image.main(input_image, xstartpix, ystartpix, xendpix, yendpix, output_image, hdu=hdu)
    return output_image


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create mask of spiral arms")
    parser.add_argument("inputImage", help="Input fits image")
    parser.add_argument("--o", help="Output fits image", type=str, default=None)
    parser.add_argument("--r", help="Galaxy region", type=str, default=None)

    parser.add_argument("--c", nargs='?', const=1, help="Optional: Galaxy center [x,y], separated by comma", type=str,
                        default=None)
    parser.add_argument("--sma", nargs='?', const=1, help="Optional: Galaxy semi-major axis", type=float, default=0.0)
    parser.add_argument("--ell", nargs='?', const=1, help="Optional: Galaxy ellipticity", type=float, default=0.0)
    parser.add_argument("--posang", nargs='?', const=1,
                        help="Optional: Galaxy position angle in degrees. Up=0, Left=90. Default is posang=0.0",
                        type=float, default=0.0)
    parser.add_argument("--method", nargs='?', const=1, help="Optional: Method: Crop or Montage", type=str,
                        default='Crop')

    args = parser.parse_args()

    input_image = args.inputImage
    output_image = args.o
    c = args.c
    sma = args.sma
    ell = args.ell
    posang = args.posang
    galaxy_region = args.r
    method = args.method

    if c is not None:
        [xcen, ycen] = c.split(',')
        xcen = float(xcen)
        ycen = float(ycen)
    else:
        xcen = None
        ycen = None

    smb = sma * (1. - ell)
    main(input_image, sma, smb, xc=xcen, yc=ycen, output_image=output_image, PA=posang, hdu=0, method=method,
         square=False, galaxy_region=galaxy_region)

'''
input_image = '/Users/mosenkov/CurrentWork/Test_Pipeline/NGC5529_SDSS_r.fits'  
output_image = 'crop.fits'
xc = 1997.
yc = 2010.
sma = 368.024235545 / 0.45
q = 1./4.97486278969
smb = sma*q
PA = 23.4680577545

main(input_image, 'crop_m.fits', xc, yc, 2.*sma, 2.*smb, PA, hdu=0, method='Montage', square=False)  
#main(input_image, 'crop_a.fits', xc, yc, 2.*sma, 2.*smb, PA, hdu=0, method='Astropy', square=False)  
'''
