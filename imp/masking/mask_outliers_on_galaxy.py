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
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from scipy.signal import find_peaks, find_peaks_cwt
from scipy.signal import argrelextrema

def outliers_iqr(ys, low_quart, up_quart):
    quartile_1, quartile_3 = np.percentile(ys, [int(low_quart), int(up_quart)])
    iqr = quartile_3 - quartile_1
    lower_bound = quartile_1 - (iqr * 1.5)
    upper_bound = quartile_3 + (iqr * 1.5)
    return np.where((ys > upper_bound) | (ys < lower_bound))

def define_backgr(intensity):
    # Function to define maximum inetnsity of background. Everything larger than that should be masked!
    hist, bins = np.histogram(intensity, bins=15) # auto?
    bincenters = 0.5*(bins[1:]+bins[:-1])
    peaks_max, _ = find_peaks(hist, prominence=1)      # BEST!
    peaks_min = argrelextrema(hist, np.less)[0]
    if len(peaks_max)>=2 and len(peaks_min)>=1 and len(peaks_max)<=2:
        I_backg_max = min(bincenters[peaks_min])
    else:
        I_backg_max = max(intensity)
    return I_backg_max


def get_inten_along_ellipse(ellipse, data, new_mask, low_quart, up_quart, mask=None): #### Added mask
    """
    Function makes an azimuthal slice along the ellipse. The smoothing is
    performed by the averaging of several ellipses with a bit different radii
    """
    exc_anomaly = np.linspace(0, 2*pi, int(10*ellipse.sma))
    intensity = np.zeros_like(exc_anomaly)
    sinpa = sin(ellipse.pa)
    cospa = cos(ellipse.pa)
    cose = np.cos(exc_anomaly)
    sine = np.sin(exc_anomaly)
    q_value = 1 - ellipse.eps
    if mask is None:
      for i in range(len(exc_anomaly)):
        fx, ix = modf(ellipse.x0 + ellipse.sma * cose[i] * cospa - ellipse.sma * q_value * sine[i] * sinpa)
        fy, iy = modf(ellipse.y0 + ellipse.sma * q_value * sine[i] * cospa + ellipse.sma * cose[i] * sinpa)
        iix = int(ix)
        iiy = int(iy)
        intensity[i] = ((1.0-fx)*(1.0-fy)*data[iiy][iix] + fx*(1.0-fy)*data[iiy][iix+1] +
                        fy*(1.0-fx)*data[iiy+1][iix] + fx*fy*data[iiy+1][iix+1])
    else:
      for i in range(len(exc_anomaly)):
        fx, ix = modf(ellipse.x0 + ellipse.sma * cose[i] * cospa - ellipse.sma * q_value * sine[i] * sinpa)
        fy, iy = modf(ellipse.y0 + ellipse.sma * q_value * sine[i] * cospa + ellipse.sma * cose[i] * sinpa)
        iix = int(ix)
        iiy = int(iy)
        if mask[iiy][iix]==0.: #### Take into account only those pixels which are not masked
            intensity[i] = ((1.0-fx)*(1.0-fy)*data[iiy][iix] + fx*(1.0-fy)*data[iiy][iix+1] +
                        fy*(1.0-fx)*data[iiy+1][iix] + fx*fy*data[iiy+1][iix+1])
    
    
    I_backg_max = define_backgr(intensity)
    
    #mean = np.mean(intensity)
    #mean, median, std = sigma_clipped_stats(intensity, sigma=3.0, maxiters=5)
    #std = np.std(intensity)
    inds = list(outliers_iqr(intensity, low_quart, up_quart)[0])
    
    for i in range(len(exc_anomaly)):
        fx, ix = modf(ellipse.x0 + ellipse.sma * cose[i] * cospa - ellipse.sma * q_value * sine[i] * sinpa)
        fy, iy = modf(ellipse.y0 + ellipse.sma * q_value * sine[i] * cospa + ellipse.sma * cose[i] * sinpa)
        iix = int(ix)
        iiy = int(iy)
        #I_res = (data[iiy][iix]-mean)/mean
        #if (I_res >=0.25 and I_res < 1.0) or i in inds: #median+1.*std:# or data[iiy][iix]<=median-1.*std:   # if i in inds:
        if i in inds:
            new_mask[iiy][iix] = 1

    if False:#ellipse.sma==44:
        print(list(outliers_iqr(intensity)[0]))
        print(median,std)
        plt.figure(0,figsize=(5, 5))
        N, bins, patches = plt.hist(intensity, 'auto', normed=False, facecolor='gray', alpha=0.75)
        peaks_max, _ = find_peaks(N, prominence=1)      # BEST!
        peaks_min = argrelextrema(N, np.less)[0]
        #print(peaks_min[0])
        #exit()
        plt.plot(bins[peaks_max], N[peaks_max], "o")
        plt.plot(bins[peaks_min], N[peaks_min], "*")        
        #print(N, bins)
        #plt.plot(bins[0:-1],N)
        plt.show()

        '''
        plt.figure(0,figsize=(5, 5))
        plt.plot(range(len(intensity)),intensity,'o')
        plt.show()        
        '''
        z=1
    return new_mask


def is_ellipse_inside_frame(ellipse, x_size, y_size, margin=5):
    """
    Function checks if the given ellipse fits inside of the image
    with a given margin
    """
    exc_anomaly = np.linspace(0, 2*pi, int(10*ellipse.sma))
    sinpa = sin(ellipse.pa)
    cospa = cos(ellipse.pa)
    q_value = 1 - ellipse.eps
    x_ell = ellipse.x0 + ellipse.sma * np.cos(exc_anomaly) * cospa - ellipse.sma * q_value * np.sin(exc_anomaly) * sinpa
    y_ell = ellipse.y0 + ellipse.sma * q_value * np.sin(exc_anomaly) * cospa + ellipse.sma * np.cos(exc_anomaly) * sinpa

    x_min = min(x_ell)
    x_max = max(x_ell)
    y_min = min(y_ell)
    y_max = max(y_ell)

    if (x_min > margin) and (y_min > margin) and (x_max < x_size-margin) and (y_max < y_size-margin):
        return True
    return False


def main(args):
    data = fits.getdata(args.image)
    if args.mask is not None:
        mask = fits.getdata(args.mask)
    else:
        mask = None
    y_size, x_size = data.shape
    image_x_center, image_y_center = x_size/2, y_size/2
    object_x_center = args.xcen if args.xcen > 0 else image_x_center
    object_y_center = args.ycen if args.ycen > 0 else image_y_center

    max_ell_sma = args.sma_max if args.sma_max > 0 else max(object_x_center, object_y_center,
                                                            x_size-object_x_center, y_size-object_y_center)

    new_mask = np.zeros(shape=(y_size, x_size))
    for sma in np.arange(args.sma_min, max_ell_sma, args.step):
        ellipse = SimpleNamespace(x0=object_x_center, y0=object_y_center, sma=sma, intens=0,
                                  eps=args.ell, pa=np.radians(90.+args.posang), grad=0, a3=0, b3=0, a4=0, b4=0)
        if not is_ellipse_inside_frame(ellipse, x_size, y_size, margin=5):
            break
        new_mask = get_inten_along_ellipse(ellipse, data, new_mask, args.low_quart, args.up_quart, mask) #### Added mask

        print('Radius: %.2f/%.2f' % (sma,max_ell_sma))
    #new_mask = ndimage.gaussian_filter(new_mask, sigma=3., order=0)
    #new_mask = new_mask[np.where(new_mask>0.3)]
    #new_mask = np.place(new_mask, new_mask<0.3, 0)
    #new_mask = np.place(new_mask, new_mask>0.3, 1)
    #new_mask[new_mask<=0.5]=0
    #new_mask[new_mask>0.5]=1
    fits.PrimaryHDU(data=new_mask).writeto(args.model, overwrite=True)
    print('Done!')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("image", help="Name of image file")
    parser.add_argument("--model", help="Name of output model",
                        type=str, default='inner_mask.fits') #'azim_model.fits'
    parser.add_argument("--azim_tab", help="Name of output model",
                        type=str, default='azim_model.txt')
    parser.add_argument("--mask", help="Mask image (zero pixels are good).",
                        type=str, default=None)
    parser.add_argument("--xcen", help="x-coordinate of the object center (image center by default).",
                        type=float, default=-1)
    parser.add_argument("--ycen", help="y-coordinate of the object center (image center by default).",
                        type=float, default=-1)
    parser.add_argument("--ell", help="Ellipticity. Default is ell=0.0", default=0.0, type=float)
    parser.add_argument("--posang", help="Position angle of the ellipse in degrees. Up=0, Left=90. Default is posang=0.0",
                        default=0.0, type=float) # Now in degrees!
    parser.add_argument("--sma-max", help="Maximal major axis of ellipse. Default: fitted to the image size.",
                        default=-1, type=float)
    parser.add_argument("--sma-min", help="Minimal major axis of ellipse. Default: 1.",
                        default=1., type=float)
    parser.add_argument("--step", help="Linear step in semi-major axis length between successive ellipses.",
                        default=1., type=float)
    parser.add_argument("--up_quart", help="Upper quartile to remove outliers.",
                        default=75, type=int)
    parser.add_argument("--low_quart", help="Lower quartile to remove outliers.",
                        default=25, type=int)
    args = parser.parse_args()
    main(args)
