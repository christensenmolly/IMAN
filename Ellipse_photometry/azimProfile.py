#!/usr/bin/env python
# EXAMPLE:

from math import sin
from math import cos
from math import pi
from math import modf
from math import sqrt
from math import ceil
import numpy as np
import argparse
from types import SimpleNamespace

from astropy.io import fits
from photutils.isophote import IsophoteList
from photutils.isophote import build_ellipse_model
from astropy.stats import sigma_clipped_stats


def ds9_to_np(x):
    '''
    Function to convert ds9 pixel coordinates to numpy format, i.e
    0.5 pix is 0, 1.5 pix is 1
    '''
    return int(round(x)) - 1

def get_inten_along_ellipse(ellipse, data, mask=None): #### Added mask
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
        try:
          if mask[iiy][iix]==0.: #### Take into account only those pixels which are not masked
            intensity[i] = ((1.0-fx)*(1.0-fy)*data[iiy][iix] + fx*(1.0-fy)*data[iiy][iix+1] +
                        fy*(1.0-fx)*data[iiy+1][iix] + fx*fy*data[iiy+1][iix+1])
          else:
            intensity[i] = float('nan')  
        except:
            intensity[i] = float('nan')
    intensity = np.array(intensity)
    intensity = intensity[~np.isnan(intensity)]
    return exc_anomaly, intensity


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


def main(input_image, output_model=None, azim_tab='azim_model.txt', mask_image=None, xcen=-1, ycen=-1, ell=0., posang=0., sma_min=-1, sma_max=-1, step=1., sigma_sky=None, sigma_cal=None, outside_frame=False, center_model=False, stop_where_negative=False, verbosity=True, linear=True):
    
    #data = fits.getdata(input_image)
    hdulist = fits.open(input_image)#, ignore_missing_end=True)
    data = hdulist[0].data
    header = hdulist[0].header
    
    if mask_image is not None:
        mask = fits.getdata(mask_image)
    else:
        mask = None
    
    y_size, x_size = data.shape
    
    if x_size%2==0:
        image_x_center = int(x_size/2. + 1)
    else:
        image_x_center = int(x_size/2. + 0.5)
    
    if y_size%2==0:
        image_y_center = int(y_size/2. + 1)
    else:
        image_y_center = int(y_size/2. + 0.5)    
    
    object_x_center = xcen if xcen > 0 else image_x_center
    object_y_center = ycen if ycen > 0 else image_y_center

    object_x_center = ds9_to_np(object_x_center)
    object_y_center = ds9_to_np(object_y_center)

    max_ell_sma = sma_max if sma_max > 0 else max(object_x_center, object_y_center,
                                                            x_size-object_x_center, y_size-object_y_center)

    min_ell_sma = sma_min if sma_min > 0 else 1.

    ellipses = IsophoteList([])
    f_res = open(azim_tab, 'w') #### Output model is saved in a text-file as well
    f_res.truncate(0)
    f_res.write("# sma[pix]\tflux[DN]\tflux_err[DN]\n")
    
    if center_model:
        NX = ceil(2.*max_ell_sma)+1
        NY = NX
        X0 = int(NX/2.+0.5)-1
        Y0 = int(NY/2.+0.5)-1
    
    if linear:
        smas = np.arange(min_ell_sma, max_ell_sma, step)
    else:
        smaa = min_ell_sma
        smas = []
        while smaa < max_ell_sma:
            smaa = smaa * (1. + step)
            smas.append(smaa)
        smas = np.array(smas)

    
    
    for k in range(len(smas)-1):
      if True:
        smaa = smas[k]
        rads = []
        ints = []
        Npix = []
        while smaa < smas[k+1]:
            rads.append(smaa)
            ellipse = SimpleNamespace(x0=object_x_center, y0=object_y_center, sma=smaa, intens=0,
                                  eps=ell, pa=np.radians(posang), grad=0, a3=0, b3=0, a4=0, b4=0, fit_ok=False)
            if outside_frame==False:
                if not is_ellipse_inside_frame(ellipse, x_size, y_size, margin=5):
                    break
            _, intensity = get_inten_along_ellipse(ellipse, data, mask) #### Added mask
            #ellipse.intens = np.mean(intensity)

            std = np.std(intensity)
            mean, median, std = sigma_clipped_stats(intensity, sigma=1.0, maxiters=5)
            if not np.isnan(median):
                ellipse.intens = median
            else:
                ellipse.intens = mean
            ints.append(ellipse.intens)
            
            Npix.append(np.shape(intensity)[0])
            if linear and step==1:
                sigma_instr = std#/sqrt(float(np.shape(intensity)[0]))
                if sigma_sky is not None:
                    Sigma_sky = sigma_sky/sqrt(float(np.shape(intensity)[0]))
                else:
                    Sigma_sky = 0.
                
                if sigma_cal is None:
                    sigma_cal = 0.
            
                sigma_sum = sqrt( sigma_instr**2 + Sigma_sky**2 + sigma_cal**2 )
            
            smaa+=1.
        
        ellipse.intens = np.mean(ints)
        ellipse.sma = np.mean(rads)
        sma = np.mean(rads)
        
        if linear and step==1:
            error = sigma_sum
        else:
            if sigma_sky is not None:
                    Sigma_sky = sigma_sky/sqrt(float(np.sum(Npix)))
            else:
                    Sigma_sky = 0.
                
            if sigma_cal is None:
                    sigma_cal = 0.
                    
            error = sqrt( (np.std(ints))**2 + Sigma_sky**2 + sigma_cal**2 )
        
        f_res.write("%7.2f\t%15.8e\t%15.8e\n" % (np.mean(rads), np.mean(ints), error))
        
        
        if center_model:
            ellipse.x0=X0
            ellipse.y0=Y0
        if stop_where_negative and  ellipse.intens<0.:
            NX = ceil(2.*sma)+1
            NY = NX
            X0 = int(NX/2.+0.5)-1
            Y0 = int(NY/2.+0.5)-1
            for i in range(len(ellipses)):
                ellipses[i].x0 =X0
                ellipses[i].y0 =Y0
                
            break
        ellipses.append(ellipse) 
        if verbosity: print('Radius: %.2f/%.2f' % (sma,max_ell_sma))
      else:
        if verbosity: print('Radius: %.2f/%.2f FAILED!' % (sma,max_ell_sma))  
        z=1
    f_res.close()
    
    
    
    
    
    
    
    if output_model is not None:
        if not center_model:
            model = build_ellipse_model(data.shape, ellipses)
        else:
            model = build_ellipse_model((NY,NX), ellipses)
        fits.PrimaryHDU(model, header).writeto(output_model, overwrite=True)
    if verbosity: print('Done!')
    return output_model, azim_tab



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input_image", help="Name of image file")
    parser.add_argument("--output_model", help="Name of output model",
                        type=str, default=None) #'azim_model.fits'
    parser.add_argument("--azim_tab", help="Name of output model",
                        type=str, default='azim_model.txt')
    parser.add_argument("--mask_image", help="Mask image (zero pixels are good).",
                        type=str, default=None)
    parser.add_argument("--xcen", help="x-coordinate of the object center (image center by default).",
                        type=float, default=-1)
    parser.add_argument("--ycen", help="y-coordinate of the object center (image center by default).",
                        type=float, default=-1)
    parser.add_argument("--ell", help="Ellipticity. Default is ell=0.0", default=0.0, type=float)
    parser.add_argument("--posang", help="Position angle of the ellipse in degrees. Up=90, Right=0, counterclockwise. Default is posang=0.0",
                        default=0.0, type=float) # Now in degrees!
    parser.add_argument("--sma_min", help="Minimal major axis of ellipse. Default: 1.",
                        default=-1, type=float)
    parser.add_argument("--sma_max", help="Maximal major axis of ellipse. Default: fitted to the image size.",
                        default=-1, type=float)
    parser.add_argument("--step", help="Linear step in semi-major axis length between successive ellipses.",
                        default=1., type=float)
    parser.add_argument("--sigma_cal", help="Calibration error (in DN).",
                        default=None, type=float)
    parser.add_argument("--sigma_sky", help="RMS of the sky (in DN).",
                        default=None, type=float)
    parser.add_argument("--outside_frame", help="Count ellipses which lie outside the frame.",
                        default=False, action="store_true")    
    parser.add_argument("--linear", help="Linear step?",
                        default=False, action="store_true")  

    args = parser.parse_args()

    input_image = args.input_image
    output_model = args.output_model
    azim_tab = args.azim_tab
    mask_image = args.mask_image
    xcen = args.xcen
    ycen = args.ycen
    ell = args.ell
    posang = args.posang
    sma_min = args.sma_min
    sma_max = args.sma_max
    step = args.step
    sigma_cal = args.sigma_cal
    sigma_sky = args.sigma_sky
    outside_frame = args.outside_frame
    linear = args.linear
    #print(outside_frame)
    #exit()
    
    
    
    main(input_image, output_model, azim_tab, mask_image, xcen, ycen, ell, posang, sma_min, sma_max, step, sigma_sky, sigma_cal, outside_frame=outside_frame, linear=linear)
