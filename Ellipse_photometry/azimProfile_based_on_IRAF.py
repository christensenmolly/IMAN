#!/usr/bin/env python

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
import os
import sys

LOCAL_DIR = "/Ellipse_photometry"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'iraf_fitting'))

import plot_iraf_results

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


def main(input_image, ellipse_file, output_model=None, azim_tab='azim_model.txt', mask_image=None, sigma_sky=None, sigma_cal=None, outside_frame=False, center_model=False, stop_where_negative=False, verbosity=True, layer='all', Ell=None, PAA=None):
    
    #data = fits.getdata(input_image)
    hdulist = fits.open(input_image)#, ignore_missing_end=True)
    

    if layer == 'all':
        layers = range(len(hdulist))
    else:
        layers = layer.split(',')
        
    layer_names = []
    for layer in layers:
        layer = int(layer)
        print('Layer %i' % (layer))
        data = hdulist[layer].data
        header = hdulist[layer].header
        
        layer_name = None
        if 'NAME_OF_LAYER' in header:
            layer_name = header['NAME_OF_LAYER']

        if layer!=0:
                if layer_name is not None:
                    layer_names.append(layer_name)
                    if layer_names.count(layer_name)>1:
                        layer_name = '%s%i' % (layer_name, layer_names.count(layer_name))
                    Azim_tab = azim_tab.split('.txt')[0]+'_%s.txt' % (layer_name)
                else:
                    Azim_tab = azim_tab.split('.txt')[0]+'_%i.txt' % (layer)
        else:
                Azim_tab = azim_tab

 
        if mask_image is not None:
            mask = fits.getdata(mask_image)
        else:
            mask = None
    
        y_size, x_size = data.shape
        
        sma, inten, inten_err, ell, errell, PA, errPA, x0, y0, B4, errB4 = plot_iraf_results.read_ell(ellipse_file)
        PA[PA<0.] = PA[PA<0.] + 180.
        PA = PA - 90.
    
        ellipses = IsophoteList([])
        f_res = open(Azim_tab, 'w') #### Output model is saved in a text-file as well
        f_res.truncate(0)
        f_res.write("# sma[pix]\tflux[DN]\tflux_err[DN]\n")
    
        if center_model:
            NX = ceil(2.*max(sma))+1
            NY = NX
            X0 = int(NX/2.+0.5)-1
            Y0 = int(NY/2.+0.5)-1

        for k in range(len(sma)):
          if sma[k]>=1.:
            try:
                if Ell is not None:
                    elll = Ell
                else:
                    elll = ell[k]

                if Ell is not None:
                    PAAA = PAA
                else:
                    PAAA = PA[k]
                    
                ellipse = SimpleNamespace(x0=x0[k], y0=y0[k], sma=sma[k], intens=0,
                                        eps=elll, pa=np.radians(PAAA), grad=0, a3=0, b3=0, a4=0, b4=0, fit_ok=False)
                if outside_frame==False:
                    if not is_ellipse_inside_frame(ellipse, x_size, y_size, margin=5):
                        break
                _, intensity = get_inten_along_ellipse(ellipse, data, mask) #### Added mask
                #ellipse.intens = np.mean(intensity)

                std = np.std(intensity)
                mean, median, std = sigma_clipped_stats(intensity, sigma=3.0, maxiters=5)
                if not np.isnan(median):
                    ellipse.intens = median
                else:
                    ellipse.intens = mean
            
                sigma_instr = std/sqrt(float(np.shape(intensity)[0]))
                if sigma_sky is not None:
                    sigma_sky = sigma_sky/sqrt(float(np.shape(intensity)[0]))
                else:
                    sigma_sky = 0.
            
                if sigma_cal is None:
                    sigma_cal = 0.
                
                sigma_sum = sqrt( sigma_instr**2 + sigma_sky**2 + sigma_cal**2 )
            

                f_res.write("%7.2f\t%15.8e\t%15.8e\n" % (sma[k], ellipse.intens, sigma_sum))
                if center_model:
                    ellipse.x0=X0
                    ellipse.y0=Y0
                if stop_where_negative and  ellipse.intens<0.:
                    NX = ceil(2.*sma[k])+1
                    NY = NX
                    X0 = int(NX/2.+0.5)-1
                    Y0 = int(NY/2.+0.5)-1
                    for i in range(len(ellipses)):
                        ellipses[i].x0 =X0
                        ellipses[i].y0 =Y0
                        
                    break
                ellipses.append(ellipse) 
                if verbosity: print('Radius: %.2f/%.2f' % (sma[k],max(sma)))
            except:
                if verbosity: print('Radius: %.2f/%.2f FAILED!' % (sma[k],max(sma)))  
                z=1
    f_res.close()
    
    
    if output_model is not None:
        if not center_model:
            model = build_ellipse_model(data.shape, ellipses)
        else:
            model = build_ellipse_model((NY,NX), ellipses)
        fits.PrimaryHDU(model, header).writeto(output_model, overwrite=True)
    if verbosity: print('Done!')
    



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input_image", help="Name of image file")
    parser.add_argument("ellipse_file", help="Name of IRAF/ELLIPSE file")
    parser.add_argument("--output_model", help="Name of output model",
                        type=str, default=None) #'azim_model.fits'
    parser.add_argument("--azim_tab", help="Name of output model",
                        type=str, default='azim_model.txt')
    parser.add_argument("--mask", help="Mask image (zero pixels are good).",
                        type=str, default=None)
    parser.add_argument("--sigma_cal", help="Calibration error (in DN).",
                        default=None, type=float)
    parser.add_argument("--sigma_sky", help="RMS of the sky (in DN).",
                        default=None, type=float)
    parser.add_argument("--outside_frame", help="Count ellipses which lie outside the frame.",
                        default=False, action="store_true")
    parser.add_argument("--hdu", help="Name of hdu layer",
                        type=str, default='all')    
    parser.add_argument("--ell", help="Ellipticity.",
                        default=None, type=float)
    parser.add_argument("--pa", help="Position angle (deg). Up=90, Right=0, counterclockwise.",
                        default=None, type=float)
    args = parser.parse_args()

    input_image = args.input_image
    ellipse_file = args.ellipse_file
    output_model = args.output_model
    azim_tab = args.azim_tab
    mask_image = args.mask

    sigma_cal = args.sigma_cal
    sigma_sky = args.sigma_sky
    outside_frame = args.outside_frame
    layer = args.hdu
    ell = args.ell
    pa = args.pa
    #print(outside_frame)
    #exit()
    
    
    
    main(input_image, ellipse_file, output_model, azim_tab, mask_image, sigma_sky, sigma_cal, outside_frame=outside_frame, layer=layer, Ell=ell, PAA = pa)
