# This research made use of Photutils, an Astropy package for
# detection and photometry of astronomical sources (Bradley et al. 2022).
#==========================================================================
# This function performs aperture photometry of the given stars

from photutils.aperture import aperture_photometry
from photutils.aperture import CircularAnnulus, CircularAperture
from photutils.aperture import ApertureStats, CircularAperture
from astropy.visualization import simple_norm
from photutils.utils import calc_total_error
from useful_functions import *
from scipy.integrate import quad
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from  scipy.stats import median_abs_deviation as mad
from astropy.stats import SigmaClip
import numpy as np
import os

import platform
opsyst = platform.system()


#def mod(x, a, c):
#    return np.array([quad(lambda z: a*np.exp(-(z*c)**2), 0, y)[0] for y in x])

#-------------------------------------------------------------------------------
def grow_curve(data, xy, gal, fwhm, k0, k2, minr=3, step=1):
    
    xs = [x[0] for x in xy]
    ys = [x[1] for x in xy]
    
    maxr = fwhm*k2
    
    sigclip = SigmaClip(sigma=3.0, maxiters=10)
    
    plt.figure(figsize=(10,10), dpi=300)
    plt.title(gal)
    mags = []
    for i, x, y  in zip(range(len(xs)), xs, ys):
        aps = np.arange(minr, maxr, step)
        apertures = [CircularAperture((x,y), r=r) for r in aps ]
        annulus_apertures =[ CircularAnnulus((x, y), r_in=2.5*fwhm, r_out=3*fwhm) for r in aps]
        fluxes = []
        for aperture, annulus_aperture in zip(apertures, annulus_apertures):
            aperstats_star = ApertureStats(data, aperture, sigma_clip=None)
            aperstats_bkg = ApertureStats(data, annulus_aperture)

            total_bkg = aperstats_bkg.median * aperstats_star.sum_aper_area.value
            flux = aperstats_star.sum - total_bkg
            fluxes.append(flux)
        f = np.array(fluxes)
        df = f
        #B = f[-1]/np.sqrt(f[-1]+total_bkg)
        #if B < 10:
        #    continue
        #f = f/f[-1]
        #phot_table = aperture_photometry(data, apertures)
        #print(phot_table)
        #ann_table = aperture_photometry(data, annulus_apertures)
        #f = [phot_table['aperture_sum_%s' %i] for i in range(len(aps))]
        #f = f/f[-1]
        #total_bkg = aperstats_bkg.median * aperstats_star.sum_aper_area.value
        #f = f - total_bkg
        #mag = -2.5*np.log10(fluxes)
        #mag = np.array([mag[i] - mag[i+1] for i in range(len(mag)-1)])
        #mag = np.append(mag, 0)
        #mag -= mag[-1]
        beta = -2.5*np.log(1 + (np.max(f) - f)/f)
        plt.plot(aps, beta, '-o')
        mags.append(beta)


    mag = np.median(mags, axis=0)
    foo = interp1d(aps, mag)
    #F1  = foo(k0*fwhm)
    #dF1 = np.max(mag) - F1
    #dm  = -2.5*np.log(1 + dF1/F1)
    dm = foo(k0*fwhm)
    #F1  = foo(k0*fwhm)
    



 #  p, covp = curve_fit(mod, xdata=aps, ydata=mag, p0=np.array([mag[-1], 1/2/fwhm]))
 #   print(p)
    #plt.plot(aps, mag, '-o', label='(num=%s)' %(i))
#    plt.plot(aps, mod(aps, p[0], p[1]), '--g')
    #plt.plot([fwhm*k0, fwhm*k0], [np.max(mag), np.min(mag)], '--r')
    #plt.xlabel('Aper, pix')
    #plt.ylabel('$F/F_{min}$')
    #plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #plt.tight_layout()
    #plt.clf()
    #plt.close()
    #plt.savefig('temp_grow_curve.pdf')
    #plt.show()
    return dm

#------------------------------------------------------------------------------
def aperture(data, apertures, annulus_apertures):
    #   def aperture(data, xy, fwhm, k0, k1, k2, apertures=None, annulus_apertures=None):
    
    fluxes = []
    fwhms  = np.array([])
    #if apertures is None:
    #    positions = xy
    #    apertures = CircularAperture(positions, r=k0*fwhm)
    #    annulus_apertures = CircularAnnulus(positions, r_in=k1*fwhm, r_out=k2*fwhm)

    from astropy.stats import SigmaClip
    sigclip = SigmaClip(sigma=3.0, maxiters=10)
    
    coord = []
    for aperture, annulus_aperture in zip(apertures, annulus_apertures):
        aperstats_star = ApertureStats(data, aperture, sigma_clip=None)
        aperstats_bkg = ApertureStats(data, annulus_aperture, sigma_clip=sigclip)
        

        total_bkg = aperstats_bkg.median * aperstats_star.sum_aper_area.value
        flux  = aperstats_star.sum - total_bkg
        coord.append([aperstats_star.xcentroid, aperstats_star.ycentroid])
        fluxes.append(flux)
        fwhms = np.append(fwhms, aperstats_star.fwhm.value)
    #echo(aperstats_star)
    mag =  -2.5*np.log10(fluxes) 
    coord = np.array(coord)
    #print(coord) 
    return mag, fwhms, coord

#-----------------------------------------------------------------------------
def plot_apers(data, apertures, annulus_apertures, out):
    
    #if apertures is None:
    #    positions = xy
    #    apertures = CircularAperture(positions, r=k0*fwhm)
    #    annulus_apertures = CircularAnnulus(positions, r_in=k1*fwhm, r_out=k2*fwhm)

    plt.figure(figsize=(10,10), dpi=300)
    norm = simple_norm(data, 'sqrt', percent=99)
    plt.imshow(data, norm=norm, interpolation='nearest', origin='lower')

    for aperture, annulus_aperture in zip(apertures, annulus_apertures):

        ap_patches = aperture.plot(color='white', lw=2,
                           label='Photometry aperture')
        ann_patches = annulus_aperture.plot(color='red', lw=2,
                                    label='Background annulus')
        handles = (ap_patches[0], ann_patches[0])

    plt.legend(loc=(0.17, 0.05), facecolor='#458989', labelcolor='white',handles=handles, prop={'weight': 'bold', 'size': 11})
    plt.tight_layout()
    plt.savefig(Path(out, 'apertures.png'))
    #plt.show()
    return

#--------------------------------------------------------------------------------------------
def create_ds9_apertures(xy, fwhm, k0, k1, k2):

    xs = [x[0] for x in xy]
    ys = [x[1] for x in xy]

    with open('aper.reg', 'w') as ds9:
        ds9.write("#Region file format: DS9 version 4.1 \n")
        ds9.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
        ds9.write("image\n")


        for i, x, y in zip(range(len(xs)),xs, ys):
            ds9.write('circle(%s, %s,%s) # num=%i \n' %(x, y, k0*fwhm, i))
            ds9.write('annulus(%s, %s, %s, %s) # color=red num=%i \n' %(x, y, k1*fwhm, k2*fwhm, i)) 
    return

#--------------------------------------------------------------------------------------------
import subprocess
def call_ds9(fname, region):
    #subprocess.run(['ds9 -scale histequ %s -region %s' %(fname, region)], shell=True)


    if opsyst=='Linux':
            ds9Proc = subprocess.Popen(["ds9", fname,
                                                    "-regions", region,
                                                    "-scale", "histequ"])
            ds9Proc.wait()
    elif opsyst=='Darwin':
            #subprocess.call(["open","-W","-n","-a","/Applications/SAOImageDS9.app",input_image,"-regions", 'general_mask.reg',"-scale", "histequ"])
            ds9Proc = subprocess.Popen(["/Applications/SAOImageDS9.app/Contents/MacOS/ds9", fname,
                                                    "-regions", region,
                                                    "-scale", "histequ"])
            ds9Proc.wait()
    return

#--------------------------------------------------------------------------------------------
def load_ds9_apertures(region): #, xy):
    positions = []
    r_c       = []
    r_in      = []
    r_out     = []
    inds = []
    with open(region, 'r') as ds9:
        for line in ds9.readlines()[3:]:
            typ, tail = line.split('(')
            if typ == 'circle':
                num = int(line.split('num=')[-1])
                inds.append(num)
                tail = tail.split(')')[0]
                x, y, r = map(float, tail.split(','))
                positions.append((x,y))
                r_c.append(r)
  
            elif typ == 'annulus':
                tail = tail.split(')')[0]
                x, y, r1, r2 = map(float, tail.split(','))
                r_in.append(r1)
                r_out.append(r2)
    aperture = [CircularAperture((p[0], p[1]), r=r) for p, r in zip(positions, r_c)]
    annulus_aperture = [CircularAnnulus((p[0], p[1]) , r_in=ri, r_out=ro) for p, ri, ro in zip(positions, r_in, r_out)] 

    return np.array(aperture), np.array(annulus_aperture), np.array(inds)
#-------------------------------------------------------------------------------------------

def update_ds9_apertures(region, inds, outdir):
    
    if region.name == 'aper_res.reg':
        out_path = Path(outdir, 'aper_res_new.reg')
    else:
        out_path = Path(outdir, 'aper_res.reg')
    out = open(out_path, 'w')

    
    with open(region, 'r') as temp:
        for line in temp.readlines():

            try:
                typ, tail = line.split('(')
            except:
                print(line, end='', file=out)
                continue
            
            if typ == 'circle':
                num = int(line.split('num=')[-1]) 
            elif typ == 'annulus':
                num = int(line.split('num=')[-1])
            else:
                continue
            
            if num in inds:
                print(line, end='', file=out)
    return



