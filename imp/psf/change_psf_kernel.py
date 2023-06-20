from astroquery.irsa_dust import IrsaDust
import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import sys
from scipy import interpolate
import math
from math import *
import pyfits

sys.path.append('/home/amosenko/MyPrograms/IMAN/Plot')
sys.path.append('/home/amosenko/MyPrograms/IMAN/FindFluxes')
import radial_profile
import ds9_contour
from Sergey_pipelinelib.rotima import crotima
import rebin_image
import shutil
import os

def main(input_psf_file, output_psf_file, fill_nan=False, Angle=None, Radius=None, Factor=1., FWHM=3., print_mes=True, hdu_level=0):
    hdulist = pyfits.open(input_psf_file, do_not_scale_image_data=True)
    scidata = hdulist[hdu_level].data
    #header = hdulist[0].header 
    ny,nx = np.shape(scidata)
    #print ny,nx
    #exit()
    if nx%2==0:
	xc = int(nx/2. + 1)
    else:
	xc = int(nx/2. + 0.5)

    if ny%2==0:
	yc = int(ny/2. + 1)
    else:
	yc = int(ny/2. + 0.5)



    if fill_nan==True:
        if print_mes==True:
            print "Interpolating the nan-values..."
        interp_data = ds9_contour.replace_nans(scidata, 5, 0.5, int(ceil(FWHM)), "localmean")
        hdu = pyfits.PrimaryHDU(interp_data)
        hdu.writeto('tmp_no_nan.fits',clobber=True)
    else:
        if hdu_level==0:
            shutil.copy(input_psf_file,'tmp_no_nan.fits')
        else:
            hdu = pyfits.PrimaryHDU(scidata)
            hdu.writeto('tmp_no_nan.fits',clobber=True)            

   
    if Angle!=None and Angle!=0.:
        if print_mes==True:
            print "Rotating the frame..."        
        xCenRot, yCenRot = crotima('tmp_no_nan.fits', 'tmp_rot.fits',xc,yc,Angle,layer=1, set_wcs=False)
    else:
        shutil.copy('tmp_no_nan.fits', 'tmp_rot.fits') 
        xCenRot, yCenRot = xc,yc

    xCenRot = int(ceil(xCenRot)) - 1
    yCenRot = int(ceil(yCenRot)) - 1
    #print xCenRot, yCenRot
    #exit()     
    if Radius!=None:
        if print_mes==True:
            print "Croping by radius..."   
        hdulist = pyfits.open('tmp_rot.fits', do_not_scale_image_data=True)
        scidata = hdulist[0].data
        #header = hdulist[0].header 
        ny,nx = np.shape(scidata)
        '''
        for k in range(ny):
            for i in range(nx):
               if (k-yCenRot)**2 + (i-xCenRot)**2 > Radius**2:                   
                scidata[k,i] = 0.
        '''
        xl = max([xCenRot-Radius,0])
        yl = max([yCenRot-Radius,0])
        xr = min([xCenRot+Radius,nx-1])
        yr = min([yCenRot+Radius,ny-1])
        hdu = pyfits.PrimaryHDU(scidata[yl:yr+1,xl:xr+1])
        hdu.writeto(output_psf_file,clobber=True) 
    else:
        shutil.copy('tmp_rot.fits', output_psf_file)
    #exit()
    if Factor!=1.:
        rebin_image.downsample(output_psf_file, Factor, 'tmp_rebin.fits', set_wcs=False,print_mes=print_mes)
        shutil.move('tmp_rebin.fits',output_psf_file)

    if Radius!=None:
        hdulist = pyfits.open(output_psf_file, do_not_scale_image_data=True)
        scidata = hdulist[0].data
        #header = hdulist[0].header 
        ny,nx = np.shape(scidata)

        if nx%2==0:
            xc = int(nx/2. + 1) - 1
        else:
            xc = int(nx/2. + 0.5) - 1

        if ny%2==0:
            yc = int(ny/2. + 1) - 1
        else:
            yc = int(ny/2. + 0.5) - 1

        for k in range(ny):
            for i in range(nx):
               if (k-yc)**2 + (i-xc)**2 > (Radius/Factor)**2:                   
                scidata[k,i] = 0.


        hdu = pyfits.PrimaryHDU(scidata/np.sum(scidata))
        hdu.writeto('tmp_radius.fits',clobber=True)
        hdulist.close()
        shutil.move('tmp_radius.fits',output_psf_file)


    
    for file in ['tmp_no_nan.fits', 'tmp_rot.fits']:
        if os.path.exists(file):
            os.remove(file)


main('/home/amosenko/Downloads/Spitzer/mips24_prf_mosaic_stinytim_2.5_4x.fits', '/home/amosenko/Downloads/Spitzer/psf.fits', fill_nan=False, Angle=30., Radius=None, Factor=1., FWHM=3., print_mes=True, hdu_level=0)

#input_psf_file = '/home/amosenko/CurrentWork/DustPedia/REDO_HERSCHEL_DECOMP/test/UGC10043_i_r_sub_rot_ext.fits'
#output_psf_file = 'test.fits'
#main(input_psf_file, output_psf_file, fill_nan=False, Angle=-10., Radius=None, Factor=1, FWHM=int(ceil(10./3.)), print_mes=True, hdu_level=0)
'''
input_psf_file = '0x5000241aL_PSW_bgmod10_1arcsec.fits'
output_psf_file = 'SPIRE_250_100.fits'
main(input_psf_file, output_psf_file, fill_nan=False, Angle=None, Radius=750, Factor=3, FWHM=int(ceil(10./3.)), print_mes=True, hdu_level=1)

output_psf_file = 'SPIRE_250_160.fits'
main(input_psf_file, output_psf_file, fill_nan=False, Angle=None, Radius=750, Factor=4, FWHM=int(ceil(10./3.)), print_mes=True, hdu_level=1)

output_psf_file = 'SPIRE_250_250.fits'
main(input_psf_file, output_psf_file, fill_nan=False, Angle=None, Radius=750, Factor=6, FWHM=int(ceil(10./3.)), print_mes=True, hdu_level=1)



input_psf_file = '0x5000241aL_PMW_bgmod10_1arcsec.fits'
output_psf_file = 'SPIRE_350_100.fits'
main(input_psf_file, output_psf_file, fill_nan=False, Angle=None, Radius=750, Factor=3, FWHM=int(ceil(10./3.)), print_mes=True, hdu_level=1)

output_psf_file = 'SPIRE_350_160.fits'
main(input_psf_file, output_psf_file, fill_nan=False, Angle=None, Radius=750, Factor=4, FWHM=int(ceil(10./3.)), print_mes=True, hdu_level=1)

output_psf_file = 'SPIRE_350_250.fits'
main(input_psf_file, output_psf_file, fill_nan=False, Angle=None, Radius=750, Factor=6, FWHM=int(ceil(10./3.)), print_mes=True, hdu_level=1)

output_psf_file = 'SPIRE_350_350.fits'
main(input_psf_file, output_psf_file, fill_nan=False, Angle=None, Radius=750, Factor=8, FWHM=int(ceil(10./3.)), print_mes=True, hdu_level=1)



input_psf_file = '0x5000241aL_PLW_bgmod10_1arcsec.fits'
output_psf_file = 'SPIRE_500_100.fits'
main(input_psf_file, output_psf_file, fill_nan=False, Angle=None, Radius=750, Factor=3, FWHM=int(ceil(10./3.)), print_mes=True, hdu_level=1)

output_psf_file = 'SPIRE_500_160.fits'
main(input_psf_file, output_psf_file, fill_nan=False, Angle=None, Radius=750, Factor=4, FWHM=int(ceil(10./3.)), print_mes=True, hdu_level=1)

output_psf_file = 'SPIRE_500_250.fits'
main(input_psf_file, output_psf_file, fill_nan=False, Angle=None, Radius=750, Factor=6, FWHM=int(ceil(10./3.)), print_mes=True, hdu_level=1)

output_psf_file = 'SPIRE_500_350.fits'
main(input_psf_file, output_psf_file, fill_nan=False, Angle=None, Radius=750, Factor=8, FWHM=int(ceil(10./3.)), print_mes=True, hdu_level=1)

output_psf_file = 'SPIRE_500_500.fits'
main(input_psf_file, output_psf_file, fill_nan=False, Angle=None, Radius=750, Factor=12, FWHM=int(ceil(10./3.)), print_mes=True, hdu_level=1)
'''