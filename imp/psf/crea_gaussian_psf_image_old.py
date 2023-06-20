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
import iraf_fit_ellipse
import pyfits

sys.path.append('/home/amosenko/MyPrograms/IMAN/Plot')
import radial_profile

def main(input_psf_file, output_psf_file):
    hdulist = pyfits.open(input_psf_file, do_not_scale_image_data=True)
    scidata = hdulist[0].data
    header = hdulist[0].header 
    ny,nx = np.shape(scidata)
    #xc = nx/2.
    #yc = ny/2.

    if nx%2==0:
	xc = int(nx/2. + 1)
    else:
	xc = int(nx/2. + 0.5)

    if ny%2==0:
	yc = int(ny/2. + 1)
    else:
	yc = int(ny/2. + 0.5)

    #print xc,yc
    #y, x = np.indices(scidata.shape)


    #center = np.array([(x.max()-x.min())/2.0, (y.max()-y.min())/2.0])
    #print center
    #exit()

    '''
    iraf_fit_ellipse.main(input_psf_file,xc,yc,20.0,1.0,minsma=0.,maxsma=ceil(math.sqrt(ny**2+nx**2)/2.)+5.,step=1.,outp_format='eps')
    

    sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = np.loadtxt(input_psf_file.split('/')[-1].split('.fits')[0]+'_iraf_ell.txt', usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 6, dtype='str')



    for k in range(len(sma)):
	    if sma[k]=='INDEF': sma[k]=0
	    if inten[k]=='INDEF': inten[k]=0
	    if inten_err[k]=='INDEF': inten_err[k]=0
	    if ell[k]=='INDEF': ell[k]=0
	    if errell[k]=='INDEF': errell[k]=0
	    if PA[k]=='INDEF': PA[k]=0
	    if errPA[k]=='INDEF': errPA[k]=0
	    if x0[k]=='INDEF': x0[k]=0
	    if y0[k]=='INDEF': y0[k]=0
	    if B4[k]=='INDEF': B4[k]=0
	    if errB4[k]=='INDEF': errB4[k]=0
    sma = np.array(sma,dtype='float')
    inten = np.array(inten,dtype='float')
    inten_err = np.array(inten_err,dtype='float')
    ell = np.array(ell,dtype='float')
    errell = np.array(errell,dtype='float')
    PA = np.array(PA,dtype='float')
    errPA = np.array(errPA,dtype='float')
    x0 = np.array(x0,dtype='float')
    y0 = np.array(y0,dtype='float')
    B4 = np.array(B4,dtype='float')
    errB4 = np.array(errB4,dtype='float')
    #rr = sma*np.sqrt(1.-ell)
    tck = interpolate.splrep(sma, inten, s=0)
    
    #ynew = interpolate.splev(4.5, tck, der=0)
    
    new_data = np.copy(scidata)
    print 'Creating new image'
    for k in range(ny):
        for i in range(nx):
            r = math.sqrt( (xc-1-i)**2 + (yc-1-k)**2 )
            new_data[k,i] = abs(interpolate.splev(r, tck, der=0) - scidata[k,i])*100./scidata[k,i] #interpolate.splev(r, tck, der=0)
            if new_data[k,i]<0. or r>nx/2.+0.5:
                new_data[k,i]=0.
    '''
    Rmax = math.sqrt(ny**2+nx**2)/2.+5.

    bin_centers, radial_prof = radial_profile.azimuthalAverage(scidata, center=[xc-1,yc-1], returnradii=True, 
    binsize=1., weights=None, interpnan=True)

    I = []
    r = []      
    for k in range(len(bin_centers)):
        if bin_centers[k]<=Rmax:
	  r.append(bin_centers[k])
	  I.append(radial_prof[k])
    r = np.array(r)
    I = np.array(I)    

    tck = interpolate.splrep(r, I, s=0)
    new_data = np.copy(scidata)
    print 'Creating new image'
    for k in range(ny):
        for i in range(nx):
            r = math.sqrt( (xc-1-i)**2 + (yc-1-k)**2 )
            new_data[k,i] = abs(interpolate.splev(r, tck, der=0) - scidata[k,i])*100./scidata[k,i]
            if new_data[k,i]<0. or r>nx/2.+0.5:
                new_data[k,i]=0.

    
    #hdu = pyfits.PrimaryHDU(new_data/np.sum(new_data),header)
    hdu = pyfits.PrimaryHDU(new_data,header)
    #hdu = pyfits.PrimaryHDU(fabs(new_data-scidata)*100./scidata,header)
    hdu.writeto(output_psf_file,clobber=True)    


input_psf_file = 'Combined_PACS100_V60p_100_100.fits'
output_psf_file = 'psf1.fits'
main(input_psf_file, output_psf_file)
