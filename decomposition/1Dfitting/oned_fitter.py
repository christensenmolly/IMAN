import numpy as np
from scipy.interpolate import interp1d
import astropy.io.fits as pyfits
import os
import shutil
import subprocess
import math

fitscopy = ''
imfit_path = ''

def np_to_ds9(x):
    '''
    Function to convert numpy coordinates to ds9 format, i.e
    0 is between 0.5 and 1.5 (not incl), 1 is between 1.5 (incl) and
    2.5 (not incl)
    '''
    return x+1. 

def psf_slice(input_psf, axis='yaxis'):
    imageHDU = pyfits.open(input_psf)[0]
    image = imageHDU.data
    ny,nx = image.shape
    if axis=='yaxis': 
        if nx%2:
            x0 = int(nx/2.+1.)
        else:
            x0 = int(nx/2.+0.5)
        array = image[:,x0-1]
        if os.path.exists('psf_yline.fits'):
            os.remove('psf_yline.fits')
        subprocess.call('%sfitscopy \'%s[%i:%i,%i:%i]\' psf_yline.fits' % (fitscopy,input_psf,x0,x0,1,ny), shell=True)

    else: 
        if ny%2:
            y0 = int(ny/2.+1.)
        else:
            y0 = int(ny/2.+0.5)
        array = image[y0-1,:]
        if os.path.exists('psf_xline.fits'):
            os.remove('psf_xline.fits')
        subprocess.call('%sfitscopy \'%s[%i:%i,%i:%i]\' psf_xline.fits' % (fitscopy,input_psf,1,nx,y0,y0), shell=True)            

    return array/sum(array)

def crea_galfit_input(m0, pix2sec, input_image, sigma_image, psf_image, mask_image, xc, yc, mu0, re, n):
    if sigma_image is None:
        sigma_image = 'none'

    if mask_image is None:
        mask_image = 'none'

    hdu = pyfits.open(input_image)
    image = hdu[0].data
    ny,nx =  np.shape(image)

    hdu_psf = pyfits.open(psf_image)
    image_psf = hdu_psf[0].data
    ny_psf,nx_psf =  np.shape(image_psf)

    if psf_image is None:
        psf_image = 'none'
        nx_psf=1
        ny_psf=1

    

    xc_ch = 0
    yc_ch = 0

       
    xc = np_to_ds9(xc)
    yc = 1.0
        
    s1="""

================================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) %s         # Input data image (FITS file)
B) model.fits          # Output data image block
C) %s          # Sigma image name (made from data if blank or "none") 
D) %s            # Input PSF image and (optional) diffusion kernel
E) 1                   # PSF fine sampling factor relative to data 
F) %s           # Bad pixel mask (FITS image or ASCII coord list)
G) none                # File with parameter constraints (ASCII file) 
H) 1  %i  1  %i  # Image region to fit (xmin xmax ymin ymax)
I) %i    %i          # Size of the convolution box (x y)
J) %f              # Magnitude photometric zeropoint 
K) %f  %f        # Plate scale (dx dy)   [arcsec per pixel]
O) regular             # Display type (regular, curses, both)
P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps

# INITIAL FITTING PARAMETERS
#
#   For component type, the allowed functions are: 
#       sersic, expdisk, edgedisk, devauc, king, nuker, psf, 
#       gaussian, moffat, ferrer, and sky. 
#  
#   Hidden parameters will only appear when they're specified:
#       Bn (n=integer, Bending Modes).
#       C0 (diskyness/boxyness), 
#       Fn (n=integer, Azimuthal Fourier Modes).
#       R0-R10 (coordinate rotation, for creating spiral structures).
#       To, Ti, T0-T10 (truncation function).
# 
# ------------------------------------------------------------------------------
#   par)    par value(s)    fit toggle(s)    # parameter description 
# ------------------------------------------------------------------------------
""" % (input_image, sigma_image, psf_image, mask_image, nx, ny, nx_psf, ny_psf, m0, pix2sec, pix2sec)# 3*n_psf, 3*n_psf)



    s2="""
# Component number: 1
 0) sersic1                 #  Component type
 1) %f %f 0 0  #  Position x, y
 3) %f     1          #  Central surface brightness
 4) %f     1          #  R_e (effective radius)   [pix]
 5) %f      1          #  Sersic index n (de Vaucouleurs n=4) 
 9) 1.0      0          #  Axis ratio (b/a)  
10) 90.0     0          #  Position angle (PA) [deg: Up=0, Left=90]
 Z) 0                      #  Skip this model in output image?  (yes=1, no=0)
""" % (xc, yc, mu0, re, n)

    s3 = "\n================================================================================"
        
    s = s1 + s2 + s3
    
    f = open("input_galf.txt", 'w')
    f.write(s)
    f.close()


def find_initial_guess(r, I):
     L = []
     L_sum = 0.
     re = max(r)/3.
     Ie = max(I)*0.187
     for k in range(len(r)):
            L_sum = L_sum + 1.*I[k] 
            L.append( L_sum)
     for k in range(len(r)):
        if L[k]>=L_sum/2.:
            re = r[k]
            Ie = I[k]
            break
     return re, Ie

def crea_imfit_input(Ie, re, xc, ron=None, gain=None, sigma_image=None, psf_image=None, mask_image=None):
    xc = np_to_ds9(xc)
    f = open('input_imf.txt','w')

    imfit_line = '# %simfit radial_line.fits -c input_imf.txt --ftol 0.00001 --max-threads 1' % (imfit_path)

    if sigma_image!=None:
        imfit_line += ' --noise=sigma_line.fits '

    if psf_image!=None:
        imfit_line += ' --psf=%s ' % (psf_image)        

    if mask_image!=None:
        imfit_line += ' --mask=mask_line.fits '  
    f.write(imfit_line)
    #exit()

    #print >>f, "%s\n" % (imfit_line)
    if ron is not None and gain is not None:
        print >>f,  "GAIN	%.2f" % (gain)
        print >>f,  "READNOISE	%.2f" % (ron)
    s = """
X0		%.2f	fixed
Y0		%.2f	fixed   
FUNCTION Sersic
PA		     90 fixed
ell		0.0 fixed
n		%f 
I_e		%f 
r_e		%f 
""" % (xc, 1.0, 1.0, Ie, re)  
    #print >>f, s
    f.write(s)
    f.close()
    #exit()


def read_ellipse_file(ellipse_file):
        sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = np.loadtxt(ellipse_file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 5, dtype='str')

        
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
        
        return sma, inten, inten_err


def interpolate_profile(sma, inten, inten_err):
    r = np.arange(math.ceil(min(sma)), math.floor(max(sma)), 1)

    f_inten = interp1d(sma, inten)
    f_inten_err = interp1d(sma, inten_err)
    
    I = f_inten(r)
    I_err = f_inten_err(r)
    return r, I, I_err

def create_slice(I_image, output_image, mode='ordinary'):
    #### Create slice images: ###
    # Data image:
    fict_data = np.random.random((2,len(I_image)))
    if mode=='equal':
        fict_data[0,:] = np.ones(len(I_image))
    else:
        fict_data[0,:] = I_image
    
    fict_data[1,:] = np.zeros(len(I_image))
    hdu = pyfits.PrimaryHDU(fict_data)
    hdu.writeto('slice.fits',overwrite=True)
    #exit()
    if os.path.exists(output_image):
        os.remove(output_image)
    subprocess.call('%sfitscopy \'%s[%i:%i,%i:%i]\' %s' % (fitscopy,'slice.fits',1,len(I_image),1,1,output_image), shell=True)
    os.remove('slice.fits')
    #input_image = 'radial_line.fits'   


def main(ellipse_file, average_psf_image, m0, pix2sec, code='imfit'):
    # Read file
    sma, inten, inten_err = read_ellipse_file(ellipse_file)

    # Interpolate the profile:
    r, I, I_err = interpolate_profile(sma, inten, inten_err)
    
    # Create 1D image:
    create_slice(I, 'radial_line.fits')
    
    # Create 1D sigma image
    create_slice(I_err, 'sigma_line.fits', 'equal')

    
    # Create 1D profile:
    psf_slice(average_psf_image, axis='xaxis')
    
    # Create input Galfit/Imfit file
    if code=='galfit':
        crea_galfit_input(m0, pix2sec, 'radial_line.fits', None, 'psf_xline.fits', None, -1, -1, m0 - 2.5*math.log10(np.max(I)/(pix2sec**2)), np.max(r)/3., 1.0)
    else:
        re, Ie = find_initial_guess(r, I)
        crea_imfit_input(Ie, re, -1, ron=None, gain=None, sigma_image='sigma_line.fits', psf_image='psf_xline.fits', mask_image=None)
    

main('ellipse.txt','extended_psf_r.fits', 28.098, 0.396)
    
    
    