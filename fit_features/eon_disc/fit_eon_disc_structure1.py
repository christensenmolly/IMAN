#! /usr/bin/env python
import numpy as np
import gzip
import shutil
from joblib import Parallel, delayed
import astropy.io.fits as pyfits
import sys
import os
import shutil
import time
import subprocess
import glob
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings("ignore")

FNULL = open(os.devnull, 'w')
sys.path.append( '/home/amosenko/MEGA/MyPrograms/IMAN/imp/rotate')
import rotate_image

imfit_path = '/home/amosenko/MEGA/MyPrograms/imfit-1.6.1/'


####
n_times = 10
####


def main(input_image,input_psf,x0,y0,x0_psf,radius_x_in,radius_x_out,radius_y,n_bins,ron=5.474,gain=4.76,z0_thin=0.,z0_thick=0.,Rbreak=0.,sigma_image=None,error_est=True):
  print bcolors.OKBLUE+'\n\n************ IMFIT analysis of vertical profiles 2015 ************' + bcolors.ENDC
  print "Analyzing..."
  PARS = {}
  PARS["Z_CUTS_IMFIT"] = [float(nan),float(nan),float(nan),float(nan),float(nan),float(nan),float(nan),float(nan)]





def crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind):
    I_image = []; I_mask = []; I_sigma = []
    ny,nx = np.shape(image_data)
    
    for i in range(ny):  
        image_array = list(image_data[i,x_min_ind:x_max_ind])
        mask_array = list(mask_data[i,x_min_ind:x_max_ind])
        sigma_array = list(sigma_data[i,x_min_ind:x_max_ind])
        
        a_image = np.ma.array(image_array, mask = mask_array)
        a_sigma = np.ma.array(sigma_array, mask = mask_array)
        try:
            I_image.append(float(np.ma.average(a_image, weights=sigma_array)))
            I_mask.append(False)
            I_sigma.append(float(np.ma.average(a_sigma, weights=sigma_array)))
        except:
            I_image.append(0.)
            I_mask.append(True)
            I_sigma.append(0.)
    return np.array(I_image), np.array(I_mask), np.array(I_sigma)
        



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
        subprocess.call('fitscopy \'%s[%i:%i,%i:%i]\' psf_yline.fits' % (input_psf,x0,x0,1,ny), shell=True)

    else: 
        if ny%2:
            y0 = int(ny/2.+1.)
        else:
            y0 = int(ny/2.+0.5)
        array = image[y0-1,:]
        if os.path.exists('psf_xline.fits'):
            os.remove('psf_xline.fits')
        subprocess.call('fitscopy \'%s[%i:%i,%i:%i]\' psf_xline.fits' % (input_psf,1,nx,y0,y0), shell=True)            

    return array/sum(array)


def disk_edge_soph(z, I0d, z0, z_c):
    #*** For edge-on disk SB in mag/arcsec^2 (along z-axis). n - index of this law. ***
    return I0d * (1.0 / np.cosh(1.*np.fabs(z-z_c)/(2.*z0)) ) **(2./1.)




def crea_imf(L0, n, z0, y0, ron=None, gain=None):
    f = open('input_imf.txt','w')
    if ron is not None and gain is not None:
        print >>f,  "GAIN	%.2f" % (gain)
        print >>f,  "READNOISE	%.2f" % (ron)
    '''
    print >>f,  "X0	0.5	fixed"
    print >>f,  "Y0	%.2f	fixed" % (y0)

    print >>f,  "FUNCTION Z_eon"
    print >>f,  "I_0	%.10f" % (I0)
    print >>f,  "n	%.2f	0.1,100" % (n)
    print >>f,  "z_0	%.2f" % (z0)
    '''
    s = """
X0		%.2f	fixed
Y0		%.2f	fixed    
FUNCTION EdgeOnDisk
PA		%f	fixed
L_0		%f
h		%f	fixed
n		%f
z_0		%f    
""" % (0.5,y0,90.,L0,10.,n,z0)  
    print >>f, s
    f.close()
    #exit()


def crea_imf_double(L0, n, z0, y0, ron=None, gain=None):
    L01 = L0
    L02 = L0/4.
    n1 = n
    n2 = n
    z01 = z0/3.
    z02 = z0*1.5
    y0 = y0

    f = open('input_imf.txt','w')
    if ron is not None and gain is not None:
        print >>f,  "GAIN	%.2f" % (gain)
        print >>f,  "READNOISE	%.2f" % (ron)
    s = """
X0		%.2f	fixed
Y0		%.2f	fixed    
FUNCTION EdgeOnDisk
PA		%f	fixed
L_0		%f
h		%f	fixed
n		%f
z_0		%f

FUNCTION EdgeOnDisk
PA		%f	fixed
L_0		%f
h		%f	fixed
n		%f
z_0		%f  
""" % (0.5,y0,90.,L01,10.,n1,z01,90.,L02,10.,n2,z02)  
    print >>f, s
    f.close()


def crea_imf_de(L0,n,z0,y0,ron,gain):
    if n>=100.:
        n=100.
    f = open('input_imf.txt','w')
    if ron is not None and gain is not None:
        print >>f,  "GAIN	%.2f" % (gain)
        print >>f,  "READNOISE	%.2f" % (ron)
    '''
    print >>f,  "X0	0.5	fixed"
    print >>f,  "Y0	%.2f	fixed" % (y0)

    print >>f,  "FUNCTION Z_eon"
    print >>f,  "I_0	%.10f" % (I0)
    print >>f,  "n	%.2f	0.1,100" % (n)
    print >>f,  "z_0	%.2f" % (z0)
    '''
    s = """
X0		%.2f	fixed
Y0		%.2f	fixed    
FUNCTION EdgeOnDisk
PA		%f	fixed
L_0		%.10f	%.10f,%.10f
h		%f	fixed
n		%.2f	0.1,100
z_0		%.2f	%.2f,%.2f    
""" % (0.5,y0,90.,L0,L0/10.,L0*10.,10.,n,z0,z0/3.,z0*3.)  
    print >>f, s
    f.close()


def crea_imf_de2(L01,n1,z01,L02,n2,z02,y0,ron,gain):
    f = open('input_imf.txt','w')
    if ron is not None and gain is not None:
        print >>f,  "GAIN	%.2f" % (gain)
        print >>f,  "READNOISE	%.2f" % (ron)
    '''
    print >>f,  "X0	0.5	fixed"
    print >>f,  "Y0	%.2f	fixed" % (y0)

    print >>f,  "FUNCTION Z_eon"
    print >>f,  "I_0	%.10f" % (I0)
    print >>f,  "n	%.2f	0.1,100" % (n)
    print >>f,  "z_0	%.2f" % (z0)
    '''
    s = """
X0		%.2f	fixed
Y0		%.2f	fixed    
FUNCTION EdgeOnDisk
PA		%f	fixed
L_0		%.10f	%.10f,%.10f
h		%f	fixed
n		%.2f	0.1,100
z_0		%.2f	%.2f,%.2f    

FUNCTION EdgeOnDisk
PA		%f	fixed
L_0		%.10f	%.10f,%.10f
h		%f	fixed
n		%.2f	0.1,100
z_0		%.2f	%.2f,%.2f  
""" % (0.5,y0,90.,L01,L01/10.,L01*10.,10.,n1,z01,z01/3.,z01*3.,L02,L02/10.,L02*10.,10.,n2,z02,z02/3.,z02*3.)  
    print >>f, s
    f.close()

def read_imf(imf_file):
  f = open(imf_file,'r')
  L0 = 99999.
  for line in f:
    if "#   Reduced value:" in line:
      chi2 = float(line.split()[3])
    if "Y0\t" in line:
        y0 = float(line.split()[1])
    if 'L_0' in line:
      L0 = float(line.split()[1])
    if 'n' in line and L0!=99999.:
      n = float(line.split()[1])
    if 'z_0' in line:
      z0 = float(line.split()[1])
  f.close()
  return y0,L0,n,z0,chi2


def read_imf_double(imf_file):
  f = open(imf_file,'r')
  func = 0
  n1 = 99999.
  for line in f:
    if "#   Reduced value:" in line:
      chi2 = float(line.split()[3])
    if "FUNCTION Z_EON" in line and func==0:
      func=1
    if 'L_0' in line and func==1:
      L01 = float(line.split()[1])      
    if 'n' in line and func==1:
      n1 = float(line.split()[1])
    if 'z_0' in line and func==1:
      z01 = float(line.split()[1])
    if "FUNCTION Z_EON" in line and func==1 and n1!=99999.:
      func=2
    if 'L_0' in line and func==2:
      L02 = float(line.split()[1])   
    if 'n' in line and func==2:
      n2 = float(line.split()[1])
    if 'z_0' in line and func==2:
      z02 = float(line.split()[1])      
  f.close()
  
  return L01,n1,z01,L02,n2,z02,chi2


def run_imf(L0, n, z0, y0, sigma_image, psf_image, mask_image, ron, gain):
    crea_imf(L0, n, z0, y0, ron, gain)

    imfit_line = '%simfit vertical_line.fits -c input_imf.txt --ftol 0.00001 --max-threads 1' % (imfit_path)
    if sigma_image!=None:
        imfit_line += ' --noise sigma_line.fits '

    if psf_image!=None:
        imfit_line += ' --psf psf_line.fits '        

    if mask_image!=None:
        imfit_line += ' --mask mask_line.fits '   


    subprocess.call(imfit_line, shell=True, stdout=FNULL)
    
    shutil.copy('bestfit_parameters_imfit.dat','bestfit_parameters_imfit_lm.dat')


def run_imf_double(L0, n, z0, y0, sigma_image, psf_image, mask_image, ron, gain):
    crea_imf_double(L0, n, z0, y0, ron=None, gain=None)

    imfit_line = '%simfit vertical_line.fits -c input_imf.txt --ftol 0.00001 --max-threads 1' % (imfit_path)
    if sigma_image!=None:
        imfit_line += ' --noise sigma_line.fits '

    if psf_image!=None:
        imfit_line += ' --psf psf_line.fits '        

    if mask_image!=None:
        imfit_line += ' --mask mask_line.fits '   


    subprocess.call(imfit_line, shell=True, stdout=FNULL)
    
    shutil.copy('bestfit_parameters_imfit.dat','bestfit_parameters_imfit_lm.dat')

    


def run_imf_de1(y0, sigma_image, psf_image, mask_image, ron, gain):
    y0,L00,N,Z0,chi2 = read_imf('bestfit_parameters_imfit.dat')
    crea_imf_de(L00,N,Z0,y0,ron,gain)
    
    Y00=[];L00=[];N=[];Z0=[];CHI2=[]
    for k in range(n_times):
        #print 'Number of GA run: %i' % (k+1)
        imfit_line = '%simfit vertical_line.fits -c input_imf.txt --ftol 1.0 --de --max-threads 1' % (imfit_path)
        if sigma_image!=None:
            imfit_line += ' --noise sigma_line.fits '

        if psf_image!=None:
            imfit_line += ' --psf psf_line.fits '        

        if mask_image!=None:
            imfit_line += ' --mask mask_line.fits '    

        
        subprocess.call(imfit_line, shell=True, stdout=FNULL)

        p0,p1,p2,p3,p4 = read_imf('bestfit_parameters_imfit.dat')
        Y00.append(p0)
        L00.append(p1)
        N.append(p2)
        Z0.append(p3)
        CHI2.append(p4)
        
    return Y00,L00,N,Z0,CHI2


def run_imf_de2(y0, sigma_image, psf_image, mask_image, ron, gain):
    L01,n1,z01,L02,n2,z02,chi2 = read_imf_double('bestfit_parameters_imfit.dat')
    crea_imf_de2(L01,n1,z01,L02,n2,z02,y0,ron,gain)
    
    L01=[];N1=[];Z01=[];L02=[];N2=[];Z02=[];CHI2=[]
    for k in range(n_times):
        #print 'Number of GA run: %i' % (k+1)
        imfit_line = '%simfit vertical_line.fits -c input_imf.txt --ftol 1.0 --de --max-threads 1' % (imfit_path)
        if sigma_image!=None:
            imfit_line += ' --noise sigma_line.fits '

        if psf_image!=None:
            imfit_line += ' --psf psf_line.fits '        

        if mask_image!=None:
            imfit_line += ' --mask mask_line.fits '    

        
        subprocess.call(imfit_line, shell=True, stdout=FNULL)

        p0,p1,p2,p3,p4,p5,p6 = read_imf_double('bestfit_parameters_imfit.dat') # L01,n1,z01,L02,n2,z02,chi2
        L01.append(p0)
        N1.append(p1)
        Z01.append(p2)
        L02.append(p3)
        N2.append(p4)
        Z02.append(p5)        
        CHI2.append(p6)
        
    return L01,N1,Z01,L02,N2,Z02,CHI2


def single_disc_fit(L0, n, z0, y0, sigma_image, psf_image, mask_image, ron, gain, error_est=True):
    # Levenberg-Marquardt:
    run_imf(L0, n, z0, y0, sigma_image, psf_image, mask_image, ron, gain)

    if error_est==True:
        # Genetic algorithms to estimate errors (optional):
        y0_de,l0_de,n_de,z0_de,CHI2 = run_imf_de1(y0, sigma_image, psf_image, mask_image, ron, gain)
        
        Y0,L00,N,Z0,chi2 = read_imf('bestfit_parameters_imfit_lm.dat')
        ind = CHI2.index(min(CHI2))
        L0_de = l0_de[ind]
        Y0_de = y0_de[ind]
        N_de = n_de[ind]
        Z0_de = z0_de[ind]
        L00_err = np.std(l0_de)
        N_err = np.std(n_de)
        Z0_err = np.std(z0_de)
        Y0_err = np.std(y0_de)
        
        print('\tResults LM: ',N,Z0)
        print('\tResults GA: ',N_de,'+/-',N_err,Z0_de,'+/-',Z0_err)
        if min(CHI2)<chi2:
            L00 = L00_err
            N = N_de
            Z0 = Z0_de
            Y0 = Y0_de
    else:
        Y0,L00,N,Z0,chi2 = read_imf('bestfit_parameters_imfit_lm.dat')
        L00_err = float('nan')
        N_err = float('nan')
        Z0_err = float('nan')
        Y0_err = float('nan')
        
    return L00,L00_err,Y0,Y0_err,N,N_err,Z0,Z0_err



def double_disc_fit(L0, n, z0, y0, sigma_image, psf_image, mask_image, ron, gain, error_est=True):
    # Levenberg-Marquardt:
    run_imf_double(L0, n, z0, y0, sigma_image, psf_image, mask_image, ron, gain)

    if error_est==True:
        # Genetic algorithms to estimate errors (optional):
        L001_de,n1_de,z01_de,L002_de,n2_de,z02_de,CHI2 = run_imf_de2(y0, sigma_image, psf_image, mask_image, ron, gain)
        
        L01,N1,Z01,L02,N2,Z02,chi2 = read_imf_double('bestfit_parameters_imfit_lm.dat')
        ind = CHI2.index(min(CHI2))

        L01_de = L001_de[ind]
        N1_de = n1_de[ind]
        Z01_de = z01_de[ind]
        L01_err = np.std(L001_de)
        N1_err = np.std(n1_de)
        Z01_err = np.std(z01_de)

        L02_de = L002_de[ind]
        N2_de = n2_de[ind]
        Z02_de = z02_de[ind]
        L02_err = np.std(L002_de)
        N2_err = np.std(n2_de)
        Z02_err = np.std(z02_de)
        
        print('\tResults LM: ',n1,z01,n2,z02)
        print('\tResults GA: ',N1_de,'+/-',N1_err,Z01_de,'+/-',Z01_err,N2_de,'+/-',N2_err,Z02_de,'+/-',Z02_err)
        if min(CHI2)<chi2:
            L01 = L01_de
            N1 = N1_de
            Z01 = Z01_de
            L02 = L02_de
            N2 = N2_de
            Z02 = Z02_de
    else:
        L01,N1,Z01,L02,N2,Z02,chi2 = read_imf_double('bestfit_parameters_imfit_lm.dat')
        L01_err = float('nan')
        L02_err = float('nan')
        N1_err = float('nan')
        Z01_err = float('nan')
        N2_err = float('nan')
        Z02_err = float('nan')
        
    return L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err
    



def do_each_bin_x(image_data, sigma_data, mask_data, psf_image, x_min_ind, x_max_ind, double_disc=True, error_est=True, ron=None, gain=None):
    print('SLICE %i:%i' % (x_min_ind, x_max_ind) )
    # Make average slices of image and sigma image:
    if x_min_ind>x_max_ind:
      x_min_ind,x_max_ind = x_max_ind,x_min_ind  
    

    I_image, I_mask, I_sigma = crea_vert_average_cut(image_data, mask_data, sigma_data, x_min_ind, x_max_ind)      
    #plt.plot(np.arange(len(list(I_image))), I_image) 
    #plt.show()
    #exit()
    I0 = max(I_image)
    y0 = float(np.where(I_image==I0)[0])
    y0_true = y0 #+ z_min_ind 
    
    

    # Create slice images:
    
    # Data image:
    fict_data = np.random.random((len(I_image),2))
    fict_data[:,0] = I_image
    fict_data[:,1] = np.zeros(len(I_image))
    hdu = pyfits.PrimaryHDU(fict_data)
    hdu.writeto('slice.fits',clobber=True)

    if os.path.exists('vertical_line.fits'):
        os.remove('vertical_line.fits')
    subprocess.call('fitscopy \'%s[%i:%i,%i:%i]\' vertical_line.fits' % ('slice.fits',1,1,1,len(I_image)), shell=True)
    os.remove('slice.fits')
    
    
    # Sigma image:
    fict_data = np.random.random((len(I_sigma),2))
    fict_data[:,0] = I_sigma
    fict_data[:,1] = np.zeros(len(I_sigma))
    hdu = pyfits.PrimaryHDU(fict_data)
    hdu.writeto('slice_sigma.fits',clobber=True)      

    if os.path.exists('sigma_line.fits'):
        os.remove('sigma_line.fits')
    subprocess.call('fitscopy \'%s[%i:%i,%i:%i]\' sigma_line.fits' % ('slice_sigma.fits',1,1,1,len(I_image)), shell=True)
    os.remove('slice_sigma.fits')
    sigma_image = 'sigma_line.fits'
    

    # Mask image:
    fict_data = np.random.random((len(I_mask),2))
    fict_data[:,0] = I_mask
    fict_data[:,1] = np.zeros(len(I_mask))
    hdu = pyfits.PrimaryHDU(fict_data)
    hdu.writeto('slice_mask.fits',clobber=True)      

    if os.path.exists('mask_line.fits'):
        os.remove('mask_line.fits')
    subprocess.call('fitscopy \'%s[%i:%i,%i:%i]\' mask_line.fits' % ('slice_mask.fits',1,1,1,len(I_image)), shell=True)
    os.remove('slice_mask.fits')
    mask_image = 'mask_line.fits'



    # Find first guess on the parameters:
    popt_soph, pcov_soph = curve_fit( disk_edge_soph, range(1, len(I_image)+1, 1), I_image, p0=(I0, len(I_image)/4., 0.) )
    L0 = popt_soph[0]
    z0 = popt_soph[1]
    y0 = popt_soph[2]


    
    # Single disc fit:
    L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err = single_disc_fit(L0, 1., z0, y0, sigma_image, psf_image, mask_image, ron, gain, error_est)


    
    if double_disc:
        # Double disc fit:
        L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err = double_disc_fit(L0, N, Z0, Y0, sigma_image, psf_image, mask_image, ron, gain, error_est)   
    else:
        L01 = float('nan')
        L01_err = float('nan')
        N1 = float('nan')
        N1_err = float('nan')
        Z01 = float('nan')
        Z01_err = float('nan')
        L02 = float('nan')
        L02_err = float('nan')
        N2 = float('nan')
        N2_err = float('nan')
        Z02 = float('nan')
        Z02_err = float('nan')
      

    for file in ['bestfit_parameters_imfit.dat',
                 'modelimage.fits',
                 'vertical_line.fits',
                 'bestfit_parameters_imfit_lm.dat',
                 'sigma_line.fits']:
      if os.path.exists(file):
        os.remove(file)

    return [L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err],[L01,L01_err,N1,N1_err,Z01,Z01_err],[L02,L02_err,N2,N2_err,Z02,Z02_err] 



























def define_xbins(rmin, rmax, Rmax, bin_type='log'):
    if rmax is None:
        rmax = Rmax

    if bin_type=='log':
        return Rmax*np.power(10, np.arange(math.log10(1./Rmax), math.log10(rmax/Rmax), 0.1))
    if bin_type=='linear':
        return Rmax*np.arange(rmin/Rmax, rmax/Rmax, 0.1)        
    if bin_type=='pix':
        return np.arange(rmin, rmax, 10)     


def define_ybins(ymin, ymax, Zmax, bin_type='log'):
    # TODO  
    z=1





def plot_single_disc_z(R, L0_single, L0_single_err, Y0_single, Y0_single_err, N_single, N_single_err, Z0_single, Z0_single_err, error_est):
    plt.figure(0, figsize=(5, 5))
    plt.xlabel(r' $R$ (arcsec) ', fontsize=15)
    plt.ylabel(r' $z_0$ (arcsec) ', fontsize=15)
    
    if error_est==True:
      plt.errorbar(R, Z0_single, yerr=Z0_single_err, fmt='o',markersize=7,color='black',markeredgecolor='black', ecolor='black', capthick=2)
    else:
      plt.plot(R, Z0_single, 'o',markersize=7,color='black',markeredgecolor='black')
    plt.ylim(0, 2.*np.median(Z0_single))
    #if Rbreak!=0. and not math.isnan(Rbreak):
    #  plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')
    plt.savefig('single_z0_r.eps', transparent = False, dpi=300)
    plt.clf()
    plt.close()  
     

    plt.figure(1, figsize=(5, 5))
    plt.xlabel(r' $R$ (arcsec) ', fontsize=15)
    plt.ylabel(r' $n$ ', fontsize=15)
    
    if error_est==True:
      plt.errorbar(R, N_single, yerr=N_single_err, fmt='o',markersize=7,color='black',markeredgecolor='black', ecolor='black', capthick=2)
    else:
      plt.plot(R, N_single, 'o',markersize=7,color='black',markeredgecolor='black')
    plt.ylim(0, 110)
    #if Rbreak!=0. and not math.isnan(Rbreak):
    #  plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')
    plt.savefig('single_n_r.eps', transparent = False, dpi=300)
    plt.clf()
    plt.close()  
    

    plt.figure(2, figsize=(5, 5))
    plt.xlabel(r' $R$ (arcsec) ', fontsize=15)
    plt.ylabel(r' $y0 (arcsec)$ ', fontsize=15)
    
    if error_est==True:
      plt.errorbar(R, Y0_single, yerr=Y0_single_err, fmt='o',markersize=7,color='black',markeredgecolor='black', ecolor='black', capthick=2)
    else:
      plt.plot(R, Y0_single, 'o',markersize=7,color='black',markeredgecolor='black')
    plt.ylim(0, 2.*np.median(Y0_single))
    #if Rbreak!=0. and not math.isnan(Rbreak):
    #  plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')
    plt.savefig('single_y0_r.eps', transparent = False, dpi=300)
    plt.clf()
    plt.close() 

    
def plot_double_disc_z(x_bins, L0_1, L0err_1, N_1, Nerr_1, Z0_1, Z0err_1, L0_2, L0err_2, N_2, Nerr_2, Z0_2, Z0err_2):
    z=1




def ds9_to_np(x):
    '''
    Function to convert ds9 pixel coordinates to numpy format, i.e
    0.5 pix is 0, 1.5 pix is 1
    '''
    return int(round(x)) - 1


def create_add_mask(nx, ny, R_in_l, R_out_l, R_in_r, R_out_r, z_in_dn, z_out_dn, z_in_up, z_out_up):
    xc = nx/2.
    yc = ny/2.
    
    # define numpy indices:
    if R_out_l is None:
        x_min_ind_l = 0
    else:
        x_min_ind_l = ds9_to_np(xc - R_out_l)
        
    x_max_ind_l = ds9_to_np(xc - R_in_l)        


    if R_out_r is None:
        x_max_ind_r = nx
    else:
        x_max_ind_r = ds9_to_np(xc + R_out_r)
        
    x_min_ind_r = ds9_to_np(xc + R_in_r)   


    if z_out_dn is None:
        z_min_ind_dn = 0
    else:
        z_min_ind_dn = ds9_to_np(yc - z_out_dn)
        
    z_max_ind_dn = ds9_to_np(yc - z_in_dn)        


    if z_out_up is None:
        z_max_ind_up = ny
    else:
        z_max_ind_up = ds9_to_np(yc + z_out_up)
        
    z_min_ind_up = ds9_to_np(yc + z_in_up)
    
    add_mask_data = np.ones((ny,nx))
    
    for k in range(ny):
        for i in range(nx):
            if (k >= z_min_ind_dn and k<=z_max_ind_dn) or (k>=z_min_ind_up and k<=z_max_ind_up):
                if (i>=x_min_ind_l and i<=x_max_ind_l) or (i>=x_min_ind_r and i<=x_max_ind_r):
                    add_mask_data[k,i]=0
    
    hdu = pyfits.PrimaryHDU(add_mask_data)
    hdu.writeto('add_mask_tmp.fits', clobber=True)
    return x_min_ind_l,x_max_ind_l,x_min_ind_r,x_max_ind_r,z_min_ind_dn,z_max_ind_dn,z_min_ind_up,z_max_ind_up
    



def main(input_image, mask_image=None, weight_image=None, psf_image=None, R_in_l=0., R_out_l=None, R_in_r=0., R_out_r=None, z_in_dn = 0., z_out_dn=None, z_in_up=0., z_out_up=None, mode='vertical', bin_type='linear'):
    # The input images should be rotated, cropped and centered so that the major axis is horizontal

    # Read in input files
    hdu = pyfits.open(input_image)
    image_data = hdu[0].data 
    ny,nx =  np.shape(image_data)
    xc = nx/2.
    yc = ny/2.



    ##### MASK ########
    # Create additional mask:
    x_min_ind_l,x_max_ind_l,x_min_ind_r,x_max_ind_r,z_min_ind_dn,z_max_ind_dn,z_min_ind_up,z_max_ind_up = create_add_mask(nx, ny, R_in_l, R_out_l, R_in_r, R_out_r, z_in_dn, z_out_dn, z_in_up, z_out_up)

    if mask_image is None:
        shutil.copy('add_mask_tmp.fits', 'mask_tmp.fits')
    else:
        # Merge all masks
        main([mask_image, 'add_mask_tmp.fits'], 'mask_tmp.fits')
    hdu = pyfits.open('mask_tmp.fits')
    mask_data = hdu[0].data
    mask_data = (mask_data>0)
    ##### MASK END ########        




    ##### PSF ########    
    # Create slice of PSF:
    if psf_image is not None:
        if mode=='vertical':
            psf_slice(psf_image, axis='yaxis')
            psf_image_cut = 'psf_yline.fits'
        else:
            psf_slice(psf_image, axis='xaxis')
            psf_image_cut = 'psf_xline.fits'            
    else:
        psf_image_cut = None
    ##### PSF END ########     





    ##### WEIGHT IMAGE ######## 
    if weight_image is not None:
        hdu = pyfits.open(weight_image)
        sigma_data = hdu[0].data
    else: 
        sigma_data = np.ones((ny,nx))
    ##### WEIGHT IMAGE END ######## 




    
    if mode=='vertical':
        # Left part:
        print('Vertical slices fitting: LEFT')
        r_bins = define_xbins(R_in_l, R_out_l, nx/2., bin_type)

        
        ind_x_bins = []
        for k in range(len(r_bins)):
            ind_x_bins.append( ds9_to_np(xc-r_bins[k]) )

    
        f = open('vertical_fits_left.dat', 'w')
        for k in range(0, len(ind_x_bins)-1):
            [L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err],[L01,L01_err,N1,N1_err,Z01,Z01_err],[L02,L02_err,N2,N2_err,Z02,Z02_err] = do_each_bin_x(image_data, sigma_data, mask_data, psf_image_cut, ind_x_bins[k], ind_x_bins[k+1], double_disc=False, error_est=False)

            R = (r_bins[k]+r_bins[k+1])/2.
            f.write('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n'   % (R,L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err,L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err ) )
        f.close()    
        #exit()
        
        # Right part:
        print('\nVertical slices fitting: RIGHT')
        r_bins = define_xbins(R_in_r, R_out_r, nx/2., bin_type)
        ind_x_bins = []
        for k in range(len(r_bins)):
            ind_x_bins.append( ds9_to_np(xc+r_bins[k]) )

        f = open('vertical_fits_right.dat', 'w')
        for k in range(0, len(ind_x_bins)-1):
            [L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err],[L01,L01_err,N1,N1_err,Z01,Z01_err],[L02,L02_err,N2,N2_err,Z02,Z02_err] = do_each_bin_x(image_data, sigma_data, mask_data, psf_image_cut, ind_x_bins[k], ind_x_bins[k+1], double_disc=False, error_est=True)

            R = (r_bins[k]+r_bins[k+1])/2.
            f.write('%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n'   % (R,L0,L0_err,Y0,Y0_err,N,N_err,Z0,Z0_err,L01,L01_err,N1,N1_err,Z01,Z01_err,L02,L02_err,N2,N2_err,Z02,Z02_err ) )    
        f.close() 






































    '''
    L0_single = []; Y0_single = []; N_single = []; Z0_single = []
    L0_single_err = []; Y0_single_err = []; N_single_err = []; Z0_single_err = []
    L0_1 = []; N_1 = []; Z0_1= [];
    L0err_1 = []; Nerr_1 = []; Z0err_1 = []
    L0_2 = []; N_2 = []; Z0_2 = []
    L0err_2 = []; Nerr_2 = []; Z0err_2 = []


    L0_single.append(L0)
    L0_single_err.append(L0_err)
    Y0_single.append(Y0)
    Y0_single_err.append(Y0_err)
    N_single.append(N)
    N_single_err.append(N_err)
    Z0_single.append(Z0)
    Z0_single_err.append(Z0_err)
        
    L0_1.append(L01)
    L0err_1.append(L01_err)
    N_1.append(N1)
    Nerr_1.append(N1_err)
    Z0_1.append(Z01)
    Z0err_1.append(Z01_err)
        
    L0_2.append(L02)
    L0err_2.append(L02_err)
    N_2.append(N2)
    Nerr_2.append(N2_err)
    Z0_2.append(Z02)
    Z0err_2.append(Z02_err)

    
    # TODO:
    plot_single_disc_z(x_bins, L0_single, L0_single_err, Y0_single, Y0_single_err, N_single, N_single_err, Z0_single, Z0_single_err)
    plot_double_disc_z(x_bins, L0_1, L0err_1, N_1, Nerr_1, Z0_1, Z0err_1, L0_2, L0err_2, N_2, Nerr_2, Z0_2, Z0err_2)
    '''


    

    
    
input_image = '/home/amosenko/MyCurrentWork/Edge_on_HERON/Test_eon/test_cut_fitting/galaxy1.fits'
main(input_image, mask_image=None, weight_image=None, psf_image=None, R_in_l=0., R_out_l=None, R_in_r=0., R_out_r=None, z_in_dn = 0., z_out_dn=None, z_in_up=0., z_out_up=None, mode='vertical', bin_type='linear')    


























'''
    plt.clf()
    plt.close()        
    os.remove('psf_line.fits')

    
    
    plt.figure(0,figsize=(5, 5))
    plt.xlabel(r' R/R$_{max}$ ', fontsize=15)
    plt.ylabel(r' z$_0$ (pixels) ', fontsize=15)
    if error_est==True:
      plt.errorbar(R/radius_x_out, Z001, yerr=Z001_err,fmt='o',markersize=9,color='blue',markeredgecolor='blue', ecolor='black', capthick=2)
      plt.errorbar(R/radius_x_out, Z002, yerr=Z002_err,fmt='^',markersize=7,color='red',markeredgecolor='red', ecolor='black', capthick=2)
    else:
      plt.plot(R/radius_x_out,Z001,'o',markersize=9,color='blue',markeredgecolor='blue')
      plt.plot(R/radius_x_out,Z002,'^',markersize=7,color='red',markeredgecolor='red')
    plt.ylim(0,2.*max([median(Z002),median(Z001)]))
    if Rbreak!=0. and not math.isnan(Rbreak):
      plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')
    plt.savefig('imf_z0_r.eps', transparent = False, dpi=300)
    plt.clf()
    plt.close()  
    
    plt.figure(1,figsize=(5.5, 5))
    plt.xlabel(r' R/R$_{max}$ ', fontsize=15)
    plt.ylabel(r' n ', fontsize=15)
    plt.yscale('log')
    if error_est==True:
      plt.errorbar(R/radius_x_out,NN1, yerr=NN1_err,fmt='o',markersize=9,color='blue',markeredgecolor='blue', ecolor='black', capthick=2)
      plt.errorbar(R/radius_x_out,NN2, yerr=NN2_err,fmt='^',markersize=7,color='red',markeredgecolor='red', ecolor='black', capthick=2)      
    else:
      plt.plot(R/radius_x_out,NN1,'o',markersize=9,color='blue',markeredgecolor='blue')
      plt.plot(R/radius_x_out,NN2,'^',markersize=7,color='red',markeredgecolor='red')
    plt.ylim(0.1,110)
    if Rbreak!=0. and not math.isnan(Rbreak):
      plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')
    plt.savefig('imf_n_r.eps', transparent = False, dpi=300)
    plt.clf()
    plt.close()      
    
    print 'Done!'
    PARS["Z_CUTS_IMFIT"] = [median(Z001),std(Z001),median(NN1),std(NN1),median(Z002),std(Z002),median(NN2),std(NN2)]
    os.remove('input_imf.txt')

    if error_est==True:
      return R, NN1, NN1_err, Z001, Z001_err, NN2, NN2_err, Z002, Z002_err, PARS   
    else:
      return R,NN1,Z001,NN2,Z002,PARS  




























































	
      os.remove('bestfit_parameters_imfit.dat')
      if os.path.exists('bestfit_parameters_imfit_lm.dat'):
	os.remove('bestfit_parameters_imfit_lm.dat')
      os.remove('modelimage.fits')
      os.remove('vertical_line.fits')
      if sigma_image!=None:
	os.remove('sigma_line.fits')








    plt.clf()
    plt.close()  
    plt.figure(0,figsize=(5, 5))
    plt.xlabel(r' R/R$_{max}$ ', fontsize=15)
    plt.ylabel(r' z$_0$ (pixels) ', fontsize=15)
    if error_est==True:
      plt.errorbar(R/radius_x_out, Z00, yerr=Z00_err,fmt='^',markersize=7,color='red',markeredgecolor='red', ecolor='black', capthick=2)
    else:
      plt.plot(R/radius_x_out,Z00,'^',markersize=7,color='red',markeredgecolor='red')
    plt.ylim(0,2.*median(Z00))
    if Rbreak!=0. and not math.isnan(Rbreak):
      plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')
    plt.savefig('imf_z0_r.eps', transparent = False, dpi=300)
    plt.clf()
    plt.close()  

    plt.clf()
    plt.close()      
    plt.figure(1,figsize=(5.5, 5))
    plt.xlabel(r' R/R$_{max}$ ', fontsize=15)
    plt.yscale('log')
    plt.ylabel(r' n ', fontsize=15)
    if error_est==True:
      plt.errorbar(R/radius_x_out, NN, yerr=NN_err,fmt='^',markersize=7,color='red',markeredgecolor='red', ecolor='black', capthick=2)
    else:
      plt.plot(R/radius_x_out,NN,'^',markersize=7,color='red',markeredgecolor='red')
    plt.ylim(0.1,110)

    if Rbreak!=0. and not math.isnan(Rbreak):
      plt.axvline(x=Rbreak/radius_x_out,color='black',linewidth=2,ls='--')

    plt.savefig('imf_n_r.eps', transparent = False, dpi=300)
    plt.clf()
    plt.close()  
      
    os.remove('psf_line.fits')
    PARS["Z_CUTS_IMFIT"] = [median(Z00),std(Z00),median(NN),std(NN),float(nan),float(nan),float(nan),float(nan)]
    if error_est==True:
      return R, np.empty(len(R)) * np.nan, np.empty(len(R)) * np.nan, np.empty(len(R)) * np.nan, np.empty(len(R)) * np.nan, NN, NN_err, Z00, Z00_err, PARS
    else:
      return R, np.empty(len(R)) * np.nan, np.empty(len(R)) * np.nan, NN, Z00, PARS
'''
        