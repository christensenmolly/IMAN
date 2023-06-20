#!/usr/bin/env python
# Module to create an initial model
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **        © Astronomical Observatory, Ghent University         **
# *****************************************************************

# Import standard modules
import pylab
import sys
import os
import shutil
import math
import numpy as np
import scipy as sp
from numpy import *
from pylab import *
import subprocess
import argparse
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import optimize
from scipy.odr.odrpack import *
from scipy import special

# Import astro modules
import pyfits
from photutils import datasets
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import EllipticalAperture
from astropy.table import hstack

# Import DECA modules
import deca_setup
import radial_profile
import do_mode

tmp_out = sys.stdout

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0].split('DECA')[0]

#sys.path.append(PATH_TO_PACKAGE+'Decomposition')

# -----------------------------------------------------------------
# OLD FUNCTION TO CREATE CONSTRAINTS FILE
def create_constrain_file(reb, model):
	f = open('constraint.txt', 'w') 
	sys.stdout = f
	
	if model == 'sersic':
            print '    1	  	   n	     0.3  to  15'
	
	#print '    1	  	   x	     %i  %i'  % (-int(max([reb/10.,2.])),int(max([reb/10.,2.])))
	#print '    1	  	   y	     %i  %i'  % (-int(max([reb/10.,2.])),int(max([reb/10.,2.])))

        if '+' in model:
         if model == 'sersic+exp_disc' or model == 'sersic+eon_disc':
            # More than one component
	    #print '    1	  	   re	     0.1  to  1500'
	    #print '    2	  	   re	     0.1  to  1500'
	    #print '    1	  	   mag	     5  to  35'
	    #print '    2	  	   mag	     5  to  35'
	    print '    1	  	   n	     0.3  to  15'


	    print '    1/2	  	   x	     1.0  to  1.0'
	    print '    1/2	  	   y	     1.0  to  1.0'

	    #print '    1-2	  	   pa	     0.0  0.0'
	    #print '    2/1	  	   re	     0.1  to  10.0'
         elif model == 'sersic+sersic':
	    print '    1	  	   n	     0.3  to  15'
	    print '    2	  	   n	     0.3  to  15'
	    print '    1-2	  	   x	     0.0  to  0.0'
	    print '    1-2	  	   y	     0.0  to  0.0'             
             
	sys.stdout = tmp_out
	f.close()
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO FIND MINIMUM X-VALUE ON THE PROFILE, WITH Y-VALUE LARGER THAN 
def find_lowest_val(x,y,y_value):
    X = []
    for k in range(len(x)):
        if y[k]<=y_value:
            X.append(x[k])
    return min(X)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# PIECEWISE LINEAR FUNCTION
def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# CHI2 FINCTION
def chi2_func(ref,mod):
  return sum( (ref-mod)**2)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
#  FUNCTION TO DEFINE BULGE-DOMINATED REGION
def define_bulge_dom_region(r,I,m0,pix2sec):
	I = m0 - 2.5*log10(np.array(I, float)) + 5.*log10(pix2sec)
	r = np.array(r)*pix2sec
	#plt.plot(r,I,'o')
	#plt.show()
	
	Res_chi2 = []
	Res_p = []

	for rr in r:
		  ind = list(r).index(rr)
		  p , e = optimize.curve_fit(piecewise_linear, r, I,[rr,I[ind],1,1])
		  chi2 = chi2_func(I,piecewise_linear(r, *p))
		  Res_chi2.append(chi2)
		  Res_p.append(p)

	results_pos = Res_p[Res_chi2.index(min(Res_chi2))]
	R_bulge_ma = fabs(results_pos[0])

        if R_bulge_ma>max(r)/2.:
            r_new = []
            I_new = []
            for k in range(len(r)):
                if r[k]<R_bulge_ma:
                    r_new.append(r[k])
                    I_new.append(I[k])
            
            r = np.array(r_new)
            I = np.array(I_new)

            Res_chi2 = []
            Res_p = []

            for rr in r:
                    ind = list(r).index(rr)
                    p , e = optimize.curve_fit(piecewise_linear, r, I,[rr,I[ind],1,1])
                    chi2 = chi2_func(I,piecewise_linear(r, *p))
                    Res_chi2.append(chi2)
                    Res_p.append(p)

            results_pos = Res_p[Res_chi2.index(min(Res_chi2))]
            R_bulge_ma = fabs(results_pos[0])

        
        return R_bulge_ma/pix2sec
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION OF A LINE
def line(B, x):
	#*** For non edge-on disk SB in mag/arcsec^2 (radial) ***
	# B[0]=h
	# B[1]=m0d
	x = np.array(x)
	return (1.0857/B[0])*fabs(x) + B[1]
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# SERSIC FUNCTION
def bulge_func(B, x):
	#*** For Sersic SB distribution ***
	# B[0]=Ieb
	# B[1]=reb 
	# B[2]=n
	return B[0]*np.exp( -(1.9992*B[2] - 0.3271)*((fabs(x)/B[1])**(1./B[2]) - 1.) )
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# SERSIC FUNCTION
def bulge_func_fix(B, x):
	#*** For Sersic SB distribution ***
	# B[0]=Ieb
	# B[1]=reb
	# n=2.5
	return B[0]*np.exp( -(1.9992*2.5 - 0.3271)*((fabs(x)/B[1])**(1./2.5) - 1.) )
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# SERSIC TOTAL MAGNITUDE
def magBulge_f(reb,meb,n,q):
	nu = 1.9987*n-0.3267
	An = 2.5*log10(2.*math.pi*n/(nu**(2.*n)) * special.gamma(2.0*n)) + 1.0857*nu
	return - 5.*log10(reb) + meb - An - 2.5*log10(q)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# SERSIC EFFECTIVE SB
def meb_bulge_f(reb,mag_bulge,n,q):
	nu = 1.9987*n-0.3267
	An = 2.5*log10(2.*math.pi*n/(nu**(2.*n)) * special.gamma(2.0*n)) + 1.0857*nu
	return mag_bulge+5.*log10(reb) + An + 2.5*log10(q)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# DISC TOTAL MAGNITUDE
def magDisk_f(m0d,h,q):
	return -2.5*log10(2.*math.pi) + m0d - 5.*log10(h) + 2.5*log10(1./q)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# DISC CENTRAL SB
def m0d_disk_f(mag_disc,h,q):
	return mag_disc + 2.5*log10(2.*math.pi) + 5.*log10(h) - 2.5*log10(1./q)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# EDGE-ON DISC TOTAL MAGNITUDE
def magDisk_edge_on_f(m0d,h,z0):
	return -2.5*log10(2.*math.pi) + m0d - 2.5*log10(z0*h)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# EDGE-ON DISC SB
def m0_disc_edge_on_f(mag_disc,h,z0):
	return +2.5*log10(2.*math.pi) + mag_disc + 2.5*log10(z0*h)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO FIT EXPONENTIAL DISC WITH A LINE IN LOG VALUES
def fit_exp_disc(r,Int,m0,pix2sec):
        r = r*pix2sec
        mag = m0-2.5*log10(Int)+5.*log10(pix2sec)

	dataToFit = RealData(r, mag)
	mo = Model(line)
	fitt = ODR(dataToFit, mo, [max(r)/4.,min(mag)])
	fitt.set_job()
	fit_result = fitt.run()
	h = fit_result.beta[0]
	m0d = fit_result.beta[1]
	return h/pix2sec,m0d   
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO FIT BULGE WITH THE SERSIC FUNCTION
def fit_ser_bulge(r,Int,m0,pix2sec):
	dataToFit = RealData(r, Int)
	mo = Model(bulge_func_fix)
	fitt = ODR(dataToFit, mo, [min(Int),max(r)/3.])
	fitt.set_job()
	fit_result = fitt.run()
	Ieb = fit_result.beta[0] / 1.5
	reb = fit_result.beta[1] / 1.5
	return reb,m0-2.5*log10(Ieb)+5.*log10(pix2sec)
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# OLD FUNCTION TO FIND INITIAL GUESS FROM SCRATCH
def initial_guess_from_scratch_old(PARS,scale,model):
 
  # From Vika et al.(2014)
  meb = PARS["First_Ser_meb"]
  Mb = PARS["First_Ser_Mb"]
  reb = PARS["First_Ser_reb"]
  n = PARS["First_Ser_n"]
  PA = PARS["First_Ser_PA"]
  q = PARS["First_Ser_q"]

  if (PARS["First_Ser_q"] > deca_setup.q_eon and model!='sersic+sersic') or model=='sersic+exp_disc':
      PARS["First_Disk_Md"] = Mb + 0.65      
      PARS["First_Disk_hd"] = reb / 1.68
      PARS["First_Disk_PA"] = PA
      PARS["First_Disk_q"] = q
      PARS["First_Disk_m0d"] = m0d_disk_f(PARS["First_Disk_Md"],PARS["First_Disk_hd"]*scale,PARS["First_Disk_q"])
    
      PARS["First_Ser_Mb"] = Mb + 0.75
      PARS["First_Ser_reb"] = 0.1 * reb  
      PARS["First_Ser_n"] = 2.
      PARS["First_Ser_q"] = 0.8
      PARS["First_Ser_PA"] = PA
      PARS["First_Ser_meb"] = meb_bulge_f(PARS["First_Ser_reb"]*scale,PARS["First_Ser_Mb"],PARS["First_Ser_n"],PARS["First_Ser_q"])

  if  model=='exp_disc':
      PARS["First_Disk_Md"] = Mb 
      PARS["First_Disk_hd"] = reb / 1.68
      PARS["First_Disk_PA"] = PA
      PARS["First_Disk_q"] = q
      PARS["First_Disk_m0d"] = m0d_disk_f(PARS["First_Disk_Md"],PARS["First_Disk_hd"]*scale,PARS["First_Disk_q"])
    
      PARS["First_Ser_meb"] = float(nan)
      PARS["First_Ser_Mb"] = float(nan)
      PARS["First_Ser_reb"] = float(nan)
      PARS["First_Ser_n"] = float(nan)
      PARS["First_Ser_PA"] = float(nan)  
      PARS["First_Ser_q"] = float(nan) 
  
  if  (PARS["First_Ser_q"] <= deca_setup.q_eon and model!='sersic+sersic') or model=='sersic+eon_disc':
      PARS["First_Disk_Md"] = Mb + 0.65
      PARS["First_Disk_hd"] = reb / 1.68
      PARS["First_Disk_PA"] = PA
      PARS["First_Disk_q"] = float(nan)
      PARS["First_Disk_z0"] = PARS["First_Disk_hd"] / 4.      
      PARS["First_Disk_m0d"] = Mb + 0.65 + 2.5*log10(2.*math.pi*PARS["First_Disk_hd"]*scale*PARS["First_Disk_z0"]*scale) + 5.*log10(scale)
    
      PARS["First_Ser_Mb"] = Mb + 0.75
      PARS["First_Ser_reb"] = 0.1 * reb  
      PARS["First_Ser_n"] = 2.
      PARS["First_Ser_q"] = 0.8
      PARS["First_Ser_PA"] = 10.        

      
  if  model=='eon_disc':
      PARS["First_Disk_Md"] = Mb
      PARS["First_Disk_hd"] = reb / 1.68
      PARS["First_Disk_PA"] = PA
      PARS["First_Disk_q"] = float(nan)
      PARS["First_Disk_z0"] = PARS["First_Disk_hd"] / 4.      
      PARS["First_Disk_m0d"] = m0_disc_edge_on_f(PARS["First_Disk_Md"],PARS["First_Disk_hd"]*scale,PARS["First_Disk_z0"]*scale)
    
      PARS["First_Ser_meb"] = float(nan)
      PARS["First_Ser_Mb"] = float(nan)
      PARS["First_Ser_reb"] = float(nan)
      PARS["First_Ser_n"] = float(nan)
      PARS["First_Ser_PA"] = float(nan)  
      PARS["First_Ser_q"] = float(nan)      

  if  model=='agn+exp_disc':
      PARS["First_Disk_Md"] = Mb 
      PARS["First_Disk_hd"] = reb / 1.68
      PARS["First_Disk_PA"] = PA
      PARS["First_Disk_q"] = q
      PARS["First_Disk_m0d"] = m0d_disk_f(PARS["First_Disk_Md"],PARS["First_Disk_hd"]*scale,PARS["First_Disk_q"])
      
      PARS["PSF_M"] = 100. - 2.5*log10( (10**(0.4*(100.-Mb))/100.) )
    
      PARS["First_Ser_meb"] = float(nan)
      PARS["First_Ser_Mb"] = float(nan)
      PARS["First_Ser_reb"] = float(nan)
      PARS["First_Ser_n"] = float(nan)
      PARS["First_Ser_PA"] = float(nan)  
      PARS["First_Ser_q"] = float(nan) 


  if  model=='agn+eon_disc':
      PARS["First_Disk_Md"] = Mb
      PARS["First_Disk_hd"] = reb / 1.68
      PARS["First_Disk_PA"] = PA
      PARS["First_Disk_q"] = float(nan)
      PARS["First_Disk_z0"] = PARS["First_Disk_hd"] / 4.      
      PARS["First_Disk_m0d"] = m0_disc_edge_on_f(PARS["First_Disk_Md"],PARS["First_Disk_hd"]*scale,PARS["First_Disk_z0"]*scale)
      
      PARS["PSF_M"] = 100. - 2.5*log10( (10**(0.4*(100.-Mb))/100.) )
    
      PARS["First_Ser_meb"] = float(nan)
      PARS["First_Ser_Mb"] = float(nan)
      PARS["First_Ser_reb"] = float(nan)
      PARS["First_Ser_n"] = float(nan)
      PARS["First_Ser_PA"] = float(nan)  
      PARS["First_Ser_q"] = float(nan)

  if  model=='sersic+sersic':   
      PARS["First_Ser_Mb"] = Mb + 0.75
      PARS["First_Ser_reb"] = 0.1 * reb  
      PARS["First_Ser_n"] = 2.
      PARS["First_Ser_q"] = 0.8
      PARS["First_Ser_PA"] = PA
      PARS["First_Ser_meb"] = meb_bulge_f(PARS["First_Ser_reb"]*scale,PARS["First_Ser_Mb"],PARS["First_Ser_n"],PARS["First_Ser_q"])      

      PARS["Second_Ser_Mb"] = Mb + 0.65
      PARS["Second_Ser_reb"] = reb
      PARS["Second_Ser_n"] = 1.
      PARS["Second_Ser_q"] = q
      PARS["Second_Ser_PA"] = PA
      PARS["Second_Ser_meb"] = meb_bulge_f(PARS["Second_Ser_reb"]*scale,PARS["Second_Ser_Mb"],PARS["Second_Ser_n"],PARS["Second_Ser_q"])  
      
  create_constrain_file(reb)
  return PARS
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO FIND INITIAL GUESS FROM SCRATCH 
def initial_guess_from_scratch(PARS,scale,model):
  
  # From Vika et al.(2014)
  
  for k in range(len(PARS["First_Ser_meb"])):
            meb = PARS["First_Ser_meb"][k]
            Mb = PARS["First_Ser_Mb"][k]
            reb = PARS["First_Ser_reb"][k]
            n = PARS["First_Ser_n"][k]
            q = PARS["First_Ser_q"][k]
            PA = PARS["First_Ser_PA"][k]
            COEFFS = deca_setup.COEFFS

            if  model[k]=='exp_disc':
                PARS["First_Disk_Md"][k] = float(eval(str(Mb)+COEFFS['exp_disc:exp_disc'][0])) 
                PARS["First_Disk_hd"][k] = float(eval(str(reb)+COEFFS['exp_disc:exp_disc'][1])) 
                PARS["First_Disk_q"][k] = float(eval(str(q)+COEFFS['exp_disc:exp_disc'][2])) 
                PARS["First_Disk_PA"][k] = float(eval(str(PA)+COEFFS['exp_disc:exp_disc'][3])) 
                PARS["First_Disk_m0d"][k] = m0d_disk_f(PARS["First_Disk_Md"][k],PARS["First_Disk_hd"][k]*scale,PARS["First_Disk_q"][k])
                
                PARS["First_Ser_meb"][k] = float(nan)
                PARS["First_Ser_Mb"][k] = float(nan)
                PARS["First_Ser_reb"][k] = float(nan)
                PARS["First_Ser_n"][k] = float(nan)
                PARS["First_Ser_q"][k] = float(nan) 
                PARS["First_Ser_PA"][k] = float(nan) 
                
                
            if  model[k]=='eon_disc':
                PARS["First_Disk_Md"][k] = float(eval(str(Mb)+COEFFS['eon_disc:eon_disc'][0])) 
                PARS["First_Disk_hd"][k] = float(eval(str(reb)+COEFFS['eon_disc:eon_disc'][1])) 
                PARS["First_Disk_q"][k] = float(nan)
                PARS["First_Disk_PA"][k] = float(eval(str(PA)+COEFFS['eon_disc:eon_disc'][2])) 
                PARS["First_Disk_z0"][k] = PARS["First_Disk_hd"] / 4.      
                PARS["First_Disk_m0d"][k] = m0_disc_edge_on_f(PARS["First_Disk_Md"][k],PARS["First_Disk_hd"][k]*scale,PARS["First_Disk_z0"][k]*scale)
                
                PARS["First_Ser_meb"][k] = float(nan)
                PARS["First_Ser_Mb"][k] = float(nan)
                PARS["First_Ser_reb"][k] = float(nan)
                PARS["First_Ser_n"][k] = float(nan)
                PARS["First_Ser_q"][k] = float(nan) 
                PARS["First_Ser_PA"][k] = float(nan)  

            
            if  model[k]=='ring':
                PARS["Ring_m0"][k] = float(eval(str(meb)+COEFFS['ring:ring'][0])) 
                PARS["Ring_sigma_r"][k] = float(eval(str(reb)+COEFFS['ring:ring'][1])) 
                PARS["Ring_R_r"][k] = float(eval(str(reb)+COEFFS['ring:ring'][2]))
                PARS["Ring_q"][k] = float(eval(str(q)+COEFFS['ring:ring'][3])) 
                PARS["Ring_PA"][k] = float(eval(str(PA)+COEFFS['ring:ring'][4]))       
                
                PARS["First_Ser_meb"][k] = float(nan)
                PARS["First_Ser_Mb"][k] = float(nan)
                PARS["First_Ser_reb"][k] = float(nan)
                PARS["First_Ser_n"][k] = float(nan)
                PARS["First_Ser_q"][k] = float(nan) 
                PARS["First_Ser_PA"][k] = float(nan)
                

            if  model[k]=='sersic+sersic':   
                PARS["First_Ser_Mb"][k] = float(eval(str(Mb)+COEFFS['sersic+sersic:sersic1'][0])) 
                PARS["First_Ser_reb"][k] = float(eval(str(reb)+COEFFS['sersic+sersic:sersic1'][1]))   
                PARS["First_Ser_n"][k] = float(COEFFS['sersic+sersic:sersic1'][2])
                PARS["First_Ser_q"][k] = float(COEFFS['sersic+sersic:sersic1'][3])   
                PARS["First_Ser_PA"][k] = float(eval(str(PA)+COEFFS['sersic+sersic:sersic1'][4]))   
                PARS["First_Ser_meb"][k] = meb_bulge_f(PARS["First_Ser_reb"][k]*scale,PARS["First_Ser_Mb"][k],PARS["First_Ser_n"][k],PARS["First_Ser_q"][k])      
                #print PARS["First_Ser_Mb"], PARS["First_Ser_meb"]
                #meb_bulge_f(reb,mag_bulge,n,q)
                #print PARS["First_Ser_reb"]*scale,PARS["First_Ser_Mb"],PARS["First_Ser_n"],PARS["First_Ser_q"]
                #exit()

                PARS["Second_Ser_Mb"][k] = float(eval(str(Mb)+COEFFS['sersic+sersic:sersic2'][0]))
                PARS["Second_Ser_reb"][k] = float(eval(str(reb)+COEFFS['sersic+sersic:sersic2'][1]))
                PARS["Second_Ser_n"][k] = float(COEFFS['sersic+sersic:sersic2'][2])
                PARS["Second_Ser_q"][k] = float(eval(str(q)+COEFFS['sersic+sersic:sersic2'][3]))
                PARS["Second_Ser_PA"][k] = float(eval(str(Mb)+COEFFS['sersic+sersic:sersic2'][4]))
                PARS["Second_Ser_meb"][k] = meb_bulge_f(PARS["Second_Ser_reb"][k]*scale,PARS["Second_Ser_Mb"][k],PARS["Second_Ser_n"][k],PARS["Second_Ser_q"][k])  
                #print PARS["First_Ser_Mb"], PARS["First_Ser_meb"], PARS["First_Ser_q"]
                #print PARS["Second_Ser_Mb"], PARS["Second_Ser_meb"], PARS["Second_Ser_q"]
                #exit()
                            

            #if (PARS["First_Ser_q"] > deca_setup.q_eon and model!='sersic+sersic') or model=='sersic+exp_disc':
            if model[k]=='sersic+exp_disc':
                PARS["First_Disk_Md"][k] = float(eval(str(Mb)+COEFFS['sersic+exp_disc:exp_disc'][0]))  
                PARS["First_Disk_hd"][k] = float(eval(str(reb)+COEFFS['sersic+exp_disc:exp_disc'][1]))
                PARS["First_Disk_q"][k] = float(eval(str(q)+COEFFS['sersic+exp_disc:exp_disc'][2]))
                PARS["First_Disk_PA"][k] = float(eval(str(PA)+COEFFS['sersic+exp_disc:exp_disc'][3]))
                PARS["First_Disk_m0d"][k] = m0d_disk_f(PARS["First_Disk_Md"][k],PARS["First_Disk_hd"][k]*scale,PARS["First_Disk_q"][k])
                
                PARS["First_Ser_Mb"][k] = float(eval(str(Mb)+COEFFS['sersic+exp_disc:sersic'][0]))
                PARS["First_Ser_reb"][k] = float(eval(str(reb)+COEFFS['sersic+exp_disc:sersic'][1]))  
                PARS["First_Ser_n"][k] = float(COEFFS['sersic+exp_disc:sersic'][2])
                PARS["First_Ser_q"][k] = float(COEFFS['sersic+exp_disc:sersic'][3])
                PARS["First_Ser_PA"][k] = float(eval(str(Mb)+COEFFS['sersic+exp_disc:sersic'][4]))
                PARS["First_Ser_meb"][k] = meb_bulge_f(PARS["First_Ser_reb"][k]*scale,PARS["First_Ser_Mb"][k],PARS["First_Ser_n"][k],PARS["First_Ser_q"][k])


                
            #if  (PARS["First_Ser_q"] <= deca_setup.q_eon and model!='sersic+sersic') or model=='sersic+eon_disc':
            if model[k]=='sersic+eon_disc':
                PARS["First_Disk_Md"][k] = float(eval(str(Mb)+COEFFS['sersic+eon_disc:eon_disc'][0])) 
                PARS["First_Disk_hd"][k] = float(eval(str(reb)+COEFFS['sersic+eon_disc:eon_disc'][1]))
                PARS["First_Disk_q"][k] = float(nan)
                PARS["First_Disk_PA"][k] = float(eval(str(PA)+COEFFS['sersic+eon_disc:eon_disc'][2]))
                PARS["First_Disk_z0"][k] = PARS["First_Disk_hd"][k] / 4.      
                PARS["First_Disk_m0d"][k] = Mb + 0.65 + 2.5*log10(2.*math.pi*PARS["First_Disk_hd"][k]*scale*PARS["First_Disk_z0"][k]*scale) + 5.*log10(scale)
                
                PARS["First_Ser_Mb"][k] = float(eval(str(Mb)+COEFFS['sersic+eon_disc:sersic'][0])) 
                PARS["First_Ser_reb"][k] = float(eval(str(reb)+COEFFS['sersic+eon_disc:sersic'][1]))   
                PARS["First_Ser_n"][k] = float(COEFFS['sersic+eon_disc:sersic'][2])
                PARS["First_Ser_q"][k] = float(COEFFS['sersic+eon_disc:sersic'][3])
                PARS["First_Ser_PA"][k] = float(eval(str(PA)+COEFFS['sersic+eon_disc:sersic'][4]))         

                

                

            if  model[k]=='agn+exp_disc':
                PARS["First_Disk_Md"][k] = float(eval(str(Mb)+COEFFS['agn+exp_disc:exp_disc'][0])) 
                PARS["First_Disk_hd"][k] = float(eval(str(reb)+COEFFS['agn+exp_disc:exp_disc'][1])) 
                PARS["First_Disk_q"][k] = float(eval(str(q)+COEFFS['agn+exp_disc:exp_disc'][2])) 
                PARS["First_Disk_PA"][k] = float(eval(str(PA)+COEFFS['agn+exp_disc:exp_disc'][3])) 
                PARS["First_Disk_m0d"][k] = m0d_disk_f(PARS["First_Disk_Md"][k],PARS["First_Disk_hd"][k]*scale,PARS["First_Disk_q"][k])
                
                PARS["PSF_M"][k] = 100. - 2.5*log10( (10**(0.4*(100.-Mb))/100.) )
                
                PARS["First_Ser_meb"][k] = float(nan)
                PARS["First_Ser_Mb"][k] = float(nan)
                PARS["First_Ser_reb"][k] = float(nan)
                PARS["First_Ser_n"][k] = float(nan)
                PARS["First_Ser_q"][k] = float(nan) 
                PARS["First_Ser_PA"][k] = float(nan)  



            if  model[k]=='agn+eon_disc':
                PARS["First_Disk_Md"][k] = float(eval(str(Mb)+COEFFS['agn+eon_disc:eon_disc'][0])) 
                PARS["First_Disk_hd"][k] = float(eval(str(reb)+COEFFS['agn+eon_disc:eon_disc'][1])) 
                PARS["First_Disk_q"][k] = float(nan)
                PARS["First_Disk_PA"][k] = float(eval(str(PA)+COEFFS['agn+eon_disc:eon_disc'][2])) 
                PARS["First_Disk_z0"][k] = PARS["First_Disk_hd"][k] / 4.      
                PARS["First_Disk_m0d"][k] = m0_disc_edge_on_f(PARS["First_Disk_Md"][k],PARS["First_Disk_hd"][k]*scale,PARS["First_Disk_z0"][k]*scale)
                
                PARS["PSF_M"][k] = 100. - 2.5*log10( (10**(0.4*(100.-Mb))/100.) )
                
                PARS["First_Ser_meb"][k] = float(nan)
                PARS["First_Ser_Mb"][k] = float(nan)
                PARS["First_Ser_reb"][k] = float(nan)
                PARS["First_Ser_n"][k] = float(nan)
                PARS["First_Ser_q"][k] = float(nan)
                PARS["First_Ser_PA"][k] = float(nan)  


            if  model[k]=='agn+ring':
                PARS["PSF_M"][k] = Mb

                PARS["Ring_m0"][k] = float(eval(str(meb)+COEFFS['agn+ring:ring'][0])) 
                PARS["Ring_sigma_r"][k] = float(eval(str(reb)+COEFFS['agn+ring:ring'][1])) 
                PARS["Ring_R_r"][k] = float(eval(str(reb)+COEFFS['agn+ring:ring'][2]))
                PARS["Ring_q"][k] = float(eval(str(q)+COEFFS['agn+ring:ring'][3])) 
                PARS["Ring_PA"][k] = float(eval(str(PA)+COEFFS['agn+ring:ring'][4]))       
                
                PARS["First_Ser_meb"][k] = float(nan)
                PARS["First_Ser_Mb"][k] = float(nan)
                PARS["First_Ser_reb"][k] = float(nan)
                PARS["First_Ser_n"][k] = float(nan)
                PARS["First_Ser_q"][k] = float(nan) 
                PARS["First_Ser_PA"][k] = float(nan)


            if  model[k]=='sersic+ring':
                PARS["First_Ser_Mb"][k] = float(eval(str(Mb)+COEFFS['sersic+ring:sersic'][0])) 
                PARS["First_Ser_reb"][k] = float(eval(str(reb)+COEFFS['sersic+ring:sersic'][1]))   
                PARS["First_Ser_n"][k] = float(COEFFS['sersic+ring:sersic'][2])
                PARS["First_Ser_q"][k] = 1.#float(COEFFS['sersic+ring:sersic'][3])
                PARS["First_Ser_PA"][k] = float(eval(str(PA)+COEFFS['sersic+ring:sersic'][4]))     
                
                PARS["Ring_m0"][k] = float(eval(str(meb)+COEFFS['sersic+ring:ring'][0])) 
                PARS["Ring_sigma_r"][k] = float(eval(str(reb)+COEFFS['sersic+ring:ring'][1])) 
                PARS["Ring_R_r"][k] = float(eval(str(reb)+COEFFS['sersic+ring:ring'][2]))
                PARS["Ring_q"][k] = 1.#float(eval(str(q)+COEFFS['sersic+ring:ring'][3])) 
                PARS["Ring_PA"][k] = float(eval(str(PA)+COEFFS['sersic+ring:ring'][4]))       


            if  model[k]=='ferrer+sersic':
                PARS["FERRER_m0"][k] = float(eval(str(Mb)+COEFFS['ferrer+sersic:ferrer'][0]))
                PARS["FERRER_Rad"][k] = float(eval(str(reb)+COEFFS['ferrer+sersic:ferrer'][1])) 
                PARS["FERRER_Alfa"][k] = float(COEFFS['ferrer+sersic:ferrer'][2])
                PARS["FERRER_Beta"][k] = float(COEFFS['ferrer+sersic:ferrer'][3])  
                PARS["FERRER_q"][k] = float(eval(str(q)+COEFFS['ferrer+sersic:ferrer'][4]))    
                PARS["FERRER_PA"][k] = float(eval(str(PA)+COEFFS['ferrer+sersic:ferrer'][5]))  

                PARS["First_Ser_Mb"][k] = float(eval(str(Mb)+COEFFS['ferrer+sersic:sersic'][0])) 
                PARS["First_Ser_reb"][k] = float(eval(str(reb)+COEFFS['ferrer+sersic:sersic'][1]))   
                PARS["First_Ser_n"][k] = float(COEFFS['ferrer+sersic:sersic'][2])
                PARS["First_Ser_q"][k] = float(COEFFS['ferrer+sersic:sersic'][3])
                PARS["First_Ser_PA"][k] = float(eval(str(PA)+COEFFS['ferrer+sersic:sersic'][4]))     
                
            
                

            if  model[k]=='sersic+sersic+sersic':   
                PARS["First_Ser_Mb"][k] = float(eval(str(Mb)+COEFFS['sersic+sersic+sersic:sersic1'][0])) 
                PARS["First_Ser_reb"][k] = float(eval(str(reb)+COEFFS['sersic+sersic+sersic:sersic1'][1]))   
                PARS["First_Ser_n"][k] = float(COEFFS['sersic+sersic+sersic:sersic1'][2])
                PARS["First_Ser_q"][k] = float(COEFFS['sersic+sersic+sersic:sersic1'][3])  
                PARS["First_Ser_PA"][k] = float(eval(str(PA)+COEFFS['sersic+sersic+sersic:sersic1'][4]))   
                PARS["First_Ser_meb"][k] = meb_bulge_f(PARS["First_Ser_reb"][k]*scale,PARS["First_Ser_Mb"][k],PARS["First_Ser_n"][k],PARS["First_Ser_q"][k])      

                PARS["Second_Ser_Mb"][k] = float(eval(str(Mb)+COEFFS['sersic+sersic+sersic:sersic2'][0]))
                PARS["Second_Ser_reb"][k] = float(eval(str(reb)+COEFFS['sersic+sersic+sersic:sersic2'][1]))
                PARS["Second_Ser_n"][k] = float(COEFFS['sersic+sersic+sersic:sersic2'][2])
                PARS["Second_Ser_q"][k] = float(eval(str(q)+COEFFS['sersic+sersic+sersic:sersic2'][3]))
                PARS["Second_Ser_PA"][k] = float(eval(str(Mb)+COEFFS['sersic+sersic+sersic:sersic2'][4]))
                PARS["Second_Ser_meb"][k] = meb_bulge_f(PARS["Second_Ser_reb"][k]*scale,PARS["Second_Ser_Mb"][k],PARS["Second_Ser_n"][k],PARS["Second_Ser_q"][k])  

                PARS["Third_Ser_Mb"][k] = float(eval(str(Mb)+COEFFS['sersic+sersic+sersic:sersic3'][0]))
                PARS["Third_Ser_reb"][k] = float(eval(str(reb)+COEFFS['sersic+sersic+sersic:sersic3'][1]))
                PARS["Third_Ser_n"][k] = float(COEFFS['sersic+sersic+sersic:sersic3'][2])
                PARS["Third_Ser_q"][k] = float(eval(str(q)+COEFFS['sersic+sersic+sersic:sersic3'][3]))
                PARS["Third_Ser_PA"][k] = float(eval(str(Mb)+COEFFS['sersic+sersic+sersic:sersic3'][4]))
                PARS["Third_Ser_meb"][k] = meb_bulge_f(PARS["Second_Ser_reb"][k]*scale,PARS["Second_Ser_Mb"][k],PARS["Second_Ser_n"][k],PARS["Second_Ser_q"][k])  


            if  model[k]=='agn+sersic+exp_disc':
                PARS["PSF_M"][k] = 100. - 2.5*log10( (10**(0.4*(100.-Mb))/100.) )
                
                PARS["First_Ser_Mb"][k] = float(eval(str(Mb)+COEFFS['agn+sersic+exp_disc:sersic'][0])) 
                PARS["First_Ser_reb"][k] = float(eval(str(reb)+COEFFS['agn+sersic+exp_disc:sersic'][1]))   
                PARS["First_Ser_n"][k] = float(COEFFS['agn+sersic+exp_disc:sersic'][2])
                PARS["First_Ser_q"][k] = float(eval(str(q)+COEFFS['agn+sersic+exp_disc:sersic'][3]))   
                PARS["First_Ser_PA"][k] = float(eval(str(PA)+COEFFS['agn+sersic+exp_disc:sersic'][4]))   
                PARS["First_Ser_meb"][k] = meb_bulge_f(PARS["First_Ser_reb"][k]*scale,PARS["First_Ser_Mb"][k],PARS["First_Ser_n"][k],PARS["First_Ser_q"][k])    
                
                PARS["First_Disk_Md"][k] = float(eval(str(Mb)+COEFFS['agn+sersic+exp_disc:exp_disc'][0])) 
                PARS["First_Disk_hd"][k] = float(eval(str(reb)+COEFFS['agn+sersic+exp_disc:exp_disc'][1])) 
                PARS["First_Disk_q"][k] = float(eval(str(q)+COEFFS['agn+sersic+exp_disc:exp_disc'][2])) 
                PARS["First_Disk_PA"][k] = float(eval(str(PA)+COEFFS['agn+sersic+exp_disc:exp_disc'][3])) 
                PARS["First_Disk_m0d"][k] = m0d_disk_f(PARS["First_Disk_Md"][k],PARS["First_Disk_hd"][k]*scale,PARS["First_Disk_q"][k])
                

            if  model[k]=='agn+sersic+eon_disc':
                PARS["PSF_M"][k] = 100. - 2.5*log10( (10**(0.4*(100.-Mb))/100.) )
                
                PARS["First_Ser_Mb"][k] = float(eval(str(Mb)+COEFFS['agn+sersic+eon_disc:sersic'][0])) 
                PARS["First_Ser_reb"][k] = float(eval(str(reb)+COEFFS['agn+sersic+eon_disc:sersic'][1]))   
                PARS["First_Ser_n"][k] = float(COEFFS['agn+sersic+eon_disc:sersic'][2])
                PARS["First_Ser_q"][k] = float(eval(str(q)+COEFFS['agn+sersic+eon_disc:sersic'][3]))   
                PARS["First_Ser_PA"][k] = float(eval(str(PA)+COEFFS['agn+sersic+eon_disc:sersic'][4]))   
                PARS["First_Ser_meb"][k] = meb_bulge_f(PARS["First_Ser_reb"][k]*scale,PARS["First_Ser_Mb"][k],PARS["First_Ser_n"][k],PARS["First_Ser_q"][k])    
                
                PARS["First_Disk_Md"][k] = float(eval(str(Mb)+COEFFS['agn+sersic+eon_disc:eon_disc'][0])) 
                PARS["First_Disk_hd"][k] = float(eval(str(reb)+COEFFS['agn+sersic+eon_disc:eon_disc'][1])) 
                PARS["First_Disk_q"][k] = float(nan)
                PARS["First_Disk_PA"][k] = float(eval(str(PA)+COEFFS['agn+sersic+eon_disc:eon_disc'][2])) 
                PARS["First_Disk_z0"][k] = PARS["First_Disk_hd"][k] / 4.      
                PARS["First_Disk_m0d"][k] = m0_disc_edge_on_f(PARS["First_Disk_Md"][k],PARS["First_Disk_hd"][k]*scale,PARS["First_Disk_z0"][k]*scale)

            if  model[k]=='sersic+ferrer+sersic':
                PARS["First_Ser_Mb"][k] = float(eval(str(Mb)+COEFFS['sersic+ferrer+sersic:sersic1'][0])) 
                PARS["First_Ser_reb"][k] = float(eval(str(reb)+COEFFS['sersic+ferrer+sersic:sersic1'][1]))   
                PARS["First_Ser_n"][k] = float(COEFFS['sersic+ferrer+sersic:sersic1'][2])
                PARS["First_Ser_q"][k] = float(COEFFS['sersic+ferrer+sersic:sersic1'][3])
                PARS["First_Ser_PA"][k] = float(eval(str(PA)+COEFFS['sersic+ferrer+sersic:sersic1'][4]))        
                
                PARS["FERRER_m0"][k] = float(eval(str(Mb)+COEFFS['sersic+ferrer+sersic:ferrer'][0]))
                PARS["FERRER_Rad"][k] = float(eval(str(reb)+COEFFS['sersic+ferrer+sersic:ferrer'][1])) 
                PARS["FERRER_Alfa"][k] = float(COEFFS['sersic+ferrer+sersic:ferrer'][2])
                PARS["FERRER_Beta"][k] = float(COEFFS['sersic+ferrer+sersic:ferrer'][3])  
                PARS["FERRER_q"][k] = float(eval(str(q)+COEFFS['sersic+ferrer+sersic:ferrer'][4]))    
                PARS["FERRER_PA"][k] = float(eval(str(PA)+COEFFS['sersic+ferrer+sersic:ferrer'][5]))  

                PARS["Second_Ser_Mb"][k] = float(eval(str(Mb)+COEFFS['sersic+ferrer+sersic:sersic2'][0])) 
                PARS["Second_Ser_reb"][k] = float(eval(str(reb)+COEFFS['sersic+ferrer+sersic:sersic2'][1]))   
                PARS["Second_Ser_n"][k] = float(COEFFS['sersic+ferrer+sersic:sersic2'][2])
                PARS["Second_Ser_q"][k] = float(COEFFS['sersic+ferrer+sersic:sersic2'][3])
                PARS["Second_Ser_PA"][k] = float(eval(str(PA)+COEFFS['sersic+ferrer+sersic:sersic2'][4]))  

  return PARS
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO FIND INITIAL GUESS BASED ON SOME INTERNAL CALCULATIONS
def initial_guess_deca(new_images,observation_info,object_info,add_info,PARS,model):
  #TODO: ADD MULTIOBJECT INITIAL GUESS!
  ## Load the input arrays 
  [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,Sky_level,SkySubtr,sampling,convolution_box,FILTER] = observation_info
  [xc,yc,name] = object_info
  [number,fitting_proc,mode,model,region_to_fit] = add_info
  [input_image,sigma_image,psf_image,mask_image] = new_images
  
  if deca_setup.code=='GALFIT':
      comp_model = './galfit_1/composed_model.fits'
  elif deca_setup.code=='IMFIT':
      comp_model = './imfit_1/composed_model.fits'

  hdu_comp_model = pyfits.open(comp_model)
  data_resid = hdu_comp_model[2].data
  

  for k in range(len(PARS["First_Ser_meb"])): 
        meb = PARS["First_Ser_meb"][k]
        Mb = PARS["First_Ser_Mb"][k]
        reb = PARS["First_Ser_reb"][k]
        n = PARS["First_Ser_n"][k]
        PA = PARS["First_Ser_PA"][k]
        q = PARS["First_Ser_q"][k]
        
        # Read in input_image
        #hdu_ima = pyfits.open(input_image)
        #data_ima = hdu_ima[0].data
        data_ima = data_resid + hdu_comp_model[4+k].data

        hdu = pyfits.PrimaryHDU(data_ima)
        hdu.writeto('one_object.fits',clobber=True)
        
        
        
            
            
        
        # If needed, cut out the galaxy image
        if region_to_fit!=None:                
                        #  x1:x2;y1:y2 format
                        xmin = int(region_to_fit.split(':')[0])  
                        xmax = int(region_to_fit.split(':')[1].split(';')[0])  
                        ymin = int(region_to_fit.split(';')[1].split(':')[0])  
                        ymax = int(region_to_fit.split(':')[2])
        else:
                        ny,nx = np.shape(data_ima)
                        xmin = 0
                        xmax = nx
                        ymin = 0 
                        ymax = ny
        
        data_ima_crop = data_ima[ymin:ymax,xmin:xmax]
        ny,nx = np.shape(data_ima_crop)

        # Mask
        if os.path.exists(mask_image):
            hdu_mask = pyfits.open(mask_image)
            data_mask = hdu_mask[0].data
            data_mask_crop = data_mask[ymin:ymax,xmin:xmax]

            # Create astropy mask
            mask_astropy = np.zeros_like(np.array(data_mask_crop,dtype=float),dtype=bool)
            
            for k in range(ny):
                for i in range(nx):
                    if data_mask_crop[k,i]!=0.:
                        mask_astropy[k,i] = True
        else:
            mask_astropy = None


        if 'eon_disc' in model[k]:
            hdu = pyfits.PrimaryHDU(data_ima_crop)
            hdu.writeto('one_object.fits',clobber=True)

            hdu = pyfits.PrimaryHDU(data_mask_crop)
            hdu.writeto('mask_one_object.fits',clobber=True)
            
            # We have to rotate the frame so that the major axis would be horizontal!
            #TODO: import rotima
            xc,yc = rotima.crotima('one_object.fits', 'one_object_rot.fits', xc-xmin, yc-ymin, 90+angle, layer=1, set_wcs=True)
            data_ima_crop = pyfits.open('one_object_rot.fits')
            rotima.crotima('mask_one_object.fits', 'mask_one_object_rot.fits', xc-xmin, yc-ymin, 90+angle, layer=1, set_wcs=True)
            data_mask_crop = pyfits.open('mask_one_object_rot.fits')            

        # Mask
        if os.path.exists(mask_image):
            # Create astropy mask
            mask_astropy = np.zeros_like(np.array(data_mask_crop,dtype=float),dtype=bool)
            
            for k in range(ny):
                for i in range(nx):
                    if data_mask_crop[k,i]!=0.:
                        mask_astropy[k,i] = True
        else:
            mask_astropy = None                        

        ## CONSIDER TWO CASES: NON-EDGE-ON GALAXY OR EDGE-ON GALAXY
        if 'eon_disc' not in model[k]:#PARS["First_Ser_q"] > deca_setup.q_eon:
                    # -------------------------NON-EDGE-ON GALAXY-------------------------
                    # Plot 1D-profile          
                    bin_centers,radial_prof = radial_profile.azimuthalAverage(data_ima_crop, returnradii=True, mask=np.invert(mask_astropy))

                    # Find std
                    I_mean, I_median, std = sigma_clipped_stats(radial_prof, sigma=3.0, iters=5)

                    # Define rmax
                    rmax = find_lowest_val(bin_centers,radial_prof,std)

                    # Crop profile at rmax:
                    r = []; I = []
                    for k in range(len(radial_prof)):
                        if bin_centers[k]<=rmax:
                            r.append(bin_centers[k])
                            I.append(radial_prof[k])
                    
                    # Find the bulge- and disc-dominated regions
                    R_bulge_dom = define_bulge_dom_region(r,I,m0,scale)

                    
                    r_disc = []; I_disc = [];r_bulge = []; I_bulge = []
                    for k in range(len(r)):
                        if r[k]>R_bulge_dom*1.2:
                            r_disc.append(r[k])
                            I_disc.append(I[k])
                        if r[k]<R_bulge_dom:
                            r_bulge.append(r[k])
                            I_bulge.append(I[k])
                    r_disc = np.array(r_disc)
                    I_disc = np.array(I_disc)
                    r_bulge = np.array(r_bulge)
                    I_bulge = np.array(I_bulge)    


                    
                    # Fit the disc
                    h_disc,m0_disc = fit_exp_disc(r_disc,I_disc,m0,scale)

        
                    mag_disc = magDisk_f(m0_disc,h_disc,q)
                    
                    # Fit the bulge
                    re_bulge,me_bulge = fit_ser_bulge(r_bulge,I_bulge,m0,scale)  
                    n_bulge = 2.5
                    q_bulge = 0.8
                    mag_bulge = magBulge_f(re_bulge,me_bulge,n_bulge,q_bulge)

                    if model[k]=='sersic+exp_disc' or model[k]=='none' or model[k]=='None' or model[k]==None:
                        PARS["First_Disk_m0d"][k] = m0_disc
                        PARS["First_Disk_Md"][k] = mag_disc
                        PARS["First_Disk_hd"][k] = h_disc
                        PARS["First_Disk_PA"][k] = PA
                        PARS["First_Disk_q"][k] = q
                        
                        PARS["First_Ser_meb"][k] = me_bulge
                        PARS["First_Ser_Mb"][k] = mag_bulge
                        PARS["First_Ser_reb"][k] = re_bulge
                        PARS["First_Ser_n"][k] = n
                        PARS["First_Ser_q"][k] = q_bulge
                        PARS["First_Ser_PA"][k] = PA
                        #create_constrain_file(re_bulge, model)


                    if model[k]=='exp_disc':
                        PARS["First_Disk_m0d"][k] = m0_disc
                        PARS["First_Disk_Md"][k] = mag_disc
                        PARS["First_Disk_hd"][k] = h_disc
                        PARS["First_Disk_PA"][k] = PA
                        PARS["First_Disk_q"][k] = q
                        
                        PARS["First_Ser_meb"][k] = float(nan)
                        PARS["First_Ser_Mb"][k] = float(nan)
                        PARS["First_Ser_reb"][k] = float(nan)
                        PARS["First_Ser_n"][k] = float(nan)
                        PARS["First_Ser_PA"][k] = float(nan)  
                        PARS["First_Ser_q"][k] = float(nan)   
                        #create_constrain_file(re_bulge, model)

                    if  model[k]=='agn+exp_disc' or model[k]=='exp_disc+agn':
                        PARS["First_Disk_m0d"][k] = m0_disc
                        PARS["First_Disk_Md"][k] = mag_disc
                        PARS["First_Disk_hd"][k] = h_disc
                        PARS["First_Disk_PA"][k] = PA
                        PARS["First_Disk_q"][k] = q
                        
                        PARS["PSF_M"][k] = 100. - 2.5*log10( (10**(0.4*(100.-mag_disc))/100.) )
                    
                        PARS["First_Ser_meb"][k] = float(nan)
                        PARS["First_Ser_Mb"][k] = float(nan)
                        PARS["First_Ser_reb"][k] = float(nan)
                        PARS["First_Ser_n"][k] = float(nan)
                        PARS["First_Ser_PA"][k] = float(nan)  
                        PARS["First_Ser_q"][k] = float(nan) 

                    if  model[k]=='sersic+sersic':                  
                        PARS["First_Ser_meb"][k] = me_bulge
                        PARS["First_Ser_Mb"][k] = mag_bulge
                        PARS["First_Ser_reb"][k] = re_bulge
                        PARS["First_Ser_n"][k] = n
                        PARS["First_Ser_q"][k]= q_bulge
                        PARS["First_Ser_PA"][k] = PA

                        PARS["Second_Ser_Mb"][k] = mag_disc
                        PARS["Second_Ser_reb"][k] = 1.68*h_disc
                        PARS["Second_Ser_n"][k] = 1.
                        PARS["Second_Ser_q"][k] = q
                        PARS["Second_Ser_PA"][k] = PA
                        PARS["Second_Ser_meb"][k] = meb_bulge_f(PARS["Second_Ser_reb"][k]*scale,PARS["Second_Ser_Mb"][k],PARS["Second_Ser_n"][k],PARS["Second_Ser_q"][k])                


        else: #PARS["First_Ser_q"] <= deca_setup.q_eon:  
                    # -------------------------EDGE-ON GALAXY-------------------------
                    # Create the 1D-profile           
                    bin_centers,radial_prof = radial_profile.azimuthalAverage(data_ima_crop,returnradii=True,mask=np.invert(mask_astropy))
                    
                    # Find std
                    I_mean, I_median, std = sigma_clipped_stats(radial_prof, sigma=3.0, iters=5)

                    # Define rmax
                    Radius = find_lowest_val(bin_centers,radial_prof,std)
                    Minor_axis = Radius * q

                    I_mean, I_median, std = sigma_clipped_stats(data_ima_crop, mask=mask_astropy, sigma=3.0, iters=5)
                    Imin = 3.*std
                    mag_level = m0 - 2.5*log10(Imin/(scale**2))



                    import edgeon_initial_guess
                    PARS_bps,file_1,file_2,file_3 = edgeon_initial_guess.main('one_object_rot.fits', 'mask_one_object_rot.fits', mag_level, Radius, Minor_axis, 0., m0, scale, xc=xc, zc=yc, inter_mode=False, min_disc_radius=0., Profile='horizon')
                    
                    [z0_disc] = PARS_bps["THICKNESS"]
                    [h_disc,m0_disc,R_disc_min] = PARS_bps["DISC_DOM"]
                    
                    if model[k]=='sersic+eon_disc' or model[k]=='none' or model[k]=='None' or model[k]==None or model[k]=='eon_disc+sersic':
                        PARS["First_Disk_m0d"][k] = m0_disc
                        PARS["First_Disk_hd"][k] = h_disc
                        PARS["First_Disk_PA"][k] = PA
                        PARS["First_Disk_q"][k] = float(nan)
                        PARS["First_Disk_z0"][k] = z0_disc
                        PARS["First_Disk_Md"][k] = magDisk_edge_on_f(PARS["First_Disk_m0d"][k],PARS["First_Disk_hd"][k]*scale,PARS["First_Disk_z0"][k]*scale)
                        
                        PARS["First_Ser_Mb"][k] = PARS["First_Disk_Md"][k] + 1.
                        PARS["First_Ser_reb"][k] = PARS["First_Disk_z0"][k]
                        PARS["First_Ser_n"][k] = 2.5
                        PARS["First_Ser_q"][k] = 0.7
                        PARS["First_Ser_PA"][k] = PA
                        PARS["First_Ser_meb"][k] = meb_bulge_f(PARS["First_Ser_reb"][k]*scale,PARS["First_Ser_Mb"][k],PARS["First_Ser_n"][k],PARS["First_Ser_q"][k])                
                        
                        #create_constrain_file(re_bulge, model)
                        
                    if model[k]=='eon_disc':
                        PARS["First_Disk_m0d"][k] = m0_disc
                        PARS["First_Disk_hd"][k] = h_disc
                        PARS["First_Disk_PA"][k] = PA
                        PARS["First_Disk_q"][k] = float(nan)
                        PARS["First_Disk_z0"][k] = z0_disc
                        PARS["First_Disk_Md"][k] = magDisk_edge_on_f(PARS["First_Disk_m0d"][k],PARS["First_Disk_hd"][k]*scale,PARS["First_Disk_z0"][k]*scale)               

                        PARS["First_Ser_meb"][k] = float(nan)
                        PARS["First_Ser_Mb"][k] = float(nan)
                        PARS["First_Ser_reb"][k] = float(nan)
                        PARS["First_Ser_n"][k] = float(nan)
                        PARS["First_Ser_PA"][k] = float(nan)  
                        PARS["First_Ser_q"][k] = float(nan) 

                    if model[k]=='agn+eon_disc' or model[k]=='eon_disc+agn':
                        PARS["First_Disk_m0d"][k] = m0_disc
                        PARS["First_Disk_hd"][k] = h_disc
                        PARS["First_Disk_PA"][k] = PA
                        PARS["First_Disk_q"][k] = float(nan)
                        PARS["First_Disk_z0"][k] = z0_disc
                        PARS["First_Disk_Md"][k] = magDisk_edge_on_f(PARS["First_Disk_m0d"][k],PARS["First_Disk_hd"][k]*scale,PARS["First_Disk_z0"][k]*scale)
                        
                        PARS["PSF_M"][k] = 100. - 2.5*log10( (10**(0.4*(100.-PARS["First_Disk_Md"][k]))/100.) )

                        PARS["First_Ser_meb"][k] = float(nan)
                        PARS["First_Ser_Mb"][k] = float(nan)
                        PARS["First_Ser_reb"][k] = float(nan)
                        PARS["First_Ser_n"][k]= float(nan)
                        PARS["First_Ser_PA"][k] = float(nan)  
                        PARS["First_Ser_q"][k] = float(nan) 
                        
                        #create_constrain_file(re_bulge, model)
                        
                    for file in [file_1,file_2,file_3]:
                        if os.path.exists(file):
                            os.remove(file)
                    
                    os.remove('one_object.fits')
                
  return PARS            
# -----------------------------------------------------------------
