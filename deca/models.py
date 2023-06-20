#!/usr/bin/python
# Module with some functions to deal with the model
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **        © Astronomical Observatory, Ghent University         **
# *****************************************************************

# Import standard modules
import random as random_number
import sys
import math
import numpy as np
from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import NullFormatter
from numpy import *
from pylab import *
import os
import shutil
import signal
import collections
import re

# Import astro modules
import pyfits

# Import DECA modules
import initial_guess

# -----------------------------------------------------------------
# FUNCTION TO DEFINE THE MODEL
def define_model_galaxies(number_of_objects):
  # The model comprises the following components:
  # First Disc, Second Disc, First Sersic, Second Sersic, Third Sersic, Ferrer, AGN (as a PSF), Ring and Sky
  # SB and magnitudes are in mag/arcsec^2 or magnitudes
  # Sizes are in pixels
  PARS = collections.OrderedDict()
  
  for k in range(number_of_objects):
        # FIRST DISK (==Thin Disc OR ==DISC in case of Bule+Disc model)
        PARS["First_Disk_m0d"].append(float(nan))
        PARS["First_Disk_Md"].append(float(nan))
        PARS["First_Disk_hd"].append(float(nan))
        PARS["First_Disk_z0"].append(float(nan))
        PARS["First_Disk_PA"].append(float(nan))  
        PARS["First_Disk_q"].append(float(nan))  
        
        # SECOND DISK (==Thick disc)
        PARS["Second_Disk_m0d"].append(float(nan))
        PARS["Second_Disk_Md"].append(float(nan))
        PARS["Second_Disk_hd"].append(float(nan))
        PARS["Second_Disk_z0"].append(float(nan))
        PARS["Second_Disk_PA"].append(float(nan)) 
        PARS["Second_Disk_q"].append(float(nan))  

        # FIRST SERSIC (==Bulge in case of Bule+Disc model)
        PARS["First_Ser_meb"].append(float(nan))
        PARS["First_Ser_Mb"].append(float(nan))
        PARS["First_Ser_reb"].append(float(nan))
        PARS["First_Ser_n"].append(float(nan))
        PARS["First_Ser_PA"].append(float(nan))  
        PARS["First_Ser_q"].append(float(nan)) 
        
        # SECOND SERSIC (==Bar or Disc)
        PARS["Second_Ser_meb"].append(float(nan))
        PARS["Second_Ser_Mb"].append(float(nan))
        PARS["Second_Ser_reb"].append(float(nan))
        PARS["Second_Ser_n"].append(float(nan))
        PARS["Second_Ser_PA"].append(float(nan))  
        PARS["Second_Ser_q"].append(float(nan)) 

        # THIRD SERSIC (==Bar or Disc)
        PARS["Third_Ser_meb"].append(float(nan))
        PARS["Third_Ser_Mb"].append(float(nan))
        PARS["Third_Ser_reb"].append(float(nan))
        PARS["Third_Ser_n"].append(float(nan))
        PARS["Third_Ser_PA"].append(float(nan))  
        PARS["Third_Ser_q"].append(float(nan)) 
        
        # FERRER (==Bar)
        PARS["FERRER_m0"].append(float(nan))
        PARS["FERRER_Rad"].append(float(nan))
        PARS["FERRER_Alfa"].append(float(nan))
        PARS["FERRER_Beta"].append(float(nan))  
        PARS["FERRER_q"].append(float(nan))    
        PARS["FERRER_PA"].append(float(nan))   

        # PSF (==AGN)
        PARS["PSF_M"].append(float(nan))
        
        # GAUSSIAN RING
        PARS["Ring_m0"].append(float(nan))
        PARS["Ring_sigma_r"].append(float(nan))
        PARS["Ring_R_r"].append(float(nan))
        PARS["Ring_q"].append(float(nan))
        PARS["Ring_PA"].append(float(nan))

  
  # Sky
  PARS["SKY"] = float(nan)
  return PARS
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO DEFINE MODEL FOR POINT SOURCES (STARS)
def define_model_stars(number_of_objects):
  # The model of stars as PSF
  PARS = collections.OrderedDict()
  
  for k in range(number_of_objects):
        # PSF (==AGN)
        PARS["PSF_M"].append(float(nan))

  return PARS
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# READ GALFIT RESULTS OF THE INITIAL SERSIC DECOMPOSITION AND PUT THEM INTO THE MODEL
def read_results_First_Ser(PARS,galfit_outp_file,pix2sec):
    with open(galfit_outp_file) as f:
	lines = f.readlines()
        comp_number = 0
	for k in range(len(lines)):    
	  if '0) sersic' in lines[k] and 'Component type' in lines[k]:
            PARS["First_Ser_reb"][comp_number] = float(lines[k+3].split()[1])
            PARS["First_Ser_n"][comp_number] = float(lines[k+4].split()[1]) 
            if '10)' in lines[k+9] and 'Position angle (PA) [deg: Up=0, Left=90]' in lines[k+9]:
                PARS["First_Ser_PA"][comp_number] = float(lines[k+9].split()[1]) 
                PARS["First_Ser_q"][comp_number] = float(lines[k+8].split()[1])  
            else:
                PARS["First_Ser_PA"][comp_number] = float(lines[k+6].split()[1]) 
                PARS["First_Ser_q"][comp_number] = float(lines[k+5].split()[1])                  
            
            PARS["First_Ser_Mb"][comp_number] = float(lines[k+2].split()[1])
            PARS["First_Ser_meb"][comp_number] = initial_guess.meb_bulge_f(PARS["First_Ser_reb"][comp_number]*pix2sec,PARS["First_Ser_Mb"][comp_number],PARS["First_Ser_n"][comp_number],PARS["First_Ser_q"][comp_number])             
            comp_number = comp_number + 1
	  if '0) sersic1' in lines[k] and 'Component type' in lines[k]:
            PARS["First_Ser_reb"][comp_number] = float(lines[k+3].split()[1])
            PARS["First_Ser_n"][comp_number] = float(lines[k+4].split()[1]) 
            if '10)' in lines[k+9] and 'Position angle (PA) [deg: Up=0, Left=90]' in lines[k+9]:
                PARS["First_Ser_PA"][comp_number] = float(lines[k+9].split()[1]) 
                PARS["First_Ser_q"][comp_number] = float(lines[k+8].split()[1])  
            else:
                PARS["First_Ser_PA"][comp_number] = float(lines[k+6].split()[1]) 
                PARS["First_Ser_q"][comp_number] = float(lines[k+5].split()[1])                  
            
            PARS["First_Ser_meb"][comp_number] = float(lines[k+2].split()[1])
            PARS["First_Ser_Mb"][comp_number] = initial_guess.magBulge_f(PARS["First_Ser_reb"][comp_number]*pix2sec,PARS["First_Ser_meb"][comp_number],PARS["First_Ser_n"][comp_number],PARS["First_Ser_q"][comp_number])
            comp_number = comp_number + 1
	return PARS  
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# READ IMFIT RESULTS OF THE INITIAL SERSIC DECOMPOSITION AND PUT THEM INTO THE MODEL
def read_results_First_Ser_IMFIT(PARS,imfit_outp_file,m0,pix2sec):
    with open(imfit_outp_file) as f:
	lines = f.readlines()
        comp_number = 0
        for k in range(len(lines)):      
	  if 'FUNCTION' in lines[k] and 'Sersic' in lines[k]:
            for i in range(5):
                if 'PA' in lines[k+i+1]:
                    PARS["First_Ser_PA"][comp_number] = float(lines[k+i+1].split('\t\t')[1].split('#')[0])
                if 'ell' in lines[k+i+1]:
                    PARS["First_Ser_q"][comp_number] = 1. - float(lines[k+i+1].split('\t\t')[1].split('#')[0])
                if lines[k+i+1].split('\t\t')[0]=='n':
                    PARS["First_Ser_n"][comp_number] = float(lines[k+i+1].split('\t\t')[1].split('#')[0])    
                if 'r_e' in lines[k+i+1]:
                    PARS["First_Ser_reb"][comp_number] = float(lines[k+i+1].split('\t\t')[1].split('#')[0])  
                if 'I_e' in lines[k+i+1]:
                    PARS["First_Ser_meb"][comp_number] = m0 - 2.5*log10(float(lines[k+i+1].split('\t\t')[1].split('#')[0])) + 5.*log10(pix2sec)
                PARS["First_Ser_Mb"] = initial_guess.magBulge_f(PARS["First_Ser_reb"][comp_number]*pix2sec,PARS["First_Ser_meb"][comp_number],PARS["First_Ser_n"][comp_number],PARS["First_Ser_q"][comp_number])
            comp_number = comp_number + 1
    return PARS
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO READ GALFIT RESULTS FOR STARS AND PUT THEM INTO THE MODEL
def read_results_stars(PARS_stars, galfit_outp_file):
    with open(galfit_outp_file) as f:
	lines = f.readlines()
        comp_number = 0
	for k in range(len(lines)):    
	  if '0) psf' in lines[k] and ' type' in lines[k]:
            PARS["PSF_M"][comp_number] = float(lines[k+2].split()[1])            
            comp_number = comp_number + 1

	return PARS  
# -----------------------------------------------------------------


# -----------------------------------------------------------------
def read_results_stars_IMFIT(PARS_stars, imfit_outp_file):
    # TODO: FUNCTION TO READ IMFIT RESULTS FOR STARS AND PUT THEM INTO THE MODEL
    return None
# -----------------------------------------------------------------
