#!/usr/bin/python
# Module to create input file for GALFIT decomposition
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **         Astronomical Observatory, Ghent University         **
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

# Import astro packages
import astropy.io.fits as pyfits

# Import DECA modules
import deca_setup

tmp_out = sys.stdout

# -----------------------------------------------------------------
# FUNCTION TO RETURN FUNCTIONS AND PARAMETERS USED IN GALFIT
def read_list_of_funcs(functions=None):
        if functions==None:
            functions = ['sersic','sersic2','nuker','devauc','expdisk','edgedisk','moffat','ferrer','ferrer2','gaussian','king','psf','sky']
            functions = sorted(functions)


        Galfit_functions = []
        for function in functions:
          if function=='sersic':
            Pars = [1,2,3,4,5,9,10]
            PARS_com = ['x_cen','y_cen','m_tot','r_e','n','q','PA']
          if function=='sersic2':
            Pars = [1,2,3,4,5,9,10]
            PARS_com = ['x_cen','y_cen','mu_re','r_e','n','q','PA']
          if function=='nuker':
            Pars = [1,2,3,4,5,6,7,9,10]
            PARS_com = ['x_cen','y_cen','mu_Rb','R_b','alpha','beta','gamma','q','PA']
          if function=='devauc':
            Pars = [1,2,3,4,9,10]
            PARS_com = ['x_cen','y_cen','m_tot','r_e','q','PA']
          if function=='expdisk':
            Pars = [1,2,3,4,9,10]
            PARS_com = ['x_cen','y_cen','m_tot','h','q','PA']
          if function=='edgedisk':
            Pars = [1,2,3,4,5,10]
            PARS_com = ['x_cen','y_cen','mu_0','z_0','h','PA']
          if function=='moffat':
            Pars = [1,2,3,4,5,9,10]
            PARS_com = ['x_cen','y_cen','m_tot','FWHM','beta','q','PA']
          if function=='ferrer':
            Pars = [1,2,3,4,5,6,9,10]
            PARS_com = ['x_cen','y_cen','mu_0','R_tr','alpha','beta','q','PA']
          if function=='ferrer2':
            Pars = [1,2,3,4,5,6,9,10]
            PARS_com = ['x_cen','y_cen','mu_re','R_tr','alpha','beta','q','PA']
          if function=='gaussian':
            Pars = [1,2,3,4,9,10]
            PARS_com = ['x_cen','y_cen','m_tot','FWHM','q','PA']
          if function=='king':
            Pars = [1,2,3,4,5,6,9,10]
            PARS_com = ['x_cen','y_cen','mu_0','R_c','R_t','alpha','q','PA']
          if function=='psf':
            Pars = [1,2,3]
            PARS_com = ['x_cen','y_cen','m_tot']
          if function=='sky':
            Pars = [1,2,3]
            PARS_com = ['sky_back','dsky/dx','dsky/dy']

          for k in range(len(Pars)):
              Par = Pars[k]
              Com = PARS_com[k]
              Galfit_functions.append(function+':'+str(Par)+':'+str(Com))
        #print(functions,Galfit_functions
        #exit()
        return functions,Galfit_functions
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO RETURN GEOMETRICAL PARAMETERS
def geom_pars():
  return ['r_e','R_b','h','z_0','FWHM','R_tr','R_c','R_t']
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO RETURN TOTAL LUMINOSITY PARAMETERS
def lum_pars():
  return ['m_tot']
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO RETURN SURFACE BRIGHTNESS PARAMETERS
def SB_pars():
  return ['mu_Rb','mu_0','mu_re','sky_back']
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO CREATE A HEADER FOR GALFIT FILE
def header(file_Image,file_Sigma,file_Psf,file_Mask,file_constraints,file_galfit_outimage,xmin,xmax,ymin,ymax,magZP,generalScaleValue,sampling,convolution_box,write_to=None):
        '''
        try:
            hdulist = pyfits.open(file_Image)
            image = hdulist[0].data    
            ySize, xSize = image.shape
            
            hdulist = pyfits.open(file_Psf)
            image = hdulist[0].data    
            ySize_psf, xSize_psf = image.shape
            
            conv_size = deca_setup.conv_size

            if conv_size==None:
              psf_box = min([ySize_psf, xSize_psf,xmax-xmin+5])
            else:
              psf_box = int(np.ceil(conv_size))
        except:
              psf_box = 100 # No PSF convolution
        '''
        if convolution_box==None:
            psf_box = deca_setup.conv_size
        else:
            psf_box = convolution_box
            
        if write_to==None:
            print("\n===============================================================================")
            print("# IMAGE and GALFIT CONTROL PARAMETERS")
            print("A) %s                  # Input data image (FITS file)" % (file_Image))
            print("B) %s                  # Output data image block" % (file_galfit_outimage))
            print("C) %s                # Sigma image name (made from data if blank or none)" % (file_Sigma))
            print("D) %s                # Input PSF image and (optional) diffusion kernel" % (file_Psf))
            print("E) %i                   # PSF fine sampling factor relative to data" % (sampling))
            print("F) %s                # Bad pixel mask (FITS image or ASCII coord list)" % (file_Mask))
            print("G) %s                # File with parameter constraints (ASCII file)" % (file_constraints))
            print("H) %i    %i   %i    %i   # Image region to fit (xmin xmax ymin ymax)" % (xmin, xmax, ymin, ymax))
            print("I) %i    %i              # Size of the convolution box (x y)" % (psf_box,psf_box))
            print("J) %.3f              # Magnitude photometric zeropoint" % (magZP))
            print("K) %.3f  %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (generalScaleValue,generalScaleValue))
            print("O) regular             # Display type (regular, curses, both)")
            print("P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n")

            print("# INITIAL FITTING PARAMETERS")
            print("#")
            print("#   For object type, the allowed functions are:" )
            print("#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat," )
            print("#       ferrer, powsersic, sky, and isophote." )
            print("#"  )
            print("#   Hidden parameters will only appear when they're specified:")
            print("#       C0 (diskyness/boxyness)," )
            print("#       Fn (n=integer, Azimuthal Fourier Modes),")
            print("#       R0-R10 (PA rotation, for creating spiral structures).")
            print("#" )
            print("# -----------------------------------------------------------------------------")
            print("#   par)    par value(s)    fit toggle(s)    # parameter description" )
            print("# -----------------------------------------------------------------------------\n")
        else:
            write_to.write( "\n===============================================================================")
            write_to.write( "# IMAGE and GALFIT CONTROL PARAMETERS")
            write_to.write( "A) %s                  # Input data image (FITS file)" % (file_Image))
            write_to.write( "B) %s                  # Output data image block" % (file_galfit_outimage))
            write_to.write( "C) %s                # Sigma image name (made from data if blank or none)" % (file_Sigma))
            write_to.write( "D) %s                # Input PSF image and (optional) diffusion kernel" % (file_Psf))
            write_to.write( "E) %i                   # PSF fine sampling factor relative to data" % (sampling))
            write_to.write( "F) %s                # Bad pixel mask (FITS image or ASCII coord list)" % (file_Mask))
            write_to.write( "G) %s                # File with parameter constraints (ASCII file)" % (file_constraints))
            write_to.write( "H) %i    %i   %i    %i   # Image region to fit (xmin xmax ymin ymax)" % (xmin, xmax, ymin, ymax))
            write_to.write( "I) %i    %i              # Size of the convolution box (x y)" % (psf_box,psf_box))
            write_to.write( "J) %.3f              # Magnitude photometric zeropoint" % (magZP))
            write_to.write( "K) %.3f  %.3f        # Plate scale (dx dy)    [arcsec per pixel]" % (generalScaleValue,generalScaleValue))
            write_to.write( "O) regular             # Display type (regular, curses, both)")
            write_to.write( "P) 0                   # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n")

            write_to.write( "# INITIAL FITTING PARAMETERS")
            write_to.write( "#")
            write_to.write( "#   For object type, the allowed functions are:" )
            write_to.write( "#       nuker, sersic, expdisk, devauc, king, psf, gaussian, moffat," )
            write_to.write( "#       ferrer, powsersic, sky, and isophote." )
            write_to.write( "#"  )
            write_to.write( "#   Hidden parameters will only appear when they're specified:")
            write_to.write( "#       C0 (diskyness/boxyness)," )
            write_to.write( "#       Fn (n=integer, Azimuthal Fourier Modes),")
            write_to.write( "#       R0-R10 (PA rotation, for creating spiral structures).")
            write_to.write( "#" )
            write_to.write( "# -----------------------------------------------------------------------------")
            write_to.write( "#   par)    par value(s)    fit toggle(s)    # parameter description" )
            write_to.write( "# -----------------------------------------------------------------------------\n"       ) 
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION FOR SERSIC 
def Sersic(component,xc,yc,meb,Mb,reb,n,q,PA):
         print("\n# Component %i: Bulge function" % (component))
         #if np.isnan(meb)==False:
         #   print(" 0) sersic2             # Object type"
         #   MM = meb
         #else:
         print(" 0) sersic              # Object type"   )
         MM = Mb
         print(" 1) %.2f  %.2f  1 1    #  position x, y        [pixel]" % (xc,yc))
         print(" 3) %.3f       1       #  total magnitude"  % (MM)   )
         print(" 4) %.3f       1       #    R_e              [Pixels]" % (reb))
         print(" 5) %.3f       1       #  Sersic exponent (deVauc=4, expdisk=1)" % (n))
         print(" 9) %.3f       1       # axis ratio (b/a)" % (q))
         print("10) %.3f       1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA))
         if deca_setup.c0==True:
           print("C0) 0.0       1       # Boxyness")
         print(" Z) 0                  #  Skip this model in output image?  (yes=1, no=0)")
# -----------------------------------------------------------------




# -----------------------------------------------------------------
# FUNCTION FOR EXPONENTIAL DISC
def ExpDisc(component,xc,yc,m0d,Md,h,q,PA):
         print("\n# Component %i: Exponential function" % (component))
         #if np.isnan(m0d)==False:
         #   print(" 0) expdisk1             # Object type"
         #   MM = m0d
         #else:
         print(" 0) expdisk              # Object type"   )
         MM = Md
         print(" 1) %.2f  %.2f  1 1    # position x, y        [pixel]" % (xc,yc))
         print(" 3) %.3f       1       # total magnitude"  % (MM)   )
         print(" 4) %.3f       1       #     R_s              [Pixels]" % (h))
         print(" 9) %.3f       1       # axis ratio (b/a)" % (q)  )
         print("10) %.3f       1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA))
         print(" Z) 0                  #  Skip this model in output image?  (yes=1, no=0)")
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION FOR EDGE-ON DISC
def EdgeDisc(component,xc,yc,m0d,z0,h,PA):
          print("\n# Component %i: Edge-on disk function\n")
          print(" 0) edgedisk           # Object type")
          print(" 1) %.2f  %.2f  1 1    # position x, y        [pixel]" % (xc,yc))
          print(" 3) %.3f       1       # central surface brightness  [mag/arcsec^2]" % (m0d))
          print(" 4) %.3f       1       # disk scale-height    [Pixels]" % (z0))
          print(" 5) %.3f       1       # disk scale-length    [Pixels]" % (h))
          print("10) %.3f       1       # position angle (PA)  [Degrees: Up=0, Left=90]" % (PA))
          print(" Z) 0                  #  Skip this model in output image?  (yes=1, no=0)\n")
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION FOR FERRERS FUNCTION
def Ferrer(component,xc,yc,m0d,rtr,alpha,beta,q,PA):  
         print("\n# Component %i: Ferrer function" % (component))
         print(" 0) ferrer             #  Component type")
         print(" 1) %.2f  %.2f  1 1    #  Position x, y" % (xc,yc) )
         print(" 3) %.3f       1       #  Central surface brghtness [mag/arcsec^2]" % (m0d))
         print(" 4) %.3f       1       #  Outer truncation radius  [pix]" % (rtr))
         print(" 5) %.3f       1       #  Alpha (outer truncation sharpness) " % (alpha))
         print(" 6) %.3f       1       #  Beta (central slope)" % (beta))
         print(" 9) %.3f       1       #  Axis ratio (b/a)  " % (q))
         print("10) %.3f       1       #  Position angle (PA) [deg: Up=0, Left=90]" % (PA))
         print(" Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"  )
# -----------------------------------------------------------------
 

# -----------------------------------------------------------------
# FUNCTION FOR AGN
def Agn(component,xc,yc,mag): 
         print("\n# Component %i: PSF fit" % (component))
         print(" 0) psf                # object type")
         print(" 1) %.2f  %.2f  1 1    # position x, y        [pixel]" % (xc,yc))
         print(" 3) %.3f       1       # total magnitude  "   % (mag))
         print(" Z) 0                  #  Skip this model in output image?  (yes=1, no=0)")
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION FOR SKY
def Sky(component,sky_level,fit_sky=1):
         if fit_sky==1:
           fit_sky = 0
         else:
           fit_sky = 1
         print("\n# Component %i: SKY" % (component))
         print(" 0) sky")
         print(" 1) %.5f       %i       # sky background       [ADU counts]" % (sky_level, fit_sky))
         print(" 2) 0.000      %i       # dsky/dx (sky gradient in x)" % (fit_sky))
         print(" 3) 0.000      %i       # dsky/dy (sky gradient in y)" % (fit_sky) )
         print(" Z) 0                  #  Skip this model in output image?  (yes=1, no=0)")
# -----------------------------------------------------------------
