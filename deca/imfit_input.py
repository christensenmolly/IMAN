#!/usr/bin/python
# -*- coding:  cp1251 -*-
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
import subprocess
import signal
import astropy.io.fits as pyfits

FNULL = open(os.devnull, 'w')

tmp_out = sys.stdout


import deca_setup


def read_list_of_funcs():
    subprocess.call("%simfit --list-parameters > imfit_funcs.txt" % (deca_setup.imfitPath), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    
    '''
    Imfit_functions = {}
    
    with open('imfit_funcs.txt') as f:
        lines = f.readlines()

    
    fill_func=False
    for line in lines:
      if line=='\n':
        fill_func = False
        
      if 'FUNCTION' in line and fill_func==False:
        name = line.split(' ')[1].split('\n')[0]
        fill_func = True
      if fill_func == True and 'FUNCTION' not in line:
        Imfit_functions.setdefault(name, []).append(line.split('\n')[0])
    '''

    #Imfit_functions = {}
    Imfit_functions = []
    Imfit_names = []
    with open('imfit_funcs.txt') as f:
        lines = f.readlines()
        
    
    fill_func=False
    for line in lines:
      if line=='\n':
        fill_func = False
        
      if 'FUNCTION' in line and fill_func==False:
        name = line.split(' ')[1].split('\n')[0]
        Imfit_names.append(name)
        fill_func = True
      if fill_func == True and 'FUNCTION' not in line:
        #Imfit_functions.setdefault(name+':'+line.split('\n')[0], []).append(float('nan'))
        Imfit_functions.append(name+':'+line.split('\n')[0]+':'+line.split('\n')[0])
    
    return Imfit_names,Imfit_functions


def geom_pars():
  return ['sigma','fwhm','r_c','r_t','h','r_e','r_b','h1','h2','r_break','R_ring','sigma_r','sigma_r_in','sigma_r_out','z_0','r','sigma_z','r_break','a_ring','h_z']

def lum_pars():
  return ['tot_mag']

def SB_pars():
  return ['I_0','I_e', 'I_sky','I_b','A','L_0','J_0']


def imfit_exec_line(input_files,observation_info,object_info,add_info,imfit_input_file):
        [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,sampling,convolution_box,FILTER] = observation_info
        [xc,yc,name,distance,Aext,Kext] = object_info
        [number,fitting_proc,mode,model,region_to_fit] = add_info
        [input_image,sigma_image,psf_image,mask_image] = input_files

        hdulist = pyfits.open(input_image)
        data = hdulist[0].data
        if region_to_fit!=None:                
                #  x1:x2;y1:y2 format
                xmin = int(region_to_fit.split(':')[0])  
                xmax = int(region_to_fit.split(':')[1].split(';')[0])  
                ymin = int(region_to_fit.split(';')[1].split(':')[0])  
                ymax = int(region_to_fit.split(':')[2])
        else:
                ny,nx = np.shape(data)
                xmin = 0
                xmax = nx
                ymin = 0 
                ymax = ny

        runString = "%simfit -c %s %s[0][%i:%i,%i:%i] " % (deca_setup.imfitPath, imfit_input_file, input_image, xmin+1,xmax+1,ymin+1,ymax+1)
        if psf_image != "none":
                runString += " --psf %s " % (psf_image)
        if mask_image != "none":
                runString += " --mask %s[0][%i:%i,%i:%i] " % (mask_image, xmin+1,xmax+1,ymin+1,ymax+1)
        if sigma_image != "none":
                runString += " --noise %s[0][%i:%i,%i:%i] " % (sigma_image, xmin+1,xmax+1,ymin+1,ymax+1)
        if read_out_noise != "none":
                runString += " --readnoise=%1.2f " % (read_out_noise)
        if gain != "none":
                runString += " --gain=%1.2f " % (gain)
        if SkySubtr==1 and deca_setup.add_sky_comp==True:
                runString += " --sky=%.5f " % (float(sky_level))
        if exptime != "none":
                runString += " --exptime=%1.2f " % (float(exptime))
        if ncombine != "none":
                runString += " --ncombined=%i " % (int(ncombine))
        runString += "--ftol 0.00001"
        runString += " --save-params imfit.01 "
        runString += " --save-model model.fits "
        #runString += " --save-residual %i_lm_residual.fits " % (ident) 

        return  runString     

def constrain_coors(xc,yc):
            if 'pix' in deca_setup.CENTRE_CONSTRS['MAX_OFFSET']:
                offset = int(math.ceil(float(deca_setup.CENTRE_CONSTRS['MAX_OFFSET'].split('pix')[0])))
            elif 'arcsec' in deca_setup.CENTRE_CONSTRS['MAX_OFFSET']:
                offset = int(math.ceil(float(deca_setup.CENTRE_CONSTRS['MAX_OFFSET'].split('arcsec')[0]) / arcsec))
            return xc-offset,xc+offset,yc-offset,yc+offset

    

def header(run_line,xc,yc,write_to=None):
    x_min,x_max,y_min,y_max = constrain_coors(xc,yc)
    
    if write_to==None:
        print('#',run_line,'\n')
        #print( 'X0                %.1f        %.1f,%.1f' % (xc,x_min,x_max)
        #print( 'Y0                %.1f        %.1f,%.1f' % (yc,y_min,y_max)  
    else:
        write_to.write('#',run_line,'\n')
        #print( >>write_to, 'X0                %.1f        %.1f,%.1f' % (xc,x_min,x_max)
        #print( >>write_to, 'Y0                %.1f        %.1f,%.1f' % (yc,y_min,y_max)  
      

def Sersic(component,xc,yc,PA,ell,n,I_e,r_e,m0,pix2sec,model,write_coords=True,comp_number=''):
        pars_value = [PA,ell,n,I_e,r_e]
        if deca_setup.set_constraints:
            try:
                [meb_lim, r_e_lim, n_lim, q_lim, PA_lim] = deca_setup.CONSTRS[model+':'+'sersic'+comp_number]
            except:
                [meb_lim, r_e_lim, n_lim, q_lim, PA_lim] = deca_setup.CONSTRS[model+':'+'sersic1'+comp_number]
            
            if meb_lim==None:
                I_e_lim = None
            else:
                I_e_lim = [10**(0.4*(m0-float(meb_lim[1])))*(pix2sec**2),10**(0.4*(m0-float(meb_lim[0])))*(pix2sec**2)]

            if q_lim==None:
                ell_lim = None
            else:
                ell_lim = [1.-float(q_lim[1]),1.-float(q_lim[0])]
            
            pars_lim = [PA_lim,ell_lim,n_lim,I_e_lim,r_e_lim]

            x_min,x_max,y_min,y_max = constrain_coors(xc,yc)
          

        if deca_setup.set_constraints:
            print( '\nX0                %.1f        %s,%s' % (xc,str(x_min),str(x_max)))
            print( 'Y0                %.1f        %s,%.1f' % (yc,str(y_min),str(y_max))  )
            print( "# Component %i: Sersic function" % (component))
            print( "FUNCTION Sersic")

            pars = ['PA','ell','n','I_e','r_e']
            for k in range(len(pars)):
                if pars_lim[k]==None:
                    print( "%s                %f        " % (pars[k],pars_value[k]))
                else:
                    print( "%s                %f        %f,%f" % (pars[k],pars_value[k],float(pars_lim[k][0]),float(pars_lim[k][1])))

            if write_coords==True:
                print( "X0                %.1f       %.1f,%.1f" % (xc,x_min,x_max))
                print( "Y0                %.1f       %.1f,%.1f" % (yc,y_min,y_max) )
        else:
            print( '\nX0                %.1f' % (xc))
            print( 'Y0                %.1f' % (yc)  )
            print( "# Component %i: Sersic function" % (component))
            print( "FUNCTION Sersic")

            pars = ['PA','ell','n','I_e','r_e']
            for k in range(len(pars)):
                print( "%s                %f        " % (pars[k],pars_value[k]))


            if write_coords==True:
                print( "X0                %.1f" % (xc))
                print( "Y0                %.1f" % (yc) )



def ExpDisc(component,xc,yc,PA,ell,I_0,h,m0,pix2sec,model,write_coords=True,comp_number=''):
        pars_value = [PA,ell,I_0,h]
        if deca_setup.set_constraints:
            [m0d_lim, h_lim, q_lim, PA_lim] = deca_setup.CONSTRS[model+':'+'exp_disc'+comp_number]
            
            if m0d_lim==None:
                I_0_lim = None
            else:
                I_0_lim = [10**(0.4*(m0-float(m0d_lim[1])))*(pix2sec**2),10**(0.4*(m0-float(m0d_lim[0])))*(pix2sec**2)]

            if q_lim==None:
                ell_lim = None
            else:
                ell_lim = [1.-float(q_lim[1]),1.-float(q_lim[0])]
            
            pars_lim = [PA_lim,ell_lim,I_0_lim,h_lim]

            x_min,x_max,y_min,y_max = constrain_coors(xc,yc)

        if deca_setup.set_constraints:            
            print( '\nX0                %.1f        %.1f,%.1f' % (xc,x_min,x_max))
            print( 'Y0                %.1f        %.1f,%.1f' % (yc,y_min,y_max)  )
            print( "# Component %i: Exponential function" % (component))
            print( "FUNCTION Exponential")

            pars = ['PA','ell','I_0','h']
            for k in range(len(pars)):
                if pars_lim[k]==None:
                    print( "%s                %f        " % (pars[k],pars_value[k]))
                else:
                    print( "%s                %f        %f,%f" % (pars[k],pars_value[k],float(pars_lim[k][0]),float(pars_lim[k][1])))

            if write_coords==True:
                print( "X0                %.1f       %.1f,%.1f" % (xc,x_min,x_max))
                print( "Y0                %.1f       %.1f,%.1f" % (yc,y_min,y_max) )
        else:
            print( '\nX0                %.1f' % (xc))
            print( 'Y0                %.1f' % (yc)  )
            print( "# Component %i: Exponential function" % (component))
            print( "FUNCTION Exponential")

            pars = ['PA','ell','I_0','h']
            for k in range(len(pars)):
                    print( "%s                %f        " % (pars[k],pars_value[k]))


            if write_coords==True:
                print( "X0                %.1f" % (xc))
                print( "Y0                %.1f" % (yc) )


#def BrokenExpDisc(component,xc,yc,m0d,Md,h,q,PA,m0,pix2sec,write_coords=True):
def BrokenExpDisc(component,xc,yc,PA,ell,I_0,h,m0,pix2sec,model,write_coords=True,comp_number=''):
        pars_value = [PA,ell,I_0,h,h,4.*h,1.]
        if deca_setup.set_constraints:    
            [m0d_lim, h_lim, q_lim, PA_lim] = deca_setup.CONSTRS[model+':'+'exp_disc'+comp_number]
            
            if m0d_lim==None:
                I_0_lim = None
            else:
                I_0_lim = [10**(0.4*(m0-float(m0d_lim[1])))*(pix2sec**2),10**(0.4*(m0-float(m0d_lim[0])))*(pix2sec**2)]

            if q_lim==None:
                ell_lim = None
            else:
                ell_lim = [1.-float(q_lim[1]),1.-float(q_lim[0])]
            
            h1_lim = h_lim
            h2_lim = h_lim
            r_break = None
            alpha_lim = None
            
            pars_lim = [PA_lim,ell_lim,I_0_lim,h1_lim,h2_lim,r_break_lim,alpha_lim]
            
            x_min,x_max,y_min,y_max = constrain_coors(xc,yc)
        
        if deca_setup.set_constraints:    
            print( '\nX0                %.1f        %.1f,%.1f' % (xc,x_min,x_max))
            print( 'Y0                %.1f        %.1f,%.1f' % (yc,y_min,y_max) ) 
            print( "# Component %i: Exponential function" % (component))
            print( "FUNCTION  BrokenExponential")

            pars = ['PA','ell','I_0','h1','h2','r_break','alpha']
            for k in range(len(pars)):
                if pars_lim[k]==None:
                    print( "%s                %f        " % (pars[k],pars_value[k]))
                else:
                    print( "%s                %f        %f,%f" % (pars[k],pars_value[k],float(pars_lim[k][0]),float(pars_lim[k][1])))
                    
            if write_coords==True:
                print( "X0                %.1f       %.1f,%.1f" % (xc,x_min,x_max))
                print( "Y0                %.1f       %.1f,%.1f" % (yc,y_min,y_max))
        else:          
            print( '\nX0                %.1f' % (xc))
            print( 'Y0                %.1f' % (yc)  )
            print( "# Component %i: Exponential function" % (component))
            print( "FUNCTION  BrokenExponential")

            pars = ['PA','ell','I_0','h1','h2','r_break','alpha']
            for k in range(len(pars)):
                    print( "%s                %f        " % (pars[k],pars_value[k]))

                    
            if write_coords==True:
                print( "X0                %.1f" % (xc))
                print( "Y0                %.1f" % (yc))


def EdgeDisc(component,xc,yc,PA,L_0,h,z_0,m0,pix2sec,model,write_coords=True,comp_number=''):
        n = deca_setup.imfit_disc_n
        
        pars_value = [PA,L_0,h,n,z_0]
        if deca_setup.set_constraints:  
            [m0d_lim, h_lim, PA_lim] = deca_setup.CONSTRS[model+':'+'eon_disc'+comp_number]
            
            if m0d_lim==None:
                L_0_lim = None
            else:
                L_0_lim = [10**(0.4*(m0-float(m0d_lim[1])))*(pix2sec**2)/(2.*h),10**(0.4*(m0-float(m0d_lim[0])))*(pix2sec**2)/(2.*h)]

            n_lim = [n,n]
            z_0_lim = None
            
            pars_lim = [PA_lim,L_0_lim,h_lim,n_lim,z_0_lim]

            x_min,x_max,y_min,y_max = constrain_coors(xc,yc)
        if deca_setup.set_constraints:              
            print( '\nX0                %.1f        %.1f,%.1f' % (xc,x_min,x_max))
            print( 'Y0                %.1f        %.1f,%.1f' % (yc,y_min,y_max)  )
            print( "# Component %i: Edge-on disk function\n")
            print( "FUNCTION EdgeOnDisk")

            pars = ['PA','L_0','h','n','z_0']
            for k in range(len(pars)):
                if pars_lim[k]==None:
                    print( "%s                %f        " % (pars[k],pars_value[k]))
                else:
                    print( "%s                %f        %f,%f" % (pars[k],pars_value[k],float(pars_lim[k][0]),float(pars_lim[k][1])))
            
            if write_coords==True:
                print( "X0                %.1f       %.1f,%.1f" % (xc,x_min,x_max))
                print( "Y0                %.1f       %.1f,%.1f" % (yc,y_min,y_max) )
        else:              
            print( '\nX0                %.1f' % (xc))
            print( 'Y0                %.1f' % (yc)  )
            print( "# Component %i: Edge-on disk function\n")
            print( "FUNCTION EdgeOnDisk")

            pars = ['PA','L_0','h','n','z_0']
            for k in range(len(pars)):
                    print( "%s                %f        " % (pars[k],pars_value[k]))
            
            if write_coords==True:
                print( "X0                %.1f" % (xc))
                print( "Y0                %.1f" % (yc) )

def GaussRing(component,xc,yc,PA,ell,A,R_ring,sigma_r,m0,pix2sec,model,write_coords=True,comp_number=''):
        
        pars_value = [PA,ell,A,R_ring,sigma_r]
        if deca_setup.set_constraints:  
            [m0d_lim, sigma_r_lim, R_ring_lim, q_lim, PA_lim] = deca_setup.CONSTRS[model+':'+'ring'+comp_number]
            
            if m0d_lim==None:
                A_lim = None
            else:
                A_lim = [10**(0.4*(m0-float(m0d_lim[1])))*(pix2sec**2),10**(0.4*(m0-float(m0d_lim[0])))*(pix2sec**2)]


            if q_lim==None:
                ell_lim = None
            else:
                ell_lim = [1.-float(q_lim[1]),1.-float(q_lim[0])]
            
            pars_lim = [PA_lim,ell_lim,A_lim,R_ring_lim,sigma_r_lim]
            
            x_min,x_max,y_min,y_max = constrain_coors(xc,yc)
        if deca_setup.set_constraints:  
            print( '\nX0                %.1f        %.1f,%.1f' % (xc,x_min,x_max))
            print( 'Y0                %.1f        %.1f,%.1f' % (yc,y_min,y_max) ) 
            print( "# Component %i: GaussianRing function\n" % (component))
            print( "FUNCTION GaussianRing")

            pars = ['PA','ell','A','R_ring','sigma_r']
            for k in range(len(pars)):
                if pars_lim[k]==None:
                    print( "%s                %f        " % (pars[k],pars_value[k]))
                else:
                    print( "%s                %f        %f,%f" % (pars[k],pars_value[k],float(pars_lim[k][0]),float(pars_lim[k][1])))


            if write_coords==True:
                print( "X0                %.1f       %.1f,%.1f" % (xc,x_min,x_max))
                print( "Y0                %.1f       %.1f,%.1f" % (yc,y_min,y_max) )
        else:  
            print( '\nX0                %.1f' % (xc))
            print( 'Y0                %.1f' % (yc)  )
            print( "# Component %i: GaussianRing function\n" % (component))
            print( "FUNCTION GaussianRing")

            pars = ['PA','ell','A','R_ring','sigma_r']
            for k in range(len(pars)):
                    print( "%s                %f        " % (pars[k],pars_value[k]))

            if write_coords==True:
                print( "X0                %.1f" % (xc))
                print( "Y0                %.1f" % (yc) )


#def Ferrer(component,xc,yc,m0d,rtr,alpha,beta,q,PA):  
#         x_min,x_max,y_min,y_max = constrain_coors(xc,yc)
#         print( "\n# Component %i: Ferrer function" % (component)
#         print( " 0) ferrer             #  Component type"
#         print( " 1) %.2f  %.2f  1 1    #  Position x, y" % (xc,yc) 
#         print( " 3) %.3f       1       #  Central surface brghtness [mag/arcsec^2]" % (m0d)
#         print( " 4) %.3f       1       #  Outer truncation radius  [pix]" % (rtr)
#         print( " 5) %.3f       1       #  Alpha (outer truncation sharpness) " % (alpha)
#         print( " 6) %.3f       1       #  Beta (central slope)" % (beta)
#         print( " 9) %.3f       1       #  Axis ratio (b/a)  " % (q)
#         print( "10) %.3f       1       #  Position angle (PA) [deg: Up=0, Left=90]" % (PA)
#         print( " Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"  
  
#def Agn(component,xc,yc,mag,model):                 
#         x_min,x_max,y_min,y_max = constrain_coors(xc,yc)
#         print( "\n# Component %i: PSF fit" % (component)
#         print( " 0) psf                # object type"
#         print( " 1) %.2f  %.2f  1 1    # position x, y        [pixel]" % (xc,yc)
#         print( " 3) %.3f       1       # total magnitude  "   % (mag)
#         print( " Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"

def Sky(component,sky_level,fit_sky=1):
        print( "\n# Component %i: SKY" % (component))
        print( "FUNCTION FlatSky")
        if fit_sky==1:
          print( "I_sky                %.5f        fixed" % (sky_level))
        else:
          print( "I_sky                %.5f             " % (sky_level))
        
#read_list_of_funcs()
