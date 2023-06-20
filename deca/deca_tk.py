#!/usr/bin/env python
# *****************************************************************
# **    DECA-TK: New Python script with interactive  interface   **
# **   Astronomical Institute, Saint Petersburge State University**
# **       Main developers: A. Mosenkov and S. Savchenko         **
# **       mosenkovAV@gmail.com, savchenko.s.s@gmail.com         **
# *****************************************************************
# BRIEF DESCRIPTION:
# This script uses the Python TKinter package to create a friendly
# graphical interface to perform multicomponent decomposition analysis
# based on galaxy images. 
# REQUIREMENTS: 
# This script should be called from the directory where one wants to
# work (later, working directory). All temporary files will be created
# in working directory. If possible, clean this direcotory beforehand
# (except for files which will be given to DECA as the input). 
# EXAMPLE: >python deca_tk.py
#
# -----------------------------------------------------------------
#*** Common modules ***
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
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from numpy import *
from pylab import *
import os
import shutil
import random
import time
import pyfits
import re
import shelve
import threading
import inspect
import subprocess as sub
from collections import OrderedDict
from collections import Counter
import glob
import subprocess

#*** TKinter modules ***
import Pmw
from Tkinter import *
import Tkinter as Tk
import tkFileDialog
import tkMessageBox



#*** IMAN modules ***
PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
sys.path.append(PATH_TO_PACKAGE.split('/DECA')[0]+'/Decomposition')
sys.path.append(PATH_TO_PACKAGE.split('/DECA')[0]+'/Plot')
sys.path.append(PATH_TO_PACKAGE.split('/DECA')[0]+'/FindFluxes')
import make_model_ima_galfit
import iraf_ell_plots
import plot_2d_profile
from misc_functions import * 

import deca_setup
galfitPath = deca_setup.galfitPath
imfitPath = deca_setup.imfitPath
show_code_run = deca_setup.show_code_run

FNULL = open(os.devnull, 'w')

Current_Dir = os.getcwd()

# We are using this directory where you've created input files
configdir = '.'
sys.path.append(configdir)
DECA_PATH = os.path.dirname(__file__)
sys.path.append(DECA_PATH+'/deca_tk_lib')


#*** DECA modules ***
import galfit_input
import imfit_input
import galfit_parser
import tex_creator

def filters(Filter):
	# from http://mips.as.arizona.edu/~cnaw/sun.html
	# B&M 1998
	# [Vega,AB]
	if Filter==365.:#U
		Msun = [5.61,6.32]
	if Filter==445.:#B
		Msun = [5.48,5.36]
	if Filter==551.:#V
		Msun = [4.83,4.79]
	if Filter==658.:#R
		Msun = [4.42,4.65]
	if Filter==806.:#I
		Msun = [4.08,4.55]
	if Filter==1220:#J
		Msun = [3.64,4.57]
	if Filter==1630.:#H
		Msun = [3.32,4.71]
	if Filter==2159.:#Ks
		Msun = [3.27,5.18]
	if Filter==2190.:#K
		Msun = [3.28,5.19]
	if Filter==3450:#L
		Msun = [3.25,5] #AB-???

	#SDSS-bands:
	if Filter==354.:#u
		Msun = [5.46,6.45]
	if Filter==475.:#g
		Msun = [5.22,5.14]
	if Filter==622:#r
		Msun = [4.50,4.65]
	if Filter==763.:#i
		Msun = [4.16,4.54]
	if Filter==905.:#z
		Msun = [4.01,4.52]

	# Spitzer
	if Filter==3600.:#3.6mkm
		Msun = [3.24,3.24] #AB-???
	if Filter==4600.:#4.6mkm
		Msun = [3.27,3.27] #AB-???
		
	if deca_setup.mag_system=='Vega':
	  Msun = Msun[0]
	else:
	  Msun = Msun[1]
	return Msun

def units_converter(value,par_units,input_geom_units,input_SB_units,input_lum_units,output_geom_units,
		    output_SB_units,output_lum_units,m0,arcsec_per_pixel,kpc_per_arcsec,Filter,Aext=0.,Kcorr=0.):
  # Function to convert between units
  #print value,par_units,input_geom_units,input_SB_units,input_lum_units,output_geom_units,output_SB_units,output_lum_units
  Distance = kpc_per_arcsec*206265./1000.
  Msun = filters(Filter)	####TODO: Make a continous function to find Msun
  
  if par_units=='geom':
	if input_geom_units=='pix' and input_geom_units!=output_geom_units:
	  if output_geom_units=='arcsec':
	    value = value * arcsec_per_pixel
	  if output_geom_units=='kpc':
	    value = value * arcsec_per_pixel * kpc_per_arcsec
	    
	if input_geom_units=='arcsec' and input_geom_units!=output_geom_units:
	  if output_geom_units=='pix':
	    value = value / arcsec_per_pixel
	  if output_geom_units=='kpc':
	    value = value * kpc_per_arcsec    
	    
	if input_geom_units=='kpc' and input_geom_units!=output_geom_units:
	  if output_geom_units=='pix':
	    value = value / kpc_per_arcsec / arcsec_per_pixel
	  if output_geom_units=='arcsec':
	    value = value / kpc_per_arcsec
	    
  elif par_units=='SB':
	if input_SB_units=='mag/arcsec2' and input_SB_units!=output_SB_units:
	  if output_SB_units=='ADU/pix2':
	    value = 10**(0.4*(m0 - value + 5.*log10(arcsec_per_pixel) + Aext + Kcorr))
	  if output_SB_units=='Lsun/pc2':
	    value = 10**(0.4*(Msun - value + 21.572 + Aext + Kcorr))
	    
	if input_SB_units=='ADU/pix2' and input_SB_units!=output_SB_units:
	  if output_SB_units=='mag/arcsec2':
	    value = m0 -2.5*log10(value) + 5.*log10(arcsec_per_pixel) - Aext - Kcorr
	  if output_SB_units=='Lsun/pc2':
	    value = 10**(0.4*(Msun+21.572-m0 -2.5*log10(value) + 5.*log10(arcsec_per_pixel) - Aext - Kcorr))
	
	if input_SB_units=='Lsun/pc2' and input_SB_units!=output_SB_units:
	  if output_SB_units=='mag/arcsec2':
	    value = Msun + 21.572 - 2.5*log10(value) - Aext - Kcorr
	  if output_SB_units=='ADU/pix2':
	    value = 10**(0.4*(m0 - Msun + 21.572 - 2.5*log10(value) - Aext - Kcorr + 5.*log10(arcsec_per_pixel)))

    
  elif par_units=='lum':
	if input_lum_units=='mag' and input_lum_units!=output_lum_units:
	  if output_lum_units=='ADU':
	    value = 10**(0.4*(m0 - value + Aext + Kcorr))
	  if output_lum_units=='Lsun':
	    value = value - 5.*log10(Distance) - 25.
	    value = 10**(0.4*(Msun - value))
	    
	if input_lum_units=='ADU' and input_lum_units!=output_lum_units:
	  if output_lum_units=='mag':
	    value = m0 -2.5*log10(value) - Aext - Kcorr
	  if output_lum_units=='Lsun':
	    value = m0 -2.5*log10(value) - Aext - Kcorr
	    value = value - 5.*log10(Distance) - 25.
	    value = 10**(0.4*(Msun - value))
	    
	if input_lum_units=='Lsun' and input_lum_units!=output_lum_units:
	  if output_lum_units=='ADU':
	    value = Msun - 2.5*log10(value)
	    value = value + 5.*log10(Distance) + 25.
	    value = 10**(0.4*(m0 - value + Aext + Kcorr))
	  if output_lum_units=='mag':    
	    value = Msun - 2.5*log10(value)
	    value = value + 5.*log10(Distance) + 25.

  else:
	value = value
  return value






# -----------------------------------------------------------------
# Class to choose components for a new model
class Checkbar(Frame):
   def __init__(self, parent=None, picks=[], load_picks=[], side=TOP, anchor=W):
      Frame.__init__(self, parent)
      self.vars = OrderedDict()
      self.cbVars = OrderedDict()
      self.vars_numbs = []
      
      Tk.Label(self, text="Component", fg="black").grid(column=1, row=0,sticky='w')
      Tk.Label(self, text="No", fg="black").grid(column=2, row=0,sticky='w')
      
      if load_picks!=[]:
	#Count duplicates:
	dupl_load_picks = dict(Counter(load_picks))
	

      for k in range(len(picks)):
	 pick = picks[k]
         var = IntVar()
         cbVar = IntVar()
         
	 if load_picks!=[]:
	   if pick in dupl_load_picks:
	     cbVar.set(1)
	     var.set(dupl_load_picks[pick])
 

         chk = Checkbutton(self, text=pick, fg='blue', variable=cbVar, offvalue=0, onvalue=1)
	 chk.grid(column=1,row=k+1, sticky=Tk.W)
         cbVar.trace("w", lambda n,m,x: self.vars[n].set(self.cbVars[n].get()))

	 NumOfCompsEntry = Tk.Entry(self, textvariable=var, width=5, state="normal", bg="white")
	 NumOfCompsEntry.grid(column=2, row=k+1, sticky=Tk.W)
         
         self.vars[cbVar._name] = var
         self.vars_numbs.append(var)
         self.cbVars[cbVar._name] = cbVar

   def state_numb(self):
       return [numb.get() for numb in self.vars.values()]
# -----------------------------------------------------------------




# -----------------------------------------------------------------
# Class to create new model from the selected components
class Parbar(Frame):
   def __init__(self, new_images, observation_info, object_info, add_info, settings, parent=None, Compnames=[], Comppars=[], side=TOP, anchor=W):
      
      Frame.__init__(self, parent)
      self.canvas = Tk.Canvas(parent, borderwidth=0)
      self.frame = Tk.Frame(self.canvas)
      self.vsb = Tk.Scrollbar(parent, orient="vertical", command=self.canvas.yview)
      self.canvas.configure(yscrollcommand=self.vsb.set)
      
      self.vsb.pack(side="right", fill="y")
      self.canvas.pack(side="left", fill="both", expand=True)
      self.canvas.create_window((6,12), window=self.frame, anchor="nw", tags="self.frame")
      
      self.frame.bind("<Configure>", self.onFrameConfigure)     
      
      self.FinalParameters = OrderedDict()
      self.Compnames = Compnames

      self.window = parent
      self.window.title('Peek model')
      self.window.geometry(("500x600"))

      self.observation_info = observation_info
      self.object_info = object_info
      self.add_info = add_info
      self.new_images = new_images
      self.settings = settings
      #[new_images,ADD_INFO,RA,DEC,galaxy_name,Sampling,Sky_level,SkySubtr,ChooseCode,ChooseGeom,ChooseSB,ChooseLum,Scale] = self.input_info

      [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,Sky_level,SkySubtr,Filter] = observation_info
      [RA,DEC,galaxy_name,distance,Aext,Kext] = object_info
      [number,fitting_proc,mode,model] = add_info
      [input_image,sigma_image,psf_image,mask_image,Sampling] = new_images
      [ChooseCode,ChooseGeom,ChooseSB,ChooseLum] = settings
      
      if ChooseCode=='GALFIT':
	  self.geom_pars = galfit_input.geom_pars()
	  self.lum_pars = galfit_input.lum_pars()
	  self.SB_pars = galfit_input.SB_pars()
      elif ChooseCode=='IMFIT':
	  self.geom_pars = imfit_input.geom_pars()
	  self.lum_pars = imfit_input.lum_pars()
	  self.SB_pars = imfit_input.SB_pars()
	  
      line = 1
      for k in range(len(self.Compnames)):
	 func_name = self.Compnames[k].split('/')[0]
	 Tk.Label(self.frame, text=func_name+' #'+self.Compnames[k].split('/')[1],fg="blue",font=("Helvetica", 13)).grid(column=0, row=line,sticky='w')
	 
	 if k==0:
		  line = line + 1
	          Tk.Label(self.frame, text='Param',fg="black").grid(column=0, row=line,sticky='w')
	          if ChooseCode=='GALFIT':
		    Tk.Label(self.frame, text='Description',fg="black").grid(column=1, row=line,sticky='w')
		  Tk.Label(self.frame, text='Units',fg="black").grid(column=2, row=line,sticky='w')
		  Tk.Label(self.frame, text='value',fg="black").grid(column=3, row=line,sticky='w')
		  Tk.Label(self.frame, text='min',fg="black").grid(column=4, row=line)
		  Tk.Label(self.frame, text='max',fg="black").grid(column=5, row=line)
		  
	 line = line + 1
	 for i in range(len(Comppars)):
	    if func_name==Comppars[i].split(':')[0] and self.Compnames[k].split('/')[1]==Comppars[i].split('/')[1].split(':')[0]:
		par_name = Comppars[i].split(':')[1]
		com_name =  Comppars[i].split(':')[2].split('/')[0]

		if (par_name=='1' or par_name=='2' or par_name=='X0' or par_name=='Y0') and float(Comppars[i].split(':')[3].split(',')[0])==0.:
		  if par_name=='1':
		    par_var = Tk.DoubleVar(value=float(RA))
		  elif par_name=='2':
		    par_var = Tk.DoubleVar(value=float(DEC))
		  elif  par_name=='X0':
		    par_var = Tk.DoubleVar(value=float(RA))
		  elif par_name=='Y0':
		    par_var = Tk.DoubleVar(value=float(DEC))
		else:  
		  par_var = Tk.DoubleVar(value=float(Comppars[i].split(':')[3].split(',')[0]))
		
		par_min_var = Tk.StringVar(value=str(Comppars[i].split(':')[3].split(',')[1]))
		par_max_var = Tk.StringVar(value=str(Comppars[i].split(':')[3].split(',')[2]))
		unit = ''

		if com_name in self.geom_pars:
		  unit = ChooseGeom
		elif com_name in self.lum_pars:
		  unit = ChooseLum
		elif com_name in self.SB_pars:
		  #unit = ChooseLum+'/'+ChooseGeom+u'\u00b2'
		  unit = ChooseSB
		elif com_name=='PA':
		  unit = 'deg'
		elif com_name=='x_cen' or com_name=='y_cen':
		  unit = 'pix'
		  
		# Name of the parameters   
		name = Tk.Label(self.frame, text=par_name,fg="green").grid(column=0, row=line,sticky='w')
		if ChooseCode=='GALFIT':
		  # Description
		  com = Tk.Label(self.frame, text=com_name,fg="green").grid(column=1, row=line,sticky='w' )

		# Units
		unitss = Tk.Label(self.frame, text=unit,fg="green").grid(column=2, row=line,sticky='w' )

		par_value = Tk.Entry(self.frame, textvariable=par_var, width=5, state="normal", bg="white")
		par_value.grid(column=3, row=line, sticky=Tk.W)
		if com_name=='PA':
		  par_min_var.set('0.')
		  par_max_var.set('360.')
		  
		par_min_value = Tk.Entry(self.frame, textvariable=par_min_var, width=5, state="normal", bg="white")
		par_min_value.grid(column=4, row=line, sticky=Tk.W)

		par_max_value = Tk.Entry(self.frame, textvariable=par_max_var, width=5, state="normal", bg="white")
		par_max_value.grid(column=5, row=line, sticky=Tk.W)

		self.FinalParameters.update({Comppars[i].split(':')[0]+':'+Comppars[i].split(':')[1]+':'+Comppars[i].split(':')[2]:[par_var,par_min_var,par_max_var]})
		line = line + 1

      CreateModel = Tk.Frame(self.frame)
      CreateModel.grid(column=0, row=line+2)
      CreateModelButton = Tk.Button(CreateModel, text="Create Input", state="normal",command=self.create_model, bg="green",font=("Helvetica", 10, "bold"), width=10)
      CreateModelButton.grid(column=0, row=line+2)   

      PlotModel = Tk.Frame(self.frame)
      PlotModel.grid(column=1, row=line+2)
      PlotModelButton = Tk.Button(PlotModel, text="Create Model", state="normal",command=self.build_model, bg="green",font=("Helvetica", 10, "bold"), width=10)
      PlotModelButton.grid(column=1, row=line+2)  
      '''
      QuitModel = Tk.Frame(self.frame)
      QuitModel.grid(column=1, row=line+3)
      QuitModelButton = Tk.Button(QuitModel, text="   Quit   ", state="normal",command=self.quit_model, bg="red",font=("Helvetica", 10, "bold"), width=10)
      QuitModelButton.grid(column=1, row=line+3)  
      '''   

   def onFrameConfigure(self, event):
	self.canvas.configure(scrollregion=self.canvas.bbox("all"))

   def quit_model(self):
	self.window.destroy()

   def build_model(self):
	#[new_images,ADD_INFO,RA,DEC,galaxy_name,Sampling,Sky_level,SkySubtr,ChooseGeom,ChooseSB,ChooseLum,Scale] = self.input_info
	[nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,Sky_level,SkySubtr,Filter] = self.observation_info
	[RA,DEC,galaxy_name,distance,Aext,Kext] = self.object_info
	[number,fitting_proc,mode,model] = self.add_info
	[input_image,sigma_image,psf_image,mask_image,Sampling] = self.new_images
	[ChooseCode,ChooseGeom,ChooseSB,ChooseLum] = self.settings

	make_model_ima_galfit.main(self.new_images[0], 'galfit.inp')
   
   def create_model(self):
     [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,Sky_level,SkySubtr,Filter] = self.observation_info
     [RA,DEC,galaxy_name,distance,Aext,Kext] = self.object_info
     [number,fitting_proc,mode,model] = self.add_info
     [input_image,sigma_image,psf_image,mask_image,Sampling] = self.new_images
     [ChooseCode,ChooseGeom,ChooseSB,ChooseLum] = self.settings
     Scale = distance * 1000. / 206265.
     
     # Save tmp file:
     TMP_PARS = []
     for key, value in self.FinalParameters.iteritems() :
       TMP_PARS.append(key+':'+str(value[0].get())+','+str(value[1].get())+','+str(value[2].get()))

     tmpfile = open('tmp_model.npy', 'w')
     np.save(tmpfile, TMP_PARS)
     tmpfile.close()
     
     if ChooseCode=='IMFIT':
	f = open('imfit.inp', "w")
	
        run_line = imfit_input.imfit_exec_line(self.new_images,self.observation_info,self.object_info,'imfit.inp') 

	imfit_input.header(run_line,RA,DEC,f)
	
	
	for k in range(len(self.Compnames)):
	  func_name = self.Compnames[k].split('/')[0]
	  numb_comp = self.Compnames[k].split('/')[1]
	  print >>f, '\n'
	  write_func_header=False
	  for key, value in self.FinalParameters.iteritems() :
	      if key.split(':')[0]==func_name and key.split('/')[1]==numb_comp:
		com_name = key.split(':')[2].split('/')[0]
		if com_name in self.geom_pars:
		  par_units = 'geom'
		elif com_name in self.lum_pars:
		  par_units = 'lum'
		elif com_name in self.SB_pars:
		  par_units = 'SB'
		else:
		  par_units = ''
		  
		if str(value[1].get())=='' or str(value[2].get())=='':
		  add_string = ''
		if str(value[1].get())=='*' or str(value[2].get())=='*':
		  add_string = '    fixed'
		if str(value[1].get())!='' and str(value[2].get())!='' and str(value[1].get())!='*' and str(value[2].get())!='*':
		  VALUE_MIN = units_converter(float(value[1].get()),par_units,ChooseGeom,ChooseSB,ChooseLum,'pix','ADU/arcsec2','ADU',m0,scale,Scale,Filter)
		  VALUE_MAX = units_converter(float(value[2].get()),par_units,ChooseGeom,ChooseSB,ChooseLum,'pix','ADU/arcsec2','ADU',m0,scale,Scale,Filter)
		  add_string = '    ' + str(VALUE_MIN) + ',' + str(VALUE_MAX)
		VALUE = units_converter(float(value[0].get()),par_units,ChooseGeom,ChooseSB,ChooseLum,'pix','ADU/arcsec2','ADU',m0,scale,Scale,Filter)
		if write_func_header==True:
		  print >>f, '# Component %s no. %s' % (func_name,numb_comp)
		  print >>f, 'Function %s' % (func_name)
		  write_func_header=False
		print >>f, key.split(':')[1]+'\t\t'+str(VALUE) + add_string
		if key.split(':')[1]=='Y0':
		  write_func_header=True
	
	if deca_setup.add_sky_comp==True and SkySubtr==0:
	    print >>f, "\n# Component %i: SKY" % (len(self.Compnames)+1)
	    print >>f, "FUNCTION FlatSky"
	    print >>f, "I_sky		%.5f        fixed" % (Sky_level)    
	  
	f.close()


     if ChooseCode=='GALFIT':
	f = open('galfit.inp', "w") 
        f_c = open('constraint.txt', "w") 

	#[input_image,sigma_image,psf_image,mask_image,Sampling] = new_images
	#[nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level] = ADD_INFO

	hdulist = pyfits.open(input_image)
	data = hdulist[0].data
	ymax,xmax = np.shape(data)
	galfit_input.header(input_image,sigma_image,psf_image,mask_image,'constraint.txt','model.fits',1,xmax,1,ymax,m0-2.5*log10(exptime),scale,Sampling,f)
	for k in range(len(self.Compnames)):
	  func_name = self.Compnames[k].split('/')[0]
	  numb_comp = self.Compnames[k].split('/')[1]
	  print >>f, '\n# Component %s no. %s' % (func_name,numb_comp)
	  print >>f, '0) %s' % (func_name)

	  
	  for key, value in self.FinalParameters.iteritems() :
	      if key.split(':')[0]==func_name and key.split('/')[1]==numb_comp:
		com_name = key.split(':')[2].split('/')[0]
		if com_name in self.geom_pars:
		  par_units = 'geom'
		elif com_name in self.lum_pars:
		  par_units = 'lum'
		elif com_name in self.SB_pars:
		  par_units = 'SB'  
		else:
		  par_units = ''
  
		if str(value[1].get())=='' or str(value[2].get())=='':
		  add_string = '    1'
		if str(value[1].get())=='*' or str(value[2].get())=='*':
		  add_string = '    0'
		if str(value[1].get())!='' and str(value[2].get())!='' and str(value[1].get())!='*' and str(value[2].get())!='*':
		  add_string = '    1'
		  VALUE_MIN = units_converter(float(value[1].get()),par_units,ChooseGeom,ChooseSB,ChooseLum,'pix','mag/arcsec2','mag',m0,scale,Scale,Filter)
		  VALUE_MAX = units_converter(float(value[2].get()),par_units,ChooseGeom,ChooseSB,ChooseLum,'pix','mag/arcsec2','mag',m0,scale,Scale,Filter)  
		  print >>f_c, '     %s              %s        %s to %s    ' % ( k+1,key.split(':')[1].split('/')[0],str(VALUE_MIN),str(VALUE_MAX) )
		if int(key.split(':')[1])==1:
		  coord_x = key.split(':')[1]+') '+str(value[0].get())
		  add_string_x = add_string
		elif int(key.split(':')[1])==2:
		  print >>f, coord_x + '   ' + str(value[0].get()) + add_string_x + add_string
		else:
		  VALUE = units_converter(float(value[0].get()),par_units,ChooseGeom,ChooseSB,ChooseLum,'pix','mag/arcsec2','mag',m0,scale,Scale,Filter)
		  print >>f, key.split(':')[1]+') '+str(value[0].get()) + add_string

	# Add sky if needed:
	if deca_setup.add_sky_comp==True:
	    print >>f, "\n# Component %i: SKY" % (len(self.Compnames)+1)
	    print >>f, "  0) sky"
	    print >>f, "  1) %.5f       0       # sky background       [ADU counts]" % (Sky_level)
	    print >>f, "  2) 0.000      0       # dsky/dx (sky gradient in x)" 
	    print >>f, "  3) 0.000      0       # dsky/dy (sky gradient in y)" 
	    print >>f, "  Z) 0                  #  Skip this model in output image?  (yes=1, no=0)"

	f.close()
	f_c.close()
# -----------------------------------------------------------------

  
  
  
# -----------------------------------------------------------------
# Class to create DECA-TK application (MAIN)
class MainApp(object):
    def __init__(self, input_par_file):
	#input_par_file = self.input_par_file
	master=self.master=Tk.Tk()
	master.title("DECA-TK")
	master.geometry(("1400x900"))
	master.configure(background='grey')
	try:
            img_icon = PhotoImage(file=DECA_PATH+'/deca_tk_lib/deca_icon.png')
            master.tk.call('wm', 'iconphoto', master._w, img_icon)
        except:
            z=1
	
	# Create layers
	Pmw.initialise()
	nb = Pmw.NoteBook(master)
	self.p1 = nb.add('Input data')
	self.p2 = nb.add('Isophotes')
	self.p3 = nb.add('1D-Profiles')
	self.p4 = nb.add('2D-Profiles')
	self.p5 = nb.add('Ellipse fitting')
	self.p6 = nb.add('Decomposition')
	nb.pack(padx=5, pady=5, fill=BOTH, expand=1)
	
	nb.tab('Input data').configure(font=("Helvetica", 12))#,fg='blue')
	nb.tab('Isophotes').configure(font=("Helvetica", 12))#,fg='blue')
	nb.tab('1D-Profiles').configure(font=("Helvetica", 12))#,fg='blue')
	nb.tab('2D-Profiles').configure(font=("Helvetica", 12))#,fg='blue')
	nb.tab('Ellipse fitting').configure(font=("Helvetica", 12))#,fg='blue')
	nb.tab('Decomposition').configure(font=("Helvetica", 12))#,fg='blue')


	############################################################
	#                Create first file                         #
	############################################################
	from deca_tk_lib import input_file
	input_file.beginning(master)


	############################################################
	#                  File 1: Input data                      #
	############################################################
	#self.nx1,self.ny1,self.file_field1,self.weight_image1,self.psf_image1,\
	#self.mask_image1,self.m01,self.pix2sec1,self.GAIN1,self.NCOMBINE1,self.EXPTIME1,\
	#self.RA1,self.Dec1,self.galaxy_name1,self.RON1,self.FWHM1,self.Sampling1,self.Sky_level1,self.CoordsCheck = input_file.inp(self.p1)
	input_files,observation_info,object_info,settings = input_file.inp(self.p1,input_par_file)
	[self.file_field1,self.weight_image1,self.psf_image1,self.mask_image1,self.Sampling1] = input_files
	[self.nx1,self.ny1,self.m01,self.pix2sec1,self.GAIN1,self.NCOMBINE1,self.EXPTIME1,self.RON1,self.FWHM1,self.Sky_level1,self.SkySubtr1,self.Filter1] = observation_info
	[self.RA1,self.Dec1,self.CoordsCheck,self.galaxy_name1,self.Dist1,self.Ext1,self.Kcorr1] = object_info
	[self.ChooseCode,self.ChooseGeom,self.ChooseSB,self.ChooseLum] = settings


	############################################################
	#                  File 2: Isophotes                       #
	############################################################ 

	############################################################
	#                  Creating right panel                    #
	############################################################

	rightPanelp2 = Tk.Frame(self.p2)
	rightPanelp2.pack(side=Tk.RIGHT ,anchor="n")
	
	# Plot isophotes
	self.min_mag_iso = Tk.DoubleVar()
	self.max_mag_iso = Tk.DoubleVar()
	self.step_iso = Tk.DoubleVar()
	self.color_iso = Tk.IntVar()
	self.units_iso = Tk.StringVar()
	
	self.min_mag_iso.set(20.0)
	self.max_mag_iso.set(25.5)
	self.step_iso.set(0.1)
	self.color_iso.set(1)
	self.units_iso.set(self.ChooseGeom.get())
	
	Tk.Label(rightPanelp2, text="Inner isophote",fg="blue").grid(column=0, row=0, sticky=Tk.W)
	InnerIsoEntry = Tk.Entry(rightPanelp2, textvariable=self.min_mag_iso, width=5, state="normal", bg="white")
	InnerIsoEntry.grid(column=1, row=0, sticky=Tk.W)
	Tk.Label(rightPanelp2, text=u"mag/\u25a1\"").grid(column=2, row=0) 

	Tk.Label(rightPanelp2, text="Outer isophote",fg="blue").grid(column=0, row=1, sticky=Tk.W)
	OuterIsoEntry = Tk.Entry(rightPanelp2, textvariable=self.max_mag_iso, width=5, state="normal", bg="white")
	OuterIsoEntry.grid(column=1, row=1, sticky=Tk.W)
	Tk.Label(rightPanelp2, text=u"mag/\u25a1\"").grid(column=2, row=1) 

	Tk.Label(rightPanelp2, text="Step",fg="blue").grid(column=0, row=2, sticky=Tk.W)
	StepIsoEntry = Tk.Entry(rightPanelp2, textvariable=self.step_iso, width=5, state="normal", bg="white")
	StepIsoEntry.grid(column=1, row=2, sticky=Tk.W)
	Tk.Label(rightPanelp2, text=u"mag/\u25a1\"").grid(column=2, row=2) 

	chk = Checkbutton(rightPanelp2, text='Color', fg='blue', variable=self.color_iso, offvalue=0, onvalue=1)
	chk.grid(column=0,row=3, sticky=Tk.W)
	
	Tk.Label(rightPanelp2, text="Units:",fg="black").grid(column=0, row=4, sticky=Tk.W)
	IsoPix = Radiobutton(rightPanelp2, text='pix', fg='blue', variable=self.units_iso, value='pix')
	IsoPix.grid(column=0,row=5, sticky=Tk.W)
	IsoArcsec = Radiobutton(rightPanelp2, text='arcsec', fg='blue', variable=self.units_iso, value='arcsec')
	IsoArcsec.grid(column=0,row=6, sticky=Tk.W)
	IsoKpc = Radiobutton(rightPanelp2, text='kpc', fg='blue', variable=self.units_iso, value='kpc')
	IsoKpc.grid(column=0,row=7, sticky=Tk.W)

	PlotIsophots = Tk.Frame(rightPanelp2)
	PlotIsophots.grid(column=0, row=8)
	PlotIsophotsButton = Tk.Button(PlotIsophots, text="Plot isophotes", state="normal",command= self.plot_isomap, fg="green", font=("Helvetica", 10, "bold"))
	PlotIsophotsButton.grid(column=0, row=8, sticky=Tk.N)

	PlotModelIsophots = Tk.Frame(rightPanelp2)
	PlotModelIsophots.grid(column=0, row=9)
	PlotModelIsophotsButton = Tk.Button(PlotModelIsophots, text="Overlay model", state="normal",command= self.plot_isomap_model, fg="green", font=("Helvetica", 10, "bold"))
	PlotModelIsophotsButton.grid(column=0, row=9, sticky=Tk.N)

	############################################################
	#                Creating graph stuff                      #
	############################################################
	self.IsoGraph = figure(1,figsize=(4, 4), dpi=100)
	self.canvas_iso = FigureCanvasTkAgg(self.IsoGraph, master=self.p2)
	#self.canvas.get_tk_widget().config(width=800, height=800)
	toolbar_iso = NavigationToolbar2TkAgg(self.canvas_iso, self.p2)
	toolbar_iso.update()
	self.canvas_iso.show()
	self.canvas_iso.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
	############################################################
	#                  End graph stuff                         #
	############################################################
	
	
	



	############################################################
	#                  File 3: 1D-Profiles                     #
	############################################################ 
	rightPanelp3 = Tk.Frame(self.p3)
	rightPanelp3.pack(side=Tk.RIGHT ,anchor="n")
	
	input_files,observation_info,object_info,settings = self.inp_f()
	RA,DEC = object_info[0],object_info[1]
	#new_images,ADD_INFO,RA,DEC,galaxy_name,Sampling,Sky_level = self.inp_f()

	# Set variables
	self.x_cen_prof1D = Tk.DoubleVar()
	self.y_cen_prof1D = Tk.DoubleVar()
	self.PA_cut = Tk.DoubleVar()
	self.Rmin_prof1D = Tk.DoubleVar()
	self.Rmax_prof1D = Tk.DoubleVar()
	self.step_prof1D = Tk.DoubleVar()
	self.zmax_prof1D = Tk.DoubleVar()
	self.zmin_prof1D = Tk.DoubleVar()
	self.type_1D = Tk.IntVar()
	
	self.x_cen_prof1D.set(RA)
	self.y_cen_prof1D.set(DEC)
	
	
	self.PA_cut.set(0.)
	self.Rmin_prof1D.set(0.)
	self.Rmax_prof1D.set(100.)
	self.step_prof1D.set(1.)
	self.zmin_prof1D.set(0.)
	self.zmax_prof1D.set(100.)
	
	self.type_1D.set(0)

	cen_x_Prof1D = Tk.Entry(rightPanelp3, 
		                             textvariable=self.x_cen_prof1D, 
		                             width=6, 
		                             bg="white")	

	cen_y_Prof1D = Tk.Entry(rightPanelp3, 
		                             textvariable=self.y_cen_prof1D, 
		                             width=6, 
		                             bg="white")

	PA_Prof1D = Tk.Entry(rightPanelp3, 
		                             textvariable=self.PA_cut, 
		                             width=6,
		                             state='disabled',
		                             bg="white")
	Rmin_Prof1D = Tk.Entry(rightPanelp3, 
		                             textvariable=self.Rmin_prof1D,
		                             state='disabled',
		                             width=6,
		                             bg="white")
	Rmax_Prof1D = Tk.Entry(rightPanelp3, 
		                             textvariable=self.Rmax_prof1D, 
		                             width=6, 
		                             bg="white")
	step_Prof1D = Tk.Entry(rightPanelp3, 
		                             textvariable=self.step_prof1D,
		                             state='normal',
		                             width=6,
		                             bg="white")
	zmin_Prof1D = Tk.Entry(rightPanelp3, 
		                             textvariable=self.zmin_prof1D,
		                             state='disabled',
		                             width=6,
		                             bg="white")

	zmax_Prof1D = Tk.Entry(rightPanelp3, 
		                             textvariable=self.zmax_prof1D,
		                             state='disabled',
		                             width=6, 
		                             bg="white")

	# Radiobuttons to choose the type of the profile
	Tk.Label(rightPanelp3, text="Profiles:").grid(column=0, row=1,sticky='w')
	type_of_prof1D_aver = Radiobutton(rightPanelp3, text='Azimuthally averaged', fg='blue', variable=self.type_1D, value=1, command=lambda e1=PA_Prof1D, e2=Rmin_Prof1D, e3=step_Prof1D, e4=zmin_Prof1D, e5=zmax_Prof1D, e6=cen_x_Prof1D, e7=cen_y_Prof1D, e8=Rmax_Prof1D, v=self.type_1D: self.dis_prof1D(e1,e2,e3,e4,e5,e6,e7,e8,v))
	type_of_prof1D_aver.grid(column=0,row=2, sticky=Tk.W)

	type_of_prof1D_cut = Radiobutton(rightPanelp3, text='Photometrical cut', fg='blue', variable=self.type_1D, value=2, command=lambda e1=PA_Prof1D, e2=Rmin_Prof1D, e3=step_Prof1D, e4=zmin_Prof1D, e5=zmax_Prof1D, e6=cen_x_Prof1D, e7=cen_y_Prof1D, e8=Rmax_Prof1D, v=self.type_1D: self.dis_prof1D(e1,e2,e3,e4,e5,e6,e7,e8,v))
	type_of_prof1D_cut.grid(column=0,row=3, sticky=Tk.W)

	type_of_prof1D_rad = Radiobutton(rightPanelp3, text='Radial (pixel)', fg='blue', variable=self.type_1D, value=3, command=lambda e1=PA_Prof1D, e2=Rmin_Prof1D, e3=step_Prof1D, e4=zmin_Prof1D, e5=zmax_Prof1D, e6=cen_x_Prof1D, e7=cen_y_Prof1D, e8=Rmax_Prof1D, v=self.type_1D: self.dis_prof1D(e1,e2,e3,e4,e5,e6,e7,e8,v))
	type_of_prof1D_rad.grid(column=0,row=4, sticky=Tk.W)

	type_of_prof1D_sumx = Radiobutton(rightPanelp3, text='Summed over x', fg='blue', variable=self.type_1D, value=4, command=lambda e1=PA_Prof1D, e2=Rmin_Prof1D, e3=step_Prof1D, e4=zmin_Prof1D, e5=zmax_Prof1D, e6=cen_x_Prof1D, e7=cen_y_Prof1D, e8=Rmax_Prof1D, v=self.type_1D: self.dis_prof1D(e1,e2,e3,e4,e5,e6,e7,e8,v))
	type_of_prof1D_sumx.grid(column=0,row=5, sticky=Tk.W)

	type_of_prof1D_sumy = Radiobutton(rightPanelp3, text='Summed over y', fg='blue', variable=self.type_1D, value=5, command=lambda e1=PA_Prof1D, e2=Rmin_Prof1D, e3=step_Prof1D, e4=zmin_Prof1D, e5=zmax_Prof1D, e6=cen_x_Prof1D, e7=cen_y_Prof1D, e8=Rmax_Prof1D, v=self.type_1D: self.dis_prof1D(e1,e2,e3,e4,e5,e6,e7,e8,v))
	type_of_prof1D_sumy.grid(column=0,row=6, sticky=Tk.W)
	

	Tk.Label(rightPanelp3, text="Center:x [pix]").grid(column=0, row=7,sticky='w')
	cen_x_Prof1D.grid(column=1, row=7, sticky=Tk.W)


	Tk.Label(rightPanelp3, text="Center:y [pix]").grid(column=0, row=8,sticky='w')
	cen_y_Prof1D.grid(column=1, row=8, sticky=Tk.W)
	

	Tk.Label(rightPanelp3, text="Position angle [deg]").grid(column=0, row=9,sticky='w')
	PA_Prof1D.grid(column=1, row=9, sticky=Tk.W)
	

	Tk.Label(rightPanelp3, text="Rmin [pix]").grid(column=0, row=10,sticky='w')
	Rmin_Prof1D.grid(column=1, row=10, sticky=Tk.W)

	Tk.Label(rightPanelp3, text="Rmax [pix]").grid(column=0, row=11,sticky='w')
	Rmax_Prof1D.grid(column=1, row=11, sticky=Tk.W)	

	Tk.Label(rightPanelp3, text="Step [pix]").grid(column=0, row=12,sticky='w')
	step_Prof1D.grid(column=1, row=12, sticky=Tk.W)	



	Tk.Label(rightPanelp3, text="zmin [pix]").grid(column=0, row=13,sticky='w')
	zmin_Prof1D.grid(column=1, row=13, sticky=Tk.W)

	Tk.Label(rightPanelp3, text="zmax [pix]").grid(column=0, row=14,sticky='w')
	zmax_Prof1D.grid(column=1, row=14, sticky=Tk.W)	

	ShowProfile = Tk.Frame(rightPanelp3)
	ShowProfile.grid(column=0, row=15)
	ShowProfileButton = Tk.Button(ShowProfile, text="Show profile alignment", state="normal",command= self.show_prof, fg="green", font=("Helvetica", 10, "bold"))
	ShowProfileButton.grid(column=0, row=15, sticky=Tk.N)

	Tk.Label(rightPanelp3, text="  ").grid(column=0, row=16,sticky='w')	# Blank line





	PlotProfiles = Tk.Frame(rightPanelp3)
	PlotProfiles.grid(column=0, row=19)
	PlotProfilesButton = Tk.Button(PlotProfiles, text="Plot profile", state="normal",command= self.plot_prof1D, width=15, fg="green", font=("Helvetica", 10, "bold"))
	PlotProfilesButton.grid(column=0, row=19, sticky=Tk.N)

	PlotModelProfiles = Tk.Frame(rightPanelp3)
	PlotModelProfiles.grid(column=0, row=20)
	PlotModelProfilesButton = Tk.Button(PlotModelProfiles, text="Overlay model", state="normal",command= self.plot_prof1Dmodel, width=15, fg="green", font=("Helvetica", 10, "bold"))
	PlotModelProfilesButton.grid(column=0, row=20, sticky=Tk.N)

	

	############################################################
	#                Creating graph stuff                      #
	############################################################
	self.Prof1DGraph = figure(2,figsize=(4, 4), dpi=100)
	self.canvas_Prof1D = FigureCanvasTkAgg(self.Prof1DGraph, master=self.p3)
	#self.canvas.get_tk_widget().config(width=800, height=800)
	toolbar_Prof1D = NavigationToolbar2TkAgg(self.canvas_Prof1D, self.p3)
	toolbar_Prof1D.update()
	self.canvas_Prof1D.show()
	self.canvas_Prof1D.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
	############################################################
	#                  End graph stuff                         #
	############################################################

	
	

	############################################################
	#                  File 4: 2D-Profiles                     #
	############################################################
	rightPanelp4 = Tk.Frame(self.p4)
	rightPanelp4.pack(side=Tk.RIGHT ,anchor="n")

	# Define variables
	self.OuterIso2DProf = Tk.DoubleVar()
	self.Color2DProf = Tk.IntVar()
	self.View2DProf = Tk.IntVar()
	#self.Arc2DProf = Tk.IntVar()
	#self.Kpc2DProf = Tk.IntVar()
	self.ScaleBarLength = Tk.StringVar()
	
	

	Tk.Label(rightPanelp4, text=u"Outer isophote [mag/\u25a1\"]").grid(column=0, row=0,sticky='w')
	self.OuterIso2DProf.set(25.5)
	OuterIso2DProfEntry = Tk.Entry(rightPanelp4, 
		                             textvariable=self.OuterIso2DProf, 
		                             width=6, 
		                             bg="white")
	OuterIso2DProfEntry.grid(column=0, row=1, sticky=Tk.W)


	Tk.Label(rightPanelp4, text="  ").grid(column=0, row=2,sticky='w')

	Tk.Label(rightPanelp4, text="Color").grid(column=0, row=3,sticky='w')
	self.Color2DProf.set(1)
	
	Gray2DProfButton = Radiobutton(rightPanelp4, text='gray', fg='blue', variable=self.Color2DProf, value=1)
	Gray2DProfButton.grid(column=0,row=4, sticky=Tk.W)
	
	Color2DProfButton = Radiobutton(rightPanelp4, text='color', fg='blue', variable=self.Color2DProf, value=2)
	Color2DProfButton.grid(column=0,row=5, sticky=Tk.W)
	
	Tk.Label(rightPanelp4, text="  ").grid(column=0, row=6,sticky='w')
	
	
	
	Tk.Label(rightPanelp4, text="View").grid(column=0, row=7,sticky='w')
	self.View2DProf.set(1)

	LineProfButton = Radiobutton(rightPanelp4, text='line', fg='blue', variable=self.View2DProf, value=1)
	LineProfButton.grid(column=0,row=8, sticky=Tk.W)
	
	ColProfButton = Radiobutton(rightPanelp4, text='column', fg='blue', variable=self.View2DProf, value=2)
	ColProfButton.grid(column=0,row=9, sticky=Tk.W)


	Tk.Label(rightPanelp4, text="  ").grid(column=0, row=10,sticky='w')

	#self.Arc2DProf.set(0)
	#self.Kpc2DProf.set(0)
	
	Tk.Label(rightPanelp4, text="Scale bar [\"]").grid(column=0, row=11,sticky='w')
	
	self.ScaleBarLength.set('')
	ScaleBarLengthEntry = Tk.Entry(rightPanelp4, 
		                             textvariable=self.ScaleBarLength, 
		                             width=6, 
		                             bg="white")
	ScaleBarLengthEntry.grid(column=1, row=11, sticky=Tk.W)
	
	#Arc2DProfButton = Checkbutton(rightPanelp4, text='arcsec', fg='blue', variable=self.Arc2DProf, offvalue=0, onvalue=1)
	#Arc2DProfButton.grid(column=0,row=12, sticky=Tk.W)

	#Kpc2DProfButton = Checkbutton(rightPanelp4, text='kpc', fg='blue', variable=self.Kpc2DProf, offvalue=0, onvalue=1)
	#Kpc2DProfButton.grid(column=0,row=13, sticky=Tk.W)



	PlotModelProfile = Tk.Frame(rightPanelp4)
	PlotModelProfile.grid(column=0, row=20)
	PlotModelProfileButton = Tk.Button(PlotModelProfile, text="Plot model", state="normal",command= self.plot_prof2Dmodel, width=15, fg="green", font=("Helvetica", 10, "bold"))
	PlotModelProfileButton.grid(column=0, row=20, sticky=Tk.N)


	############################################################
	#                Creating graph stuff                      #
	############################################################
	self.Prof2DGraph = figure(33,figsize=(4, 4), dpi=100)
	self.canvas_Prof2D = FigureCanvasTkAgg(self.Prof2DGraph, master=self.p4)
	#self.canvas.get_tk_widget().config(width=800, height=800)
	toolbar_Prof2D = NavigationToolbar2TkAgg(self.canvas_Prof2D, self.p4)
	toolbar_Prof2D.update()
	self.canvas_Prof2D.show()
	self.canvas_Prof2D.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
	############################################################
	#                  End graph stuff                         #
	############################################################




	############################################################
	#                  File 5: Ellipse fitting                 #
	############################################################ 
	rightPanelp5 = Tk.Frame(self.p5)
	rightPanelp5.pack(side=Tk.RIGHT ,anchor="n")
	
	self.min_rad_ell = Tk.DoubleVar()
	self.max_rad_ell = Tk.DoubleVar()
	self.step_rad_ell = Tk.DoubleVar()
	self.units_ell = Tk.IntVar()
	self.Ell_files = Tk.StringVar()
	self.Labels = Tk.StringVar()


	self.min_rad_ell.set(0.)
	self.max_rad_ell.set(100.)
	self.step_rad_ell.set(1.0)
	self.units_ell.set(1)
	
	Tk.Label(rightPanelp5, text="Inner radius",fg="blue").grid(column=0, row=0, sticky=Tk.W)
	InnerRadEllEntry = Tk.Entry(rightPanelp5, textvariable=self.min_rad_ell, width=5, state="normal", bg="white")
	InnerRadEllEntry.grid(column=1, row=0, sticky=Tk.W)


	Tk.Label(rightPanelp5, text="Outer radius",fg="blue").grid(column=0, row=1, sticky=Tk.W)
	OuterIsoEntry = Tk.Entry(rightPanelp5, textvariable=self.max_rad_ell, width=5, state="normal", bg="white")
	OuterIsoEntry.grid(column=1, row=1, sticky=Tk.W)


	Tk.Label(rightPanelp5, text="Step",fg="blue").grid(column=0, row=2, sticky=Tk.W)
	StepIsoEntry = Tk.Entry(rightPanelp5, textvariable=self.step_rad_ell, width=5, state="normal", bg="white")
	StepIsoEntry.grid(column=1, row=2, sticky=Tk.W)


	#radio_ell_units_pix = Radiobutton(rightPanelp5, text='pixel', fg='blue', variable=self.units_ell, value=1)
	#radio_ell_units_pix.grid(column=0,row=3, sticky=Tk.W)

	#radio_ell_units_arcsec = Radiobutton(rightPanelp5, text='arcsec', fg='blue', variable=self.units_ell, value=2)
	#radio_ell_units_arcsec.grid(column=1,row=3, sticky=Tk.W)

	LaunchEllipse = Tk.Frame(rightPanelp5)
	LaunchEllipse.grid(column=0, row=4)
	LaunchEllipseButton = Tk.Button(LaunchEllipse, text="Galaxy fitting", state="normal",command= self.launch_ellipse, width=10, fg="green", font=("Helvetica", 10, "bold"))
	LaunchEllipseButton.grid(column=0, row=4, sticky=Tk.N)


	PlotEllipse = Tk.Frame(rightPanelp5)
	PlotEllipse.grid(column=0, row=5)
	PlotEllipseButton = Tk.Button(PlotEllipse, text="Plot", state="normal",command= self.plot_ellipse, width=10, fg="green", font=("Helvetica", 10, "bold"))
	PlotEllipseButton.grid(column=0, row=5, sticky=Tk.N)

	LaunchEllipseModel = Tk.Frame(rightPanelp5)
	LaunchEllipseModel.grid(column=0, row=6)
	LaunchEllipseModelButton = Tk.Button(LaunchEllipseModel, text="Model fitting", state="normal",command= self.launch_ellipse_model, width=10, fg="green", font=("Helvetica", 10, "bold"))
	LaunchEllipseModelButton.grid(column=0, row=6, sticky=Tk.N)

	PlotModelEllipse = Tk.Frame(rightPanelp5)
	PlotModelEllipse.grid(column=0, row=7)
	PlotModelEllipseButton = Tk.Button(PlotModelEllipse , text="Overlay model", state="normal",command= self.plot_ellipse_model, width=10, fg="green", font=("Helvetica", 10, "bold"))
	PlotModelEllipseButton.grid(column=0, row=7, sticky=Tk.N)


	############################################################
	#                Creating graph stuff                      #
	############################################################
	self.EllipseGraph = figure(44,figsize=(4, 4), dpi=100)
	self.canvas_ell = FigureCanvasTkAgg(self.EllipseGraph, master=self.p5)
	#self.canvas_ell.get_tk_widget().config(width=800, height=800)
	toolbar_ell = NavigationToolbar2TkAgg(self.canvas_ell, self.p5)
	toolbar_ell.update()
	self.canvas_ell.show()
	self.canvas_ell.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
	############################################################




	

	############################################################
	#                  File 6: Decomposition                   #
	############################################################ 
	DecompPanel = Tk.Frame(self.p6)
	DecompPanel.pack(side=Tk.LEFT ,anchor="n")
	self.DecompPanel = DecompPanel
	
	self.Model_file = Tk.StringVar()
	self.Model_file.set('composed_model.fits')

	Tk.Label(DecompPanel, text="               ", font=("Helvetica", 10, "bold"), fg='red', width=25).grid(column=0, row=0,sticky='e')
	Tk.Label(DecompPanel, text="               ", font=("Helvetica", 10, "bold"), fg='red', width=25).grid(column=0, row=1,sticky='e')

	Tk.Label(DecompPanel, text="STEP I:", font=("Helvetica", 10, "bold"), fg='red',width=25).grid(column=0, row=2,sticky=Tk.W)
	Tk.Label(DecompPanel, text="FIND INITIAL GUESS   ", font=("Helvetica", 10, "bold"), width=25).grid(column=1, row=2,sticky='e')


	# Create a new object
	CreateNewObject = Tk.Frame(DecompPanel)
	CreateNewObject.grid(column=1, row=3)
	CreateNewObjectButton = Tk.Button(CreateNewObject, text="Create initial model", state="normal",command=lambda: self.create_window([]), width=25, bg="green", font=("Helvetica", 10, "bold"))
	CreateNewObjectButton.grid(column=1, row=3)


	# Load object
	LoadObject = Tk.Frame(DecompPanel)
	LoadObject.grid(column=2, row=3)
	LoadObjectButton = Tk.Button(LoadObject, text="Load model", state="normal",command=self.load_window, width=25, fg="green", font=("Helvetica", 10, "bold"))
	LoadObjectButton.grid(column=2, row=3)

	Tk.Label(DecompPanel, text="               ", font=("Helvetica", 10, "bold"), fg='red', width=25).grid(column=0, row=4,sticky='e')
	Tk.Label(DecompPanel, text="               ", font=("Helvetica", 10, "bold"), fg='red', width=25).grid(column=0, row=5,sticky='e')

	Tk.Label(DecompPanel, text="STEP II:", font=("Helvetica", 10, "bold"), fg='red',width=25).grid(column=0, row=6,sticky=Tk.W)
	Tk.Label(DecompPanel, text="CHECK INPUT FILE(S)   ", font=("Helvetica", 10, "bold"), width=25).grid(column=1, row=6,sticky='e')


	# Edit input files
	EditFiles = Tk.Frame(DecompPanel)
	EditFiles.grid(column=1, row=7)
	EditFilesButton = Tk.Button(EditFiles, text="Edit input files", state="normal",command= self.edit_files, width=25, fg="green", font=("Helvetica", 10, "bold"))
	EditFilesButton.grid(column=1, row=7)


	Tk.Label(DecompPanel, text="               ", font=("Helvetica", 10, "bold"), fg='red', width=25).grid(column=0, row=8,sticky='e')
	Tk.Label(DecompPanel, text="               ", font=("Helvetica", 10, "bold"), fg='red', width=25).grid(column=0, row=9,sticky='e')


	Tk.Label(DecompPanel, text="STEP III:", font=("Helvetica", 10, "bold"), fg='red',width=25).grid(column=0, row=10,sticky=Tk.W)
	Tk.Label(DecompPanel, text="APPROVE MASK FILE   ", font=("Helvetica", 10, "bold"), width=25).grid(column=1, row=10,sticky='e')

	# Additional mask
	self.FinalMask = Tk.StringVar()
	self.FinalMask.set(self.mask_image1.get())
	
	AddMask = Tk.Frame(DecompPanel)
	AddMask.grid(column=1, row=11)
	AddMaskButton = Tk.Button(AddMask, text="Create additional mask", state="normal",command= self.add_mask, width=25, fg="green", font=("Helvetica", 10, "bold"))
	AddMaskButton.grid(column=1, row=11)

	# Choose the mask file
	self.ChooseMaskFileEntry = Tk.Entry(DecompPanel, 
		                             textvariable=self.FinalMask, 
		                             width=25, 
		                             bg="white")

	ChooseMaskFile = Tk.Frame(DecompPanel)
	ChooseMaskFile.grid(column=2, row=11)
	ChooseMaskFileButton = Tk.Button(ChooseMaskFile, text="Choose final mask", state="normal",command=self.choose_mask, width=25, fg="green", font=("Helvetica", 10, "bold"))
	ChooseMaskFileButton.grid(column=2, row=11)
	
	self.ChooseMaskFileEntry.grid(column=3, row=11, sticky=Tk.W)





	Tk.Label(DecompPanel, text="               ", font=("Helvetica", 10, "bold"), fg='red', width=25).grid(column=0, row=12,sticky='e')
	Tk.Label(DecompPanel, text="               ", font=("Helvetica", 10, "bold"), fg='red', width=25).grid(column=0, row=13,sticky='e')


	Tk.Label(DecompPanel, text="STEP IV:", font=("Helvetica", 10, "bold"), fg='red',width=25).grid(column=0, row=14,sticky=Tk.W)
	Tk.Label(DecompPanel, text="LAUNCH DECOMPOSITION ", font=("Helvetica", 10, "bold"), fg='red', width=25).grid(column=1, row=14,sticky='e')


	# Print the line to launch code:
	self.LaunchLine = Tk.StringVar()

	#######################################
	#new_images,ADD_INFO,RA,DEC,galaxy_name,Sampling,Sky_level = self.inp_f()
	input_files,observation_info,object_info,settings = self.inp_f()
	input_files[3] = self.FinalMask.get()
      
	if os.path.exists('model.fits'):
	   os.remove('model.fits')

	if self.ChooseCode.get()=='GALFIT':
	  self.LaunchLine.set(galfitPath+'galfit ' + 'galfit.inp')
	elif self.ChooseCode.get()=='IMFIT':
	  self.LaunchLine.set(imfit_input.imfit_exec_line(input_files,observation_info,object_info,'imfit.inp'))  
	#######################################


	Tk.Label(DecompPanel, text="Line to launch the code:", font=("Helvetica", 10, "bold")).grid(column=2, row=15,sticky='e')
	#self.LaunchLine.set("")
	LaunchLineEntry = Tk.Entry(DecompPanel, 
		                             textvariable=self.LaunchLine, 
		                             width=100, 
		                             bg="white")
	LaunchLineEntry.grid(column=3, row=15, sticky=Tk.W)

	Tk.Label(DecompPanel, text="(change it if you wish)  ").grid(column=2, row=16,sticky='e')


	# Auto decomposition
	AutoDec = Tk.Frame(DecompPanel)
	AutoDec.grid(column=1, row=17)
	AutoDecButton = Tk.Button(AutoDec, text="Auto decomposition", state="normal",command= self.auto_dec, bg="yellow", width=25, font=("Helvetica", 10, "bold"))
	AutoDecButton.grid(column=1, row=17)


	# Launch decomposition
	LaunchDec = Tk.Frame(DecompPanel)
	LaunchDec.grid(column=2, row=17)
	LaunchDecButton = Tk.Button(LaunchDec, text="Final decomposition", state="normal",command= self.launch_dec, bg="green", width=25, font=("Helvetica", 10, "bold"))
	LaunchDecButton.grid(column=2, row=17)


	Tk.Label(DecompPanel, text="               ", font=("Helvetica", 10, "bold"), fg='red', width=25).grid(column=0, row=18,sticky='e')
	Tk.Label(DecompPanel, text="               ", font=("Helvetica", 10, "bold"), fg='red', width=25).grid(column=0, row=19,sticky='e')


	Tk.Label(DecompPanel, text="STEP V:", font=("Helvetica", 10, "bold"), fg='red',width=25).grid(column=0, row=20,sticky=Tk.W)
	Tk.Label(DecompPanel, text="RESULTS    ", font=("Helvetica", 10, "bold"), width=25).grid(column=1, row=20,sticky='e')



	# Choose the model to be plotted:
	ChooseModel = Tk.Frame(DecompPanel)
	ChooseModel.grid(column=1, row=21)
	ChooseModelButton = Tk.Button(ChooseModel, text="Choose model to be plotted", state="normal",command= self.choose_plot, width=25, font=("Helvetica", 10, "bold"))
	ChooseModelButton.grid(column=1, row=21)
	

	# Create output file
	CreateOutput = Tk.Frame(DecompPanel)
	CreateOutput.grid(column=1, row=22)
	CreateOutputButton = Tk.Button(CreateOutput, text="Create output", state="normal",command= self.create_output, bg="green", width=25, font=("Helvetica", 10, "bold"))
	CreateOutputButton.grid(column=1, row=22)


	# Save the model to:
	self.SaveModel = Tk.StringVar()
	Tk.Label(DecompPanel, text="Save output to (folder name):", font=("Helvetica", 10, "bold")).grid(column=2, row=22,sticky='e')
	self.SaveModel.set("")
	SaveModelEntry = Tk.Entry(DecompPanel, 
		                             textvariable=self.SaveModel, 
		                             width=100, 
		                             bg="white")
	SaveModelEntry.grid(column=3, row=22, sticky=Tk.E)


















    ####********************** FUNCTIONS **********************###

    def convert_coords(self):
	  from astropy import wcs
	  
	  # Convert RA,DEC to x,y  
	  hdulist = pyfits.open(self.file_field1.get())
	  w = wcs.WCS(hdulist[0].header)
	  [[XX,YY]] = w.wcs_world2pix([[self.RA1.get(),self.Dec1.get()]],1)
	  return XX,YY
    
    #-------------------------------------------------------------
    # ELLIPSE FITTING FUNCTIONS
    def launch_ellipse(self):
	    print 'HERE'
	    iraf_ell_plots.crea_ell(self.file_field1.get(),float(self.RA1.get()),self.Dec1.get(),self.step_rad_ell.get(),self.min_rad_ell.get(),self.max_rad_ell.get(),'ellipse.txt')
	    os.chmod(r"ell.cl",0777)
	    subprocess.call("cl < ell.cl -o", shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)

    def launch_ellipse_model(self):
	    hdu = pyfits.open(self.Model_file.get())
	    number_of_layers = len(hdu)
	    labels = []
	    ell_files = []
	    for layer in range(1,number_of_layers):
	      prihdr = hdu[layer].header
	      Label = prihdr['NAME_OF_LAYER']

	      ell_file = 'ellipse_%s.' % (Label)
	      iraf_ell_plots.crea_ell(self.file_field1.get(),float(self.RA1.get()),self.Dec1.get(),self.step_rad_ell.get(),self.min_rad_ell.get(),self.max_rad_ell.get(),ell_file)
	      os.chmod(r"ell.cl",0777)
	      subprocess.call("cl < ell.cl -o", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
	      ell_files.append(ell_file)
	      labels.append(Label)
	    self.Ell_files.set(ell_files)
	    self.Labels.set(labels)
    
    def plot_ellipse(self):
	import graph_creator
	input_files,observation_info,object_info,settings = self.inp_f()
	[input_image,sigma_image,psf_image,mask_image,sampling] = input_files
	[nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,Filter] = observation_info
	[RA,DEC,galaxy_name,Dist,Ext,Kcorr] = object_info
	[ChooseCode,ChooseGeom,ChooseSB,ChooseLum] = settings
	Scale = Dist * 1000. / 206265.   
	graph_creator.Plot_IRAF(['ellipse.txt'],['galaxy'], self.EllipseGraph, self.canvas_ell, m0,scale,Scale,geom_units=ChooseGeom,SB_units=ChooseSB)

    
    def plot_ellipse_model(self):
	import graph_creator
	input_files,observation_info,object_info,settings = self.inp_f()
	[input_image,sigma_image,psf_image,mask_image,sampling] = input_files
	[nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,Filter] = observation_info
	[RA,DEC,galaxy_name,Dist,Ext,Kcorr] = object_info
	[ChooseCode,ChooseGeom,ChooseSB,ChooseLum] = settings
	Scale = Dist * 1000. / 206265.   
	graph_creator.Plot_IRAF(self.Ell_files.get(),self.Labels.get(), self.EllipseGraph, self.canvas_ell, m0,scale,Scale,geom_units=ChooseGeom,SB_units=ChooseSB)
    #-------------------------------------------------------------
    

    #-------------------------------------------------------------
    # ISOPHOTES FUNCTIONS
    def plot_isomap(self):
      import graph_creator
      input_files,observation_info,object_info,settings = self.inp_f()
      [input_image,sigma_image,psf_image,mask_image,sampling] = input_files
      [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,Filter] = observation_info
      [RA,DEC,galaxy_name,Dist,Ext,Kcorr] = object_info
      Scale = Dist * 1000. / 206265.

      if self.color_iso.get()==1:
	color = None
      else:
	color = 'white'
      graph_creator.Plot_isomap(self.file_field1.get(), self.mask_image1.get(), self.IsoGraph, self.canvas_iso, float(self.m01.get()), float(self.pix2sec1.get()),RA,DEC, self.min_mag_iso.get(), self.max_mag_iso.get(), self.step_iso.get(), color, self.units_iso.get(),Scale) 

    def plot_isomap_model(self):
      import graph_creator
      input_files,observation_info,object_info,settings = self.inp_f()
      [input_image,sigma_image,psf_image,mask_image,sampling] = input_files
      [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,Filter] = observation_info
      [RA,DEC,galaxy_name,Dist,Ext,Kcorr] = object_info
      Scale = Dist * 1000. / 206265.
      
      if self.color_iso.get()==1:
	color = None
      else:
	color = 'white'
      graph_creator.Plot_isomap(self.Model_file.get(), self.mask_image1.get(), self.IsoGraph, self.canvas_iso, float(self.m01.get()), float(self.pix2sec1.get()),RA,DEC, self.min_mag_iso.get(), self.max_mag_iso.get(), self.step_iso.get(), color, self.units_iso.get(),Scale)     
      
    #-------------------------------------------------------------


    #-------------------------------------------------------------
    # 1D-PROFILE FUNCTIONS
    def plot_prof1D(self):
      input_files,observation_info,object_info,settings = self.inp_f()
      [input_image,sigma_image,psf_image,mask_image,sampling] = input_files
      [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,Filter] = observation_info
      [RA,DEC,galaxy_name,Dist,Ext,Kcorr] = object_info
      [ChooseCode,ChooseGeom,ChooseSB,ChooseLum] = settings
      Scale = Dist * 1000. / 206265.   
      import graph_creator
      if self.type_1D.get()==1:
	type_prof = 'azim'
      elif self.type_1D.get()==2:
	type_prof = 'cut'
      elif self.type_1D.get()==3:
	type_prof = 'rad'
      elif self.type_1D.get()==4:
	type_prof = 'vert'
      elif self.type_1D.get()==5:
	type_prof = 'summed'
      else:
	tkMessageBox.showerror('Warning','Please choose the profile!')
	return 1
      graph_creator.Plot_Prof1D(self.file_field1.get(), self.mask_image1.get(), self.Prof1DGraph, self.canvas_Prof1D,
				float(self.m01.get()), float(self.pix2sec1.get()),type_prof,
				xc=self.x_cen_prof1D.get(), yc=self.y_cen_prof1D.get(),
				PA=self.PA_cut.get(), Rmin=self.Rmin_prof1D.get(),
				Rmax=self.Rmax_prof1D.get(), step=self.step_prof1D.get(),
				zmin=self.zmin_prof1D.get(),zmax=self.zmax_prof1D.get(),
				geom_units=ChooseGeom,SB_units=ChooseSB,Scale=Scale)

    def plot_prof1Dmodel(self):
      input_files,observation_info,object_info,settings = self.inp_f()
      [input_image,sigma_image,psf_image,mask_image,sampling] = input_files
      [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,Filter] = observation_info
      [RA,DEC,galaxy_name,Dist,Ext,Kcorr] = object_info
      [ChooseCode,ChooseGeom,ChooseSB,ChooseLum] = settings
      Scale = Dist * 1000. / 206265.      
      
      import graph_creator
      if self.type_1D.get()==1:
	type_prof = 'azim'
      elif self.type_1D.get()==2:
	type_prof = 'cut'
      elif self.type_1D.get()==3:
	type_prof = 'rad'
      elif self.type_1D.get()==4:
	type_prof = 'vert'
      elif self.type_1D.get()==5:
	type_prof = 'summed'
      else:
	tkMessageBox.showerror('Warning','Please choose the profile!')
	return 1
      graph_creator.Plot_Prof1D(self.Model_file.get(), self.mask_image1.get(), self.Prof1DGraph, self.canvas_Prof1D,
				float(self.m01.get()), float(self.pix2sec1.get()),type_prof,
				xc=self.x_cen_prof1D.get(), yc=self.y_cen_prof1D.get(),
				PA=self.PA_cut.get(), Rmin=self.Rmin_prof1D.get(),
				Rmax=self.Rmax_prof1D.get(), step=self.step_prof1D.get(),
			        zmin=self.zmin_prof1D.get(), zmax=self.zmax_prof1D.get(),
			        geom_units=ChooseGeom,SB_units=ChooseSB,Scale=Scale)

    def dis_prof1D(self,e1,e2,e3,e4,e5,e6,e7,e8,var):
      #e1=PA_Prof1D, e2=Rmin_Prof1D, e3=step_Prof1D, e4=zmin_Prof1D, e5=zmax_Prof1D,e8=Rmax_Prof1D
      #new_images,ADD_INFO,RA,DEC,galaxy_name,Sampling,Sky_level = self.inp_f()
      input_files,observation_info,object_info,settings = self.inp_f()
      RA,DEC = object_info[0],object_info[1]
      
      e6.configure(textvariable=self.x_cen_prof1D.set(RA))
      e7.configure(textvariable=self.y_cen_prof1D.set(DEC))
      e8.configure(textvariable=self.Rmax_prof1D.set(max([RA,self.nx1-RA,DEC,self.ny1-DEC])))
      
      if var.get() == 1:
	e1.configure(state='disabled')
	e2.configure(state='disabled')
	e3.configure(state='normal')
	e4.configure(state='disabled')
	e5.configure(state='disabled')
	e1.update();e2.update();e3.update();e4.update();e5.update()
      elif var.get() == 2:
	e1.configure(state='normal')
	e2.configure(state='disabled')
	e3.configure(state='disabled')
	e4.configure(state='disabled')
	e5.configure(state='disabled')
	e1.update();e2.update();e3.update();e4.update();e5.update()     
      elif var.get() == 3:
	e1.configure(state='disabled')
	e2.configure(state='disabled')
	e3.configure(state='disabled')
	e4.configure(state='disabled')
	e5.configure(state='disabled')
	e1.update();e2.update();e3.update();e4.update();e5.update()   
      elif var.get() == 4:
	e1.configure(state='disabled')
	e2.configure(state='normal')
	e3.configure(state='disabled')
	e4.configure(state='disabled')
	e5.configure(state='normal')
	e1.update();e2.update();e3.update();e4.update();e5.update() 
      elif var.get() == 5:
	e1.configure(state='disabled')
	e2.configure(state='disabled')
	e3.configure(state='disabled')
	e4.configure(state='normal')
	e5.configure(state='normal')
	e1.update();e2.update();e3.update();e4.update();e5.update() 

    def show_prof(self):
	if self.type_1D.get()==1:
	  type_prof = 'azim'
	elif self.type_1D.get()==2:
	  type_prof = 'cut'
	elif self.type_1D.get()==3:
	  type_prof = 'rad'
	elif self.type_1D.get()==4:
	  type_prof = 'vert'
	elif self.type_1D.get()==5:
	  type_prof = 'summed'
	else:
	  tkMessageBox.showerror('Warning','Please choose the profile!')
	  return 1
	  
	window = Tk.Toplevel(self.master)
	window.title('Alignment')

	Prof1DAlign = figure(22,figsize=(4, 4), dpi=100)
	canvas_Prof1DAlign = FigureCanvasTkAgg(Prof1DAlign, master=window)
	#self.canvas.get_tk_widget().config(width=800, height=800)
	toolbar_Prof1DAlign = NavigationToolbar2TkAgg(canvas_Prof1DAlign, window)
	toolbar_Prof1DAlign.update()
	canvas_Prof1DAlign.show()
	canvas_Prof1DAlign.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

	import graph_creator
	graph_creator.Plot_Alignment(self.file_field1.get(), Prof1DAlign, canvas_Prof1DAlign,
			      type_prof, self.x_cen_prof1D.get(), self.y_cen_prof1D.get(),
			      self.Rmax_prof1D.get(), self.PA_cut.get(), self.Rmin_prof1D.get(),
			      self.zmin_prof1D.get(), self.zmax_prof1D.get())
    #-------------------------------------------------------------


    #-------------------------------------------------------------
    # OUTPUT FUNCTIONS
    def plot_prof2Dmodel(self):
      if self.Color2DProf.get()==1:
	color='gray_r'
      elif self.Color2DProf.get()==2:
	color='nipy_spectral_r'
      if self.ScaleBarLength.get()=='':
	scale_bar2D = None
      else:
	scale_bar2D = self.ScaleBarLength.get()
      import graph_creator
      graph_creator.Plot_Prof2D(self.Model_file.get(),self.Prof2DGraph, self.canvas_Prof2D,self.OuterIso2DProf.get(),
				float(self.pix2sec1.get()),float(self.m01.get()),self.mask_image1.get(),
				self.View2DProf.get(),color,scale_bar2D)
    #-------------------------------------------------------------



    

    #-------------------------------------------------------------
    # OUTPUT FUNCTIONS
    def create_output(self):
      window_out = Tk.Toplevel(self.master)
      window_out.title("Settings")
      window_out.geometry(("600x300"))

      input_files,observation_info,object_info,settings = self.inp_f()
      [input_image,sigma_image,psf_image,mask_image,sampling] = input_files
      [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,Filter] = observation_info
      [RA,DEC,galaxy_name,Dist,Ext,Kcorr] = object_info
      [ChooseCode,ChooseGeom,ChooseSB,ChooseLum] = settings

      self.ChooseGeomOut = Tk.StringVar()
      self.ChooseLumOut = Tk.StringVar()
      self.ChooseSBOut = Tk.StringVar()
      self.ChooseGeomOut.set(ChooseGeom)
      self.ChooseLumOut.set(ChooseLum)
      self.ChooseSBOut.set(ChooseSB)

      Tk.Label(window_out, text="Choose units for",anchor='w', font=("Helvetica", 10, "bold"), width=35).grid(column=0, row=5,sticky='w')
      
      Tk.Label(window_out, text="     - Geometrical parameters:",anchor='w', font=("Helvetica", 10, "bold"), width=35).grid(column=0, row=6,sticky='w')

      PixButton = Radiobutton(window_out, text='pix', fg='blue', variable=self.ChooseGeomOut, value='pix')
      PixButton.grid(column=1,row=6, sticky=Tk.W)        

      ArcButton = Radiobutton(window_out, text='arcsec', fg='blue', variable=self.ChooseGeomOut, value='arcsec')
      ArcButton.grid(column=2,row=6, sticky=Tk.W)   

      kpcButton = Radiobutton(window_out, text='kpc', fg='blue', variable=self.ChooseGeomOut, value='kpc')
      kpcButton.grid(column=3,row=6, sticky=Tk.W)    
      
      
      
      
      Tk.Label(window_out, text="     - Surface brightness:",anchor='w', font=("Helvetica", 10, "bold"), width=35).grid(column=0, row=7,sticky='w')

      SBADUButton = Radiobutton(window_out, text=u'ADU/pix\u00b2', fg='blue', variable=self.ChooseSBOut, value='ADU/pix2')
      SBADUButton.grid(column=1,row=7, sticky=Tk.W)        

      SBMagButton = Radiobutton(window_out, text=u'mag/\u25a1\"', fg='blue', variable=self.ChooseSBOut, value='mag/arcsec2')
      SBMagButton.grid(column=2,row=7, sticky=Tk.W)   

      #LumButton = Radiobutton(window_out, text=u'L\u2609', fg='blue', variable=self.ChooseSBOut, value='Lsun')
      #LumButton.grid(column=3,row=7, sticky=Tk.W)   




      Tk.Label(window_out, text="     - Luminosity:",anchor='w', font=("Helvetica", 10, "bold"), width=35).grid(column=0, row=8,sticky='w')

      ADUButton = Radiobutton(window_out, text='ADU', fg='blue', variable=self.ChooseLumOut, value='ADU')
      ADUButton.grid(column=1,row=8, sticky=Tk.W)        

      MagButton = Radiobutton(window_out, text='mag', fg='blue', variable=self.ChooseLumOut, value='mag')
      MagButton.grid(column=2,row=8, sticky=Tk.W)   

      #LumButton = Radiobutton(window_out, text=u'L\u2609', fg='blue', variable=self.ChooseLumOut, value='Lsun')
      #LumButton.grid(column=3,row=8, sticky=Tk.W)   
      

      DirFrame = Tk.Frame(window_out)
      DirFrame.grid(column=1, row=9)
      DirFrameButton = Tk.Button(window_out, text="Choose directory", state="normal",command=self.choose_dir, bg="green",font=("Helvetica", 10, "bold"), width=15)
      DirFrameButton.grid(column=1, row=9)
      
      
      TexFrame = Tk.Frame(window_out)
      TexFrame.grid(column=1, row=10)
      TexFrameButton = Tk.Button(window_out, text="Create Tex file", state="normal",command=self.create_tex, bg="green",font=("Helvetica", 10, "bold"), width=15)
      TexFrameButton.grid(column=1, row=10)


      EditTexFrame = Tk.Frame(window_out)
      EditTexFrame.grid(column=1, row=11)
      EditTexFrameButton = Tk.Button(window_out, text="Edit Tex file", state="normal",command=self.edit_tex, bg="green",font=("Helvetica", 10, "bold"), width=15)
      EditTexFrameButton.grid(column=1, row=11)      


      CompileTexFrame = Tk.Frame(window_out)
      CompileTexFrame.grid(column=1, row=12)
      CompileTexFrameButton = Tk.Button(window_out, text="Compile", state="normal",command=self.compile_tex, bg="green",font=("Helvetica", 10, "bold"), width=15)
      CompileTexFrameButton.grid(column=1, row=12)


      ShowTexRes = Tk.Frame(window_out)
      ShowTexRes.grid(column=1, row=13)
      ShowTexResButton = Tk.Button(window_out, text="Show ps", state="normal",command=self.show_tex_res, bg="green",font=("Helvetica", 10, "bold"), width=15)
      ShowTexResButton.grid(column=1, row=13)
      
      
    def choose_dir(self):
      Res_Directory = tkFileDialog.askdirectory()
      self.SaveModel.set(Res_Directory)
    
    
    def create_tex(self):
      input_files,observation_info,object_info,settings = self.inp_f()
      [input_image,sigma_image,psf_image,mask_image,sampling] = input_files
      [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,Filter] = observation_info
      [RA,DEC,galaxy_name,Dist,Ext,Kcorr] = object_info
      [ChooseCode,ChooseGeom,ChooseSB,ChooseLum] = settings
      Scale = Dist * 1000. / 206265.
      
      # Define results file and composed image
      composed_image = self.SaveModel.get()+'/'+'composed_model.fits'
      if ChooseCode=='IMFIT':
	results_file = self.SaveModel.get()+'/'+'imfit.01'
      elif ChooseCode=='GALFIT':
	results_file = self.SaveModel.get()+'/'+'galfit.01'
      
      # Find luminosities
      '''
      if ChooseCode=='IMFIT':
	import crea_ski
	flux,mag,frac,Mag,Flux,magg,MAG = crea_ski.find_lum(results_file,input_image,np.max([nx,ny]),m0,Dist,Ext,z=0.,filter_name='none',colour_name='none',colour_value='none')
	luminosities = []
	for i in range(len(flux)):
	  luminosities.append(flux[i]/Flux)
	luminosities.append(Flux)
	
      elif ChooseCode=='GALFIT':
      '''
      hdu = pyfits.open(composed_image)
      number_of_layers = len(hdu)
      fluxes = []
      for k in range(number_of_layers):
	  if k==1:
	    Flux = np.sum(hdu[k].data)
	  if k>=4:
	    fluxes.append(np.sum(hdu[k].data))
      luminosities = []
      for k in range(len(fluxes)):
	luminosities.append(fluxes[k]/Flux)
      luminosities.append(Flux)
      
      # Create pictures
      pictures = []
      pictures.append(plot_2d_profile.main(composed_image,25.5,scale,m0,mask_file=mask_image))
      print pictures
      tex_creator.main(galaxy_name,results_file,ChooseCode,luminosities,pictures,m0,scale,Scale,Filter,Ext,Kcorr,self.ChooseGeomOut.get(),self.ChooseLumOut.get(),self.ChooseSBOut.get())
      
    def edit_tex(self):
	from text_editor import App
	app = App()
	app.root.mainloop()      
      
    def compile_tex(self):
	subprocess.call('latex results.tex', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
	subprocess.call('latex results.tex', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
	subprocess.call('dvips results.dvi -o', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
	subprocess.call('rm -f *.aux *.dvi *.log', shell=True)

    def show_tex_res(self):
	subprocess.call('%s results.ps' % (deca_setup.image_viewer), shell=True)

    #------------------------------------------------------------- 
 
 
    #-------------------------------------------------------------
    # MASK FUNCTIONS
    def add_mask(self):
	# Check reg. files before the DS9 launch
	arr_before = glob.glob('*.reg')
	
	ds9Proc = subprocess.Popen(["ds9", self.file_field1.get(), "-scale", "log"])
	raw_input('Press Enter to continue...')

	arr_after = glob.glob('*.reg')

	
	new_reg_file = list( set(arr_after) - set(arr_before) )
	
	import reprocess_mask
	new_mask_file = reprocess_mask.add_mask_to_mask(self.mask_image1.get(),new_reg_file[0])
	self.FinalMask.set(new_mask_file)
	#return new_mask_file

    def choose_mask(self):
	fileName = tkFileDialog.askopenfilename(parent=self.DecompPanel,
				  filetypes=[("Fits files", "*.fit*"),("All files", "*.*")],
				  title="Open data file")
	if fileName!='':
	  self.ChooseMaskFileEntry.configure(textvariable=self.FinalMask.set(fileName.split('/')[-1]))
	else:
	  self.ChooseMaskFileEntry.configure(textvariable=self.FinalMask.set(self.mask_image1.get()))
    #-------------------------------------------------------------

    def choose_plot(self):
	fileName = tkFileDialog.askopenfilename(parent=self.DecompPanel,
				  filetypes=[("Fits files", "*.fit*"),("All files", "*.*")],
				  title="Open data file")
	
	if os.path.exists(fileName):
	  try:
	    shutil.copy(fileName,fileName.split('/')[-1])
	    self.Model_file.set(fileName.split('/')[-1])
	  except:
	    self.Model_file.set(fileName)

	    
	    
 
 
    #-------------------------------------------------------------
    # EDIT FILE FUNCTION
    def edit_files(self):
	from text_editor import App
	app = App()
	app.root.mainloop()

    #-------------------------------------------------------------

    #-------------------------------------------------------------
    # AUTO DECOMPOSITION FUNCTION
    def auto_dec(self):
      self.ChooseModel = IntVar()
      self.ChooseModel.set(1)

      
      # Create window where you can choose the model
      window_auto = Tk.Toplevel(self.master)

      Tk.Label(window_auto, text="Choose the model:", font=("Helvetica", 10, "bold"), width=25).grid(column=0, row=1,sticky='e')
      
      SersicButton = Radiobutton(window_auto, text='Sersic', fg='blue', variable=self.ChooseModel, value=1)
      SersicButton.grid(column=0,row=2, sticky=Tk.W)

      ExpDiscButton = Radiobutton(window_auto, text='Exponential disc', fg='blue', variable=self.ChooseModel, value=2)
      ExpDiscButton.grid(column=0,row=2, sticky=Tk.W)

      EdgeDiscButton = Radiobutton(window_auto, text='Edge-on disc', fg='blue', variable=self.ChooseModel, value=3)
      EdgeDiscButton.grid(column=0,row=3, sticky=Tk.W)      
      
      ExpDiscSersicButton = Radiobutton(window_auto, text='Exponential disc + Sersic', fg='blue', variable=self.ChooseModel, value=4)
      ExpDiscSersicButton.grid(column=0,row=4, sticky=Tk.W)         
      

      EdgeDiscSersicButton = Radiobutton(window_auto, text='Edge-on disc + Sersic', fg='blue', variable=self.ChooseModel, value=5)
      EdgeDiscSersicButton.grid(column=0,row=5, sticky=Tk.W)     

      ExpDiscAGNButton = Radiobutton(window_auto, text='Exponential disc + AGN', fg='blue', variable=self.ChooseModel, value=6)
      ExpDiscAGNButton.grid(column=0,row=6, sticky=Tk.W)      

      EdgeDiscAGNButton = Radiobutton(window_auto, text='Edge-on disc + AGN', fg='blue', variable=self.ChooseModel, value=7)
      EdgeDiscAGNButton.grid(column=0,row=7, sticky=Tk.W)    

      SersicSersicButton = Radiobutton(window_auto, text='Sersic + Sersic', fg='blue', variable=self.ChooseModel, value=8)
      SersicSersicButton.grid(column=0,row=8, sticky=Tk.W)    
      
      LaunchModel = Tk.Frame(window_auto)
      LaunchModel.grid(column=0, row=9)
      LaunchModelButton = Tk.Button(window_auto, text="Launch", state="normal",command=self.launch_auto, bg="green",font=("Helvetica", 10, "bold"), width=10)
      LaunchModelButton.grid(column=0, row=9)   


    #-------------------------------------------------------------


    #-------------------------------------------------------------
    # FINAL DECOMPOSITION FUNCTION
    def launch_dec(self):     
	if os.path.exists('model.fits'):
	   os.remove('model.fits')
  
	self.run_code()
	
	if not os.path.exists('model.fits'):
	  tkMessageBox.showerror('Warning','The code crashed! Please verify the input files and the command line!')
	  return 1







    #-------------------------------------------------------------

    def launch_auto(self):
      import do_mode
      if self.ChooseModel.get()==1:
	model='sersic'
      if self.ChooseModel.get()==2:
	model='exp_disc'
      if self.ChooseModel.get()==3:
	model='eon_disc'
      if self.ChooseModel.get()==4:
	model='sersic+exp_disc'
      if self.ChooseModel.get()==5:
	model='sersic+eon_disc'
      if self.ChooseModel.get()==6:
	model='agn+exp_disc'
      if self.ChooseModel.get()==7:
	model='agn+eon_disc'
      if self.ChooseModel.get()==8:
	model='sersic+sersic'
	
      #new_images,ADD_INFO,RA,DEC,galaxy_name,Sampling,Sky_level = self.inp_f()
      input_files,observation_info,object_info,settings = self.inp_f()
      [input_image,sigma_image,psf_image,mask_image,sampling] = input_files
      [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,Filter] = observation_info
      [RA,DEC,galaxy_name,Dist,Ext,Kcorr] = object_info

      add_info = [1,'unknown','usual',model]

      do_mode.main(input_files,observation_info,object_info,add_info,del_files=False)

      for file in ['fit.log','galfit.01','imfit.01']:
	if os.path.exists(file): 
	  os.remove(file)

    def inp_f(self):
		# Read the input
		file_field = self.file_field1.get()
		weight_image = self.weight_image1.get()
		psf_image = self.psf_image1.get()
		mask_image = self.mask_image1.get()
		Sampling = int(self.Sampling1.get())
		
		m0 = float(self.m01.get())
		pix2sec = float(self.pix2sec1.get())
		GAIN = float(self.GAIN1.get())
		NCOMBINE = float(self.NCOMBINE1.get())
		EXPTIME = float(self.EXPTIME1.get())
		RON = float(self.RON1.get())
		FWHM = float(self.FWHM1.get())
		Sky_level = float(self.Sky_level1.get())
		SkySubtr = self.SkySubtr1.get()
		Filter = float(self.Filter1.get().split(', ')[-1])
		
		RA = float(self.RA1.get())
		DEC = float(self.Dec1.get())
		if RA==0. and DEC==0.: ### NOTE: If a user didn't change the default values of RA=DEC=0 then RA=nx/2 and DEC=ny/2
		  RA = float(self.nx1/2.)
		  DEC = float(self.ny1/2.)
		else:
		  if self.CoordsCheck.get()==2:
		    RA,DEC = self.convert_coords()
		
		galaxy_name = str(self.galaxy_name1.get())
		Dist = self.Dist1.get()
		Ext = self.Ext1.get()
		Kcorr = self.Kcorr1.get()
		
		

		
		input_files = [file_field,weight_image,psf_image,mask_image,Sampling]
		observation_info = [self.nx1,self.ny1,m0,pix2sec,GAIN,NCOMBINE,EXPTIME,RON,FWHM,Sky_level,SkySubtr,Filter]
		object_info = [RA,DEC,galaxy_name,Dist,Ext,Kcorr]
		settings = [self.ChooseCode.get(),self.ChooseGeom.get(),self.ChooseSB.get(),self.ChooseLum.get()]
		
		return input_files,observation_info,object_info,settings
		#return new_images,ADD_INFO,RA,DEC,galaxy_name,Sampling,Sky_level

    def create_window(self, load_model_funcs):
      # Create window where you can choose the parameters for the new model
      window = Tk.Toplevel(self.master)
      if self.ChooseCode.get()=='GALFIT':      
	window.title("GALFIT FUNCTIONS")
      elif self.ChooseCode.get()=='IMFIT':      
	window.title("IMFIT FUNCTIONS") 

      if self.ChooseCode.get()=='GALFIT':
	      CompNames,ComponentFunctions = galfit_input.read_list_of_funcs()
      if self.ChooseCode.get()=='IMFIT':
	      CompNames,ComponentFunctions = imfit_input.read_list_of_funcs()
      
      
      Comps = Checkbar(window, CompNames, load_model_funcs)
      Comps.pack(side=TOP,  fill=X)

      Comps.config(relief=GROOVE, bd=2)

      Tk.Button(window, text='Quit', command=lambda: self.destroy_window(window)).pack(side=RIGHT)
      Tk.Button(window, text='Peek', command=lambda: self.allstates(Comps)).pack(side=RIGHT)

	
    def load_window(self):
      # With this function the user can upload the input galfit or imfit file and then its conent should be converted to the dictionary
      # Ask for the input file
      fileName = tkFileDialog.askopenfilename(parent=self.master,
				  filetypes=[("All files", "*")],
				  title="Open data file")
		  
      fileName = self.copy_file(fileName)  

      if self.ChooseCode.get()=='GALFIT':
	  CompNames,ComponentFunctions = galfit_input.read_list_of_funcs()
      if self.ChooseCode.get()=='IMFIT':
	  CompNames,ComponentFunctions = imfit_input.read_list_of_funcs()
      if fileName!=None:
	  if self.ChooseCode.get()=='GALFIT': 
		TMP_PARS,func_names = galfit_parser.to_deca_format(fileName, CompNames, ComponentFunctions)
	  
	  elif self.ChooseCode.get()=='IMFIT':
		load_imfit_model = imfit_parser.ImfitModel(fileName)
		TMP_PARS,func_names = load_imfit_model.to_deca_format()

	  # Save tmp file:
	  tmpfile = open('tmp_model.npy', 'w')
	  np.save(tmpfile, TMP_PARS)
	  tmpfile.close()
	
	  #print TMP_PARS
	  #exit()
	  self.create_window(func_names)
      
    def copy_file(self, fileName):
	    # Copy the file if it is not in the current directory:
	    if os.path.exists(fileName):
	      if os.path.dirname(fileName)==Current_Dir:
		fileName = fileName.split('/')[-1]
	      else:
		shutil.copy(fileName,fileName.split('/')[-1])
		fileName = fileName.split('/')[-1]
	    else:
	      fileName = None
	    return fileName      

    def destroy_window(self,w):
	w.destroy()

    def define_component(self,Component,TMP_COMPS,TMP_PARS):
	input_files,observation_info,object_info,settings = self.inp_f()
	[input_image,sigma_image,psf_image,mask_image,sampling] = input_files
	[nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,Filter] = observation_info
	[RA,DEC,galaxy_name,Dist,Ext,Kcorr] = object_info
	Scale = Dist * 1000. / 206265.

	if self.ChooseCode.get()=='GALFIT':
	    CompNames,ComponentFunctions = galfit_input.read_list_of_funcs()
	if self.ChooseCode.get()=='IMFIT':
	    CompNames,ComponentFunctions = imfit_input.read_list_of_funcs()


	if self.ChooseCode.get()=='GALFIT':
	    geom_pars = galfit_input.geom_pars()
	    lum_pars = galfit_input.lum_pars()
	    SB_pars = galfit_input.SB_pars()
	elif self.ChooseCode.get()=='IMFIT':
	    geom_pars = imfit_input.geom_pars()
	    lum_pars = imfit_input.geom_pars()
	    SB_pars = imfit_input.geom_pars()

	COMP = []
	for i in range(len(ComponentFunctions)):
	  if Component.split('/')[0] in ComponentFunctions[i]:
		func_name = ComponentFunctions[i].split(':')[0]
		func_par = ComponentFunctions[i].split(':')[1]
		func_comment = ComponentFunctions[i].split(':')[2]
		func_numb = Component.split('/')[1]

		if func_name+':'+func_par+':'+func_comment+'/'+func_numb not in TMP_COMPS: 
		  COMP.append(func_name+':'+func_par+':'+func_comment+'/'+func_numb+':'+'0'+','+''+','+'')
		else:
		  lll = TMP_COMPS.index(func_name+':'+func_par+':'+func_comment+'/'+func_numb)
		  tmp_arr = TMP_PARS[lll].split(':')[3].split(',')

		  if func_comment in geom_pars:
		    par_units = 'geom'
		  elif func_comment in lum_pars:
		    par_units = 'lum'
		  elif func_comment in SB_pars:
		    par_units = 'SB'  
		  else:
		    par_units = ''  

		  tmp_val = tmp_arr[0]
		  min_tmp_val = tmp_arr[1]
		  max_tmp_val = tmp_arr[2]

		  if self.ChooseCode.get()=='IMFIT':
		    VALUE_MIN = str(units_converter(min_tmp_val,par_units,'pix','ADU/arcsec2','ADU',self.ChooseGeom,self.ChooseSB,self.ChooseLum,m0,scale,Scale,Filter))
		    VALUE_MAX = str(units_converter(max_tmp_val,par_units,'pix','ADU/arcsec2','ADU',self.ChooseGeom,self.ChooseSB,self.ChooseLum,m0,scale,Scale,Filter))
		    VALUE = str(units_converter(tmp_val,par_units,'pix','ADU/arcsec2','ADU',self.ChooseGeom,self.ChooseSB,self.ChooseLum,m0,scale,Scale,Filter))
		  elif self.ChooseCode.get()=='GALFIT':
		    VALUE_MIN = str(units_converter(min_tmp_val,par_units,'pix','mag/arcsec2','mag',self.ChooseGeom,self.ChooseSB,self.ChooseLum,m0,scale,Scale,Filter))
		    VALUE_MAX = str(units_converter(max_tmp_val,par_units,'pix','mag/arcsec2','mag',self.ChooseGeom,self.ChooseSB,self.ChooseLum,m0,scale,Scale,Filter))
		    VALUE = str(units_converter(tmp_val,par_units,'pix','mag/arcsec2','mag',self.ChooseGeom,self.ChooseSB,self.ChooseLum,m0,scale,Scale,Filter))

		  COMP.append(func_name+':'+func_par+':'+func_comment+'/'+func_numb+':'+VALUE+','+VALUE_MIN+','+VALUE_MAX)

	return COMP

    def allstates(self,Comps):
	  #new_images,ADD_INFO,RA,DEC,galaxy_name,Sampling,Sky_level = self.inp_f()
	  #new_images[1] = self.FinalMask.get()
	  input_files,observation_info,object_info,settings = self.inp_f()
	  [input_image,sigma_image,psf_image,mask_image,sampling] = input_files
	  [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,Filter] = observation_info
	  [RA,DEC,galaxy_name,Dist,Ext,Kcorr] = object_info
	  [ChooseCode,ChooseGeom,ChooseSB,ChooseLum] = settings
	  Scale = Dist * 1000. / 206265.

	  input_info = [input_files,observation_info,RA,DEC,galaxy_name,sampling,sky_level,SkySubtr,ChooseCode,ChooseGeom,ChooseSB,ChooseLum,Scale]

	  if self.ChooseCode.get()=='GALFIT':
	      CompNames,ComponentFunctions = galfit_input.read_list_of_funcs()
	  if self.ChooseCode.get()=='IMFIT':
	      CompNames,ComponentFunctions = imfit_input.read_list_of_funcs()

      
	  try:
	    destroy_window(self.window_mod)
	  except:
	    zz = 1
	  list_of_comps = list(Comps.state_numb())
	  Components = []
	  for k in range(len(CompNames)):
	    if list_of_comps[k]>0:
	      for i in range(list_of_comps[k]):
		Components.append(CompNames[k]+'/'+str(i+1))

	  if os.path.exists('tmp_model.npy'):
	    # Load the temporary file
	    TMP_PARS = np.load('tmp_model.npy')
	    TMP_COMPS = []
	    for k in range(len(TMP_PARS)):
	      TMP_COMPS.append(TMP_PARS[k].split(':')[0]+':'+TMP_PARS[k].split(':')[1]+':'+TMP_PARS[k].split(':')[2])
	  else:
	    TMP_PARS = []
	    TMP_COMPS = []
	  
	  COMPS = []
	  for k in range(len(Components)):
		COMPS.append(self.define_component(Components[k],TMP_COMPS,TMP_PARS))
	  COMPS = list(flatten(COMPS))


	  if os.path.exists('tmp_model.npy'):
	    os.remove('tmp_model.npy')

	  self.window_mod = Tk.Toplevel(self.master)

	  add_info = [1,'unknown','usual','none']
	  CompsWindow = Parbar(input_files, observation_info, object_info, add_info, settings, self.window_mod, Components, COMPS)
	  CompsWindow.pack(side=TOP,  fill=X)

	  CompsWindow.config(relief=GROOVE, bd=2)
	  return self.window_mod

    def run_code(self):
      input_files,observation_info,object_info,settings = self.inp_f()
      [input_image,sigma_image,psf_image,mask_image,sampling] = input_files
      [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,sky_level,SkySubtr,Filter] = observation_info
      [RA,DEC,galaxy_name,Dist,Ext,Kcorr] = object_info
      [ChooseCode,ChooseGeom,ChooseSB,ChooseLum] = settings

      if show_code_run==True:
	subprocess.call(self.LaunchLine.get(), shell=True)
      else:
	subprocess.call(self.LaunchLine.get(), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
      
      if os.path.exists('model.fits'):
	  if ChooseCode=='GALFIT':
		import make_model_ima_galfit
		#1. Create new folder for this model
		arr = glob.glob('galfit_*')
		dir_numbers = []
		for ar in arr:
		  dir_numbers.append(int(ar.split('_')[-1]))
		if dir_numbers==[]:
		  new_dir = './galfit_1'
		else:
		  new_dir = './galfit_' + str(max(dir_numbers)+1)
		os.makedirs(new_dir)

		#2. Create log file if it hasn't been created
		if not os.path.exists("galfit.log"):
		  shutil.move("fit.log","galfit.log")
		else:
		  with open("fit.log") as f:
		    with open("galfit.log", "a") as f1:
			for line in f:
				f1.write(line) 
		
		#3. Create composed model
		make_model_ima_galfit.main(input_image, 'galfit.01')
		COPY_IMAGES(['galfit.01','composed_model.fits'],new_dir + '/')
		os.remove('galfit.01')
		
	  elif ChooseCode=='IMFIT':
		import make_model_ima_imfit
		#1. Create new folder for this model
		arr = glob.glob('imfit_*')
		dir_numbers = []
		for ar in arr:
		  dir_numbers.append(int(ar.split('_')[-1]))
		if dir_numbers==[]:
		  new_dir = './imfit_1'
		else:
		  new_dir = './imfit_' + str(max(dir_numbers)+1)
		os.makedirs(new_dir)

		#2. Create log file if it hasn't been created
		if not os.path.exists("imfit.log"):
		  shutil.copy("imfit.01","imfit.log")
		else:
		  with open("imfit.01") as f:
		    with open("imfit.log", "a") as f1:
			for line in f:
				f1.write(line) 
		
		#3. Create composed model
		make_model_ima_imfit.main(input_image, 'imfit.01', psf_image)
		COPY_IMAGES(['imfit.01','composed_model.fits'],new_dir + '/')
		os.remove('galfit.01')
	
  
def main(input_file=None):
	m = MainApp(input_file)
	m.master.mainloop()
	m.master.quit()

  
if __name__ == "__main__":
  if len(sys.argv)>1:
    input_file = str(sys.argv[1])
  else:
    input_file = None
  sys.exit(main(input_file))
  
