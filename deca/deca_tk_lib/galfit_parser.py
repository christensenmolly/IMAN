#! /usr/bin/env python

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
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from numpy import *
from pylab import *
import os
import shutil
import random
import time
import astropy.io.fits as pyfits
import re
import shelve
import threading
import inspect
import subprocess as sub
from collections import OrderedDict
from collections import Counter
import glob
import subprocess

def read_list_of_funcs():
        functions = ['sersic','sersic2','nuker','devauc','expdisk','edgedisk','moffat','ferrer','gaussian','king','psf','sky']
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
        #print functions,Galfit_functions
        #exit()
        return functions,Galfit_functions
    
def to_deca_format(fileName, CompNames, ComponentFunctions, galfitm=False):
            with open(fileName) as f:
                  lines = f.readlines()

            FinalParameters_load = OrderedDict()
            
            # Read componenet blocks
            Component_blocks = []
            fill_func=False

            
            for k in range(len(lines)):

              if ('0)' in lines[k]) and (lines[k].split()[1] in CompNames):
                #dupl_functions.append(lines[k].split()[1])
                if fill_func==True:
                  Component_blocks.append(Component_block)
                Component_block = []
                Component_block.append(lines[k])
                fill_func = True
                continue
              if fill_func==True and lines[k].startswith('#')==False and lines[k]!='\n':
                Component_block.append(lines[k])

            Component_blocks.append(Component_block)


            Dupl_Functions = {}
            
            # Define objects
            func_names = []
            for Component_block in Component_blocks:
              func_name = Component_block[0].split()[1]
              func_names.append(func_name)
              comp_name = func_names.count(func_name)
              if func_name not in Dupl_Functions:
                Dupl_Functions[func_name] = 1
              else:
                Dupl_Functions[func_name] = Dupl_Functions[func_name] + 1
              for par in ComponentFunctions:
                if par.split(':')[0]==func_name:
                    par_num = par.split(':')[1]
                    com_name = par.split(':')[2]
                    for i in range(len(Component_block)):
                      if par_num+')' in Component_block[i] and 'Sersic' not in Component_block[i]: #and re.match(r'\w[)]',Component_block[i].lstrip()):
                        if Component_block[i].lstrip().startswith('1)') and str(func_name)!='sky':
                          stri = Component_block[i].lstrip().split()
                          par_var_x = float(stri[1])
                          par_var_y = float(stri[2])
                          par_min_var = ''
                          par_max_var = ''
                          FinalParameters_load.update({str(func_name)+':'+str(1)+':'+'x_cen'+'/'+str(Dupl_Functions[func_name]):[par_var_x,par_min_var,par_max_var]})
                          FinalParameters_load.update({str(func_name)+':'+str(2)+':'+'y_cen'+'/'+str(Dupl_Functions[func_name]):[par_var_y,par_min_var,par_max_var]})      
                        else:
                          stri = Component_block[i].lstrip().split()
                          par_var = float(stri[1])
                          par_min_var = ''
                          par_max_var = ''
                          FinalParameters_load.update({str(func_name)+':'+str(par_num)+':'+str(com_name)+'/'+str(Dupl_Functions[func_name]):[par_var,par_min_var,par_max_var]})
                      elif par_num+')' in Component_block[i] and par_num=='5':
                          #print 'Sersic n line'
                          stri = Component_block[i].lstrip().split()
                          par_var = float(stri[1])
                          par_min_var = ''
                          par_max_var = ''
                          FinalParameters_load.update({str(func_name)+':'+str(par_num)+':'+str(com_name)+'/'+str(Dupl_Functions[func_name]):[par_var,par_min_var,par_max_var]}) 
            TMP_PARS = []
            for key, value in FinalParameters_load.items() :
              TMP_PARS.append(key+':'+str(value[0])+','+str(value[1])+','+str(value[2]))
              
            return TMP_PARS,func_names

#CompNames, ComponentFunctions = read_list_of_funcs()
#TMP_PARS,func_names = to_deca_format('galfit.01', CompNames, ComponentFunctions)
#print TMP_PARS,func_names
