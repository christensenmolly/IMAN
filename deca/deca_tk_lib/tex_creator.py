#! /usr/bin/env python

import pylab
#import Pmw
#from Tkinter import *
#import Tkinter as Tk
#import tkFileDialog
#import tkMessageBox
import astropy.io.fits as pyfits
import os
import shutil
import shelve
import numpy as np
import subprocess
import math
from math import *
import sys
import pickle
import collections

FNULL = open(os.devnull, 'w')

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
sys.path.append(PATH_TO_PACKAGE)

import imfit_input
import galfit_input
#import deca_tk
from misc_functions import units_converter
import imfit_parser
import galfit_parser
import galfitm_parser
import deca_setup

def read_pickle():
    with open('results.pkl','rb') as file_object:
        raw_data = file_object.read()

    results = pickle.loads(raw_data)
    return results  

def write_pickle(results=None):
    output = open('results.pkl', 'wb')
    if results is None:
        results = collections.OrderedDict()
        results['Collected parameters']=''
   
    pickle.dump(results, output)
    output.close()  

def header():
  s = r"""
\documentclass[12pt,english]{article}
\usepackage[a4paper,bindingoffset=0.2in,%
            left=1in,right=1in,top=1in,bottom=1in,%
            footskip=.25in]{geometry}
\usepackage{amsmath}    % need for subequations
\usepackage{graphicx}   % need for figures
\usepackage{verbatim}   % useful for program listings
\usepackage{color}      % use if color is used in text
\usepackage{subfigure}  % use for side-by-side figures
\usepackage{hyperref}   % use for hypertext links, including those to external documents and URLs
\usepackage{fancyhdr}
\usepackage{longtable}
\usepackage{siunitx}
\usepackage{datetime}
\usepackage{listings}
\lstset{
   breaklines=true,
   basicstyle=\ttfamily}


\sisetup{round-mode=places,round-precision=2}



\begin{comment}
\pagestyle{empty} % use if page numbers not wanted
\end{comment}


% above is the preamble

\begin{document}

\begin{center}"""
  return s


def find_lums(composed_model, comps):
    hdulist = pyfits.open(composed_model)
    model_lum = np.nansum(hdulist[1].data)
    lums = []

    for k in range(len(comps)):
      if comps[k]!='sky':
        lums.append(np.nansum(hdulist[4+k].data))
    return model_lum,lums
    
    

def split_lines(L):
  List = []
  for k in range(len(L)):
    if L[k].split(':')[0]!='sky':
      List.append(L[k].split(':')[0]+'/'+L[k].split('/')[1].split(':')[0])
  #print List
  
  name = List[0]
  first_lines = []
  cont_lines = []
  last_lines = []
  first_lines.append(0)
  
  for k in range(1,len(List)-1):
    if List[k]==name and List[k+1]==List[k]:
      cont_lines.append(k)
    elif List[k]==name and List[k+1]!=List[k]:
      last_lines.append(k)
    elif List[k]!=List[k-1]:
      first_lines.append(k)
      name = List[k]
      
  last_lines.append(len(List)-1) 
  return first_lines,cont_lines,last_lines

def res_table(object_name,file,composed_model,code,m0,scale,Distance,Scale,Filter,Aext=0.,Kcorr=0.,geom_units='kpc',lum_units='mag',SB_units='mag/arcsec2'):

  
  tex_line = []

  #tex_line.append('{\\large {\color{red} Object %s}} \\\\' % (object_name))

  tex_line.append('\\section{%s}\n' % (object_name))
  tex_line.append('\\label{%s}\n' % (object_name))

  #tex_line.append('\\today, \\currenttime')
  #tex_line.append('\\begin{verbatim} ' + 'Output directory:' + os.getcwd() + '\\end{verbatim}\\\\')
  tex_line.append('\\begin{lstlisting}\n')
  tex_line.append('Output directory:\n' + os.path.dirname(file))
  tex_line.append('\\end{lstlisting}\n')
  tex_line.append('\\begin{verbatim} ' + 'Code:' + code + '\\end{verbatim}\\\\\n')
  tex_line.append('Zero-point=%.5f, scale=%.3f~arcsec/pix, Distance=%.3f~Mpc, Scale=%.3f~kpc/arcsec, Extinction(at %s~nm)=%.3f~mag, K-correction=%.3f~mag.\\\\\n' % (m0,scale,Distance,Scale,str(Filter),Aext,Kcorr))

  #tex_line.append('\\begin{center}\n')
  tex_line.append('\\begin{longtable}{|c|c|c|c|}\n')
  tex_line.append('\\caption{Results of decomposition for %s.} \label{Table%s} \\\\\n' % (object_name,object_name))
  tex_line.append('\\hline\n')
  #\\multicolumn{1}{|c|}{\\textbf{Component}} & \\multicolumn{1}{c|}{\\textbf{Parameter}} & \\multicolumn{1}{c|}{\\textbf{Value}} & \\multicolumn{1}{c|}{\\textbf{Units}} \\\\ )
  tex_line.append('\\textbf{Component} & \\textbf{Parameter} & \\textbf{Value} & \\textbf{Units} \\\\\n')
  tex_line.append('\\hline\n')
  tex_line.append('\\endfirsthead\n')
  tex_line.append('\\multicolumn{4}{c}\n')
  tex_line.append('{{\\bfseries \\tablename\\ \\thetable{} -- continued from previous page}} \\\\\n')
  tex_line.append('\\hline\n')
  #tex_line.append('\\hline \\multicolumn{1}{|c|}{\\textbf{Component}} & \\multicolumn{1}{c|}{\\textbf{Parameter}} & \\multicolumn{1}{c|}{\\textbf{Value}} & \\multicolumn{1}{c|}{\\textbf{Units}} \\\\ \hline')
  tex_line.append('\\textbf{Component} & \\textbf{Parameter} & \\textbf{Value} & \\textbf{Units} \\\\\n')
  tex_line.append('\\hline\n')
  tex_line.append('\\endhead\n')
  tex_line.append('\\hline \\multicolumn{4}{r}{{Continued on next page}} \\\\ \\hline\n')  
  tex_line.append('\\endfoot\n')
  tex_line.append('\\hline\n')
  tex_line.append('\\endlastfoot\n')    


  num_func = -1
  if code=='GALFIT':
    CompNames,ComponentFunctions = galfit_input.read_list_of_funcs()
    geom_pars = galfit_input.geom_pars()
    lum_pars = galfit_input.lum_pars()
    SB_pars = galfit_input.SB_pars()
    decaList,func_names = galfit_parser.to_deca_format(file, CompNames, ComponentFunctions)
  elif code=='IMFIT':
    geom_pars = imfit_input.geom_pars()
    lum_pars = imfit_input.lum_pars()
    SB_pars = imfit_input.SB_pars()
    model = imfit_parser.ImfitModel(file)
    decaList,func_names = model.to_deca_format()
  if code=='GALFITM':
    CompNames,ComponentFunctions = galfit_input.read_list_of_funcs()
    geom_pars = galfit_input.geom_pars()
    lum_pars = galfit_input.lum_pars()
    SB_pars = galfit_input.SB_pars()
    decaList,func_names = galfitm_parser.to_deca_format(file, CompNames, ComponentFunctions)
    
  first_lines,cont_lines,last_lines = split_lines(decaList)
  #print first_lines,cont_lines,last_lines
  #exit()

  model_lum,lums = find_lums(composed_model, func_names)
  
  new_func=True
  last_func = ''
  N_comp = 0
  
  write_pickle()
  results = read_pickle()
  
  for k in range(len(decaList)):
      
      line = decaList[k].split(':')
      func = line[0]

      if func=='sky' or func=='FlatSky':
        continue
      par = line[1]
      com_name = line[2].split('/')[0]
      if code!='GALFITM':
        val = float(line[3].split(',')[0])
      else:
        val = list(np.array(line[3].split(',')[0].split(), float))        
      
      
      
      units = ''
      unit_cat = ''
      if com_name in geom_pars:
        units = geom_units
        unit_cat = 'geom'
      elif com_name in lum_pars:
        units = lum_units
        unit_cat = 'lum'
      elif com_name in SB_pars:
        units = SB_units
        if SB_units=='mag/arcsec2':
          units = 'mag/arcsec$^2$'
        unit_cat = 'SB'
      elif com_name=='PA' or com_name=='inc':
        units = 'deg'
      elif com_name=='X0' or com_name=='Y0' or com_name=='x_cen' or com_name=='y_cen':
        units = 'pix'
        if com_name=='X0' or com_name=='x_cen':
          com_name='x_0'
        if com_name=='Y0'or com_name=='y_cen':
          com_name='y_0'
      if 'mu_' in com_name:
          com_name = '\mu_\mathrm{' + com_name.split('_')[-1] + '}'
      
      
      if code=='IMFIT':
        value = units_converter(val,unit_cat,'pix','ADU/pix2','ADU',geom_units,
                    SB_units,lum_units,m0,scale,Distance,Scale,Filter,Aext,Kcorr)
        #apa_mag = m0 - 2.5*log10(lums[k]) - Aext - Kcorr
        #abs_mag = apa_mag - 5.*log10(Distance) - 25.
      elif code=='GALFIT':
        #print(val,unit_cat,'pix','mag/arcsec2','mag',geom_units,
        #            SB_units,lum_units,m0,scale,Distance,Scale,Filter,Aext,Kcorr)
        #exit()
        value = units_converter(val,unit_cat,'pix','mag/arcsec2','mag',geom_units,
                    SB_units,lum_units,m0,scale,Distance,Scale,Filter,Aext,Kcorr)
        #apa_mag = m0 - 2.5*log10(lums[k]) - Aext - Kcorr
        #abs_mag = apa_mag - 5.*log10(Distance) - 25.
      elif code=='GALFITM':
        value = ''
        for ii in range(len(val)):  
            valu = units_converter(val[ii],unit_cat,'pix','mag/arcsec2','mag',geom_units,
                        SB_units,lum_units,m0,scale,Distance,Scale,Filter,Aext,Kcorr)  
            if ii!=len(val):
                value = value + '\\num{'+str(valu)+'}, '
            else:
                value = value + '\\num{'+str(valu)+'}'     

            #apa_mag = m0 - 2.5*log10(lums[k][ii]) - Aext[ii] - Kcorr[ii]
            #abs_mag = apa_mag - 5.*log10(Distance) - 25.


      if '_' in com_name:
        com_name = com_name.split('_')[0]+'_'+'{'+com_name.split('_')[1]+'}'

      '''
      if last_func!=func:
        new_func=True
      else:
        new_func=False

      if new_func==True:
        print func
        num_func = num_func + 1
        if last_func!='':
          tex_line.append('      & $frac$ & %s &   \\\\[+0.5ex]' % (luminosities[num_func]))
        tex_line.append('    %i. %s & $%s$ & %s & %s \\\\[+0.5ex]' % (num_func+1, func, com_name, value, units))
      else:
        tex_line.append('        & $%s$ & %s & %s \\\\[+0.5ex]' % (com_name, value, units))
      last_func = func
      '''
      if k in first_lines:
        num_func = num_func + 1
        tex_line.append('    %i. %s & $%s$ & %s & %s \\\\[+0.5ex]\n' % (num_func+1, func, com_name, value, units))
      elif k in cont_lines:
        tex_line.append('        & $%s$ & %s & %s \\\\[+0.5ex]\n' % (com_name, value, units))
      elif k in last_lines:
        tex_line.append('        & $%s$ & %s & %s \\\\[+0.5ex]\n' % (com_name, value, units))
        tex_line.append('      & $frac$ & %s &   \\\\[+0.5ex]\n' % (round(lums[N_comp]/model_lum,3)))
        N_comp = N_comp + 1
      if unit_cat=='geom':
        results['%s: %s [%s]' % (func, com_name, 'arcsec')] = val*scale
      results['%s: %s [%s]' % (func, com_name, units)] = value

  apa_mag_total = m0 - 2.5*log10(model_lum) - Aext - Kcorr
  abs_mag_total = apa_mag_total - 5.*log10(Distance) - 25.
      
  tex_line.append('{\\bf Total} & $m_\mathrm{tot}$ & %.2f & AB-mag \\\\[+0.5ex]\n' % (float(apa_mag_total)) )   
  tex_line.append('             & $L_\mathrm{tot}$ & %.2f & AB-mag \\\\[+0.5ex]\n' % (float(abs_mag_total)) )    
  tex_line.append('\\end{longtable}\n')
  #tex_line.append('\\end{center}\n')
  #print(results)
  write_pickle(results)
  
  return tex_line

def add_figure(Number,object_name,image_file,caption,width=15.0):
  line= '\n%%Figure %i\n\\begin{figure}[h!]\n\\centering\n\\includegraphics[width=%fcm, angle=0, clip=]{%s}\n\\caption{%s}\n\\label{Figure_%s}\n\\end{figure}\n' % (Number,width,image_file,caption,object_name)
  return line


















def end():
  line = '\\end{center}\n\\end{document}'
  return line

def define_figure(picture, object_name, plot_2d_allignment='line'):
  if 'isomap' in picture:
    caption = 'Isophote map of the object %s.' % (object_name)
    width = 10.0
  elif '_prof_ver.' in picture:
    caption = 'Summed vertical profiles of the object %s.' % (object_name)
    width = 10.0
  elif '_cut.' in picture:
    caption = 'Photometric cut of the object %s.' % (object_name)
    width = 10.0
  elif '_prof_rad.' in picture:
    caption = 'Pixel distribution along the radius for the object %s.' % (object_name)
    width = 10.0
  elif '_prof_sum.' in picture:
    caption = 'Summed horizontal profiles of the object %s.' % (object_name)
    width = 10.0
  elif '_prof_azim.' in picture:
    caption = 'Azimuthally averaged profile of the object %s.' % (object_name)
    width = 10.0
  elif 'plot_2d.' in picture:
    if plot_2d_allignment == 'line':
        caption = '%s. The object image (left-hand), the best fitting image (middle), and the residual image which indicates the relative deviation between the fit and the image in absolute values (right-hand). The right-hand color bar shows the scaling of the residual maps.' % (object_name)    
    else:
        caption = '%s. The object image (top), the best fitting image (middle), and the residual image which indicates the relative deviation between the fit and the image in absolute values (bottom). The bottom color bar shows the scaling of the residual maps.' % (object_name)
    width = 15.
  #print caption
  return caption,width

  

def main(object_name,output_file,composed_model,code,pictures,m0,scale,Distance,Scale,Filter,Aext,Kcorr,geom_units,lum_units,SB_units, continue_file=False, last_object = True, plot_2d_allignment='line'):
  if continue_file == False:
    ff = open('results.tex','w')  
    ff.write(header())
  else:
    ff = open('results.tex','a')  
    
  lines_table = res_table(object_name,output_file,composed_model,code,m0,scale,Distance,Scale,Filter,Aext,Kcorr,geom_units,lum_units,SB_units)
  
  # FINAL TABLE
  for k in range(len(lines_table)):
    ff.write(lines_table[k])

  if last_object == True:
    # Figures
    for k in range(len(pictures)):
        caption,width = define_figure(pictures[k], object_name, plot_2d_allignment=plot_2d_allignment)
        figure = add_figure(k+1, object_name, pictures[k],caption, width=width)
        ff.write(figure)
    ff.write('\\clearpage\n')
    end_line = end()
    ff.write(end_line)
    ff.close()
  subprocess.call("pdflatex results.tex", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
  subprocess.call("pdflatex results.tex", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
  os.remove('results.aux')
  os.remove('results.log')
  os.remove('results.out')
    
#main('PGC10198','galfit.01','composed_model.fits','GALFIT',['plot_2d.png','major_axis_cut.png'],28.0,0.75,30.,0.1,3600.,0.01,0.,'kpc','mag','mag/arcsec2')





