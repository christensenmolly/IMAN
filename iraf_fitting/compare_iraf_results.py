#!/usr/bin/python
# DESCRIPTION:
# Compare results from two output files from IRAF/ELLIPSE

import sys
import math
import numpy as np
#from scipy import stats
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy import *
from pylab import *
import os
import shutil
import subprocess
import astropy.io.fits as pyfits
from astropy.stats import sigma_clipped_stats
import argparse
import glob


import iraf_ellipse

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

tmp_out = sys.stdout

DPI=150   #### TODO: CHANGE!!!!!!!!!!!!!!!!!!!!!!

fsize=16
FNULL = open(os.devnull, 'w')


def plot_iso(ax1,ax2,ax3,output_image,inten,inten_err,PA,errPA,ell,errell,B4,errB4,a,rmax,pix2sec,m0,overplot=False,comp=False,PA_ylim=0.,ell_ylim=0.,B4_ylim=0.,mag_ylim=0.,label=None):
        # Function to create output pictures
        mag = m0 - 2.5*log10(inten)+ 5.*log10(pix2sec)
        mag_err = fabs((2.5/log(10.0)) * inten_err/inten)


        if overplot==True:
          color = 'gray'
          color_bar = 'black'
          symb = 'o'
        else:
          color = 'red'
          color_bar = 'red'
          symb = '*'
          
        if comp==False:
          color = 'gray'
          color_bar = 'black'
          symb = 'o'  

        ax1.errorbar(a*pix2sec, fabs(PA),yerr=errPA,fmt=symb,color=color_bar, markersize=3)
        ax1.plot(a*pix2sec, fabs(PA),symb,color=color,markeredgecolor=color, markersize=3)
        ax1.set_ylabel(r'PA (deg)', fontsize=fsize)
        if overplot==False:
          if comp==False:
            ax1.set_ylim(np.min(fabs(PA)-errPA)*0.99,np.max(fabs(PA)+errPA)*1.01)
          else:
            [min_PA,max_PA] = PA_ylim

            ax1.set_ylim(min([np.min(fabs(PA)-errPA)*0.99,min_PA]),max([max_PA,np.max(fabs(PA)+errPA)*1.01]))
        else:
            PA_ylim = [np.min(fabs(PA)-errPA)*0.99,np.max(fabs(PA)+errPA)*1.01]

        ax2.errorbar(a*pix2sec, fabs(ell),yerr=errell,fmt=symb,color=color_bar, markersize=3)
        ax2.plot(a*pix2sec, fabs(ell),symb,color=color,markeredgecolor=color, markersize=3)
        ax2.set_ylabel(r'$\epsilon$ ', fontsize=fsize)
        if overplot==False:
          if comp==False:
            ax2.set_ylim(np.min(fabs(ell)-errell)*0.99,np.max(fabs(ell)+errell)*1.01)
          else:
            [min_ell,max_ell] = ell_ylim
            ax2.set_ylim(min([np.min(fabs(ell)-errell)*0.99,min_ell]),max([max_ell,np.max(fabs(ell)+errell)*1.01]))    
        else:
            ell_ylim = [np.min(fabs(ell)-errell)*0.99,np.max(fabs(ell)+errell)*1.01]

        ax3.errorbar(a*pix2sec, B4,yerr=errB4,fmt=symb,color=color_bar, markersize=3)
        mean_B4, median_B4, std_B4 = sigma_clipped_stats(B4, sigma=3.0, iters=5)
        ax3.plot(a*pix2sec, B4, symb,color=color,markeredgecolor=color, markersize=3)
        ax3.set_xlabel(r'r (arcsec)', fontsize=fsize)
        ax3.set_ylabel(r'B$_4$ ', fontsize=fsize)
        if overplot==False:
          if comp==False:
            if fabs(median_B4-5*std_B4)<2. and fabs(median_B4+5*std_B4)<2.:
              ax3.set_ylim(median_B4-5*std_B4,median_B4+5*std_B4)
            else:
              ax3.set_ylim(-2.,2.)
          else:
            [min_B4,max_B4] = B4_ylim
            if fabs(median_B4-5*std_B4)<2. and fabs(median_B4+5*std_B4)<2.:
              ax3.set_ylim(min([min_B4,median_B4-5*std_B4]),max([median_B4+5*std_B4,max_B4]))
            else:
              ax3.set_ylim(min([-2.,min_B4]),max([2.,max_B4]))    
        else:
            if fabs(median_B4-5*std_B4)<2. and fabs(median_B4+5*std_B4)<2.:
              B4_ylim=[median_B4-5*std_B4,median_B4+5*std_B4]
            else:
              B4_ylim=[-2,2]

  
        xticklabels = ax1.get_xticklabels()+ax2.get_xticklabels()
        setp(xticklabels, visible=False)

        for ax in [ax1,ax2,ax3]:
          for tick in ax.xaxis.get_major_ticks():
              tick.label.set_fontsize(fsize-6)
              if overplot==False:  
                    ax.set_xlim(0.,rmax*pix2sec)
        
        if overplot==False:
            plt.savefig(output_image, transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
            plt.clf()
            plt.close() 

        if overplot==False:
          return 0,0,0,0
        else:
          return PA_ylim,ell_ylim,B4_ylim,mag_ylim



def plot_azim(ax,a,inten,inten_err,pix2sec,m0,rmax,overplot,output_file,comp,xmax=None,ymin=None,ymax=None, label=None, SBlim=None, Rlim=None):
        # Function to plot azimuthally averaged profile
        mag = m0 - 2.5*log10(inten)+ 5.*log10(pix2sec)
        mag_err = fabs((2.5/log(10.0)) * inten_err/inten)


        if overplot==True:
          color = 'gray'
          color_bar = 'black'
          symb = 'o'
        else:
          color = 'red'
          color_bar = 'red'
          symb = '*'
          
        if comp==False:
          color = 'gray'
          color_bar = 'black'
          symb = 'o'  
        '''
        if comp==True:
          if overplot==True:
            label='galaxy'
          else:
            label='model'
        else:
            label=''
        '''



        '''
        if overplot==False:  
            # Overplot PSF 
            #file_ext_psf = '/home/amosenko/MyCurrentWork/Edge_on_HERON/PSF_HERON/NEW_PSF/psf_extended_for_paper.txt'
            #sma_psf, inten_psf = np.loadtxt(file_ext_psf, usecols=[0,1], unpack=True, skiprows = 1, dtype=float)
            
            sma_psf,inten_psf = np.loadtxt('/home/amosenko/MyCurrentWork/Rich_new/Test_PSF_impact/models/n_1/ellipse_psf_0.834.txt', usecols=[1,2], unpack=True, skiprows = 6, dtype=float)
            #print(m0 - 2.5*log10(inten_psf*inten[0]/inten_psf[0])+ 5.*log10(pix2sec))
            #exit()
            
            ax.plot(sma_psf*pix2sec,  m0 - 2.5*log10(inten_psf*inten[0]/inten_psf[0])+ 5.*log10(pix2sec), '^', color='lime',markeredgecolor='lime', markersize=5, zorder=0)
            #plt.show()
        '''

        ax.errorbar(a*pix2sec,mag,yerr=mag_err,fmt=symb,color=color,ecolor=color, markersize=5,label=label, zorder=0)
        ax.plot(a*pix2sec, mag, symb,color=color,markeredgecolor=color, markersize=3, zorder=1)
        mean_mag, median_mag, std_mag = sigma_clipped_stats(mag, sigma=3.0, iters=5)

        if overplot==True:  
            ax.set_ylim(max(mag)+std_mag,min(mag)-std_mag)
            ax.set_xlim(0.,rmax*pix2sec)
        else:
            if xmax==None:
                ax.set_xlim(0.,rmax*pix2sec)
            else:
                ax.set_xlim(0.,xmax)
            if ymax==None:
                ax.set_ylim(max(mag)+std_mag,min(mag)-std_mag)
            else:
                ax.set_ylim(32.,ymin)
                
            #mag_ylim = [max(mag)+std_mag,min(mag)-std_mag]
        
        if SBlim is not None and Rlim is not None:
            for k in range(len(SBlim)):
                ax.axhline(y=SBlim[k], linestyle='--', zorder=2, color='black', lw=1)
                #plt.axvline(x=Rlim[k]*pix2sec, linestyle='--')
                ax.plot([Rlim[k]*pix2sec],[SBlim[k]],'o', markersize=6, color=color, markeredgecolor='black', zorder=3)


        


        ax.set_xlabel(r'$r$ (arcsec)', fontsize=fsize)
        ax.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
        
        if comp==True:
          plt.legend(loc=0, borderaxespad=0.,fontsize=10,numpoints=1,markerscale=2)
        
        
        if overplot==False:  
            plt.savefig(output_file, transparent = False, dpi=DPI, bbox_inches='tight', pad_inches=0.05)
            plt.clf()
            plt.close()
            return None,None,None
        else:
            return rmax*pix2sec, min(mag)-std_mag,max(mag)+std_mag


def main(ellipse_file_1, ellipse_file_2, m01, m02, pix2sec1, pix2sec2, label1=None, label2=None, Rmax=None, output_files = None, Rlims1=None, Rlims2=None, SBlims1 = None, SBlims2 = None):
          #f= figure(0,figsize=(4,10))
          
          if output_files is not None:
              output_file1 = output_files[0]
              output_file2 = output_files[1]
          else:
              output_file1 = 'azim_aver.png'
              output_file2 = 'iraf.png'    

          fig = figure(1,figsize=(5,5))
          ax4 = fig.add_subplot(111)

          #gs = gridspec.GridSpec(3, 1)
          #gs.update(left=0.25, right=3.00, hspace=0.0)

          #ax1 = plt.subplot(gs[0])
          #ax2 = plt.subplot(gs[1])
          #ax3 = plt.subplot(gs[2])

          N=0
          PA_ylim=0.;ell_ylim=0.;B4_ylim=0.;mag_ylim=0.
          mag_min = None
          mag_max = None
          xmax = None

          overplot = True
          comp = True
          pix2sec = [pix2sec1, pix2sec2]
          m0 = [m01, m02]
          ell_files = [ellipse_file_1, ellipse_file_2]
          labels = [label1, label2]
          SBlims = [SBlims1, SBlims2]
          Rlims = [Rlims1, Rlims2]
          
          for kk in range(len(ell_files)):
            ell_file = ell_files[kk]
            sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = np.loadtxt(ell_file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 6, dtype='str')

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
            
            if Rmax is None:
                Rmax = float(np.max(sma))

            if N>0:
                overplot = False
            #PA_ylim,ell_ylim,B4_ylim,mag_ylim = iraf_ellipse.plot_iso(ax1,ax2,ax3,'iraf.%s' % (outp_format),inten,inten_err,PA,errPA,ell,errell,B4,errB4,sma,float(np.max(sma)),pix2sec[kk],m0[kk],outp_format,overplot,comp=comp,\
            #  PA_ylim=PA_ylim,ell_ylim=ell_ylim,B4_ylim=B4_ylim,mag_ylim=mag_ylim)
            
            xmax,mag_min,mag_max = plot_azim(ax4,sma,inten,inten_err,pix2sec[kk],m0[kk],Rmax,overplot,output_file1,comp,xmax,mag_min,mag_max,label=labels[kk], SBlim=SBlims[kk], Rlim=Rlims[kk])
            
            #plot_azim(ax4,sma,inten,inten_err,pix2sec,m0,maxsma,overplot,outp_format,comp)
            
            N = N + 1









          f= figure(0,figsize=(4,10))

          #fig = figure(1,figsize=(5,5))
          #ax4 = fig.add_subplot(111)

          gs = gridspec.GridSpec(3, 1)
          gs.update(left=0.25, right=3.00, hspace=0.0)

          ax1 = plt.subplot(gs[0])
          ax2 = plt.subplot(gs[1])
          ax3 = plt.subplot(gs[2])

          N=0
          PA_ylim=0.;ell_ylim=0.;B4_ylim=0.;mag_ylim=0.
          mag_min = None
          mag_max = None
          xmax = None

          overplot = True
          comp = True
          pix2sec = [pix2sec1, pix2sec2]
          m0 = [m01, m02]
          ell_files = [ellipse_file_1, ellipse_file_2]

          
          SMA = []; INTEN = []
          for kk in range(len(ell_files)):
            ell_file = ell_files[kk]
            sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = np.loadtxt(ell_file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 6, dtype='str')

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
            
            if Rmax is None:
                Rmax = float(np.max(sma))

            if N>0:
                overplot = False
            PA_ylim,ell_ylim,B4_ylim,mag_ylim = plot_iso(ax1,ax2,ax3,output_file2,inten,inten_err,PA,errPA,ell,errell,B4,errB4,sma,Rmax,pix2sec[kk],m0[kk],overplot,comp=comp,\
              PA_ylim=PA_ylim,ell_ylim=ell_ylim,B4_ylim=B4_ylim,mag_ylim=mag_ylim,label=labels[kk])
            
            #xmax,mag_min,mag_max = iraf_ellipse.plot_azim(ax4,sma,inten,inten_err,pix2sec[kk],m0[kk],float(np.max(sma)),overplot,outp_format,comp,xmax,mag_min,mag_max)
            
            #plot_azim(ax4,sma,inten,inten_err,pix2sec,m0,maxsma,overplot,outp_format,comp)
            
            N = N + 1
            SMA.append(sma)
            INTEN.append(inten)

          return SMA, INTEN        


'''
m0 = 28.49836973762612
pix2sec = 1.6674448657806298    

main('ellipse.txt', 'ellipse_surf.txt', m0, m0, pix2sec, pix2sec, label1='Mos-original', label2='Mos-surf',Rmax=400./pix2sec, output_files=['Mos_azim_aver.png','Mos_iraf.png'])
main('Ellipse.txt', 'Isofit.txt', m0, m0, pix2sec, pix2sec, label1='Fusco-original', label2='Fusco-surf',Rmax=400./pix2sec, output_files=['Fusco_azim_aver.png','Fusco_iraf.png'])
main('ellipse_surf.txt', 'Isofit.txt', m0, m0, pix2sec, pix2sec, label1='Mos-surf', label2='Fusco-surf',Rmax=400./pix2sec, output_files=['Mos_Fusco_azim_aver.png','Mos_Fusco_iraf.png'])
'''
