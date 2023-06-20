#!/usr/bin/python
# DESCRIPTION:
# STSDAS IRAF/ELLIPSE fitting. IRAF is required. Do not create your own login.cl where you launch this script. 
# MINIMAL USAGE: python  iraf_ellipse.py [input_image]
# EXAMPLE:  python iraf_ellipse.py galaxy_clean_galf.fits --xc 1122 --yc 1021 --maxsma 850 --ZP 29.5962 --pix2sec 0.833 --mask mask.fits

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
import warnings
warnings.filterwarnings("ignore")


matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

tmp_out = sys.stdout


fsize=16
FNULL = open(os.devnull, 'w')


def create_reg_with_ellipses(x, y, sma, ell, PA):
    # Function to create region file with IRAF/ELLIPSE ellipses
    f = open('iraf_ellipses.reg', 'w')
    f.write('%s\n' % ('image') )
    for k in range(len(x)):
        f.write('ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # color=red\n' % (x[k], y[k], sma[k], sma[k]*(1.-ell[k]), PA[k]+90.))
    f.close()    


def crea_ell(file_in, xc, yc, ellip=0.2, pa=20., sma=10., step=0.03, minsma=0.0, maxsma=-1., ell_file='ellipse.txt', fits_mask="", fflag=1.0, olthresh=0.0, linear='no', hcenter='no'):
        if os.path.exists('gal.tab'):
            os.remove('gal.tab')

        
        
        # Function to create cl script to run IRAF/ELLIPSE
        f = open("ell.cl", "w") 
        sys.stdout = f

        print("# Script: ellipse")
        print("!rm -f ellipse.txt")
        print("stsdas")
        print("analysis")
        print("isophote")
        
        print("geompar.linear=%s" % (linear))
        print("geompar.step=%.2f" % (step))
        print("geompar.minsma=%.2f" % (minsma))
        if maxsma==-1.:
            print("geompar.maxsma=INDEF")
        else:
            print("geompar.maxsma=%.2f" % (maxsma))
        
        print("geompar.x0=%.3f" % (xc))
        print("geompar.y0=%.3f" % (yc))
        print("geompar.ellip0=%.3f" % (ellip))
        print("geompar.pa0=%.3f" % (pa))
        print("geompar.sma0=%.3f" % (sma))
        
        print("samplepar.fflag=%.1f" % (fflag)) # Acceptable fraction of flagged data points in intensity sample.
        
        print("controlpar.hcenter=%s" % (hcenter)) 
        print("controlpar.olthresh=%.1f" % (olthresh)) # Threshold for the object locator algorithm. By lowering this value the object locator becomes less strict, in the sense that it will accept lower signal-to-noise data. If set to zero, the x0, y0 values found in the geompar pset are used without questioning.
        #if fits_mask!="no":
        #    print("geompar.dqf=%s" % (fits_mask))
        print("controlpar.minit=10")
        print("controlpar.maxit=100")
        
        print("ellipse %s gal.tab dqf=%s" % (file_in, fits_mask))
        print("tprint gal.tab pwidth=600 plength=5000 > %s" % (ell_file))
        print("!rm -f gal.tab")
        print("logout")

        sys.stdout = tmp_out
        f.close()


def plot_iso(ax1,ax2,ax3,output_image,inten,inten_err,PA,errPA,ell,errell,B4,errB4,a,rmax,pix2sec,m0,outp_format='eps',overplot=False,comp=False,PA_ylim=0.,ell_ylim=0.,B4_ylim=0.,mag_ylim=0.):
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
        mean_B4, median_B4, std_B4 = sigma_clipped_stats(B4, sigma=3.0)
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
          ax.get_xaxis().set_tick_params(direction='in', width=1)
          ax.get_yaxis().set_tick_params(direction='in', width=1)
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



def plot_azim(ax,a,inten,inten_err,pix2sec,m0,rmax,overplot,outp_format,comp,xmax=None,ymin=None,ymax=None):
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

        if comp==True:
          if overplot==True:
            label='galaxy'
          else:
            label='model'
        else:
            label=''

        ax.errorbar(a*pix2sec,mag,yerr=mag_err,fmt=symb,color=color,ecolor=color, markersize=5,label=label)
        ax.plot(a*pix2sec, mag, symb,color=color,markeredgecolor=color, markersize=3)
        mean_mag, median_mag, std_mag = sigma_clipped_stats(mag, sigma=3.0)

        if overplot==True:  
            ax.set_ylim(max(mag)+std_mag,min(mag)-std_mag)
            ax.set_xlim(0.,rmax*pix2sec)
        else:
            if xmax==None:
                ax.set_xlim(0.,rmax*pix2sec)
            else:
                ax.set_xlim(0.,xmax)
            if ymax==None:
                ax.set_ylim(np.nanmax(mag)+std_mag,np.nanmin(mag)-std_mag)
            else:
                ax.set_ylim(ymax,ymin)
                
            #mag_ylim = [max(mag)+std_mag,min(mag)-std_mag]

        ax.set_xlabel(r'$r$ (arcsec)', fontsize=fsize)
        ax.set_ylabel(r'$\mu$ (mag arcsec$^{-2}$)', fontsize=fsize)
        
        if comp==True:
          plt.legend(loc=0, borderaxespad=0.,fontsize=10,numpoints=1,markerscale=2)
        
        ax.get_xaxis().set_tick_params(direction='in', width=1)
        ax.get_yaxis().set_tick_params(direction='in', width=1)
        if overplot==False:  
            plt.savefig("azim_aver."+outp_format, transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)
            plt.clf()
            plt.close()
            return None,None,None
        else:
            return rmax*pix2sec, min(mag)-std_mag,max(mag)+std_mag


def main_ell(file_in, xc, yc, ellip=0.2, pa=20., sma0=10., m0=28., pix2sec=1.0, step=0.03, minsma=0., maxsma=None, outp_format='png', ell_file='ellipse.txt', fits_mask="", fflag=1.0, olthresh=0.0, linear='no', hcenter='no'):
        sma0 = math.ceil(sma0)
        if fits_mask=="":
            fits_mask = "no"

        f= figure(0,figsize=(4,10))

        gs = gridspec.GridSpec(3, 1)
        gs.update(left=0.25, right=3.00, hspace=0.0)

        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        ax3 = plt.subplot(gs[2])

        hdulist = pyfits.open(file_in)
        comp = False

        if len(hdulist)>1:
          comp = True
          data_1 = hdulist[0].data
          header = hdulist[0].header
          ny,nx = np.shape(data_1)
          if xc==None and yc==None:
              xc = nx/2.
              yc = ny/2.
          else:
              xc = float(xc)
              yc = float(yc)
          if maxsma==None:
              maxsma = math.sqrt( (nx/2.)**2 + (ny/2.)**2 )
          else:
              maxsma = float(maxsma)
              
          data_2 = hdulist[1].data
          #hdu = pyfits.PrimaryHDU(data_1)
          #hdu.writeto('1.fits',clobber=True)  

          outHDU = pyfits.PrimaryHDU(data_1, header=header)
          outHDU.writeto('1.fits', overwrite=True)  

          #hdu = pyfits.PrimaryHDU(data_2)
          #hdu.writeto('2.fits',clobber=True)  

          outHDU = pyfits.PrimaryHDU(data_2, header=header)
          outHDU.writeto('2.fits', overwrite=True)  


          N=0
          PA_ylim=0.;ell_ylim=0.;B4_ylim=0.;mag_ylim=0.
          overplot = True
          for file in ['1.fits','2.fits']:
            if file=='2.fits':
                fits_maskk = "no"
            else:
                fits_maskk = fits_mask
            crea_ell(file, xc, yc, ellip=ellip, pa=pa, sma=sma0, step=step, minsma=minsma, maxsma=float(maxsma), ell_file=ell_file, fits_mask=fits_maskk, fflag=fflag, olthresh=olthresh, linear=linear, hcenter=hcenter)
            os.chmod(r"ell.cl",0o777)
            subprocess.call("cl < ell.cl -o", shell=True)
            #exit()
            sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = loadtxt(ell_file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 5, dtype='str')

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

            if N>0:
                overplot = False
            PA_ylim,ell_ylim,B4_ylim,mag_ylim = plot_iso(ax1,ax2,ax3,'iraf.%s' % (outp_format),inten,inten_err,PA,errPA,ell,errell,B4,errB4,sma,float(maxsma),pix2sec,m0,outp_format,overplot,comp=comp,\
              PA_ylim=PA_ylim,ell_ylim=ell_ylim,B4_ylim=B4_ylim,mag_ylim=mag_ylim)
            
            #plot_azim(ax4,sma,inten,inten_err,pix2sec,m0,maxsma,overplot,outp_format,comp)
            
            N = N + 1


          fig = figure(1,figsize=(5,5))
          ax4 = fig.add_subplot(111)


          N=0
          PA_ylim=0.;ell_ylim=0.;B4_ylim=0.;mag_ylim=0.
          overplot = True
          mag_min = None
          mag_max = None
          xmax = None
          for file in ['1.fits','2.fits']:
            #crea_ell(file, xc, yc, ellip=ellip, pa=pa, sma=sma0, step=step, minsma=minsma, maxsma=float(maxsma), ell_file=ell_file, fits_mask=fits_mask, fflag=fflag, olthresh=olthresh, linear=linear, hcenter=hcenter)
            #os.chmod(r"ell.cl",0o777)
            #subprocess.call("cl < ell.cl -o", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

            if file=='2.fits':
                fits_maskk = "no"
            else:
                fits_maskk = fits_mask
                
            hdulist = pyfits.open(file)
            data_1 = hdulist[0].data
            ny,nx = np.shape(data_1)
            if xc==None and yc==None:
                xc = nx/2.
                yc = ny/2.
            else:
                xc = float(xc)
                yc = float(yc)
            if maxsma==None:
                maxsma = math.sqrt( (nx/2.)**2 + (ny/2.)**2 )
            else:
                maxsma = float(maxsma)

            crea_ell(file, xc, yc, ellip=ellip, pa=pa, sma=sma0, step=step, minsma=minsma, maxsma=float(maxsma), ell_file=ell_file, fits_mask=fits_maskk, fflag=fflag, olthresh=olthresh, linear=linear, hcenter=hcenter)
            os.chmod(r"ell.cl",0o777)
            subprocess.call("cl < ell.cl -o", shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)

            sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = loadtxt(ell_file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 5, dtype='str')

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

            if N>0:
                overplot = False
            #PA_ylim,ell_ylim,B4_ylim,mag_ylim = plot_iso(ax1,ax2,ax3,output_image,inten,inten_err,PA,errPA,ell,errell,B4,errB4,sma,maxsma,pix2sec,m0,outp_format,overplot,comp=comp,\
            #PA_ylim=PA_ylim,ell_ylim=ell_ylim,B4_ylim=B4_ylim,mag_ylim=mag_ylim)
    
            xmax,mag_min,mag_max = plot_azim(ax4,sma,inten,inten_err,pix2sec,m0,float(maxsma),overplot,outp_format,comp,xmax,mag_min,mag_max)
            #print xmax,mag_min,mag_max
            N = N + 1

        else:
            hdulist = pyfits.open(file_in)
            data_1 = hdulist[0].data
            nx,ny = np.shape(data_1)
            if xc==None and yc==None:
                xc = nx/2.
                yc = ny/2.
            else:
                xc = float(xc)
                yc = float(yc)
            if maxsma==None:
                maxsma = math.sqrt( (nx/2.)**2 + (ny/2.)**2 )
            else:
                maxsma = float(maxsma)

            if False: #### CHANGE HERE!!!
                from pyraf import iraf
                iraf.stsdas()
                iraf.analysis()
                iraf.isophote()
                iraf.tables()
                iraf.ttools()
                iraf.ellipse(input=file_in,output='ell.tab',dqf=fits_mask,x0=xc,y0=yc,hcenter='no',recenter=0,sma0=5.0,minsma=minsma,maxsma=float(maxsma),pa=0.,hpa='no',ellip=0.05,hellip='no',interactive='no',step=step,linear='yes',verbose='no')
                iraf.tprint(table='ell.tab',pwidth='INDEF',plength='INDEF',showhdr='yes', Stdout=ell_file)

            else:                
                crea_ell(file_in, xc, yc, ellip=ellip, pa=pa, sma=sma0, step=step, minsma=minsma, maxsma=float(maxsma), ell_file=ell_file, fits_mask=fits_mask, fflag=fflag, olthresh=olthresh, linear=linear, hcenter=hcenter)
                os.chmod(r"ell.cl",0o777)
                subprocess.call("cl < ell.cl -o", shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)

            sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = loadtxt(ell_file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 5, dtype='str')

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

            f= figure(0,figsize=(4,10))

            gs = gridspec.GridSpec(3, 1)#,width_ratios=[1,3])
            gs.update(left=0.25, right=3.00, hspace=0.0)

            ax1 = plt.subplot(gs[0])
            ax2 = plt.subplot(gs[1])
            ax3 = plt.subplot(gs[2])
            overplot = False
            comp = False
            PA_ylim=0.;ell_ylim=0.;B4_ylim=0.;mag_ylim=0.
            PA_ylim,ell_ylim,B4_ylim,mag_ylim = plot_iso(ax1,ax2,ax3,'iraf.%s' % (outp_format),inten,inten_err,PA,errPA,ell,errell,B4,errB4,sma,float(maxsma),pix2sec,m0,outp_format,overplot,comp=comp,\
              PA_ylim=PA_ylim,ell_ylim=ell_ylim,B4_ylim=B4_ylim,mag_ylim=mag_ylim)
            
            fig = figure(1,figsize=(5,5))
            ax4 = fig.add_subplot(111)    
            plot_azim(ax4,sma,inten,inten_err,pix2sec,m0,float(maxsma),overplot,outp_format,comp)
        files = glob.glob('upar*.par')
        for file in files:
            os.remove(file)
        
        create_reg_with_ellipses(x0, y0, sma, ell, PA)
        os.remove('ell.cl')
        print('Done!')
        




def read_ell(ell_file,radius,step):
        sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = loadtxt(ell_file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 5, dtype='str')

        for k in range(len(sma)):
                if sma[k]=='INDEF': sma[k]=99999.
                if inten[k]=='INDEF': inten[k]=99999.
                if inten_err[k]=='INDEF': inten_err[k]=99999.
                if ell[k]=='INDEF': ell[k]=99999.
                if errell[k]=='INDEF': errell[k]=99999.
                if PA[k]=='INDEF': PA[k]=99999.
                if errPA[k]=='INDEF': errPA[k]=99999.
                if x0[k]=='INDEF': x0[k]=99999.
                if y0[k]=='INDEF': y0[k]=99999.
                if B4[k]=='INDEF': B4[k]=99999.
                if errB4[k]=='INDEF': errB4[k]=99999.
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



        BB4 = []
        for k in range(len(sma)):
                if B4[k]!=99999. and sma[k]<radius*2.:
                        BB4.append(B4[k])

        ellip = 0.
        return ellip,median(BB4)


def rmax_bar(rmax,ell_file='ellipse.txt'):
        sma,inten,inten_err,ell,errell,PA,errPA,x0,y0,B4,errB4 = loadtxt(ell_file, usecols=[1,2,3,6,7,8,9,10,12,33,34], unpack=True, skiprows = 5, dtype='str')

        for k in range(len(sma)):
                if sma[k]=='INDEF': sma[k]=99999.
                if inten[k]=='INDEF': inten[k]=99999.
                if inten_err[k]=='INDEF': inten_err[k]=99999.
                if ell[k]=='INDEF': ell[k]=99999.
                if errell[k]=='INDEF': errell[k]=99999.
                if PA[k]=='INDEF': PA[k]=99999.
                if errPA[k]=='INDEF': errPA[k]=99999.
                if x0[k]=='INDEF': x0[k]=99999.
                if y0[k]=='INDEF': y0[k]=99999.
                if B4[k]=='INDEF': B4[k]=99999.
                if errB4[k]=='INDEF': errB4[k]=99999.

        ELL = []
        PAA = []
        R = []

        for k in range(len(sma)):
                if float(ell[k])!=99999. and float(PA[k])!=99999. and float(sma[k])<rmax/3. and float(sma[k])>4.:
                        ELL.append(float(ell[k]))
                        PAA.append(float(PA[k]))
                        R.append(float(sma[k]))


        rmaxBAR = min([R[ELL.index(max(ELL))],R[PAA.index(min(PAA))]])
        return float(rmaxBAR)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="IRAF/ELLIPSE fitting")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("--outp_format", nargs='?', const=1, help="Optional: Input the name of the format for the output picture ", type=str, default='eps') 
    parser.add_argument("--mask", nargs='?', const=1, help="Optional: Input the name of the mask file",type=str,default='')

    parser.add_argument("--xc", nargs='?', const=1, help="Optional: Input the x-coordinate of the centre",type=str,default=None) 
    parser.add_argument("--yc", nargs='?', const=1, help="Optional: Input the y-coordinate of the centre",type=str,default=None) 

    parser.add_argument("--ellip", nargs='?', const=1, help="Optional: Input initial ellipticity",type=float, default=0.2)
    parser.add_argument("--pa", nargs='?', const=1, help="Optional: Input initial PA",type=float, default=20.)
    parser.add_argument("--sma", nargs='?', const=1, help="Optional: Input initial sma",type=float, default=10.)
    
    parser.add_argument("--maxsma", nargs='?', const=1, help="Optional: Input the maximum radius",type=float, default=None)
    parser.add_argument("--minsma", nargs='?', const=1, help="Optional: Input the minimum radius",type=float, default=1.)
    parser.add_argument("--step", nargs='?', const=1, help="Optional: Input the step",type=float, default=0.03)    
    parser.add_argument("--ZP", nargs='?', const=1, help="Optional: Input the y-coordinate of the centre",type=float,default=28.) 
    parser.add_argument("--pix2sec", nargs='?', const=1, help="Optional: Input the y-coordinate of the centre",type=float,default=1.) 
    parser.add_argument("--fflag", nargs='?', const=1, help="Optional: Acceptable fraction of flagged data points in intensity sample.", type=float, default=1.0) 
    parser.add_argument("--olthresh", nargs='?', const=1, help="Optional: Threshold for the object locator algorithm. See http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?controlpar.hlp", type=float, default=0.0)
    
    
  
    parser.add_argument("--linear", action="store_true", default=False,
                        help="Linear geometric step? Default False.")  


    args = parser.parse_args()

    input_image = args.inputImage
    outp_format = args.outp_format
    xc = args.xc
    yc = args.yc
    ellip = args.ellip
    pa = args.pa
    sma = args.sma
    maxsma = args.maxsma
    minsma = args.minsma
    ZP = args.ZP
    pix2sec = args.pix2sec
    mask = args.mask
    step = args.step
    fflag = args.fflag
    olthresh = args.olthresh
    linear = args.linear
    if linear:
        linear='yes'
    else:
        linear='no'

    main_ell(input_image, xc, yc, ellip=ellip, pa=pa, sma0=sma, m0=ZP, pix2sec=pix2sec, step=step, minsma=minsma, maxsma=maxsma, outp_format=outp_format, ell_file='ellipse.txt', fits_mask=mask, fflag=fflag, olthresh=olthresh, linear=linear)
