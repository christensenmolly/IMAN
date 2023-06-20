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
from scipy.interpolate import interp1d
warnings.filterwarnings("ignore")


LOCAL_DIR = "/iraf_fitting"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'plotting/1dprofile'))

import grid_plots

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams["mathtext.fontset"] = "cm"

tmp_out = sys.stdout

axes_fontsize = 20
axes_label_fontsize = 0.8*axes_fontsize
#fsize=30
FNULL = open(os.devnull, 'w')


def read_ell(ell_file):
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
        
        return sma, inten, inten_err, ell, errell, PA, errPA, x0, y0, B4, errB4


def set_line_props(number):
    # The first symbol is the data, the second is the total model,
    # the rest are individual components.
    COLORS = ['gray','red','blue','green','magenta','orange','c','lime','orange']
    marker_line_styles = ['o','-','-.','--',':','-.','-.',':','--']
    marker_line_size = [5,3,2,2,2,2,2,2,2] 
    return COLORS[number],marker_line_styles[number],marker_line_size[number]
    


def plot_SB_profile(ax, sma, inten, inten_err, m0, pix2sec, R_units='arcsec', SB_units='mag/arcsec2', R_scale='linear', show_error_bars=True, label='', number=0, SB_lim=None, R_lim=None, text=None, final_plot=False):
    
         color, symb, size = set_line_props(number)
    
         if R_units=='arcsec':
             sma = sma*pix2sec
         
         if SB_units=='mag/arcsec2':
            inten_err = fabs((2.5/log(10.0)) * inten_err/inten)
            inten = m0 - 2.5*log10(inten)+ 5.*log10(pix2sec)
            

         if R_scale=='1/4':
                sma = sma**(0.25)
            
         if R_scale=='log':
                sma = np.log10(sma)
     
         
         if show_error_bars and number==0:
                ax.errorbar(sma, inten, yerr=inten_err, fmt=symb, color=color, ecolor=color, markersize=size, label=label, zorder=number)
         else:
                if number==0:
                    ax.plot(sma, inten, symb, color=color, markeredgecolor=color, markersize=size, zorder=number, label=label)
                else:
                    ax.plot(sma, inten, ls=symb, color=color, markeredgecolor=color, lw=size, zorder=number, label=label) 
 
 
         if number==0:
            # Define axes ranges:
            if SB_lim is not None:
                [ymin,ymax] = SB_lim
                ymin = float(ymin)
                ymax = float(ymax)
                if ymin>ymax: ymin,ymax = ymax,ymin
                if SB_units=='mag/arcsec2':
                    ax.set_ylim(ymax, ymin)
                else:
                    ax.set_ylim(ymin, ymax)
            else:
                mean_inten, median_inten, std_inten = sigma_clipped_stats(inten, sigma=3.0)
                ymax = np.nanmax(inten[inten != np.inf]) + std_inten
                ymin = np.nanmin(inten[inten != np.inf]) - std_inten
                if SB_units=='mag/arcsec2':
                    ax.set_ylim(ymax, ymin)
                else:
                    ax.set_ylim(ymin, ymax)

            if R_lim is not None:
                [xmin,xmax] = R_lim
                xmin = float(xmin)
                xmax = float(xmax)
                if xmin>xmax: xmin,xmax = xmax,xmin
                ax.set_xlim(xmin, xmax)
            else:
                ax.set_xlim(np.nanmin(sma), np.nanmax(sma))
                R_lim = [np.nanmin(sma),np.nanmax(sma)]
            
            # Define axes labels:
            if R_scale=='linear':
                ax.set_xlabel(r' $r$ (%s) ' % (R_units), fontsize=axes_fontsize)
            elif R_scale=='1/4':
                ax.set_xlabel(r' $r^{1/4}$ (%s$^{1/4}$) ' % (R_units), fontsize=axes_fontsize)
            elif R_scale=='log':
                ax.set_xlabel(r' $\log\,r$ (%s) ' % (R_units), fontsize=axes_fontsize)
            
            if SB_units=='mag/arcsec2':
                ax.set_ylabel(r' $\mu$ (mag arcsec$^{-2}$) ', fontsize=axes_fontsize)
            elif SB_units=='ADU/pix2':
                ax.set_ylabel(r' Intensity (ADU) ', fontsize=axes_fontsize)
                
            ax.minorticks_on()
            ax.tick_params(direction='in', length=3, which='minor', top=True, right=True, labelsize=axes_label_fontsize)
            ax.tick_params(direction='in', length=6, which='major', top=True, right=True, labelsize=axes_label_fontsize)
            
            
         if final_plot:
            if text is not None:
                ax.text(0.3, 0.9, text, fontsize=axes_fontsize*1.5, color='black',transform=ax.transAxes, horizontalalignment='center',verticalalignment='baseline')
             
            ax.legend(loc=0, borderaxespad=0.5, fontsize=10, numpoints=1, markerscale=1)

         return sma, inten, R_lim  

def plot_SB_residual(ax, sma_data, mag_data, sma_model, mag_model, R_scale='linear', R_units='arcsec', R_lim=None):
        f_data = interp1d(sma_data, mag_data)
        f_model = interp1d(sma_model, mag_model)
        
        if np.max(sma_data)>np.max(sma_model):
            sma = sma_model
        else:
            sma = sma_data
            
        delta = f_data(sma) - f_model(sma)
        
        rms = sigma_clipped_stats(delta, sigma=3.0)[2]
        
        ax.plot(sma, delta, 'o', color='gray', markeredgecolor='gray', markersize=5, zorder=0)
        ax.axhline(y=0)

        # Define axes labels:
        if R_scale=='linear':
                ax.set_xlabel(r' $r$ (%s) ' % (R_units), fontsize=axes_fontsize)
        elif R_scale=='1/4':
                ax.set_xlabel(r' $r^{1/4}$ (%s$^{1/4}$) ' % (R_units), fontsize=axes_fontsize)
        elif R_scale=='log':
                ax.set_xlabel(r' $\log\,r$ (%s) ' % (R_units), fontsize=axes_fontsize)        
        
        ax.set_ylabel(r'$\Delta \mu$', fontsize=axes_fontsize)

        ax.minorticks_on()
        ax.tick_params(direction='in', length=3, which='minor', top=True, right=True, labelsize=axes_label_fontsize)
        ax.tick_params(direction='in', length=6, which='major', top=True, right=True, labelsize=axes_label_fontsize)

        [xmin,xmax] = R_lim
        if xmin>xmax: xmin,xmax = xmax,xmin
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(-5*rms,5*rms)
     

def plot_ell_profile(ax, sma, ell, ell_err, pix2sec, R_units='arcsec', R_scale='linear', show_error_bars=True, label='', number=0, R_lim=None, ell_lim=None):
 
         color, symb, size = set_line_props(number)
    
         if R_units=='arcsec':
             sma = sma*pix2sec

         if ell_lim is not None:
            [min_ell,max_ell] = ell_lim
            min_ell = float(min_ell)
            max_ell = float(max_ell)
         else:
            ell_mean,ell_median,ell_rms = sigma_clipped_stats(ell, sigma=3.0)
            #min_ell = np.min(fabs(ell)-errell)*0.99
            #max_ell = np.max(fabs(ell)+errell)*1.01 
            min_ell = max(0.001, ell_median - 5.*ell_rms)
            max_ell = min(0.999, ell_median + 5.*ell_rms)

         if R_scale=='1/4':
                sma = sma**(0.25)
            
         if R_scale=='log':
                sma = np.log10(sma) 

         if show_error_bars and number==0:
                ax.errorbar(sma, ell, yerr=ell_err, fmt=symb, color=color, ecolor=color, markersize=size, label=label, zorder=number)
         else:
                if number==0:
                    ax.plot(sma, ell, symb, color=color, markeredgecolor=color, markersize=size, zorder=number, label=label)
                else:
                    ax.plot(sma, ell, ls=symb, color=color, markeredgecolor=color, lw=size, zorder=number, label=label) 


         # Define axes labels:
         if R_scale=='linear':
                ax.set_xlabel(r' $r$ (%s) ' % (R_units), fontsize=axes_fontsize)
         elif R_scale=='1/4':
                ax.set_xlabel(r' $r^{1/4}$ (%s$^{1/4}$) ' % (R_units), fontsize=axes_fontsize)
         elif R_scale=='log':
                ax.set_xlabel(r' $\log\,r$ (%s) ' % (R_units), fontsize=axes_fontsize)        
        
         ax.set_ylabel(r'$\epsilon$', fontsize=axes_fontsize)
        
         ax.minorticks_on()
         ax.tick_params(direction='in', length=3, which='minor', top=True, right=True, labelsize=axes_label_fontsize)
         ax.tick_params(direction='in', length=6, which='major', top=True, right=True, labelsize=axes_label_fontsize)

         [xmin,xmax] = R_lim
         if xmin>xmax: xmin,xmax = xmax,xmin
         ax.set_xlim(xmin, xmax)
         
         if number==0:
            ax.set_ylim(float(min_ell), float(max_ell))


def plot_PA_profile(ax, sma, PA, PA_err, pix2sec, R_units='arcsec', R_scale='linear', show_error_bars=True, label='', number=0, R_lim=None, PA_lim=None, show_neg_PA=False):
         if not show_neg_PA:
            PA[PA<0.] = PA[PA<0.] + 180.
         color, symb, size = set_line_props(number)
    
         if R_units=='arcsec':
             sma = sma*pix2sec

         if PA_lim is not None:
            [min_PA,max_PA] = PA_lim
            min_PA = float(min_PA)
            max_PA = float(max_PA)
         else:
            PA_mean,PA_median,PA_rms = sigma_clipped_stats(PA, sigma=3.0)
            #min_PA = np.min(fabs(PA)-PA_err)*0.99
            #max_PA = np.max(fabs(PA)+PA_err)*1.01
            min_PA = PA_median - 5.*PA_rms
            max_PA = PA_median + 5.*PA_rms

         if R_scale=='1/4':
                sma = sma**(0.25)
            
         if R_scale=='log':
                sma = np.log10(sma) 

         if show_error_bars and number==0:
                ax.errorbar(sma, PA, yerr=PA_err, fmt=symb, color=color, ecolor=color, markersize=size, label=label, zorder=number)
         else:
                if number==0:
                    ax.plot(sma, PA, symb, color=color, markeredgecolor=color, markersize=size, zorder=number, label=label)
                else:
                    ax.plot(sma, PA, ls=symb, color=color, markeredgecolor=color, lw=size, zorder=number, label=label) 


         # Define axes labels:
         if R_scale=='linear':
                ax.set_xlabel(r' $r$ (%s) ' % (R_units), fontsize=axes_fontsize)
         elif R_scale=='1/4':
                ax.set_xlabel(r' $r^{1/4}$ (%s$^{1/4}$) ' % (R_units), fontsize=axes_fontsize)
         elif R_scale=='log':
                ax.set_xlabel(r' $\log\,r$ (%s) ' % (R_units), fontsize=axes_fontsize)        
        
         ax.set_ylabel(r'PA (deg)', fontsize=axes_fontsize)
        
         ax.minorticks_on()
         ax.tick_params(direction='in', length=3, which='minor', top=True, right=True, labelsize=axes_label_fontsize)
         ax.tick_params(direction='in', length=6, which='major', top=True, right=True, labelsize=axes_label_fontsize)

         [xmin,xmax] = R_lim
         if xmin>xmax: xmin,xmax = xmax,xmin
         ax.set_xlim(xmin, xmax)
         
         if number==0:
            ax.set_ylim(float(min_PA), float(max_PA))
        

def plot_B4_profile(ax, sma, B4, B4_err, pix2sec, R_units='arcsec', R_scale='linear', show_error_bars=True, label='', number=0, R_lim=None, B4_lim=None):
 
         color, symb, size = set_line_props(number)
    
         if R_units=='arcsec':
             sma = sma*pix2sec

         if B4_lim is not None:
            [min_B4,max_B4] = B4_lim
            min_B4 = float(min_B4)
            max_B4 = float(max_B4)
         else:
            mean_B4, median_B4, std_B4 = sigma_clipped_stats(B4, sigma=3.0)
            if fabs(median_B4-5*std_B4)<2. and fabs(median_B4+5*std_B4)<2.:
              min_B4 = median_B4 - 5*std_B4
              max_B4 = median_B4 + 5*std_B4
            else:
              min_B4 = -2
              max_B4 = 2

         if R_scale=='1/4':
                sma = sma**(0.25)
            
         if R_scale=='log':
                sma = np.log10(sma) 

         if show_error_bars and number==0:
                ax.errorbar(sma, B4, yerr=B4_err, fmt=symb, color=color, ecolor=color, markersize=size, label=label, zorder=number)
         else:
                if number==0:
                    ax.plot(sma, B4, symb, color=color, markeredgecolor=color, markersize=size, zorder=number, label=label)
                else:
                    ax.plot(sma, B4, ls=symb, color=color, markeredgecolor=color, lw=size, zorder=number, label=label) 


         # Define axes labels:
         if R_scale=='linear':
                ax.set_xlabel(r' $r$ (%s) ' % (R_units), fontsize=axes_fontsize)
         elif R_scale=='1/4':
                ax.set_xlabel(r' $r^{1/4}$ (%s$^{1/4}$) ' % (R_units), fontsize=axes_fontsize)
         elif R_scale=='log':
                ax.set_xlabel(r' $\log\,r$ (%s) ' % (R_units), fontsize=axes_fontsize)        
         
         ax.axhline(y=0., ls='--')
         ax.set_ylabel(r'B$_4$', fontsize=axes_fontsize)
        
         ax.minorticks_on()
         ax.tick_params(direction='in', length=3, which='minor', top=True, right=True, labelsize=axes_label_fontsize)
         ax.tick_params(direction='in', length=6, which='major', top=True, right=True, labelsize=axes_label_fontsize)

         [xmin,xmax] = R_lim
         if xmin>xmax: xmin,xmax = xmax,xmin
         ax.set_xlim(xmin, xmax)
         
         if number==0:
            ax.set_ylim(float(min_B4), float(max_B4))
        

def main(ellipse_files, m0, pix2sec,
         names_of_layers=None,
         PA_lim=None,
         ell_lim=None,
         B4_lim=None,
         R_lim=None,
         SB_lim=None,
         R_units='arcsec',
         R_scale='linear',
         SB_units='mag/arcsec2',
         show_error_bars=True,
         text=None,
         verbosity=True,
         output_format='png',
         PA_offset=0.):
        

        
        # initialize plots
        f_SB, ax_SB = grid_plots.crea_one_plot(fig_plot=1)
        f_SBres, ax_SBres_SB, ax_SBres_res = grid_plots.crea_two_row_grids(fig_plot=2, xsize=2, ysize=6.7, left=0.25, right=3.0, hspace=0.0)
        f_SBell, ax_SBell_ell, ax_SBell_SB = grid_plots.crea_two_row_grids(fig_plot=3, xsize=2, ysize=6.7, left=0.25, right=3.0, hspace=0.0, ratios=[1,3])
        f_iso, ax_iso_PA, ax_iso_ell, ax_iso_B4 = grid_plots.crea_three_row_grids(fig_plot=4, xsize=2, ysize=5, left=0.25, right=3.0, hspace=0.0)






        # Plotting results
        if verbosity:
            print('\n\nPlotting results:')        

        
        
        
        for k in range(len(ellipse_files)):
            ellipse_file = ellipse_files[k]
            if os.stat(ellipse_file).st_size == 0:
                return 1
            
            if verbosity:
                print('\t File %s ...' % (ellipse_files[k]))        
            
            try:
                sma, inten, inten_err, ell, errell, PA, errPA, x0, y0, B4, errB4 = read_ell(ellipse_file) # TODO
                PA = PA + PA_offset
            except:
                sma, inten, inten_err = np.loadtxt(ellipse_file, usecols=[0,1,2], dtype=float, unpack=True, skiprows = 1, delimiter='\t')
                ell = np.ones(len(sma))
                errell = np.ones(len(sma))
                PA = np.ones(len(sma))
                errPA = np.ones(len(sma))
                x0 = np.ones(len(sma))
                y0 = np.ones(len(sma))
                B4 = np.ones(len(sma))
                errB4 = np.ones(len(sma))
                
                
            
            if k==len(ellipse_files)-1:
                final_plot = True
            else:
                final_plot = False
            
            if len(ellipse_files)==1:
                label = ''
                final_plot = True
            else:
                if names_of_layers is not None:
                    label = names_of_layers[k]
                else:
                    if k==0:
                        label='data'
                    else:
                        label = ellipse_file.split('_')[-1].split('.txt')[0]
            
            SMA, INTEN, R_lim = plot_SB_profile(ax_SB, sma, inten, inten_err, m0, pix2sec, R_units=R_units, SB_units=SB_units, R_scale=R_scale, show_error_bars=show_error_bars, label=label, number=k, SB_lim=SB_lim, R_lim=R_lim, text=text, final_plot=final_plot)
            
            if len(ellipse_files)!=1:
                SMA, INTEN, R_lim = plot_SB_profile(ax_SBres_SB, sma, inten, inten_err, m0, pix2sec, R_units=R_units, SB_units=SB_units, R_scale=R_scale, show_error_bars=show_error_bars, label=label, number=k, SB_lim=SB_lim, R_lim=R_lim, text=text, final_plot=final_plot)

            SMA, INTEN, R_lim = plot_SB_profile(ax_SBell_SB, sma, inten, inten_err, m0, pix2sec, R_units=R_units, SB_units=SB_units, R_scale=R_scale, show_error_bars=show_error_bars, label=label, number=k, SB_lim=SB_lim, R_lim=R_lim, text=text, final_plot=final_plot)
            
            if k>=0:
                plot_ell_profile(ax_SBell_ell, sma, ell, errell, pix2sec, R_units=R_units, R_scale=R_scale, show_error_bars=show_error_bars, label=label, number=k, R_lim=R_lim, ell_lim=ell_lim)
                
                plot_PA_profile(ax_iso_PA, sma, PA, errPA, pix2sec, R_units=R_units, R_scale=R_scale, show_error_bars=show_error_bars, label=label, number=k, R_lim=R_lim, PA_lim=PA_lim)                

                plot_ell_profile(ax_iso_ell, sma, ell, errell, pix2sec, R_units=R_units, R_scale=R_scale, show_error_bars=show_error_bars, label=label, number=k, R_lim=R_lim, ell_lim=ell_lim)
                
                plot_B4_profile(ax_iso_B4, sma, B4, errB4, pix2sec, R_units=R_units, R_scale=R_scale, show_error_bars=show_error_bars, label=label, number=k, R_lim=R_lim, B4_lim=B4_lim)
                
            
            if k==0:
                inten_data = INTEN
                sma_data = SMA
            if k==1:
                inten_model = INTEN
                sma_model = SMA
       
        if len(ellipse_files)!=1:
            plot_SB_residual(ax_SBres_res, sma_data, inten_data, sma_model, inten_model, R_scale=R_scale, R_units=R_units, R_lim=R_lim)
            grid_plots.merge_two_row_grids(ax_SBres_SB, ax_SBres_res)
            
        grid_plots.merge_two_row_grids(ax_SBell_ell, ax_SBell_SB)
        grid_plots.merge_three_row_grids(ax_iso_PA,ax_iso_ell, ax_iso_B4)
        
        # Save plots:
        grid_plots.save_plot(f_SB, 'profile.%s' % (output_format), DPI=300)
        grid_plots.save_plot(f_SBres, 'profile_res.%s' % (output_format), DPI=300)
        grid_plots.save_plot(f_SBell, 'profile_ell.%s' % (output_format), DPI=300)
        grid_plots.save_plot(f_iso, 'profile_iso.%s' % (output_format), DPI=300)
            
            
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="IRAF/ELLIPSE plot results")
    parser.add_argument("inputFile", help="File with IRAF/ELLIPSE results (txt)")
    parser.add_argument("--ZP", nargs='?', const=1, help="Optional: Zero point",type=float,default=28.) 
    parser.add_argument("--pix2sec", nargs='?', const=1, help="Optional: Pixel scale (arcsec/pix)",type=float,default=1.)     
    
 
    parser.add_argument("--comps", nargs='?', const=1, help="Optional: Names of components (fro composed_model.fits) ", type=str, default=None) 
    parser.add_argument("--outp_format", nargs='?', const=1, help="Optional: Input the name of the format for the output pictures ", type=str, default='png') 
    parser.add_argument("--ell_lim", nargs='?', const=1, help="Optional: Ellipticity limits (separated by comma)",type=str, default=None)
    parser.add_argument("--PA_lim", nargs='?', const=1, help="Optional: PA limits (separated by comma)",type=str, default=None)
    parser.add_argument("--B4_lim", nargs='?', const=1, help="Optional: B4 limits (separated by comma)",type=str, default=None)
    parser.add_argument("--R_lim", nargs='?', const=1, help="Optional: Radius limits (separated by comma)",type=str, default=None)
    parser.add_argument("--SB_lim", nargs='?', const=1, help="Optional: Surface brightness limits (separated by comma)",type=str, default=None)  
    parser.add_argument("--R_units", nargs='?', const=1, help="Optional: Surface brightness limits (separated by comma)",type=str, default='arcsec')  
    parser.add_argument("--R_scale", nargs='?', const=1, help="Optional: Surface brightness limits (separated by comma)",type=str, default='linear')  
    parser.add_argument("--SB_units", nargs='?', const=1, help="Optional: Surface brightness limits (separated by comma)",type=str, default='mag/arcsec2') 
    parser.add_argument("--text", nargs='?', const=1, help="Optional: Surface brightness limits (separated by comma)",type=str, default=None)  
        
    parser.add_argument("--errors", action="store_true", default=True,
                        help="Show error bars (False).")  

    parser.add_argument("--PA_offset", nargs='?', const=1, help="Optional: PA offset to add (0.)",type=float,default=0.) 


    args = parser.parse_args()


    input_files = args.inputFile.split(',')    
    ZP = args.ZP
    pix2sec = args.pix2sec
    output_format = args.outp_format
    
    ell_lim = args.ell_lim
    PA_lim = args.PA_lim
    B4_lim = args.B4_lim
    R_lim = args.R_lim
    SB_lim = args.SB_lim
    
    if ell_lim is not None:
        ell_lim = ell_lim.split(',')
 
    if PA_lim is not None:
        PA_lim = PA_lim.split(',')

    if B4_lim is not None:
        B4_lim = B4_lim.split(',') 
 
    if R_lim is not None:
        R_lim = R_lim.split(',') 
 
    if SB_lim is not None:
        SB_lim = SB_lim.split(',') 
 
    
    R_units = args.R_units
    R_scale = args.R_scale
    SB_units = args.SB_units
    show_error_bars = args.errors
    text = args.text
    if args.comps is not None:
        comps = args.comps.split(',')
    else:
        comps = args.comps
    
    PA_offset = args.PA_offset
    
    #PA_lim = [-29.,65.] # SPRC185
    main(input_files, ZP, pix2sec,
         names_of_layers=comps,
         PA_lim=PA_lim,
         ell_lim=ell_lim,
         B4_lim=B4_lim,
         R_lim=R_lim,
         SB_lim=SB_lim,
         R_units=R_units,
         R_scale=R_scale,
         SB_units=SB_units,
         show_error_bars=show_error_bars,
         text=text,
         verbosity=True,
         output_format=output_format,
         PA_offset=PA_offset)


             

            
            
            
