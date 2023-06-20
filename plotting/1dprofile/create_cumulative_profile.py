#! /usr/bin/env python
# Function to plot profiles for ordinary- and cube-fits images, the last one will be interpreted as the galfit file with the model.
# You can also specify several files which will be plotted together
# Examples:
#  python3 /home/amosenko/MyGit/IMAN/plotting/1dprofile/create_cumulative_profile.py composed_model.fits 20.0 1.0 --mask_image mask.fits  --profile summed --zmax 80 --Rmax 160 --method mean --SBmin 30 --SBmax 40
#  python3 /home/amosenko/MyGit/IMAN/plotting/1dprofile/create_cumulative_profile.py composed_model.fits 20.0 1.0 --mask_image mask.fits  --profile vert --zmax 80 --Rmax 160 --method mean

import os
import sys

from pylab import *
import astropy.io.fits as pyfits
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import glob
from matplotlib.pyplot import cm
import re
from scipy.odr.odrpack import *
import argparse
from scipy.ndimage import rotate
from astropy.stats import sigma_clipped_stats
import warnings
import radial_profile

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True



warnings.filterwarnings("ignore")

COLORS = ['b','g','m','y','c','lime','orange']
line_styles = ['-.','--',':','-.','-.',':','--']

axes_fontsize = 20
axes_label_fontsize = 0.8*axes_fontsize

def main(input_image,m0,pix2sec,mask_image=None,profile = 'azim',xc=0.,yc=0.,PA=0.,Rmin=0.,Rmax=0.,step=1.,zmin=0.,zmax=0.,output_file=None,AX=None, geom_units='arcsec',SB_units='mag/arcsec2',Scale=0.1, legend_size=6, interp=False, FWHM=3.,max_SB=None,min_SB=None, do_not_show_full_model=False, plot_symbs='o',text=None, comps=None, method='mean', start_from_zero=True,add_legend=True):

      hdu = pyfits.open(input_image)
      number_of_layers = len(hdu)    


      if mask_image is not None:
        # Read in the mask file
        hdu_mask = pyfits.open(mask_image)
        mask = hdu_mask[0].data
      else:
        mask = np.zeros_like(np.array(hdu[0].data,dtype=float),dtype=bool)
      
      
      
      mask = np.ma.make_mask(mask, shrink=False)

      if AX==None:
          save_out_ima = True

          fig = plt.figure(0,figsize=(5, 5))
          ax = fig.add_subplot(111)
          if profile == 'vert':
            plt.xlabel(r' $z$ (%s) ' % (geom_units), fontsize=15)
          else:
            plt.xlabel(r' $r$ (%s) ' % (geom_units), fontsize=15)
          if SB_units=='mag/arcsec2':
            plt.ylabel(r' $\mu$ (mag arcsec$^{-2}$) ', fontsize=15)
          elif SB_units=='ADU/pix2':
            plt.ylabel(r' Intensity (ADU) ', fontsize=15)
          plt.gca().invert_yaxis()
      else:
          save_out_ima = False
          if profile == 'vert':
            AX.set_xlabel(r' $z$ (%s) ' % (geom_units), fontsize=15)
          else:
            AX.set_xlabel(r' $r$ (%s) ' % (geom_units), fontsize=15)  
          if SB_units=='mag/arcsec2':
            AX.set_ylabel(r' $\mu$ (mag arcsec$^{-2}$) ', fontsize=15)
          elif SB_units=='ADU/pix2':
            AX.set_ylabel(r' Intensity (ADU) ', fontsize=15)
          AX.invert_yaxis()        ####???


      layers = range(number_of_layers)
      
            
      RADIUS = []
      MAG = []
      for ik in range(len(layers)):
            if len(layers)==5 and ik==1:
              continue
            layer = layers[ik]
            #print(layer)
            prihdr = hdu[layer].header
            try:
              Label = prihdr['NAME_OF_LAYER']
              if number_of_layers==5 and do_not_show_full_model==False and layer==1:
                 Label = hdu[4].header['NAME_OF_LAYER'] 
                  
            except:
              Label = 'galaxy'
              
            if '_' in Label:
                Label = Label.replace('_', ' ')
            
            if comps is not None:
                Label = comps[ik]
            
            data = hdu[layer].data
            mask[np.isnan(data)]=True

                
            
            if PA!=0.:
                dataa = np.copy(data)
                dataa[np.isnan(data)]=0.
                data = rotate(dataa, PA, order=5)
                mask = rotate(mask, PA, order=0)

            ySize, xSize = data.shape

            if xc==0. and yc==0.:
                  xc = xSize/2.
                  yc = ySize/2.
                 
                 
            if zmax==0.:
                  zmax = min([ySize-yc,yc])
                
            if Rmax==0.:
                  Rmax = max([xSize-xc,xc])
    
            
            if profile=='summed':
                I = []
                r = []

                for x in range(0,xSize,1):
                    if fabs(x-xc)<=Rmax:  
                      II = []
                      mmask = []
                      yranges = list(range(int(yc-zmax),int(yc-zmin),1)) + list(range(int(yc+zmin),int(yc+zmax)))
                      for y in yranges:
                        II.append(data[y,x])
                        mmask.append(mask[y,x])


                      if method=='sum': # Do not use it!
                            data[np.isnan(data)]=0.
                            I.append(sum(II))
                      if method=='sigma_clipping':
                            mean, median, std = sigma_clipped_stats(II, sigma=5.0, maxiters=5, mask=mmask)
                            I.append(median)
                      if method=='median':  
                            median = np.ma.median(np.ma.array(np.array(II), mask=np.array(mmask)))
                            I.append(median)
                      if method=='mean': # Gives better results!
                            mean = np.ma.mean(np.ma.array(np.array(II), mask=np.array(mmask)))
                            I.append(mean)
                      r.append(x-xc)
                r = np.array(r)
                I = np.array(I)
  
 
            elif profile=='vert':
                I = []
                I_err = []
                r = []
                
                if zmax==0.:
                    yranges = range(0,ySize,1)
                else:
                    if zmax>yc or zmax>ySize-yc:
                        Rmax = max([ySize-yc,yc])
                    yranges = list(range(int(yc-zmax),int(yc-zmin),1))+list(range(int(yc+zmin),int(yc+zmax),1))

                for y in yranges:
                  II = []
                  mmask = []
                  xranges = list(range(int(xc-Rmax),int(xc-Rmin),1))+list(range(int(xc+Rmin),int(xc+Rmax),1))
                  for x in xranges:
                    II.append(data[y,x])
                    mmask.append(mask[y,x])
                  if method=='sum': # Do not use it!
                    data[np.isnan(data)]=0.
                    I.append(sum(II))
                  if method=='sigma_clipping':
                    mean, median, std = sigma_clipped_stats(II, sigma=5.0, maxiters=5, mask=mmask)
                    I.append(median)
                  if method=='median':  
                    median = np.ma.median(np.ma.array(np.array(II), mask=np.array(mmask)))
                    I.append(median)
                  if method=='mean': # Gives better results!
                    mean = np.ma.mean(np.ma.array(np.array(II), mask=np.array(mmask)))
                    I.append(mean)
                  r.append(y-yc)
                r = np.array(r)
                I = np.array(I)


              
            if layer==0:
              color1 = 'gray'
              color2 = 'white'
              mecolor = 'gray'
              msize=5
            if layer==1:
              color1 = 'red'
              color2 = 'red'
              mecolor = 'red'
              msize=3
              lyne_style = '-'
            if layer==2 or layer==3:
              color1 = 'c'
              color2 = 'c'
              mecolor = 'c'
              msize=3
              lyne_style = '.-'
            if layer>3:
              color1 = COLORS[layer-4]
              color2 = COLORS[layer-4]
              mecolor = COLORS[layer-4]
              lyne_style = line_styles[layer-4]
              msize=3




            # ********************** Plotting ********************** 
            if profile=='summed':
                if SB_units=='mag/arcsec2':
                  mag = m0 - 2.5*log10(I)+ 5.*log10(pix2sec)
                  min_mag = np.nanmin(mag[~np.isnan(mag)])
                  max_mag = np.nanmax(mag[~np.isnan(mag) & ~np.isinf(mag)])
                  if max_SB!=None:
                        max_mag = min_SB
                        min_mag = max_SB
                elif SB_units=='ADU/pix2':
                  mag = I
                  max_mag = np.nanmin(mag[~np.isnan(mag)])
                  min_mag = np.nanmax(mag[~np.isnan(mag) & ~np.isinf(mag)])  
                  if max_SB!=None:
                        max_mag = max_SB
                        min_mag = min_SB


                if geom_units=='pix':
                  Radius = r
                elif geom_units=='arcsec':
                  Radius = r*pix2sec
                elif geom_units=='kpc':
                  Radius = r*pix2sec*Scale

                if layer==0:
                  if AX==None:
                    plt.plot(Radius, mag,'o',color='gray',markeredgecolor=mecolor,markersize=msize,label = Label)
                    ylim(max_mag,min_mag-0.5)
                  else:
                    AX.plot(Radius, mag,'o',color='gray',markeredgecolor=mecolor,markersize=msize,label = Label)
                    AX.set_ylim(max_mag,min_mag-0.5)  
                
                if layer!=0 and layer!=2 and layer!=3:
                  if AX==None:
                    plt.plot(Radius, mag,lyne_style,color=color2,lw=2,label = Label)
                  else:
                    AX.plot(Radius, mag,lyne_style,color=color2,lw=2,label = Label)
                if add_legend:
                    if AX==None:
                        if start_from_zero:
                            plt.legend(loc=3, borderaxespad=0.,fontsize=legend_size+2,numpoints=1)
                        else:
                            plt.legend(loc=2, borderaxespad=0.,fontsize=legend_size+2,numpoints=1)
                    else:
                        if start_from_zero:
                            AX.legend(loc=3, borderaxespad=0.,fontsize=legend_size+2,numpoints=1)
                        else:
                            AX.legend(loc=2, borderaxespad=0.,fontsize=legend_size+2,numpoints=1)



                if AX==None:
                    if start_from_zero:
                            xlim(0,max(Radius))
                    else:
                            xlim(-max(Radius),max(Radius))
                else:
                    if start_from_zero:
                        AX.set_xlim(0.,max(Radius))
                    else:
                        AX.set_xlim(-max(Radius),max(Radius))


            if profile=='vert':
                if SB_units=='mag/arcsec2':
                  mag = m0 - 2.5*log10(I)+ 5.*log10(pix2sec)
                  #mag_err = 2.5*I_err/(I*2.303)
                  min_mag = np.nanmin(mag[~np.isnan(mag)])
                  max_mag = np.nanmax(mag[~np.isnan(mag) & ~np.isinf(mag)])
                  if max_SB!=None:
                        max_mag = min_SB
                        min_mag = max_SB
                        
                elif SB_units=='ADU/pix2':
                  mag = I
                  #mag_err = I_err
                  max_mag = np.nanmin(mag[~np.isnan(mag)])
                  min_mag = np.nanmax(mag[~np.isnan(mag) & ~np.isinf(mag)])  
                
                  if max_SB!=None:
                        max_mag = max_SB
                        min_mag = min_SB
                
                if geom_units=='pix':
                  Radius = r
                elif geom_units=='arcsec':
                  Radius = r*pix2sec
                elif geom_units=='kpc':
                  Radius = r*pix2sec*Scale

                if layer==0:
                  if AX==None:
                    #plt.plot(Radius, mag,'o',color='gray',markeredgecolor=mecolor,markersize=msize,label = Label)
                    plt.errorbar(Radius, mag, yerr=None, fmt='o',color='gray',markeredgecolor=mecolor,markersize=msize,label = Label, zorder=layer) # markersize='4', markeredgecolor='lime', markerfacecolor='lime', ecolor='lime', zorder=overlayed_number+2)
                    
                    ylim(max_mag,min_mag-0.5)
                  else:
                    #AX.plot(Radius, mag,'o',color='gray',markeredgecolor=mecolor,markersize=msize,label = Label)
                    AX.errorbar(Radius, mag, yerr=None, fmt='o',color='gray',markeredgecolor=mecolor,markersize=msize,label = Label, zorder=layer)
                    AX.set_ylim(max_mag,min_mag-0.5)  
                
                if layer!=0 and layer!=2 and layer!=3:
                  if AX==None:
                    plt.plot(Radius, mag,lyne_style,color=color1,lw=2,label = Label)
                    if start_from_zero:
                        xlim(0,max(Radius))
                    else:
                        xlim(-max(Radius),max(Radius))
                  else:
                    AX.plot(Radius, mag,lyne_style,color=color1,lw=2,label = Label)
                    if start_from_zero:
                        AX.set_xlim(0.,max(Radius))
                    else:
                        AX.set_xlim(-max(Radius),max(Radius))
                    
                if add_legend:
                    if AX==None:
                        if start_from_zero:
                            plt.legend(loc=3, borderaxespad=0.,fontsize=legend_size+2,numpoints=1)
                        else:
                            plt.legend(loc=2, borderaxespad=0.,fontsize=legend_size+2,numpoints=1)
                    else:
                        if start_from_zero:
                            AX.legend(loc=3, borderaxespad=0.,fontsize=legend_size+2,numpoints=1)
                        else:
                            AX.legend(loc=2, borderaxespad=0.,fontsize=legend_size+2,numpoints=1)


            RADIUS.append(Radius)
            MAG.append(mag)
    
      if text is not None:
           ax.text(0.7, 0.9, text, fontsize=15, color='black',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline')
           
      if output_file==None and save_out_ima == True:
        

        if profile=='vert':
          output_file = input_image.split('.fits')[0] + '_prof_ver.eps'
        elif profile=='cut':
          output_file = input_image.split('.fits')[0] + '_prof_cut.eps'
        elif profile=='rad':
          output_file = input_image.split('.fits')[0] + '_prof_rad.eps'
        elif profile=='summed':
          output_file = input_image.split('.fits')[0] + '_prof_sum.eps'
        elif profile=='azim': 
          output_file = input_image.split('.fits')[0] + '_prof_azim.eps'



      ax.minorticks_on()
      ax.tick_params(direction='in', length=3, which='minor', top=True, right=True, labelsize=axes_label_fontsize)
      ax.tick_params(direction='in', length=6, which='major', top=True, right=True, labelsize=axes_label_fontsize)
      
      if save_out_ima == True:
        plt.savefig(output_file, transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.05)    
        plt.clf()
        plt.close()
        plt.close('all')
      return RADIUS,MAG
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Plotting profiles")
    parser.add_argument("input_image", help="Input fits image with the centered object")
    parser.add_argument("ZeroPoint", help="Input Zero Point in [mag/arcsec^2]") 
    parser.add_argument("Scale", help="Input scale in [arcsec/pix]")
    parser.add_argument("--mask_image", nargs='?', help="Input mask image (where I>0. is for masked pixels)",default=None)
    parser.add_argument("--profile", nargs='?', const=1, help="Input profile: summed for only summed horizontal profile, all - for both summed and cut, vert - for the vertical summed profile.",type=str,default='azim')
    parser.add_argument("--xc", nargs='?', const=1, help="Input x-coordinate of the center [pixels]",type=float,default=0.)
    parser.add_argument("--yc", nargs='?', const=1, help="Input y-coordinate of the center [pixels]",type=float,default=0.)
    parser.add_argument("--PA", nargs='?', const=1, help="Input position angle of the major axis [deg]; counterclockwise: up = 90 deg, right=0.",type=float,default=0.)     
    parser.add_argument("--Rmin", nargs='?', const=1, help="Input minimal radius of the disc [pixels]",type=float,default=0.) 
    parser.add_argument("--Rmax", nargs='?', const=1, help="Input maximal radius of the disc [pixels]",type=float,default=0.)
    parser.add_argument("--SBmin", nargs='?', const=1, help="Input minimal SB",type=float,default=None) 
    parser.add_argument("--SBmax", nargs='?', const=1, help="Input maximal SB",type=float,default=None)
    parser.add_argument("--zmin", nargs='?', const=1, help="Input minimal z-extension of the disc [pixels]",type=float,default=0.) 
    parser.add_argument("--zmax", nargs='?', const=1, help="Input maximal z-extension of the disc [pixels]",type=float,default=0.)    
    parser.add_argument("--output_file", nargs='?', help="Output file (with .png, .eps or .pdf extension)",default=None)
    parser.add_argument("--geom_units", nargs='?', const=1, help="Input units for X-axis [pixels, arcsec,kpc]",type=str,default='arcsec')    
    parser.add_argument("--SB_units", nargs='?', const=1, help="Input units for Y-axis [ADU/pixel2, mag/arcsec2]",type=str,default='mag/arcsec2')
    parser.add_argument("--Scale", nargs='?', help="Input Scale in [kpc/arcsec]",type=float,default=0.1)
    parser.add_argument("--legend_fsize", nargs='?', help="Input fontsize for legend",type=float,default=9)
    parser.add_argument("--comps", nargs='?', const=1, help="Optional: Names of components (from composed_model.fits) ", type=str, default=None) 
    parser.add_argument("--method", nargs='?', const=1, help="Optional: Method to use: sum, mean, median, sigma_clipping", type=str, default='mean') 

    args = parser.parse_args()

    input_image = args.input_image
    m0 = float(args.ZeroPoint)
    pix2sec = float(args.Scale)
    mask_image = args.mask_image
    profile = str(args.profile)
    xc = float(args.xc)
    yc = float(args.yc)
    PA = float(args.PA)
    Rmin = float(args.Rmin)
    Rmax = float(args.Rmax)
    zmin = float(args.zmin)
    zmax = float(args.zmax)    
    output_file = args.output_file
    geom_units = args.geom_units
    SB_units = args.SB_units
    Scale = float(args.Scale)
    legend_fsize = float(args.legend_fsize)

    SBmin = args.SBmin
    SBmax = args.SBmax

    if args.comps is not None:
        comps = args.comps.split(',')
    else:
        comps = args.comps
    
    method = args.method
    
    main(input_image,m0,pix2sec,mask_image=mask_image,profile=profile,xc=xc,yc=yc,PA=PA,Rmin=Rmin,Rmax=Rmax,zmin=zmin,zmax=zmax,output_file=output_file,AX=None,geom_units=geom_units,SB_units=SB_units,Scale=Scale,legend_size=legend_fsize,min_SB=SBmin,max_SB=SBmax, comps=comps, method=method)
