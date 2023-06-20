#! /usr/bin/env python
# Function to plot profiles for ordinary- and cube-fits images, the last one will be interpreted as the galfit file with the model.
# You can also specify several files which will be plotted together
# EXAMPLE: python ~/CurrentWork/ImaPrep/IMAN/Plot/plot_profile.py composed_model.fits 27. 0.396 summed 177 173 --o png
# python ~/MEGA/MyPrograms/IMAN/Plot/plot_profile.py composed_model.fits 21.581 0.6 summed
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
#matplotlib.use('TkAgg')
import warnings
#PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
#PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0]
#sys.path.append(PATH_TO_PACKAGE+'/FindFluxes')

#import iraf_fit_ellipse
import radial_profile
#import ds9_contour

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

#matplotlib.use('agg',warn=False, force=True)   ####!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

warnings.filterwarnings("ignore")

COLORS = ['b','g','m','y','c','lime','orange']
line_styles = ['-.','--',':','-.','-.',':','--']

axes_fontsize = 20
axes_label_fontsize = 0.8*axes_fontsize

def get_slice(fitsName, xOrig, yOrig, posang, Rmax, layer=0):
    """
    Function gets slice in fits file along specified line.
    Parameters:
        fitsName -- name of fits file.
        xOrig, yOrig -- coordinates of reference point
        posang [degrees] -- position angle of line. posang=0 -> vertical line (North-South slice).
                            positive posang values are for counterclockwise rotation, i.e. slice
                            with posang=45 is from upper-left corner till buttom right
        Rmax -- maximum radial distance
        layer -- number of hdu in multilayered images
    """
    xOrig, yOrig = yOrig, xOrig
    hdu = pyfits.open(fitsName)
    data = hdu[layer].data
    xSize, ySize = data.shape
    if xOrig==0:
      xOrig = xSize/2.
    if yOrig==0:
      yOrig = ySize/2.      
    rArray = []
    iArray = []
    posang += 90
    # tan of 90 and 270 degrees is infinity, so we have to avoid this nombers
    if (not (89.5 <= posang <= 90.5)) and (not (269.5 <= posang <= 270.5)):
        m = -tan(radians(posang))
        xbegin = xOrig - yOrig/m
        xend = min((ySize-yOrig)/m + xOrig, xSize)
        if xend < xbegin:
            xend, xbegin = xbegin, xend
        if xend > xSize:
            xend = xSize
        if xbegin < 0:
            xbegin = 0.0
        for x in linspace(xbegin, xend, xSize):
            y = m * (x-xOrig)+yOrig
            if (y<0) or (y>ySize-2) or (x>xSize-2) or (x<0):
                continue
            fx, ix = modf(x)
            fy, iy = modf(y)
            ix = int(ix)
            iy = int(iy)
            I = ((1.0-fx)*(1.0-fy)*data[ix, iy] + fx*(1.0-fy)*data[ix+1, iy] +
                fy*(1.0-fx)*data[ix, iy+1] + fx*fy*data[ix+1, iy+1])
            rArray.append(copysign(hypot(x-xOrig, y-yOrig), xOrig-x))
            iArray.append(I)
    else: # if posang is near 90 or 270 degrees then x is constant
        for y in arange(0, ySize-1):
            fx, ix = modf(xOrig)
            fy, iy = modf(y)
            ix = int(ix)
            iy = int(iy)
            I = ((1.0-fx)*(1.0-fy)*data[ix, iy] + fx*(1.0-fy)*data[ix+1, iy] +
                fy*(1.0-fx)*data[ix, iy+1] + fx*fy*data[ix+1, iy+1])
            rArray.append(y-yOrig)
            iArray.append(I)
    rArray = np.array(rArray,float)
    iArray = np.array(iArray,float)
    if Rmax!=0.:
      r = []; I = []
      for k in range(len(rArray)):
        if fabs(rArray[k])<=Rmax:
          r.append(rArray[k])
          I.append(iArray[k])
      r = np.array(r,float)
      I = np.array(I,float)
      return r,I
    else:
      return rArray, iArray

def line(B, x):
    return (1.0857/B[0])*fabs(x-B[2]) + B[1]


def azim_aver_profile_iraf(input_image,xc,yc,maxsma,m0,pix2sec):
    iraf_fit_ellipse.main(input_image,xc,yc,m0,pix2sec,minsma=0.,maxsma=maxsma,step=1.)



def main(input_image,m0,pix2sec,mask_image=None,profile = 'azim',xc=0.,yc=0.,PA=0.,Rmin=0.,Rmax=0.,step=1.,zmin=0.,zmax=0.,output_file=None,AX=None, geom_units='arcsec',SB_units='mag/arcsec2',Scale=0.1, legend_size=6, interp=False, FWHM=3.,max_SB=None,min_SB=None, do_not_show_full_model=False, plot_symbs='o',text=None, comps=None):
      # As an output one will receive the radius and mag lists with observed and modelled profiles

      #azim_aver_profile_iraf(input_image,xc,yc,Rmax,m0,pix2sec)
      #return 'azim_aver.eps'
        
      # Find number of layers:
      hdu = pyfits.open(input_image)
      number_of_layers = len(hdu)    
      
      #number_of_layers = 5
      
      if number_of_layers>=5:
        model_to_fil = hdu[1].data
      else:
        model_to_fil = hdu[0].data

      if mask_image is not None:
        # Read in the mask file
        hdu_mask = pyfits.open(mask_image)
        mask = hdu_mask[0].data
      else:
        mask = np.zeros_like(np.array(model_to_fil,dtype=float),dtype=bool)


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


      if number_of_layers<5:
        layers = [0]
      elif number_of_layers==5:
        layers = [0,1]
        if do_not_show_full_model==True:
            layers = [0,4]
      elif number_of_layers>5:
        layers = [0]+list(range(4,number_of_layers))+[1]
        if do_not_show_full_model==True:
            layers = [0]+list(range(4,number_of_layers))
            
      RADIUS = []
      MAG = []
      for ik in range(len(layers)):
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
            '''
            ySize, xSize = data.shape
            for k in range(ySize):
                for i in range(xSize):
                    if data[k,i]<0:
                        mask[k,i]=1
            '''                    
            
            if profile=='summed' or profile=='vert':
                #matplotlib.use('TkAgg')
                data = rotate(data, PA, order=5)
                mask = rotate(mask, PA, order=0)
                model_to_fil = rotate(model_to_fil, PA, order=5)
                #im = plt.imshow(mask, cmap='hot')
                #plt.show()
                #exit()
                #plt.savefig('test.png', bbox_inches='tight', pad_inches=0.02, dpi = 300)
                #exit()
            ySize, xSize = data.shape

            if xc==0. and yc==0.:
                  xc = xSize/2.
                  yc = ySize/2.
    
            ySize, xSize = data.shape
            '''
            for k in range(ySize):
                for i in range(xSize):
                    if data[k,i]<0:
                        mask[k,i]=1        
            '''

            if layer==0 and (profile=='summed' or profile=='vert'):
              # Fill masked areas with the model

              for k in range(ySize):
                for i in range(xSize):
                  if data[k,i]<=0. or mask[k,i]>0.9:
                    data[k,i]=model_to_fil[k,i]

              
              '''
              interp=True
              if interp==True:
                    nan_frame = np.copy(data)
                    for k in range(ySize):
                        for i in range(xSize):
                            if mask[k,i]>0:
                                nan_frame[k,i] = float('nan')       
                    
                    interpolated = ds9_contour.replace_nans(nan_frame, 5, 0.5, int(ceil(FWHM)), "localmean")
                    data = interpolated
                    #im = plt.imshow(data, cmap='hot')
                    #plt.savefig('test.png', bbox_inches='tight', pad_inches=0.02, dpi = 300)
                    #exit()
              '''
              zzz=1
            
            if profile=='summed':
                I = []
                r = []
                if zmax==0.:
                  zmax = min([ySize-yc,yc])
                  
                if Rmax!=0.:
                  #zmin = int(zmin)
                  #zmax = int(zmax)
                  for x in range(0,xSize,1):
                    if fabs(x-xc)<=Rmax:
                      #
                      #I.append(sum(data[int(yc+zmin):int(yc+zmax),x])+sum(data[int(yc-zmax):int(yc-zmin),x]))
                      #'''
                      II = []
                      yranges = range(int(yc-zmax),int(yc-zmin),1) + range(int(yc+zmin),int(yc+zmax))
                      for y in yranges:
                        #if mask[y,x]<0.9:
                        II.append(data[y,x])                          
                      mean, median, std = sigma_clipped_stats(II, sigma=5.0, iters=5)
                      I.append(median)
                      #'''
                      r.append(x-xc)
                else:
                  for x in range(0,xSize,1):
                    I.append(sum(data[int(yc+zmin):int(yc+zmax),x])+sum(data[int(yc-zmax):int(yc-zmin),x]))
                    r.append(x-xc)  
                r = np.array(r)
                I = np.array(I)
            
            elif profile=='vert_old':
              if xc==0. and Rmin==0. and Rmax==0.:
                I = []
                r = []
                for y in range(0,ySize,1):
                  II = []
                  for x in range(0,xSize,1):
                      #if hdu[0].data[y,x]>0.:
                      II.append(data[y,x])
                  I.append(sum(II))
                  r.append(y-yc)
                r = np.array(r)
                I = np.array(I)   
              else:
                if Rmax==0.:
                  Rmax = max([xSize-xc,xc])
                I = []
                r = []
                for y in range(0,ySize,1):
                  II = []
                  ranges = range(int(xc-Rmax),int(xc-Rmin),1)+range(int(xc+Rmin),int(xc+Rmax),1)
                  for x in ranges:
                      #if hdu[0].data[y,x]>0.:
                      II.append(data[y,x])
                  I.append(sum(II))
                  r.append(y-yc)
                r = np.array(r)
                I = np.array(I)                   

            elif profile=='vert':
              if xc==0. and Rmin==0. and Rmax==0.:
                I = []
                r = []
                for y in range(0,ySize,1):
                  II = []
                  for x in range(0,xSize,1):
                      #if hdu[0].data[y,x]>0.:
                      II.append(data[y,x])
                  I.append(sum(II))
                  r.append(y-yc)
                r = np.array(r)
                I = np.array(I)   
              else:
                if Rmax==0. or Rmax>xc or Rmax>xSize-xc:
                  Rmax = max([xSize-xc,xc])
                I = []
                I_err = []
                r = []
                
                if zmin==0. and zmax==0.:
                    yranges = range(0,ySize,1)
                else:
                    if zmax>yc or zmax>ySize-yc:
                        Rmax = max([ySize-yc,yc])
                    yranges = range(int(yc-zmax),int(yc-zmin),1)+range(int(yc+zmin),int(yc+zmax),1)
                
                for y in yranges:
                  II = []
                  xranges = list(range(int(xc-Rmax),int(xc-Rmin),1))+list(range(int(xc+Rmin),int(xc+Rmax),1))
                  for x in xranges:
                    #if hdu[0].data[y,x]>0.:
                    #if mask[y,x]<0.9:
                    II.append(data[y,x])
                  #I.append(sum(II))
                  mean, median, std = sigma_clipped_stats(II, sigma=5.0, maxiters=5)
                  #mean = np.median(II)
                  I.append(median)
                  I_err.append(std)
                  r.append(y-yc)
                r = np.array(r)
                I = np.array(I)
                I_err = np.array(I_err)



            elif profile=='rad':
                I = []
                r = []
                from random import randrange as rd

                NN = 0

                Npix=30000
                if Rmax==0.:
                  while NN<Npix:
                          i = rd(0,xSize,1)
                          k = rd(0,ySize,1)
                          RR = sqrt( (i-xc+0.5)**2 + (k-yc+0.5)**2 )
                          if mask[k,i]==0.:
                            I.append(data[k, i])
                            r.append( RR )
                            NN = NN + 1
                else:
                  while NN<Npix:
                          i = rd(0,xSize,1)
                          k = rd(0,ySize,1)
                          RR = sqrt( (i-xc+0.5)**2 + (k-yc+0.5)**2 )
                          if RR<=Rmax and mask[k,i]==0.:
                            I.append(data[k, i])
                            r.append( RR )
                            NN = NN + 1

                r = np.array(r)
                I = np.array(I)

                maxI = np.max(data[int(yc)-10:int(yc)+10,int(xc)-10:int(xc)+10])
                
            elif profile=='azim':
              if mask_image!=None:
                mask_astropy = np.ones_like(np.array(mask,dtype=float),dtype=bool)
              else:
                mask_astropy = None
              
              mask_astropy = np.ones(data.shape,dtype='bool')
              mask_astropy[mask>0] = False

                
              if Rmax==0.:
                Rmax = min([xSize-xc,xc,ySize-yc,yc])
              if layer==0:
                  data=data ####################### WARNING data=data+1. for UGC4999
              bin_centers, radial_prof = radial_profile.azimuthalAverage(data, center=[xc,yc], returnradii=True, 
              binsize=step, weights=None, interpnan=False, mask=mask_astropy )

              I = []
              r = []      
              for k in range(len(bin_centers)):
                if bin_centers[k]<=Rmax:
                  r.append(bin_centers[k])
                  I.append(radial_prof[k])
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
                  min_mag = np.min(mag[~np.isnan(mag)])
                  max_mag = np.max(mag[~np.isnan(mag) & ~np.isinf(mag)])
                elif SB_units=='ADU/pix2':
                  mag = I
                  max_mag = np.min(mag[~np.isnan(mag)])
                  min_mag = np.max(mag[~np.isnan(mag) & ~np.isinf(mag)])  

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

                if AX==None:
                  plt.legend(loc=2, borderaxespad=0.,fontsize=legend_size,numpoints=1)
                else:
                  AX.legend(loc=2, borderaxespad=0.,fontsize=legend_size,numpoints=1)



            if profile=='cut':
                r, I = get_slice(input_image, xc, yc, PA, Rmax, layer=layer)
                #print(input_image, xc, yc, PA, Rmax)
                #exit()
                if SB_units=='mag/arcsec2':
                  mag = m0 - 2.5*log10(I)+ 5.*log10(pix2sec)
                  min_mag = np.min(mag[~np.isnan(mag)])
                  max_mag = np.max(mag[~np.isnan(mag) & ~np.isinf(mag)])
                elif SB_units=='ADU/pix2':
                  mag = I
                  max_mag = np.min(mag[~np.isnan(mag)])
                  min_mag = np.max(mag[~np.isnan(mag) & ~np.isinf(mag)])  
                if min_SB is not None:
                    min_mag = min_SB

                if max_SB is not None:
                    max_mag = max_SB                
                
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
                  r, I = get_slice(input_image, xc, yc, PA, Rmax, layer=layer)
                  if SB_units=='mag/arcsec2':
                    mag = m0 - 2.5*log10(I)+ 5.*log10(pix2sec)
                  elif SB_units=='ADU/pix2':
                    mag = I

                  if geom_units=='pix':
                    Radius = r
                  elif geom_units=='arcsec':
                    Radius = r*pix2sec
                  elif geom_units=='kpc':
                    Radius = r*pix2sec*Scale
                    
                  if AX==None:
                    plt.plot(Radius, mag,lyne_style,color=color2,lw=2,label = Label) 
                  else:
                    AX.plot(Radius, mag,lyne_style,color=color2,lw=2,label = Label)
                  #plt.ylim(27,18.)  
                if AX==None:
                  plt.legend(loc=2, borderaxespad=0.,fontsize=legend_size,numpoints=1)
                else:
                  AX.legend(loc=2, borderaxespad=0.,fontsize=legend_size,numpoints=1)



            if profile=='vert':
                if SB_units=='mag/arcsec2':
                  mag = m0 - 2.5*log10(I)+ 5.*log10(pix2sec)
                  #mag_err = 2.5*I_err/(I*2.303)
                  min_mag = np.min(mag[~np.isnan(mag)])
                  max_mag = np.max(mag[~np.isnan(mag) & ~np.isinf(mag)])
                  if max_SB!=None:
                        max_mag = min_SB
                        min_mag = max_SB
                        
                elif SB_units=='ADU/pix2':
                  mag = I
                  #mag_err = I_err
                  max_mag = np.min(mag[~np.isnan(mag)])
                  min_mag = np.max(mag[~np.isnan(mag) & ~np.isinf(mag)])  
                
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
                    plt.errorbar(Radius, mag, yerr=None, fmt='o',color='gray',markeredgecolor=mecolor,markersize=msize,label = Label) # markersize='4', markeredgecolor='lime', markerfacecolor='lime', ecolor='lime', zorder=overlayed_number+2)
                    
                    ylim(max_mag,min_mag-0.5)
                  else:
                    #AX.plot(Radius, mag,'o',color='gray',markeredgecolor=mecolor,markersize=msize,label = Label)
                    AX.errorbar(Radius, mag, yerr=None, fmt='o',color='gray',markeredgecolor=mecolor,markersize=msize,label = Label)
                    AX.set_ylim(max_mag,min_mag-0.5)  
                
                if layer!=0 and layer!=2 and layer!=3:
                  if AX==None:
                    plt.plot(Radius, mag,lyne_style,color=color1,lw=2,label = Label)
                    xlim(-max(Radius),max(Radius))
                  else:
                    AX.plot(Radius, mag,lyne_style,color=color1,lw=2,label = Label)
                    AX.set_xlim(-max(Radius),max(Radius))
                    
                
                if AX==None:
                  plt.legend(loc=2, borderaxespad=0.,fontsize=legend_size,numpoints=1)
                else:
                  AX.legend(loc=2, borderaxespad=0.,fontsize=legend_size,numpoints=1)


    
            if profile=='rad':
                  xsize = 9
                  ysize = 9
                  fsize = 9

                  xticks_length = 3.
                  xticks_width = 0.5

                  if SB_units=='mag/arcsec2':
                    mag = m0 - 2.5*log10(I)+ 5.*log10(pix2sec)
                    min_mag = np.min(mag[~np.isnan(mag)])
                    max_mag = np.max(mag[~np.isnan(mag) & ~np.isinf(mag)])
                  elif SB_units=='ADU/pix2':
                    mag = I
                    max_mag = np.min(mag[~np.isnan(mag)])
                    min_mag = np.max(mag[~np.isnan(mag) & ~np.isinf(mag)])  

                  if geom_units=='pix':
                    Radius = r
                  elif geom_units=='arcsec':
                    Radius = r*pix2sec
                  elif geom_units=='kpc':
                    Radius = r*pix2sec*Scale  
                  
                  if layer==0:
                    if AX==None:
                      plt.plot(Radius,mag,'o',markeredgecolor=color1,markerfacecolor = mecolor,markersize=1.5, label = Label, alpha=0.2)
                      if SB_units=='mag/arcsec2':
                        ylim(33., m0 - 2.5*log10(maxI)+ 5.*log10(pix2sec) - 0.5)
                    else:
                      AX.plot(Radius,mag,'o',markeredgecolor=color1,markerfacecolor = mecolor,markersize=1.5, label = Label, alpha=0.2)
                      if SB_units=='mag/arcsec2':
                        AX.set_ylim(33., m0 - 2.5*log10(maxI)+ 5.*log10(pix2sec) - 0.5)
                  
                  if layer==1:
                    if AX==None:
                      plt.plot(Radius,mag,'o',markeredgecolor=mecolor,markerfacecolor = mecolor,markersize=0.5, label = Label, alpha=0.7)
                    else:
                      AX.plot(Radius,mag,'o',markeredgecolor=mecolor,markerfacecolor = mecolor,markersize=0.5, label = Label, alpha=0.7)    

                  elif layer>1:
                    if AX==None:
                      plt.plot(Radius,mag,'o',markeredgecolor=color1,markerfacecolor = mecolor,markersize=0.5, label = Label, alpha=0.7)
                    else:
                      AX.plot(Radius,mag,'o',markeredgecolor=color1,markerfacecolor = mecolor,markersize=0.5, label = Label, alpha=0.7)
                    
                  if Rmax==0.:
                    Rmax = min([xSize-xc,xc,ySize-yc,yc])
                  
                  if AX==None:
                    if geom_units=='arcsec':
                      xlim(0,Rmax*pix2sec)
                    elif geom_units=='pix':
                      xlim(0,Rmax)
                    elif geom_units=='kpc':
                      xlim(0,Rmax*pix2sec*Scale)
                  else:
                    if geom_units=='arcsec':
                      AX.set_xlim(0,Rmax*pix2sec)
                    elif geom_units=='pix':
                      AX.set_xlim(0,Rmax) 
                    elif geom_units=='kpc':
                      AX.set_xlim(0,Rmax*pix2sec*Scale)
                      
                  if AX==None:
                    plt.yticks(size=ysize)
                    plt.xticks(size=xsize)
                    plt.tick_params(length=xticks_length, width=xticks_width)
                    plt.minorticks_on()
                    plt.tick_params(length=xticks_length-5, width=xticks_width,which='minor',direction='out')
                    plt.tight_layout(pad=1.0, h_pad=None, w_pad=None, rect=[0.03,0.03,1.01, 1.04])

                    plt.legend(loc=1, numpoints=1,prop={'size':legend_size},markerscale=5)
                  
                  else:
                    #AX.set_yticks(size=ysize)
                    #AX.set_xticks(size=xsize)
                    AX.tick_params(axis='both',color='black',which='major',length=7)
                    AX.tick_params(axis='both',color='black',which='minor',length=4)
                    #plt.tight_layout(pad=1.0, h_pad=None, w_pad=None, rect=[0.03,0.03,1.01, 1.04])

                    AX.legend(loc=1, numpoints=1,prop={'size':legend_size},markerscale=5)
  
 
 
 
            if profile=='azim':
                if SB_units=='mag/arcsec2':
                  mag = m0 - 2.5*log10(I)+ 5.*log10(pix2sec)
                  min_mag = np.min(mag[~np.isnan(mag)])
                  max_mag = np.max(mag[~np.isnan(mag) & ~np.isinf(mag)])
                elif SB_units=='ADU/pix2':
                  mag = I
                  max_mag = np.min(mag[~np.isnan(mag)])
                  min_mag = np.max(mag[~np.isnan(mag) & ~np.isinf(mag)])  
                if max_SB is not None:
                    max_mag = max_SB
                if min_SB is not None:
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
                    
                if AX==None:
                  plt.legend(loc=1, borderaxespad=0.,fontsize=legend_size,numpoints=1)
                else:
                  AX.legend(loc=1, borderaxespad=0.,fontsize=legend_size,numpoints=1)
                
                plt.xlim(1,max(Radius))

            RADIUS.append(Radius)
            MAG.append(mag)
    
      if text is not None:
           plt.text(0.8, 0.9, text, fontsize=9., color='black',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline')
           
      if output_file==None and save_out_ima == True:
        

        print(output_file)
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
      #plt.show()
      
      
      #plt.axvline(x=146., ls='--', color='blue') # NGC4302
      #plt.axvline(x=-146., ls='--', color='blue') # NGC4302

      #plt.axvline(x=269.7, ls=':', color='green') # NGC4302
      #plt.axvline(x=-269.7, ls=':', color='green') # NGC4302
      
      
      
      
      #plt.axvline(x=240., ls='--', color='blue') # NGC891
      #plt.axvline(x=-240., ls='--', color='blue') # NGC891      
      #plt.axvline(x=446.3, ls=':', color='green') # NGC891
      #plt.axvline(x=-446.3, ls=':', color='green') # NGC891    
      
      
      
      
      
      
      #plt.axvline(x=185., ls='-.', color='red') # NGC3628
      #plt.axvline(x=-185., ls='-.', color='red') # NGC3628  
      
      #plt.axvline(x=366., ls='--', color='blue') # NGC3628
      #plt.axvline(x=-366., ls='--', color='blue') # NGC3628    

      #plt.axvline(x=486.0, ls=':', color='green') # NGC3628
      #plt.axvline(x=-486.0, ls=':', color='green') # NGC3628    

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
    parser.add_argument("--comps", nargs='?', const=1, help="Optional: Names of components (fro composed_model.fits) ", type=str, default=None) 

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
    zmin = float(args.Rmin)
    zmax = float(args.Rmax)    
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
    
    main(input_image,m0,pix2sec,mask_image=mask_image,profile=profile,xc=xc,yc=yc,PA=PA,Rmin=Rmin,Rmax=Rmax,zmin=zmin,zmax=zmax,output_file=output_file,AX=None,geom_units=geom_units,SB_units=SB_units,Scale=Scale,legend_size=legend_fsize,min_SB=SBmin,max_SB=SBmax, comps=comps)
