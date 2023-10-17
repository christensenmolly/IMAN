#!/usr/bin/python
import sys
import collections
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
from matplotlib.ticker import NullFormatter
import numpy as np
from math import *
from scipy import stats
import warnings
import matplotlib
import types
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FormatStrFormatter
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True
from matplotlib.ticker import MaxNLocator
warnings.filterwarnings("ignore")
import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')
mpl.rcParams['errorbar.capsize'] = 3

def find_regression(xx, yy):
        #k1, b1, r_value, p_value, err, err1 = stats.linregress(xx,yy)
        res = stats.linregress(xx,yy)
        k1 = res.slope
        b1 = res.intercept
        r_value = res.rvalue
        p_value = res.pvalue
        slope_err = res.stderr
        intercept_err = res.intercept_stderr

        k2, b2, r_value, p_value, std_err = stats.linregress(yy,xx)

        t =  (1.0-k1*k2)/(k1+k2)
        t2 =  t/(1.0+sqrt(1.0+t*t))
        k4 = (t2+k1)/(1.0-t2*k1)
        xc = ((b2+b1*k2)/(1.0-k1*k2))
        yc = k1*xc + b1
        b4 = yc - k4*xc
        return k1,b1,k2,b2,k4,b4,r_value,slope_err,intercept_err
    
def RegrLine(xxx, yyy, xx1, xx2, yy1, yy2, i, sigma_outliers=None):
        xx = []
        yy = []
        if xx1>xx2:
            xx1,xx2 = xx2,xx1
        if yy1>yy2:
            yy1,yy2 = yy2,yy1


        for k in range(0,len(xxx)):
            if xxx[k]>xx1 and xxx[k]<xx2 and yyy[k]>yy1 and yyy[k]<yy2:
                xx.append(xxx[k])
                yy.append(yyy[k])

        k1,b1,k2,b2,k4,b4,r_value,slope_err,intercept_err = find_regression(xx,yy)

        if sigma_outliers is not None:
            XX = []; YY = []; Dist = []
            for k in range(len(xx)):
                Dist.append( abs(k4*xx[k]-yy[k]+b4) / ( sqrt(k4*k4+1.) ) )
                
            std_dist = np.std(Dist)    

            for k in range(len(xx)):                
                if Dist[k]<=sigma_outliers*std_dist:
                    XX.append(xx[k])
                    YY.append(yy[k])
            k1,b1,k2,b2,k4,b4,r_value,slope_err,intercept_err = find_regression(XX,YY)    

        if i==1:
                return k1,b1,r_value,slope_err,intercept_err 
        if i==2:
                return 1./k2,-b2/k2,r_value,slope_err,intercept_err
        if i==3:
                return k4,b4,r_value,slope_err,intercept_err 

def xticks_restrict_to_integer(ax):
    aa_x = ax.get_xticks().tolist()

    xint = []
    for each in aa_x:
        xint.append(int(each))
    return xint


def yticks_restrict_to_integer(ax):
    aa_y = ax.get_yticks().tolist()
    yint = []
    for each in aa_y:
        yint.append(int(each))
    return yint

def myLogFormat(y,pos):
    # Find the number of decimal places required
    decimalplaces = int(np.maximum(-np.log10(y),0))     # =0 for numbers >=1
    # Insert that number into a format string
    formatstring = '{{:.{:1d}f}}'.format(decimalplaces)
    # Return the formatted tick label
    return formatstring.format(y)

def multicolor_ylabel(ax,list_of_strings,list_of_colors,axis='x',anchorpad=0,**kw):
    """this function creates axes labels with multiple colors
    ax specifies the axes object where the labels should be drawn
    list_of_strings is a list of all of the text items
    list_if_colors is a corresponding list of colors for the strings
    axis='x', 'y', or 'both' and specifies which label(s) should be drawn"""
    from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker

    # x-axis label
    if axis=='x' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',**kw)) 
                    for text,color in zip(list_of_strings,list_of_colors) ]
        xbox = HPacker(children=boxes,align="center",pad=0, sep=5)
        anchored_xbox = AnchoredOffsetbox(loc=3, child=xbox, pad=anchorpad,frameon=False,bbox_to_anchor=(0.2, -0.09),
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_xbox)

    # y-axis label
    if axis=='y' or axis=='both':
        boxes = [TextArea(text, textprops=dict(color=color, ha='left',va='bottom',rotation=90,**kw)) 
                     for text,color in zip(list_of_strings[::-1],list_of_colors) ]
        ybox = VPacker(children=boxes,align="center", pad=0, sep=5)
        anchored_ybox = AnchoredOffsetbox(loc=3, child=ybox, pad=anchorpad, frameon=False, bbox_to_anchor=(-0.10, 0.2), 
                                          bbox_transform=ax.transAxes, borderpad=0.)
        ax.add_artist(anchored_ybox)


def main(x,y,x_err=None,y_err=None,x_dim=None,y_dim=None,
	 figure=[None,0,None,0,[5,5],True],
	 text=None,
	 regression_line=None,
	 marker=['o',5,'black','black'],
	 line=[None,None,None],
	 labels=[None,None,13,'black'],
	 legend=None,
	 diag=False,
	 xticks=None,
	 yticks=None,
	 add_letter=None,
	 add_line=None,
	 add_points=None,
	 x_logscale=False,
	 y_logscale=False,
	 capsize=0.,
	 verbose=True):
  #print xticks, 'HERE'
  # Define all graph parameters:
  [ax, overlayed_number, output_image, plot_number, image_dims, clean_image] = figure
  [marker_symbol, marker_size, marker_face_color, marker_edge_color] = marker
  [line_type, line_width, line_color] = line
  
  [label_x, label_y, label_font_size, label_color] = labels
  if verbose:
    print('Plotting: %s VS. %s' % (label_x, label_y))
  
  if regression_line!=None:
    if len(regression_line)==4:
        [regr_line_type, regr_line_width, regr_line_color, sel_galaxies] = regression_line
    else:
        [regr_line_type, regr_line_width, regr_line_color] = regression_line
        sel_galaxies = None
  if legend!=None:
    [legend_text,legend_location, legend_font_size, legend_color] = legend
  else:
     legend_text = None 
  if x_dim==None:
    xmin = np.min(x)
    xmax = np.max(x)
  else:
    [xmin,xmax] = x_dim

  if y_dim==None:
    ymin = np.min(y)
    ymax = np.max(y)
  else:
    [ymin,ymax] = y_dim
  
  #If overlay the data using existing ax:
  if ax==None:
    fig = plt.figure(plot_number, figsize=(image_dims[0], image_dims[1]))
    ax = fig.add_subplot(111)
    #ax=plt.subplot()
    
    
    #fig, ax= plt.subplots(figsize=(image_dims[0], image_dims[1]))
    
  # Plot the data:
  if line_width is None:
    ax.errorbar(np.array(x,float), np.array(y,float), xerr=x_err, yerr=y_err, fmt=marker_symbol, markersize=marker_size, markeredgecolor=marker_edge_color, markerfacecolor=marker_face_color, ecolor=marker_edge_color, color=marker_face_color, zorder=overlayed_number, capsize=capsize, label=legend_text)
  else:
    ax.errorbar(np.array(x,float), np.array(y,float), xerr=x_err, yerr=y_err, fmt=marker_symbol, markersize=marker_size, markeredgecolor=marker_edge_color, markerfacecolor=marker_face_color, ecolor=marker_edge_color, color=line_color, zorder=overlayed_number, capsize=capsize, label=legend_text, lw=line_width, ls=line_type)
  overlayed_number = overlayed_number + 1
  
  if add_points!=None:
      for point in add_points:
          [xx,yy,p_marker_symbol,p_marker_size,p_marker_color] = point
          ax.plot(xx,yy, p_marker_symbol, markeredgecolor=p_marker_color, markerfacecolor=p_marker_color, markersize=p_marker_size)
          
  # Add a regression line:
  plot_all_regr_line = 'yes'
  if regression_line!=None:
    
    if isinstance(regr_line_type, list):
        [regr_line_type,coords] = regr_line_type
        if len(coords)==4:
            [minx,maxx,miny,maxy] = coords
            plot_all_regr_line = 'no'
        else:
            [minx,maxx,miny,maxy,plot_all_regr_line] = coords
    else:        
        minx = xmin
        maxx = xmax
        miny = ymin
        maxy = ymax
    if sel_galaxies!=None:
        k,b,r_value,slope_err,intercept_err = RegrLine(np.array(x)[sel_galaxies],np.array(y)[sel_galaxies],minx,maxx,miny,maxy,3)
    else:
        k,b,r_value,slope_err,intercept_err = RegrLine(x,y,minx,maxx,miny,maxy,3)
    if b>=0.:
        print('\t Regression line: Y = %.3f * X + %.3f (r=%.2f)\n' % (k,b,r_value))
    else:
        print('\t Regression line: Y = %.3f * X - %.3f (r=%.2f)\n' % (k,abs(b),r_value))
        
    if plot_all_regr_line == 'yes':
        x_regr = np.arange(xmin,xmax,(xmax-xmin)/100.)
    else:
        x_regr = np.arange(minx,maxx,(maxx-minx)/100.)
    y_regr = k*x_regr + b
    ax.plot(x_regr, y_regr, regr_line_type, lw=regr_line_width, color=regr_line_color, zorder=overlayed_number)
    overlayed_number = overlayed_number + 1

  if add_line!=None: 
    [k_line, k_line_min, k_line_max, b_line, b_line_min, b_line_max] = add_line
    x_regr = np.arange(xmin,xmax,(xmax-xmin)/100.)
    y_regr = k_line*x_regr + b_line
    ax.plot(x_regr, y_regr, line[0], lw=line[1], color=line[2], zorder=overlayed_number)
    
    if k_line_max!=None:
        y1 = k_line_max*x_regr + b_line_max
        y2 = k_line_min*x_regr + b_line_min
        ax.fill_between(x_regr, y1, y2, color='paleturquoise', alpha='0.5')
    
    overlayed_number = overlayed_number + 1
    
  if diag==True:
    xx = np.arange(xmin,xmax,(xmax-xmin)/100.)
    yy = xx
    ax.plot(xx, yy, '-.', lw=1, color='black', zorder=overlayed_number)
    overlayed_number = overlayed_number + 1
    



  #print ymin, ymax,xmin, xmax
  #plt.xticks(fontsize=int(label_font_size/1.3))#, rotation=90)
  aa_x = ax.get_xticks().tolist()


  if fabs(xmax - xmin)>5:
    aa_x = xticks_restrict_to_integer(ax)
  ax.set_xticklabels(aa_x, fontsize=int(label_font_size/1.3))
  #plt.yticks(fontsize=int(label_font_size/1.3))#, rotation=90)
  aa_y = ax.get_yticks().tolist()
  if fabs(ymax - ymin)>5:
    aa_y = yticks_restrict_to_integer(ax)

  ax.set_yticklabels(aa_y, fontsize=int(label_font_size/1.3))
  
  
  
  
  '''
  if abs(xmax - xmin)>2:
    ticks_restrict_to_integer(ax.xaxis)

  if abs(ymax - ymin)>2:
    ticks_restrict_to_integer(ax.yaxis)
  '''  

  
  # Set axes labels:
  if isinstance(label_color, str):
      label_color_x = label_color
      label_color_y = label_color
  else:
      label_color_x = label_color[0]
      label_color_y = label_color[1]
      
  if label_x!=None:
    ax.set_xlabel(label_x, fontsize=label_font_size, color='black')
  if label_y!=None:  
    ax.set_ylabel(label_y, fontsize=label_font_size, color=label_color_y)
  
  # Set legend:
  if legend!=None:
    legend = ax.legend(loc=legend_location, shadow=False, fontsize=legend_font_size, facecolor=legend_color, numpoints=1)

  if text is not None:
    if len(text)==4:
        [Text, text_loc_x, text_loc_y, color] = text
        plt.figtext(text_loc_x, text_loc_y, Text, color=color, fontweight='bold', fontsize=label_font_size+2, horizontalalignment='center',verticalalignment='center',backgroundcolor='white')
    else:
        plt.figtext(0.5, 0.82, text, color='black', fontweight='bold', fontsize=label_font_size+2, horizontalalignment='center',verticalalignment='center',backgroundcolor='white')
  
  if add_letter!=None:
      [letter, hor_pos, vert_pos, text_size] = add_letter
      plt.text(hor_pos, vert_pos, letter, fontsize=text_size, color='black',transform=ax.transAxes, horizontalalignment='left',verticalalignment='baseline', bbox=dict(facecolor='white', edgecolor='black'))
  
  if x_logscale==True:
      ax.set_xscale('log')
      #ax.tick_params(direction='in', length=4, width=1, colors='black')
      #ax.xaxis.set_major_formatter(ticker.FuncFormatter(myLogFormat))

      plt.tick_params(axis='x', which='minor')
      ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))
      

  if y_logscale==True:
      ax.set_yscale('log')
      ax.yaxis.set_major_formatter(ticker.FuncFormatter(myLogFormat))

  for axis in [ax.xaxis, ax.yaxis]:
    axis.set_major_formatter(ScalarFormatter())

  # Set axes borders:
  ax.set_xlim(xmin, xmax)
  ax.set_ylim(ymin, ymax)

  if xticks is not None:
    if type(xticks[0])!=float and type(xticks[0])!=int:

        ax.set_xticks(list(x))#, list(xticks))
        ax.set_xticklabels(list(xticks))
    else:
        ax.set_xticks(xticks)

  if yticks is not None:
    ax.set_yticks(yticks) 
 
  if output_image!=None:
    plt.savefig(output_image, transparent = False, dpi=200, bbox_inches='tight', pad_inches=0.05)
  else:
    if clean_image==True:
      plt.show()
      
  if clean_image==True:
    plt.clf()
    plt.close()
    
  return ax,overlayed_number
