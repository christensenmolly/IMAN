#! /usr/bin/env python
# python ~/MEGA/MyPrograms/IMAN/Functions/correlator.py final_results.dat
import sys
from os.path import exists
from os import remove
from scipy.ndimage.filters import convolve
from scipy import signal
import numpy as np
from numpy import *
from numpy import max, sum
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.patches as patches
import matplotlib.path as path
import subprocess
import os
import argparse
import collections
from scipy import stats
matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True




def RegrLine(xxx,yyy,xx1,xx2,yy1,yy2,i,typ):
    xx = []
    yy = []
    if xx1>xx2:
        xxx1=xx2
        xxx2=xx1
    else:
        xxx1=xx1
        xxx2=xx2
    if yy1>yy2:
        yyy1=yy2
        yyy2=yy1
    else:
        yyy1=yy1
        yyy2=yy2
    for k in range(0,len(xxx)):
        if xxx[k]>xxx1 and xxx[k]<xxx2 and yyy[k]>yyy1 and yyy[k]<yyy2:
            xx.append(xxx[k])
            yy.append(yyy[k])

    k1, b1, r_value, p_value, std_err = stats.linregress(xx,yy)
    x = arange(xx1,xx2,0.01)
    y = k1*x + b1

    k2, b2, r_value, p_value, std_err = stats.linregress(yy,xx)
    x = arange(yy1,yy2,0.01)
    y1 = k2*x + b2
    x = arange(xx1,xx2,0.01)
    y2 = x/k2 - b2/k2

    t =  (1.0-k1*k2)/(k1+k2)
    t2 =  t/(1.0+sqrt(1.0+t*t))
    k4 = (t2+k1)/(1.0-t2*k1)
    xc = ((b2+b1*k2)/(1.0-k1*k2))
    yc = k1*xc + b1
    b4 = yc - k4*xc
    y4 = k4*x + b4
    if typ==1:    typ = 'r--'
    if typ==2:    typ = 'r-'
    if i==1:
        #plt.plot(x,y,typ,lw=4,color='black')
        if typ=='r-':    plt.plot(x,y,typ,lw=4,color='black')
        else:    plt.plot(x,y,'r-',lw=4,color='red')
        return k1,b1,r_value
    if i==2:
        #plt.plot(x,y2,typ,lw=4,color='black')
        if typ=='r-':    plt.plot(x,y2,typ,lw=4,color='black')
        else:
            y2 = -0.116*x + 1.384
            plt.plot(x,y2,'r-',lw=4,color='red')
        return 1./k2,-b2/k2,r_value
    if i==3:
        if typ=='r--':    plt.plot(x,y4,typ,lw=2,color='red')
        else:    plt.plot(x,y4,'r--',lw=2,color='red')
        return k4,b4,r_value













def main(Input_dict, Items=None, min_rc=0.5):
  if not os.path.exists('./found_cors'):
      os.makedirs('./found_cors')
  f_res = open('./found_cors/results.dat', 'w')


  #OBJ_NAMES = COLUMNS['ID/NAME']

  if Items is None:
    Items = []
    for key, value in results_all.items():
      Items.append(key)


  for k in range(len(Items)):
    print(Items[k])
    name_var1 = Items[k]
    var1 = Input_dict[Items[k]]
    
    try:
      var1 = np.array(var1, float)
    except:
      continue
    
    for i in range(k+1, len(Items)):
      name_var2 = Items[i]
      print('\t',Items[i])
      var2 = Input_dict[Items[i]]
      
      try:
            var2 = np.array(var2, float)
      except:
            continue
      
      min_var1 = np.min(var1[~np.isnan(var1)])
      max_var1 = np.max(var1[~np.isnan(var1)])

      min_var2 = np.min(var2[~np.isnan(var2)])
      max_var2 = np.max(var2[~np.isnan(var2)])

      try:
            kk,b,r = RegrLine(var1,var2,min_var1,max_var1,min_var2,max_var2,3,1)
      except:
            continue

      if fabs(r)<min_rc:
            continue
      else:
            f_res.write('%s\t%s\t%.5f\t%.5f\t%.3f\n' % (name_var1,name_var2, kk, b, r))



      figure = plt.figure(figsize=(5,5))
      plt.errorbar(var1, var2, xerr=None, yerr=None, fmt='o')

      delta_x = ( max(var1) - min(var1) ) / 50.
      delta_y = ( max(var2) - min(var2) ) / 50.
      if np.isnan(delta_x)==True or np.isnan(delta_y)==True or np.isinf(delta_x)==True or np.isinf(delta_y)==True:
        delta_x = 0.5
        delta_y = 0.5

      #for ii in range(len(var1)):
      #  plt.annotate(OBJ_NAMES[ii], xy = (var1[ii],var2[ii]), xytext = (float(var1[ii])+delta_x,float(var2[ii])+delta_y))


      plt.xlabel(name_var1, fontsize=15)
      plt.ylabel(name_var2, fontsize=15)


      filename = name_var1 + '_' + name_var2 + '.png'


      plt.savefig('./found_cors/'+filename, bbox_inches='tight', pad_inches=0.25)



  f_res.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Correlator")
    parser.add_argument("input_file", help="Input data file")
    parser.add_argument("--columns", nargs='?', const=1, help="Optional: Input the the columns separated by comma starting from 1, e.g. 1,2,8,10",type=str,default='all') 
    parser.add_argument("--separator", nargs='?', const=1, help="Optional: Input the separator between the columns",type=str,default='\t')
    parser.add_argument("--min_rc", nargs='?', const=1, help="Optional: Input the minimum regression coefficient to select correlations",type=float,default=0.5)




    args = parser.parse_args()

    input_file = args.input_file
    columns = args.columns
    separator = args.separator
    min_rc = args.min_rc

    f = open(input_file, 'r')
    lines = f.readlines()

    # Read the header lines
    col_names = lines[0].split(separator)
    col_descr = lines[1].split(separator)
    Total_col_number = len(col_names)
    #print Total_col_number, len(col_descr)
    #exit()

    COLUMNS = collections.OrderedDict()

    COLUMNS['ID/NAME'] = []

    for i in range(2,len(lines)):
        COLUMNS['ID/NAME'].append(lines[i].split(separator)[0])
        for k in range(1,Total_col_number):
          COLUMNS[col_names[k]+'/'+col_descr[k]] = []#'A'#.append(lines[i].split(separator))


    for i in range(2,len(lines)):
        for k in range(1,Total_col_number):
          COLUMNS[col_names[k]+'/'+col_descr[k]].append(lines[i].split(separator)[k])


    main(COLUMNS)
