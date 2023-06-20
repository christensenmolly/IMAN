#! /usr/bin/env python

import sys
from os.path import exists
from os import remove
from scipy.ndimage.filters import convolve
from numpy import max, sum
import pyfits
import numpy as np
import collections
if True:
    # 1.Read catalogue of galaxies
    f = open('fullsqlmean.cgi', 'r')
    lines = f.readlines()
    header = lines[1].split('|')
    
    for k in range(len(header)):
      keyword = header[k]
      header[k] = keyword.replace(" ", "")
    #print header
    #exit()
    data = collections.OrderedDict()
    for i in range(len(header)):
      data[header[i]] = []

    for k in range(3,len(lines)):
      line = lines[k].split('|')
      
      for i in range(len(line)):
	  try:
	    if header[i]=='pgc':
	      line[i] = str(line[i].replace(" ", ""))
	    else:
	      line[i] = float(line[i])
	  except:
	    line[i] = str(line[i])
	  data[header[i]].append(line[i])

    f.close()
    

    
f = open('IRAS_table.dat', 'w')
f.write('|      ra    |       dec       |\n')
f.write('|  double    |    double       |\n')
f.write('|     deg    |       deg       |\n')

for k in range(len(data['al2000'])):
 if k<1000:
  f.write('  %.6f     %.6f\n' % (float(data['al2000'][k])*15.,float(data['de2000'][k])))

f.close()
