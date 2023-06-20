#! /usr/bin/env python
"""
Genetic algorithms applied to GALFIT optimization problem
"""

from random import random
from math import sqrt

import time
import sys
import os
import pyfits
import subprocess
import shutil
import glob
import numpy as np

try:
  inputGalfFile = str(sys.argv[1])
except:
  print "No galfit file is specified. No comments are added in the outptut file."
  inputGalfFile = 'none'

# Find all generations_res.txt files:
directories = [x[0] for x in os.walk('./')]
files = []
for directory in directories:
  for file in glob.glob(directory+"/"+"generations_res.txt"):
      files.append(file)
      
      
# Read the header:
Header = []
for file in files:
  with open(file) as f:
    Header.append(f.readline())
    
s = set([x for x in Header if Header.count(x) >= len(files)])

if len(s)!=1:
  print "ERROR: Headers in outptut files are not identical!"
  exit()

f_gen = open('Results_GA.txt', "w")

print >>f_gen, "Number of files is: %i\n\n" % (len(files))
print >>f_gen, "Par\tBest\tMean\tstd\tComments"

# Read all last lines and stack the results:
Pars = []
for file in files:
  myvar = open(file, 'r')
  lst = myvar.readlines()
  last_line = lst[len(lst)-1]
  myvar.close()
  
  pars = last_line.split('\t')
  Pars.append(pars)

Pars = np.array(Pars, float)

for k in range(len(Pars[0])):
  vector1 = Pars[:,k]
  #print np.mean(vector1),np.std(vector1)

min_Chi2 = min(vector1)
ind_best = list(vector1).index(min_Chi2)
best_par = Pars[ind_best]


header = Header[0].split('\t')
header[-1] = header[-1].split('\n')[0]

comments = []

try:
  # read galfit input file for comments:
  fp = open(inputGalfFile,'r')
  lines_inp = fp.readlines()
  comments.append('#')
  for k in range(1,len(Pars[0])):
      for i in range(len(lines_inp)):
	  if header[k].split('/')[0]!='Chi2':
	    if int(header[k].split('/')[0]) == i + 1:
	      comments.append(lines_inp[i].split('#')[1].split('|')[1])
	  if header[k].split('/')[0]=='Chi2' and 'chi2' not in comments:
	    comments.append('Chi^2/nu')
	  if header[k].split('/')[0]=='#' and 'generation' not in comments:
	    comments.append('generation')  

  fp.close()


  for k in range(1,len(Pars[0])):
      vector1 = Pars[:,k]
      print >>f_gen, '%s\t%.3f\t%.3f\t%.3f\t%s' % (header[k],best_par[k],np.mean(vector1),np.std(vector1),comments[k].split('\n')[0])
   
except:
  for k in range(1,len(Pars[0])):
    vector1 = Pars[:,k]
    print >>f_gen, '%s\t%.3f\t%.3f\t%.3f' % (header[k],best_par[k],np.mean(vector1),np.std(vector1))




print >>f_gen, '\n\nMinimum Chi2 is: %.3f' % (min_Chi2)
print >>f_gen, 'Best model is in: %s' % (files[ind_best])
f_gen.close()



