# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import math
import sys
import shutil
import argparse
import os
import glob
import collections
import warnings
warnings.filterwarnings("ignore")


#Func
def read_decription(module):
    f = open(module, 'r')
    lines = f.readlines()
    desc_lines = ''
    note_lines = ''
    usage_lines = ''
    for k in range(len(lines)):
      if k<20:
        line = lines[k]
        if '# DESCRIPTION:' in line:
            for i in range(0,20):
              try:
                if '#' in lines[k+i] and '# NOTE:' not in lines[k+i] and '# MINIMAL USAGE:' not in lines[k+i]:
                   desc_lines = lines[k+i][2:]
                else:
                    break
              except:
                  z=1

        if '# NOTE:' in line:
            for i in range(0,20):
              try:
                if '#' in lines[k+i] and '# MINIMAL USAGE:' not in lines[k+i]:
                   note_lines = lines[k+i]
                else:
                    break
              except:
                  z=1


        if '# MINIMAL USAGE:' in line:
            for i in range(0,20):
              try:
                if '#' in lines[k+i] and '# MINIMAL USAGE:' in lines[k+i]:
                   usage_lines = lines[k+i][2:]
                else:
                    break
              except:
                  z=1

                
    return desc_lines,note_lines,usage_lines


def fast_scandir(dir):
    subfolders= [f.path for f in os.scandir(dir) if f.is_dir()]
    for dir in list(subfolders):
        if '.git' not in dir:
            subfolders.extend(fast_scandir(dir))
    return subfolders

subfolders = fast_scandir('.')


packages = collections.OrderedDict()
for folder in subfolders:
    packages[folder] = glob.glob('%s/*.py' % (folder))


f = open('functions.txt','w')

for k, v in packages.items():
    print('Package %s' % (k))
    #f.write('Package %s\n' % (k))
    Package_written = False
    for file in v:
      if '__init__.py' not in file and '.git' not in file:
        print('\t%s' % (file))
        #f.write('\t%s\n' % (file))
        desc_lines,note_lines,usage_lines = read_decription(file)
        
        if (desc_lines!='' or note_lines!='' or usage_lines!='') and Package_written!=True:
            f.write('%s\n\n' % (k))
            Package_written = True
        
        if desc_lines!='' or note_lines!='' or usage_lines!='':
            f.write('\t%s\n' % (file))
        
        if desc_lines!='':
            f.write('\t\tDESCRIPTION: %s' % (desc_lines))
            
        if note_lines!='':# and ' # NOTE' in note_lines:
            f.write('\t\t%s' % (note_lines))
        
        if usage_lines!='':
            f.write('\t\t%s' % (usage_lines))
        
        if desc_lines!='' or note_lines!='' or usage_lines!='':
            f.write('\n\n')
f.close()
