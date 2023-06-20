#!/usr/bin/env python
# -*- coding: utf8 -*-

# Import the necessary modules
import pyfits
import numpy as np
import itertools
import matplotlib.pyplot as plt
import sys
import math

def main(first_file,second_file,operation,output_file):
    hdulist1 = pyfits.open(first_file)
    frame1 = hdulist1[0].data

    try:
      frame2 = float(second_file)
    except:
      hdulist2 = pyfits.open(second_file)
      frame2 = hdulist2[0].data

    sframe = np.copy(frame1)

    if operation=='add':
      sframe = frame1 + frame2
    elif operation=='sub':
      sframe = frame1 - frame2  
    elif operation=='mul':
      sframe = frame1 * frame2    
    elif operation=='div':
      sframe = frame1 / frame2  


    hdu = pyfits.PrimaryHDU(sframe)
    hdu.writeto(output_file,clobber=True)

