#! /usr/bin/env python
# Inspired by http://www.bdnyc.org/2012/10/decimal-deg-to-hms/
import numpy as np
from numpy import *
from math import *

'''
def sexa_to_deci(ra,dec):
    if ':' in ra:
        ra = ra.split(':')
        hours = float(ra[0])
        mins = float(ra[1])
        secs = float(ra[2])
    elif 'h' in ra:
        hours = float(ra.split('h')[0])
        mins = float(ra.split('h')[-1].split('m')[0])
        secs = float(ra.split('m')[-1].split('s')[0])
        
    ra = (hours*15.*3600. + mins*15.*60. + secs*15.) / 3600.
    

    if ':' in dec:
        dec = dec.split(':')
        degs = float(dec[0])
        mins = float(dec[1])
        secs = float(dec[2])
    elif 'd' in dec:
        degs = float(dec.split('d')[0])
        mins = float(dec.split('d')[-1].split('m')[0])
        secs = float(dec.split('m')[-1].split('s')[0])    

    dec = np.sign(degs)*(abs(degs*3600.) + mins*60. + secs) / 3600.
    
    return ra,dec
'''




def deci_to_sexa(ra='', dec='', format_out='space'):
  RA, DEC, rs, ds = '', '', '', ''
  if dec:
    if str(dec)[0] == '-':
      ds, dec = '-', abs(dec)
    deg = int(dec)
    decM = abs(int((dec-deg)*60))
    decS = (abs((dec-deg)*60)-decM)*60

    if deg<10:
        if format_out=='space':
            DEC = '{0}0{1} {2} {3}'.format(ds, deg, decM, round(decS, 2))
        elif format_out==':':    
            DEC = '{0}0{1}:{2}:{3}'.format(ds, deg, decM, round(decS, 2))
        else:
            DEC = '{0}0{1}d{2}m{3}s'.format(ds, deg, decM, round(decS, 2))        
    else:
        if format_out=='space':
            DEC = '{0}{1} {2} {3}'.format(ds, deg, decM, round(decS, 2))
        elif format_out==':':    
            DEC = '{0}{1}:{2}:{3}'.format(ds, deg, decM, round(decS, 2))
        else:
            DEC = '{0}{1}d{2}m{3}s'.format(ds, deg, decM, round(decS, 2))
  
  if ra:
    if str(ra)[0] == '-':
      rs, ra = '-', abs(ra)
    raH = int(ra/15)
    raM = int(((ra/15)-raH)*60)
    raS = ((((ra/15)-raH)*60)-raM)*60

    if raH<10:
        if format_out=='space':
            RA = '{0}0{1} {2} {3}'.format(rs, raH, raM, round(raS, 2))
        elif format_out==':':
            RA = '{0}0{1}:{2}:{3}'.format(rs, raH, raM, round(raS, 2)) 
        else:
            RA = '{0}0{1}h{2}m{3}s'.format(rs, raH, raM, round(raS, 2))        
    else:
        if format_out=='space':
            RA = '{0}{1} {2} {3}'.format(rs, raH, raM, round(raS, 2))
        elif format_out==':':
            RA = '{0}{1}:{2}:{3}'.format(rs, raH, raM, round(raS, 2)) 
        else:
            RA = '{0}{1}h{2}m{3}s'.format(rs, raH, raM, round(raS, 2))
  if ra and dec:
    return (RA, DEC)
  else:
    return RA or DEC



def sexa_to_deci(ra='', dec='', format_in='space'):
  RA, DEC, rs, ds = '', '', 1, 1
  if dec:
    if format_in=='space':
        D, M, S = [float(i) for i in dec.split()]
    elif format_in==':':
        D, M, S = [float(i) for i in dec.split(':')]        
    else:
        D = float(dec.split('d')[0])
        M = float(dec.split('d')[-1].split('m')[0])
        S = float(dec.split('m')[-1].split('s')[0])
        
    
    if str(D)[0] == '-':
      ds, D = -1, abs(D)
    deg = D + (M/60) + (S/3600)
    DEC = '{0}'.format(deg*ds)
  
  if ra:
    if format_in=='space':
        H, M, S = [float(i) for i in ra.split()]
    elif format_in==':':    
        H, M, S = [float(i) for i in ra.split(':')]
    else:
        H = float(ra.split('h')[0])
        M = float(ra.split('h')[-1].split('m')[0])
        S = float(ra.split('m')[-1].split('s')[0])
        
    if str(H)[0] == '-':
      rs, H = -1, abs(H)
    deg = (H*15) + (M/4) + (S/240)
    RA = '{0}'.format(deg*rs)
  
  if ra and dec:
    return (float(RA), float(DEC))
  else:
    return float(RA) or float(DEC)

'''
ra,dec = sexa_to_deci(ra='00h41m03.44s',dec='-09d56m28.1s',format_in='h')
print ra,dec
ra,dec = deci_to_sexa(ra=ra, dec=dec,format_out='h')
print ra,dec
'''