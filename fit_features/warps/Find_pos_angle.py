from skimage.morphology import skeletonize
from skimage import draw
import numpy as np
import matplotlib.pyplot as plt
import pyfits
from math import *
import sys
import os
import subprocess
# an empty image


import crea_mask
import crea_iso

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0].split('Warps')[0]
sys.path.append(PATH_TO_PACKAGE + '/IMP_NEW/Sergey_pipelinelib')
sys.path.append(PATH_TO_PACKAGE+'/FindFluxes')

from rotima import crotima
import superell_fit
import ds9_contour

def crop_image(input_image,line):
  xx = []
  yy = []
  if 'polygon' in line:
	    coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
	    for k in range(0,len(coords)-1,2):
	      xx.append(float(coords[k]))
	      yy.append(float(coords[k+1]))  
  xl = int(ceil(min(xx)))
  xr = int(floor(max(xx)))
  yl = int(ceil(min(yy)))
  yr = int(floor(max(yy)))
  output_image = 'crop.fits'
  if os.path.exists(output_image):
    os.remove('crop.fits')
  
  command = 'fitscopy \'' + str(input_image)+'['+str(xl)+':'+str(xr)+','+str(yl)+':'+str(yr)+']\' '+str(output_image)
  subprocess.call(command, shell=True)
  
  line = 'polygon('
  for k in range(len(xx)):
    xx[k] = xx[k]-xl
    yy[k] = yy[k]-yl
    if k!=len(xx)-1:
      line = line + str(xx[k])+','+str(yy[k])+','
    else:
      line = line + str(xx[k])+','+str(yy[k])
  line = line + ') # color=black'
  #print line
  return line,xl,yl

def find_PA(input_image,line,std_delta_y_mean):
  new_line,xl,yl = crop_image(input_image,line)
  image_new = crea_mask.mask('crop.fits',new_line)
  hdulist = pyfits.open('crop.fits')
  image = hdulist[0].data
  nx, ny =image.shape[1], image.shape[0]


  # perform skeletonization
  skeleton = skeletonize(image_new)

  x = []
  y = []

  for k in range(ny):
    for i in range(nx):
	  if skeleton[k,i]==True:
	    x.append(i)
	    y.append(k)

  ZZ = np.polyfit(x, y, 3)
  ff = np.poly1d(ZZ)
  #plt.plot(x,ff(x),'o',color='green')
  plt.plot(x,y,'*',color='green')
  
  x_new = []
  y_new = []
  y = np.array(y)
  x = np.array(x)
  delta_y = np.fabs(y-ff(x))
  std_delta_y = np.std(delta_y)
  if std_delta_y_mean!=0.:
    Std_delta_y = std_delta_y_mean
  else:
    Std_delta_y = std_delta_y
  for k in range(len(x)):
    if delta_y[k]<3.*Std_delta_y:
	    x_new.append(x[k])
	    y_new.append(y[k])

  ZZ = np.polyfit(x_new, y_new, 1)
  ff = np.poly1d(ZZ)
  
  PA = degrees(atan(ZZ[0]))
  

  '''
  xx = []
  yy = []
  if 'polygon' in new_line:
	    coords = new_line.replace('(',',').replace(')',',').split(',')[1:-1]
	    for k in range(0,len(coords)-1,2):
	      xx.append(float(coords[k]))
	      yy.append(float(coords[k+1]))  
  
  plt.plot(x_new,y_new,'s',markersize=8,color='blue')
  plt.plot(x,ff(x),'o',color='cyan')
  plt.plot(xx,yy,'o',color='red')
  plt.show()
  '''
  os.remove('crop.fits')  
  return PA


def main(input_image,m0,pix2sec,xc,yc,inner_level,outer_level):
  isophote_file = crea_iso.main(input_image,m0,pix2sec,xc,yc,inner_level,outer_level)
  
  reg = open(isophote_file,'r')
  for line in reg:
    if 'polygon' in line:
      PA = find_PA(input_image,line,0.)
      print 'PA=%.1f deg' % (PA)
      break


  
  # Retrieve outer isophote
  with open(isophote_file) as f:
    last = None
    for line in (line for line in f if line.rstrip('\n')):
      last = line

  x = []
  y = []
  if 'polygon' in last:
    coords = last.replace('(',',').replace(')',',').split(',')[1:-1]
    for k in range(0,len(coords)-1,2):
      x.append(float(coords[k]))
      y.append(float(coords[k+1]))
      
  x = np.array(x)
  y = np.array(y)
  
  a = superell_fit.fitEllipse(x,y)
  center = ds9_contour.ellipse_center(a)
  axes = ds9_contour.ellipse_axis_length(a)
  phi = degrees(ds9_contour.ellipse_angle_of_rotation(a))
  xc = center[0]
  yc = center[1]
  sma = axes[0]
  smb = axes[1]
  phi = 90. + phi

  Sma = sma
  if sma<smb:
    sma = smb
    smb = Sma
  
  #x_l = xc - sma
  #y_l = yc - smb
  #x_r = xc + sma
  #y_r = yc + smb
  #outFitsName = input_image.split('/')[-1].split('.fits')[0] + '_rot.fits'
  #xCenRot, yCenRot = crotima(input_image,outFitsName,xc,yc,PA)
  
  
  return PA,xc,yc,sma,smb

#main('frame-g-003813-1-0201.fits',22.5,0.396,1207.967800,1163.702300,23.5,24.5)
#input_image = '/home/amosenko/CurrentWork/S4G_eon/sample/data/92_411/filled_new.fits'
#main(input_image,21.097,0.75,224,73,22.5,22.5)