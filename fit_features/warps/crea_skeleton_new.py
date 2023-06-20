from skimage.morphology import skeletonize
from skimage import draw
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from math import *
import sys
import os
import subprocess
# an empty image
from PIL import Image
import cv2
import crea_mask


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


def convert_to_png(img_old):
    from skimage import img_as_bool, io, color, morphology
    img = np.zeros((img_old.shape[0], img_old.shape[1], 3), dtype=float)
    img[:,:,0] = img_old#img_scale.sqrt(img_old, scale_min=0, scale_max=10000)
    img[:,:,1] = img_old#img_scale.sqrt(img_old, scale_min=0, scale_max=10000)
    img[:,:,2] = img_old#img_scale.sqrt(img_old, scale_min=0, scale_max=10000)
    plt.clf()
    plt.imshow(img, aspect='equal')
    plt.savefig('test.png')
    img = img_as_bool(color.rgb2gray(io.imread('test.png')))
    #img = img//255
    return img
    
    

def main(input_image,line,std_delta_y_mean):
  from skimage import img_as_bool, io, color, morphology
  new_line,xl,yl = crop_image(input_image,line)
  image_new = crea_mask.mask('crop.fits',new_line)
  hdulist = pyfits.open('crop.fits')
  image = hdulist[0].data
  nx, ny =image.shape[1], image.shape[0]
  
  #image_new = convert_to_png(image_new)



  # perform skeletonization
  skeleton = skeletonize(image_new)
  #''' 
  x = []
  y = []

  for k in range(ny):
    for i in range(nx):
        if skeleton[k,i]==True:
            x.append(i)
            y.append(k)
  #print(x,y)
  #exit()
  #'''

  ZZ = np.polyfit(x, y, 3)
  ff = np.poly1d(ZZ)
  plt.plot(x,ff(x),'o',color='green')
  #plt.plot(x,y,'*',color='green')
  
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

  ZZ = np.polyfit(x_new, y_new, 3)
  ff = np.poly1d(ZZ)

  #'''
  xx = []
  yy = []
  if 'polygon' in new_line:
            coords = new_line.replace('(',',').replace(')',',').split(',')[1:-1]
            for k in range(0,len(coords)-1,2):
                xx.append(float(coords[k]))
                yy.append(float(coords[k+1]))  
  
  plt.plot(x,y,'s',markersize=8,color='blue')
  #plt.plot(x,ff(x),'o',color='cyan')
  plt.plot(xx,yy,'o',color='red')
  plt.show()
  #'''
  os.remove('crop.fits')
  x = np.array(x)
  #return x +float(xl),ff(x)+float(yl),std_delta_y
  return x +float(xl),y+float(yl),std_delta_y,x +float(xl),ff(x)+float(yl)


