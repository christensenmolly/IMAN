#!/usr/bin/python
# DESCRIPTION:
# Script to convert region mask to fits mask.
# Pixels with counts >0 are a mask in the output mask image.
# MINIMAL USAGE: python convert_reg_to_mask.py [input_image] [region_file]

# Import the necessary modules
from astropy.io import fits as pyfits
#import pyfits
import numpy as np
import math
import itertools
from scipy import ndimage
import sys
from itertools import product
from matplotlib.path import Path
from math import hypot, cos, sin, radians, pi
from numpy import linspace, sign, zeros_like
import shutil
import argparse
import os
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


def show_complete(i,N):
  percent = 100. * float(i) / float(N)
  sys.stdout.write("\r%2d%%" % percent)
  sys.stdout.flush()

class PPoint:
    def __init__(self, x, y):
        self.x = x
        self.y = y


def rot_point(p, orig, angle):
    """ performs rotation of a point (x,y) around point (x0,y0)
        by angle radians"""
    x1 = cos(angle)*(p.x-orig.x) - sin(angle)*(p.y-orig.y) + orig.x
    y1 = sin(angle)*(p.x-orig.x) + cos(angle)*(p.y-orig.y) + orig.y
    return PPoint(x1, y1)

def listRightIndex(alist, value):
    return len(alist) - alist[-1::-1].index(value) -1


def ellipse_mask(cen, ellA, ellB, ellPA, xSize, ySize, img_inp=None, img_mask=None, outer = False, mask_value=1, value=0.):
        cospa = cos(radians(ellPA))
        sinpa = sin(radians(ellPA))
        # Check if whole ellipse is inside of the image
        # and obtain size of the ellipse in xy plane
        xMax = 0.0
        yMax = 0.0
        xMin = 1e10
        yMin = 1e10
        for e in linspace(0, 4*pi, 1000):
            cose = cos(e)
            sine = sin(e)
            x = cen.x + ellA * cose * cospa - ellB * sine * sinpa
            y = cen.y + ellB * sine * cospa + ellA * cose * sinpa
            if x > xMax:
                xMax = x
            if y > yMax:
                yMax = y
            if x < xMin:
                xMin = x
            if y < yMin:
                yMin = y
        xMin = max(0, int(round(xMin)))
        xMax = min(xSize, int(round(xMax)))
        yMin = max(0, int(round(yMin)))
        yMax = min(ySize, int(round(yMax)))
        focusR = (ellA ** 2.0 - ellB ** 2.0) ** 0.5
        focus10 = PPoint(cen.x + focusR, cen.y)  # Unrotated
        focus20 = PPoint(cen.x - focusR, cen.y)  #
        focus1 = rot_point(focus10, cen, radians(ellPA))
        focus2 = rot_point(focus20, cen, radians(ellPA))
        # Find pixels inside of the ellipse
        dEll = 2 * ellA
        
        if outer==True:
          xMin = 0
          xMax = xSize -1
          yMin = 0
          yMax = ySize -1
  
          if img_inp==None:
            for x in range(xMin, xMax+1):
               for y in range(yMin, yMax+1):
                    dFocus1 = hypot(x-focus1.x, y-focus1.y)
                    dFocus2 = hypot(x-focus2.x, y-focus2.y)
                    dPoint = dFocus1 + dFocus2
                    if dPoint >= dEll:
                      try:
                        #inframe[y-1,x-1] = value
                        img_mask[y-1,x-1] = mask_value
                      except:
                        zz=1
                    else:
                      try:
                        #inframe[y-1,x-1] = inframe[y-1,x-1]
                        img_mask[y-1,x-1] = 0
                      except:
                        zz=1      
          else:
            for x in range(xMin, xMax+1):
               for y in range(yMin, yMax+1):
                    dFocus1 = hypot(x-focus1.x, y-focus1.y)
                    dFocus2 = hypot(x-focus2.x, y-focus2.y)
                    dPoint = dFocus1 + dFocus2
                    if dPoint >= dEll:
                      try:
                        #inframe[y-1,x-1] = img_inp[y-1,x-1]
                        z=1
                      except:
                        zz=1

        else:
          if img_inp==None:
            for x in range(xMin, xMax+1):
               for y in range(yMin, yMax+1):
                    dFocus1 = hypot(x-focus1.x, y-focus1.y)
                    dFocus2 = hypot(x-focus2.x, y-focus2.y)
                    dPoint = dFocus1 + dFocus2
                    if dPoint < dEll:
                      try:
                         #inframe[y-1,x-1] = value
                         img_mask[y-1,x-1] = mask_value
                      except:
                         zz=1
          else:
            for x in range(xMin, xMax+1):
               for y in range(yMin, yMax+1):
                    dFocus1 = hypot(x-focus1.x, y-focus1.y)
                    dFocus2 = hypot(x-focus2.x, y-focus2.y)
                    dPoint = dFocus1 + dFocus2
                    if dPoint < dEll:
                      try:
                        #inframe[y-1,x-1] = img_inp[y-1,x-1]
                        z=1
                      except:
                        zz=1
        return img_mask

def ellipse_annulus_mask(inner_ellipse, outer_ellipse, mask):
        ySize,xSize = np.shape(mask)
        MASK = np.ones(shape=(ySize,xSize))
        def define_ellipse(cen, ellA, ellB, ellPA):
            cospa = cos(radians(ellPA))
            sinpa = sin(radians(ellPA))
            # Check if whole ellipse is inside of the image
            # and obtain size of the ellipse in xy plane
            xMax = 0.0
            yMax = 0.0
            xMin = 1e10
            yMin = 1e10
            for e in linspace(0, 4*pi, 1000):
                cose = cos(e)
                sine = sin(e)
                x = cen.x + ellA * cose * cospa - ellB * sine * sinpa
                y = cen.y + ellB * sine * cospa + ellA * cose * sinpa
                if x > xMax:
                   xMax = x
                if y > yMax:
                   yMax = y
                if x < xMin:
                   xMin = x
                if y < yMin:
                   yMin = y
            xMin = max(0, int(round(xMin)))
            xMax = min(xSize, int(round(xMax)))
            yMin = max(0, int(round(yMin)))
            yMax = min(ySize, int(round(yMax)))

            focusR = (ellA ** 2.0 - ellB ** 2.0) ** 0.5
            focus10 = Point(cen.x + focusR, cen.y)  # Unrotated
            focus20 = Point(cen.x - focusR, cen.y)  #
            focus1 = rot_point(focus10, cen, radians(ellPA))
            focus2 = rot_point(focus20, cen, radians(ellPA))
            return xMin,xMax,yMin,yMax,focus1,focus2

        cen_in, ellA_in, ellB_in, ellPA_in = inner_ellipse
        cen_out, ellA_out, ellB_out, ellPA_out = outer_ellipse

        ## Inner ellipse
        xMin_in,xMax_in,yMin_in,yMax_in,focus1_in,focus2_in = define_ellipse(cen_in, ellA_in, ellB_in, ellPA_in)

        ## Outer ellipse
        xMin_out,xMax_out,yMin_out,yMax_out,focus1_out,focus2_out = define_ellipse(cen_out, ellA_out, ellB_out, ellPA_out)

        xMin = min([xMin_in,xMin_out])
        yMin = min([yMin_in,yMin_out])

        xMax = max([xMax_in,xMax_out])
        yMax = max([yMax_in,yMax_out])

        # Find pixels inside of the ellipse annulus
        dEll_in = 2 * ellA_in
        dEll_out = 2 * ellA_out
        

        for x in range(xMin, xMax):
            for y in range(yMin, yMax):
                dFocus1_in = hypot(x-focus1_in.x, y-focus1_in.y)
                dFocus2_in = hypot(x-focus2_in.x, y-focus2_in.y)
                dPoint_in = dFocus1_in + dFocus2_in

                dFocus1_out = hypot(x-focus1_out.x, y-focus1_out.y)
                dFocus2_out = hypot(x-focus2_out.x, y-focus2_out.y)
                dPoint_out = dFocus1_out + dFocus2_out                
                if dPoint_out < dEll_out and dPoint_in > dEll_in and mask[y,x]==0.:
                  try:
                    MASK[y,x] = 0.
                  except:
                    zz=1
        return MASK



def mask1(input_image, reg_file, output_image=None, output_mask=None, mask_value=1, show_running=True, mask_DN=None):
                if show_running==True:
                    print('Masking contaminants...')
                if output_image is None:
                    output_image = input_image.split('.fits')[0]+'_clean.fits'  
                    
                if output_mask is None:
                    output_image = reg_file.split('.reg')[0]+'.fits'  

                shutil.copy(input_image,output_image)



                if output_mask!=None:
                    shutil.copy(input_image,output_mask)
                    hdulist2 = pyfits.open(output_mask, do_not_scale_image_data=True, mode='update')
                    img2 = hdulist2[0].data
  
                hdulist3 = pyfits.open(output_image, do_not_scale_image_data=True, mode='update')
                img3 = hdulist3[0].data


                (dimy,dimx) = img3.shape
                value = 0#float('nan')
                paths=[]
                borders=[]
                f = open(reg_file, "r")
                all_lines = f.readlines()
                Number_of_lines = len(all_lines)
                f.close()
                f = open(reg_file, "r")

                coun = 0
                if True:
                    for line in f:
                            # Detect circle regions and fill them with the mask value
                            if 'circle(' in line:
                                params = line.split("(")[-1].split(",")[0]
                                cen = PPoint(float(line.split("(")[-1].split(",")[0]),
						  float(line.split("(")[-1].split(",")[1]))
                                ellA = float(line.split("(")[-1].split(",")[2].split(')')[0])
                                ellB = ellA
                                ellPA = 0.
                                ellipse_mask(cen, ellA, ellB, ellPA, img3, dimx, dimy, mask_value=mask_value)
    
                            # Detect polygon regions and add them to the path list
                            elif 'polygon(' in line:
                                param = map(float, ((line.split('polygon(')[1]).split(')')[0]).split(','))
                                param = list(param)
                                param2 = [None]*(int(len(param)/2))
                                for i in range(0,len(param2)):  param2[i] = (int(param[2*i]),int(param[2*i+1])) 
                                param2.append(param2[0])
                                codes = []
                                codes.append(Path.MOVETO)
                                for i in range(1,len(param2)-1):  codes.append(Path.LINETO)
                                codes.append(Path.CLOSEPOLY)
                                path = Path(param2, codes)
                                paths.append(path)
                                coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
                                X = []; Y = []
                                for k in range(0,len(coords)-1,2):
                                   X.append(int(float(coords[k])))
                                   Y.append(int(float(coords[k+1])))

                                borders.append([min(X),min(Y),max(X),max(Y)])
    
                            elif "ellipse" in line:
                                params = line.split(",")
                                cen = PPoint(float(params[0].split('(')[1]),
						  float(params[1]))
                                ellA = float(params[2])
                                ellB = float(params[3])
                                ellPA = float(params[4].split(')')[0])
                                if ellA < ellB:
                                    ellA, ellB = ellB, ellA
                                    ellPA += 90
   
                                if 'text={outer}' in line:
                                    outer = True
                                else:
                                    outer = False
                                ellipse_mask(cen, ellA, ellB, ellPA, img3, dimx, dimy, outer=outer, mask_value=mask_value)
                            coun = coun + 1

                            if show_running==True:
                                show_complete(coun,Number_of_lines)


                if True:
                      # Loop over the image and all polygons and fill them with the mask value	
                      nx, ny =img3.shape[1], img3.shape[0]
   
                      for kk in range(len(paths)):
                        path = paths[kk]
                        border = borders[kk]
                        xmin = border[0]
                        ymin = border[1]
                        xmax = border[2]
                        ymax = border[3]
                        for i, j in product(range(xmin,xmax), range(ymin,ymax)):
                          try:
                            img3[j-1][i-1]=value if path.contains_point((i,j)) else img3[j-1][i-1]
                          except:
                            zz=1
                      hdulist3.flush()

                      if output_mask!=None:
                        if mask_DN==None:
                          for i in range(dimy):
                                for k in range(dimx):
                                    if img3[i,k] == value:
                                       img2[i,k] = mask_value
                                    else:
                                       img2[i,k] = 0
                        else:
                          for i in range(dimy):
                             for k in range(dimx):
                                    if img3[i,k] == value:
                                       img2[i,k] = mask_value
                                    else:
                                       img2[i,k] = 0
                                    
                                    if img3[i,k] == mask_DN and mask_DN!=value:
                                       img2[i,k] = mask_value
                                       img3[i,k] = value

                if show_running==True:
                    print('\nDone!')
                return output_image,output_mask


def mask(input_image, reg_file, output_image=None, output_mask=None, mask_value=1, show_running=True, mask_DN=None):

                if show_running==True:
                    print('Masking contaminants...')
                
                #if output_image is None:
                #    remove_clean_image = True
                #else:
                #    remove_clean_image = False
                
                #if output_image is None:
                #    output_image = input_image.split('.fits')[0]+'_clean.fits'  
                    
                if output_mask is None:
                    output_mask = reg_file.split('.reg')[0]+'.fits'  

                #shutil.copy(input_image, output_image)
                
                hdulist = pyfits.open(input_image)
                img = hdulist[0].data
                header = hdulist[0].header
                hdulist.close()
                
 
  
                #hdulist3 = pyfits.open(output_image, do_not_scale_image_data=True, ignore_missing_end=True)
                #img3 = hdulist3[0].data
                #header = hdulist3[0].header
                


                (dimy,dimx) = img.shape
                img2 = np.zeros(np.shape(img))

                value = 0#float('nan')
                paths=[]
                borders=[]
                f = open(reg_file, "r")
                all_lines = f.readlines()
                Number_of_lines = len(all_lines)
                f.close()
                f = open(reg_file, "r")

                coun = 0
                print('here1')
                if True:
                    for line in f:
                            # Detect circle regions and fill them with the mask value
                            if 'circle(' in line:
                                params = line.split("(")[-1].split(",")[0]
                                cen = PPoint(float(line.split("(")[-1].split(",")[0]),
						  float(line.split("(")[-1].split(",")[1]))
                                ellA = float(line.split("(")[-1].split(",")[2].split(')')[0])
                                ellB = ellA
                                ellPA = 0.
                                img2 = ellipse_mask(cen, ellA, ellB, ellPA, dimx, dimy, img_mask=img2, mask_value=mask_value, value=value)
    
                            # Detect polygon regions and add them to the path list
                            elif 'polygon(' in line:
                                param = map(float, ((line.split('polygon(')[1]).split(')')[0]).split(','))
                                param = list(param)
                                #print(len(param))
                                #exit()
                                param2 = [None]*(int(len(param)/2))
                                for i in range(0,len(param2)):  param2[i] = (int(param[2*i]),int(param[2*i+1])) 
                                param2.append(param2[0])
                                codes = []
                                codes.append(Path.MOVETO)
                                for i in range(1,len(param2)-1):  codes.append(Path.LINETO)
                                codes.append(Path.CLOSEPOLY)
                                path = Path(param2, codes)
                                paths.append(path)
                                coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
                                X = []; Y = []
                                for k in range(0,len(coords)-1,2):
                                   X.append(int(float(coords[k])))
                                   Y.append(int(float(coords[k+1])))

                                borders.append([min(X),min(Y),max(X),max(Y)])
    
                            elif "ellipse" in line:
                                params = line.split(",")
                                cen = PPoint(float(params[0].split('(')[1]),
						  float(params[1]))
                                ellA = float(params[2])
                                ellB = float(params[3])
                                ellPA = float(params[4].split(')')[0])
                                if ellA < ellB:
                                    ellA, ellB = ellB, ellA
                                    ellPA += 90
   
                                if 'text={outer}' in line:
                                    outer = True
                                else:
                                    outer = False
                                img2 = ellipse_mask(cen, ellA, ellB, ellPA, dimx, dimy, outer=outer, img_mask=img2, mask_value=mask_value, value=value)
                            coun = coun + 1

                            if show_running==True:
                                show_complete(coun,Number_of_lines)

                print('here2')
                Number_of_lines = len(paths)
                coun = 0
                #exit()
                if True:
                      # Loop over the image and all polygons and fill them with the mask value	 
                      for kk in range(len(paths)):
                        path = paths[kk]
                        border = borders[kk]
                        xmin = border[0]
                        ymin = border[1]
                        xmax = border[2]
                        ymax = border[3]
                        for i, j in product(range(xmin,xmax), range(ymin,ymax)):
                          try:
                            if path.contains_point((i,j)):
                                #img3[j-1][i-1] = value
                                img2[j-1][i-1] = mask_value
                          except:
                            zz=1
                        coun = coun + 1
                        if show_running==True:
                                show_complete(coun, Number_of_lines)
                      '''
                      if mask_DN is not None:
                          mask_DN = float(mask_DN)
                          if np.isnan(mask_DN): 
                            for i in range(dimy):
                                for k in range(dimx):                                   
                                        if np.isnan(img[i,k]):
                                            img3[i,k] = value  
                                            img2[i,k] = mask_value
                          else:
                            for i in range(dimy):
                                for k in range(dimx):                                   
                                        if img[i,k] == mask_DN:
                                            img3[i,k] = value  
                                            img2[i,k] = mask_value                     
                      '''
                      
                      
                      #hdulist3.flush()
                      #hdulist3.close()
                      outHDU = pyfits.PrimaryHDU(img2, header=header)
                      outHDU.writeto(output_mask,overwrite=True)                      

                      #outHDU = pyfits.PrimaryHDU(img3, header=header)
                      #outHDU.writeto(output_image,overwrite=True)      
                
                
                #if remove_clean_image:
                #    os.remove(output_image)
                    
                if show_running==True:
                    print('\nDone!')
                    f.close()
                print('here3')
                return output_image,output_mask

            
            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert region file with a mask to a fits segmentation image")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("reg_file", help="Input region file with a mask")
    parser.add_argument("--output_image", nargs='?', const=1, help="Optional: Output image with empty (0) masked pixels", type=str, default=None) 
    parser.add_argument("--output_mask", nargs='?', const=1, help="Optional: Output mask with masked values >0", type=str, default=None) 
    parser.add_argument("--mask_value", nargs='?', const=1, help="Optional: Masking value in output_mask", type=int, default=1) 
    parser.add_argument("--mask_DN", nargs='?', const=1, help="Optional: Masked value in input_image", type=float, default=None) 
    args = parser.parse_args()

    input_image = args.input_image
    reg_file = args.reg_file
    output_image = args.output_image
    output_mask = args.output_mask   
    mask_value = args.mask_value
    mask_DN = args.mask_DN

    mask(input_image, reg_file, output_image=output_image, output_mask=output_mask, mask_value=mask_value, show_running=True, mask_DN=mask_DN)
