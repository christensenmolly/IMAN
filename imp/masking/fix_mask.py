# Import the necessary modules
from astropy.io import fits as pyfits
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
import subprocess

import merge_masks
import convert_reg_to_mask
import convert_segm_to_region
import read_mto_output
import auto_masking
import shapely
r_target = 1  #### Change this if the target is masked out

def ds9_to_np(x):
    '''
    Function to convert ds9 pixel coordinates to numpy format, i.e
    0.5 pix is 0, 1.5 pix is 1
    '''
    return int(round(x)) - 1


def check_target_object(xc, yc, segm, mto_csv):
    ID,X,Y,A,B,theta,total_flux,mu_max,mu_median,mu_mean,R_fwhm,R_e,R10,R90,RA,DEC = read_mto_output.main(mto_csv, input_image=None)
    ID = np.array(ID, dtype=int)
    ID = list(ID)
    
    ind = ID.index(int(segm[yc,xc]))
    
    print(segm[yc,xc], xc, yc, ID[ind], X[ind], Y[ind])
    return ID[ind],X[ind],Y[ind]



def convert_polygonline_to_shapleyregion(line):
    coords = line.replace('(',',').replace(')',',').split(',')[1:-1]
    pol = []
    for kk in range(0,len(coords)-1,2):
        pol.append((int(float(coords[kk])),int(float(coords[kk+1]))))
    polygon = Polygon(pol)
    return polygon


def fix_galaxy_region1(input_image, galaxy_region_file, mask_region_file):
    f = open(galaxy_region_file, 'r')
    lines = f.readlines()
    new_lines = []
    point_numbers = []
    for line in lines:
        if 'polygon(' in line:
            point_numbers.append(len(line))
            new_lines.append(line)
    f.close()
    
    ff = open('region_tmp.reg', 'w')
    ff.write('image\n')
    ff.write(new_lines[point_numbers.index(max(point_numbers))])
    ff.close()
    shutil.move('region_tmp.reg',galaxy_region_file)
    
    print(galaxy_region_file)
    subprocess.call("ds9 %s -scale histequ -regions %s" % (input_image, galaxy_region_file), shell=True)
    
    
    
    galaxy_polygon = convert_polygonline_to_shapleyregion(new_lines[point_numbers.index(max(point_numbers))])
    
    f = open(mask_region_file, "r")
    ff = open('inner_mask_mto.reg', 'w')
    ff.write('image\n')
    for line in f:
            if 'polygon(' in line:
                try:
                    polygon = convert_polygonline_to_shapleyregion(line)
                    
                    if galaxy_polygon.intersects(polygon):
                        ff.write(line)
                except:
                    zz=1
    ff.close()
    
    
    
    convert_segm_to_region.do_offseting('inner_mask_mto.reg', 'inner_mask_mto_tmp.reg', offset_size=1.3, offset_pix=0., xc=None, yc=None, system='image')
    shutil.move('inner_mask_mto_tmp.reg', 'inner_mask_mto.reg')               
    

    
def fix_galaxy_region2(input_image, galaxy_region_file, mask_region_file):
    subprocess.call("ds9 %s -scale histequ -regions %s" % (input_image, galaxy_region_file), shell=True)
    
    f = open(galaxy_region_file, 'r')
    f_reg = open(mask_region_file, "r")

    ff = open('inner_mask_mto.reg', 'w')
    ff.write('image\n')    
    
    lines = f.readlines()
    new_lines = []
    point_numbers = []
    for line in lines:
        if 'polygon(' in line:
            point_numbers.append(len(line))
            try:
                galaxy_polygon = convert_polygonline_to_shapleyregion(line)
                for line1 in f_reg:
                        if 'polygon(' in line1:
                            try:
                                polygon = convert_polygonline_to_shapleyregion(line1)
                                
                                if galaxy_polygon.intersects(polygon):
                                    ff.write(line1)
                            except:
                                zz=1
                            
            except:
                zz=1
    

    ff.close()
    f.close()
    f_reg.close()
    
    
    
    
    convert_segm_to_region.do_offseting('inner_mask_mto.reg', 'inner_mask_mto_tmp.reg', offset_size=1.3, offset_pix=0., xc=None, yc=None, system='image')
    shutil.move('inner_mask_mto_tmp.reg', 'inner_mask_mto.reg')          
    

def fix_galaxy_region(input_image, galaxy_region_file, mask_region_file, output_region='inner_mask_mto.reg'):
    if True:
            xc,yc,SMA,SMB,PA = auto_masking.read_region(galaxy_region_file)

            circle = Point(xc,yc).buffer(1)
            if SMA is not None:
                if SMB is not None and PA is not None:
                    ellipse = shapely.affinity.scale(circle,SMA,SMB)
                    rot_ellipse = shapely.affinity.rotate(ellipse, PA, origin='center', use_radians=False)
                else:
                    ellipse = shapely.affinity.scale(circle,SMA,SMA)



    #subprocess.call("ds9 %s -scale histequ -regions %s" % (input_image, galaxy_region_file), shell=True)
    
    f_reg = open(mask_region_file, "r")

    ff = open('inner_mask_mto.reg', 'w')
    ff.write('image\n')    
    

    for line1 in f_reg:
                        if 'polygon(' in line1:
                            try:
                                polygon = convert_polygonline_to_shapleyregion(line1)
                                
                                if rot_ellipse.intersects(polygon) and not polygon.contains(rot_ellipse):
                                    ff.write(line1)
                            except:
                                zz=1


    ff.close()
    f_reg.close()
    
    
    
    
    convert_segm_to_region.do_offseting(output_region, 'inner_mask_mto_tmp.reg', offset_size=1.2, offset_pix=0., xc=None, yc=None, system='image')
    shutil.move('inner_mask_mto_tmp.reg', output_region)              
    
    


def main_repeat(input_image, segm_mto_image, fix_region_file, new_segm_image='new_segm.fits', clean_image=None, mto_csv='parameters.csv', unmask_value=-1, verbosity=True):
    # Open input image
    hdulist = pyfits.open(input_image)
    img = hdulist[0].data 
    header = hdulist[0].header 
    
    # Open MTO segmentation image

    hdulist1 = pyfits.open(segm_mto_image)
    segm_mto = hdulist1[0].data
    
    new_segm_mto = np.copy(segm_mto)

    # Open region file with marked objects to unmask (as 'X') and regions to mask (polygons, circles and ellipses)
    # Here we only read in X-regions:
    f = open(fix_region_file, "r")
    xc = []; yc = []
    for line in f:
        if '# point=x' in line:
            xc.append(ds9_to_np(float(line.split('(')[1].split(',')[0])))
            yc.append(ds9_to_np(float(line.split(',')[1].split(')')[0])))
    f.close()
    

    

    if xc!=[]:
        # Here we get IDs of these X-marked objects and unmask them
        for k in range(len(xc)):
            try:
                if yc[k] > -1 and xc[k] > -1:
                    new_segm_mto[segm_mto == segm_mto[yc[k],xc[k]]] = unmask_value
                else:
                    zz = 1
            except IndexError as e:
                print(e)
                print("k == {}".format(k))
                print("yc[k] == {}".format(yc[k]))
                print("xc[k] == {}".format(xc[k]))


    
    # Convert -1 to 0 (usual unmasked value):
    new_segm_mto[new_segm_mto == unmask_value] = 0

    convert_reg_to_mask.mask(input_image, fix_region_file, output_image=None, output_mask='tmp_mask.fits', mask_value=1000.,mask_DN=None, verbosity=verbosity) # Check this line!!!


    if os.path.exists('tmp_mask.fits'):
        hdulist2 = pyfits.open('tmp_mask.fits')
        manual_mask = hdulist2[0].data    
        final_mask = new_segm_mto + manual_mask
        hdulist2.close()
    else:
        final_mask = new_segm_mto
    

    outHDU = pyfits.PrimaryHDU(final_mask, header=header)
    outHDU.writeto(new_segm_image, clobber=True)   
    
    if os.path.exists('tmp_mask.fits'):
        os.remove('tmp_mask.fits')
    

    
    if clean_image is not None:
        ny,nx = np.shape(final_mask)
        for k in range(ny):
            for i in range(nx):
                if final_mask[k,i]>0:
                    img[k,i]=0.

        outHDU = pyfits.PrimaryHDU(img, header=header)
        outHDU.writeto(clean_image, clobber=True)     
    
    

    hdulist1.close()    
    hdulist.close()    
    





def main(input_image, segm_mto, fix_region_file, galaxy_region, new_segm_image='new_segm.fits', clean_image=None, mto_csv='parameters.csv', unmask_value=-1, verbosity=True):
    # Open input image
    hdulist = pyfits.open(input_image)
    img = hdulist[0].data 
    header = hdulist[0].header 
    
    # Open MTO segmentation image
    hdulist1 = pyfits.open(segm_mto)
    segm_mto = hdulist1[0].data
    
    new_segm_mto = np.copy(segm_mto)

    # Open region file with marked objects to unmask (as 'X') and regions to mask (polygons, circles and ellipses)
    # Here we only read in X-regions:
    f = open(fix_region_file, "r")
    xc = []; yc = []
    for line in f:
        if '# point=x' in line:
            xc.append(ds9_to_np(float(line.split('(')[1].split(',')[0])))
            yc.append(ds9_to_np(float(line.split(',')[1].split(')')[0])))
    f.close()
    

    
    if xc!=[]:
        # Here we get IDs of these X-marked objects and unmask them
        for k in range(len(xc)):
            if True:
                new_segm_mto[segm_mto == segm_mto[yc[k],xc[k]]] = unmask_value
            else:
                zz=1

    
    # Convert -1 to 0 (usual unmasked value):
    new_segm_mto[new_segm_mto == unmask_value] = 0

    
    # Create galaxy region:
    if True:
        maskHDU = pyfits.PrimaryHDU(new_segm_mto, header=header)
        maskHDU.writeto('mask_segm_mto.fits', clobber=True)         

        
        convert_segm_to_region.main('mask_segm_mto.fits', 'mask_segm_mto.reg', output_mask_image=None, fits_slice = 0, offset_size=1.0, offset_pix=0., xc=None, yc=None, system='image', ignore_value=None, verbosity=verbosity)        
       
        fix_galaxy_region(input_image, galaxy_region, 'mask_segm_mto.reg', output_region='inner_mask_mto.reg')

        
    
    merge_masks.main([fix_region_file,'inner_mask_mto.reg'], 'tmp.reg', verbosity=verbosity)
    shutil.move('tmp.reg',fix_region_file)
    
    convert_reg_to_mask.mask(input_image, fix_region_file, output_image=None, output_mask='tmp_mask.fits', mask_value=1000., mask_DN=None, verbosity=verbosity)

    if os.path.exists('tmp_mask.fits'):
        hdulist2 = pyfits.open('tmp_mask.fits')
        manual_mask = hdulist2[0].data    
        final_mask = new_segm_mto + manual_mask
        hdulist2.close()
    else:
        final_mask = new_segm_mto
    

    outHDU = pyfits.PrimaryHDU(final_mask, header=header)
    outHDU.writeto(new_segm_image, clobber=True)   
    
    if os.path.exists('tmp_mask.fits'):
        os.remove('tmp_mask.fits')
    

    
    if clean_image is not None:
        ny,nx = np.shape(final_mask)
        for k in range(ny):
            for i in range(nx):
                if final_mask[k,i]>0:
                    img[k,i]=0.

        outHDU = pyfits.PrimaryHDU(img, header=header)
        outHDU.writeto(clean_image, clobber=True)     
    
    

    hdulist1.close()    
    hdulist.close()        
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky background estimation")
    parser.add_argument("inputImage", help="Input fits image with the object")
    parser.add_argument("segmImage", help="Input segmenation image")
    parser.add_argument("region", help="Input region file with a new mask (x mark those regions which should be removed!")
    parser.add_argument("--new_segm", nargs='?', const=1, help="Optional: New segmentation image", type=str, default='new_segm.fits') 

    args = parser.parse_args()

    input_image = args.inputImage
    segm_image = args.segmImage
    fix_region_file = args.region 
    new_segm_image = args.new_segm

    main(input_image, segm_image, fix_region_file, new_segm_image=new_segm_image)
        
        
        
        
