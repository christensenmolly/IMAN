import os
import shutil
import subprocess
import numpy as np
import argparse
import urllib.request
import sys
import glob 
import math 

LOCAL_DIR = "/misc_funcs"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]


sys.path.append(os.path.join(IMAN_DIR, 'detect_objects'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))

import crop_image
import get_galaxy_center_ned


#USAGE: python3 ~/MyGit/IMAN/misc_funcs/download_Stripe82.py 20.807111 -0.385269 --name UGC929 --bands g,r --frames rec,psf


if True:
    # Check what name has SWarp package on this system
    rCode = subprocess.call("which swarp >/dev/null", shell=True)
    if rCode == 0:
        swarpName = "swarp"
    else:
        rCode = subprocess.call("which SWarp >/dev/null", shell=True)
        if rCode == 0:
            swarpName = "SWarp"
        else:
            print('Error: SWarp was not found on your system.')
            print('The command has to be either swarp or SWarp')
            print('Install SWarp package or try to run this script without -s option')
            exit(1)





def find_scanline_number(RA, DEC):
    scanline = None
    number = None
    
    if DEC>-1.26 and DEC<-1.05:
        scanline = '01'
    if DEC>=-1.05 and DEC<-0.84:
        scanline = '07'

    if DEC>=-0.84 and DEC<-0.63:
        scanline = '02'
    if DEC>=-0.63 and DEC<-0.42:
        scanline = '08'        

    if DEC>=-0.42 and DEC<-0.21:
        scanline = '03'
    if DEC>=-0.21 and DEC<0.0:
        scanline = '09' 
        
    if DEC>=0.0 and DEC<0.21:
        scanline = '04'
    if DEC>=0.21 and DEC<0.42:
        scanline = '10' 

    if DEC>=0.42 and DEC<0.63:
        scanline = '05'
    if DEC>=0.63 and DEC<0.84:
        scanline = '11' 

    if DEC>=0.84 and DEC<1.05:
        scanline = '06'
    if DEC>=1.05 and DEC<1.26:
        scanline = '12'
    

    #stripe82_ra_range = np.linspace(299.992, 420.074, 401)
    
    stripe82_ra_range = np.arange(300, 420., 0.29942)


    fields = []
    for k in range(len(stripe82_ra_range)-1):
        ra_min = stripe82_ra_range[k] - 65.*0.396/3600.
        ra_max = stripe82_ra_range[k+1] + 63.*0.396/3600.
        fields.append([ra_min,ra_max]) 


    if RA>=0 and RA<=60.:
        RA = RA + 360.

    number = []
    for k in range(len(fields)):
        if RA>=fields[k][0] and RA<=fields[k][1]:
            number.append(str(k+1).zfill(3))
    
    return scanline,number
        
def download(number, type_of_file, name=None):
        if type_of_file=='ps.cat' or type_of_file=='gal.cat':
            path = 'ftp://stripero:s0l0le0@ftp.iac.es/catalogs/'
        else:
            path = 'ftp://stripero:s0l0le0@ftp.iac.es/coadds/'
        print(type_of_file)

        file = 'f%s_' % (number) + type_of_file 
        print('%s%s' % (path, file))
        urllib.request.urlretrieve('%s%s' % (path, file), file)
        
def download_jiang(scanline, number, band, log, file):
    try:
        path = 'http://das.sdss.org/ge/sample/stripe82/sdss/col%s/%s/' % (scanline, band)

        url = path + file

        print('Downloading %s...' % (url))
        urllib.request.urlretrieve(url, file)

        log.write('%s Jiang DONE\n' % (file))

    except:
        log.write('%s Jiang FAILED\n' % (file))



def swarping(galName, band, files):
    print("Running SWarp for %s" % (galName))
    callSt = "%s -verbose_type quiet -BACK_TYPE MANUAL " % (swarpName)
    callSt += " ".join(["%s" % (s) for s in files])
    subprocess.call(callSt, shell="True")
    shutil.move("./coadd.fits", "./%s_%s.fits" % (galName, band))
    os.remove("coadd.weight.fits")
    os.remove("swarp.xml")    


def main(ra=[None], dec=[None], name=[None], radius=[None], bands=['u','g','r','i','z'], frames=['rec','weight','psf'], pix2sec=0.396): # type_of_file=['gal.cat']):

  
  DECS = np.arange(-1.0,1.5,0.5)
  RAS = np.arange(310.25,420.25,0.5)
  
    
  for l in range(len(ra)):
    log = open('log_%s.txt' % (name[l]), 'w')  
    RA = ra[l]
    DEC = dec[l]
    NAME = name[l]
    RADIUS = radius[l] # In arcmin

    if NAME is None:
        NAME = '%.5f_%.5f' % (RA,DEC)    
        

    if ((RA>=0. and RA<60.) or (RA>309 and RA<=360.)) and DEC>-1.26 and DEC<1.26:
        pass
    else:
        print('The coordinates %f,%f are beyond the limits of the Stripe82' % (RA,DEC))
        log.write('%s NO DATA FOUND\n' % (NAME))
        continue
    

    if not os.path.exists(NAME):
        os.makedirs(NAME)
    
    os.chdir(NAME)

    ra_min = RA - RADIUS/60.
    dec_min = DEC - RADIUS/60.
    ra_max = RA + RADIUS/60.
    dec_max = DEC + RADIUS/60.
    coords_possible = [[ra_min,dec_min],[ra_min,dec_max],[ra_max,dec_max],[ra_max,dec_min]]
      
    # Determine the numbers of patches we need taking into account the radius of the galaxy:
    NUMBERS = []
    for k in range(len(RAS)):
      if RAS[k]>360.0:
        RAS[k] = RAS[k] - 360.
      
      xxx = k + 1
      for i in range(len(DECS)):
        y = i+1
        number = "%03d%i" % (xxx,y)
        
        for coords in coords_possible:
            if coords[0]>=RAS[k]-0.25 and coords[0]<=RAS[k]+0.25 and coords[1]>=DECS[i]-0.25 and coords[1]<=DECS[i]+0.25:
                if number not in NUMBERS:
                    NUMBERS.append(number)

    

    #'''
    # Download specific patches
    stop_loop = False
    for number in NUMBERS:
                print("Patch number for %f,%f is %s" % (RA,DEC,number))   
                for band in bands:
                    if not os.path.exists(band):
                        os.makedirs(band)
                    os.chdir(band)
                
                    for frame in frames:
                        if frame!='ps' and frame!='gal':
                            type_of_file = '%s.%s.fits' % (band,frame)
                        else:
                            type_of_file = '%s.cat' % (frame)
                        
                        try:
                            download(number, type_of_file, name=NAME)
                        except:
                            log.write('%s %s %s %s FAILED\n' % (NAME, number, band, frame))
                            stop_loop = True
                            break
                    os.chdir('..')
                    if stop_loop:
                        break
    
    
    
    if stop_loop:
        scanline,number = find_scanline_number(RA, DEC)  # Jiang
        for band in bands:
                    if not os.path.exists(band):
                        os.makedirs(band)
                    os.chdir(band)
                
                    for frame in frames:
                        for num in number:
                            if frame=='rec':
                                file = 'S82_%s%s_%s.fits' % (scanline,band,num)
                            elif frame=='weight':
                                file = 'S82_%s%s_%s.wht.fits' % (scanline,band,num)
                            #else:
                            #    file = 'S82_%s%s_%s.cat.gz' % (scanline,band,num)
                            download_jiang(scanline, num, band, log, file)
                    os.chdir('..')
      
    

    
    
    #'''
    #if len(NUMBERS)>1:
    #    return 0
        
    if len(NUMBERS)>1:
        # SWarping:
        for band in bands:
            os.chdir(band)
            for frame in frames:
                if frame=='rec':
                    files = glob.glob('*_%s.%s.fits' % (band,frame))
                    try:
                        swarping(NAME, band, files)
                    except:
                        log.write('%s %s %s FAILED\n' % (NAME, band, frame))
                        continue
            os.chdir('..')
    else:
        for band in bands:
            os.chdir(band)
            for frame in frames:
                if frame=='rec':
                    type_of_file = '%s.%s.fits' % (band,frame)
                    try:
                        shutil.move('./f%s_' % (NUMBERS[0]) + type_of_file, "./%s_%s.fits" % (NAME, band))
                    except:
                        log.write('%s %s %s FAILED\n' % (NAME, band, frame))
                        continue
            os.chdir('..')

    # Cropping:
    if RADIUS is not None:
        for band in bands:
            os.chdir(band)
            for frame in frames:
                if frame=='rec':
                    if os.path.exists("./%s_%s.fits" % (NAME, band)):
                        xc,yc = get_galaxy_center_ned.convert_radec_to_image("./%s_%s.fits" % (NAME, band), RA, DEC)
                        x_max = xc + int(math.ceil(RADIUS*60./pix2sec))
                        x_min = xc - int(math.ceil(RADIUS*60./pix2sec))
                        y_max = yc + int(math.ceil(RADIUS*60./pix2sec))
                        y_min = yc - int(math.ceil(RADIUS*60./pix2sec))
                        crop_image.main("./%s_%s.fits" % (NAME, band), x_min, y_min, x_max, y_max, output_image="./%s_%s_cropped.fits" % (NAME, band), hdu=0)
                        os.remove("./%s_%s.fits" % (NAME, band))
                    else:
                        z=1
            os.chdir('..')
        
    
    os.chdir('..')
    
  log.close()



'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Download from Stripe82") 
    parser.add_argument("ra", help="Right ascension [deg]", type=float)
    parser.add_argument("dec", help="Declination [deg]", type=float)
    parser.add_argument("--name", help="Optional: Name of the object", type=str, default=None)
    parser.add_argument("--bands", help="Optional: Bands to be downloaded, e.g. g,r,i", type=str, default='u,g,r,i,z') 
    parser.add_argument("--frames", help="Optional: Types of frames to be downloaded separated by comma [rec,weight,ps,gal,psf], e.g. rec,weight", type=str, default='rec')     

    args = parser.parse_args()

    ra = args.ra
    dec = args.dec
    
    name = args.name
    bands = args.bands
    frames = args.frames
    
    bands = bands.split(',')
    frames = frames.split(',')
    
    main(ra=[ra], dec=[dec], name=[name], bands=bands, frames=frames)
'''

#main(ra=[318.718], dec=[-0.97777], name=['1'], radius=[223.*0.396/60.], bands=['g','r','i'], frames=['rec'], pix2sec=0.396)
