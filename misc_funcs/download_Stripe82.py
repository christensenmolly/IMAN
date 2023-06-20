import os
import shutil
import subprocess
import numpy as np
import argparse
import urllib.request
#USAGE: python3 ~/MyGit/IMAN/misc_funcs/download_Stripe82.py 20.807111 -0.385269 --name UGC929 --bands g,r --frames rec,psf

def download(number,log, type_of_file, name=None):
    try:
        if type_of_file=='ps.cat' or type_of_file=='gal.cat':
            path = 'ftp://stripero:s0l0le0@ftp.iac.es/catalogs/'
        else:
            path = 'ftp://stripero:s0l0le0@ftp.iac.es/coadds/'
        print(type_of_file)

        file = 'f%s_' % (number) + type_of_file 
        print('%s%s' % (path, file))
        urllib.request.urlretrieve('%s%s' % (path, file), file)

        log.write('%s DONE\n' % (type_of_file))

    except:
        log.write('%s FAILED\n' % (type_of_file))


def main(ra=[None], dec=[None], name=[None],  bands=['u','g','r','i','z'], frames=['rec','weight','psf']): # type_of_file=['gal.cat']):
  log = open('log_cats.txt', 'w')
  DECS = np.arange(-1.0,1.5,0.5)
  RAS = np.arange(310.25,420.25,0.5)
  
    
  for l in range(len(ra)):
    RA = ra[l]
    DEC = dec[l]
    NAME = name[l]
    log.write('%s\n' % (NAME))

    if ((RA>=0. and RA<60.) or (RA>310 and RA<=360.)) and DEC>-1.25 and DEC<1.25:
        z=1
    else:
        print('The coordinates %f,%f are beyond the limits of the Stripe82' % (RA,DEC))
        continue
    
    if NAME is None:
        NAME = '%.5f_%.5f' % (RA,DEC)


    if not os.path.exists(NAME):
        os.makedirs(NAME)
    os.chdir(NAME)
    
    
    
    for k in range(len(RAS)):
      if RAS[k]>360.0:
        RAS[k] = RAS[k] - 360.
      xxx = k + 1
      for i in range(len(DECS)):
        y = i+1
        number = "%03d%i" % (xxx,y)
        if RA is not None and DEC is not None:
            if RA>=RAS[k]-0.25 and RA<=RAS[k]+0.25 and DEC>=DECS[i]-0.25 and DEC<=DECS[i]+0.25:
                # Download specific patches
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
                        download(number, log, type_of_file, name=NAME)
                    os.chdir('..')
                break; break
        else:
                # Download all patches
                print("Patch number for %f,%f is %s" % (RAS[k],DECS[i],number))
                download(number,log, type_of_file, name=NAME)
    os.chdir('..')
    
  log.close()





#main(ra=[358.546], dec=[0.38], name=['HCG98'],  bands=['u','g'], frames=['psf','ps'])

#'''
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
#'''

#main(ra=[358.546], dec=[0.38], name=['HCG98'],  bands=['rdeep'], frames=['rec','psf']) # 'weight' does not work! Why?
