import os
import shutil
import subprocess
import numpy as np
import argparse
import urllib.request

def download(scanline, number, band, log, file):
    try:
        path = 'http://das.sdss.org/ge/sample/stripe82/sdss/col%s/%s/' % (scanline, band)

        url = path + file

        print('Downloading %s...' % (url))
        urllib.request.urlretrieve(url, file)


        log.write('%s DONE\n' % (file))

    except:
        log.write('%s FAILED\n' % (file))


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
        
        

def main(ra=[None], dec=[None], name=[None],  bands=['u','g','r','i','z'], frames=['rec','weight','psf']): # type_of_file=['gal.cat']):
  log = open('log_cats.txt', 'w')

    
  for l in range(len(ra)):
    RA = ra[l]
    DEC = dec[l]
    NAME = name[l]
    log.write('%s\n' % (NAME))

    if ((RA>=0. and RA<60.) or (RA>300 and RA<=360.)) and DEC>-1.26 and DEC<1.26:
        z=1
    else:
        print('The coordinates %f,%f are beyond the limits of the Stripe82' % (RA,DEC))
        continue
    
    if NAME is None:
        NAME = '%.5f_%.5f' % (RA,DEC)


    if not os.path.exists(NAME):
        os.makedirs(NAME)
    os.chdir(NAME)
    
    if True:
                scanline,number = find_scanline_number(RA, DEC)
                print("Patch scanline and number for %f,%f are %s and %s" % (RA,DEC,scanline,str(number)[1:-1]))   
                for band in bands:
                    if not os.path.exists(band):
                        os.makedirs(band)
                    os.chdir(band)
                
                    for frame in frames:
                        for num in number:
                            if frame=='ima':
                                file = 'S82_%s%s_%s.fits' % (scanline,band,num)
                            elif frame=='weight':
                                file = 'S82_%s%s_%s.wht.fits' % (scanline,band,num)
                            else:
                                file = 'S82_%s%s_%s.cat.gz' % (scanline,band,num)
                            download(scanline, num, band, log, file)
                    os.chdir('..')
                break; break
    else:
        z=1
        
  os.chdir('..')
    
  log.close()
#c_EON_344.174_0.162_r.fits	344.17380	0.16169	47.62	13.60	9.2
#main(ra=[344.17380], dec=[0.16169], name=['EON_344.174_0.162'],  bands=['r'], frames=['ima'])
main(ra=[320.827], dec=[1.25502], name=['sdss_13078'],  bands=['g','r','i'], frames=['ima'])

'''
log = open('log_cats.txt', 'w')
scanline = '01'
num = '001'
band = 'r'
file = 'S82_%s%s_%s.fits' % (scanline,band,num)
download(scanline, num, band, log, file)

exit()
'''

#main(ra=[2.29825], dec=[-0.6151944], name=['SPRC1'],  bands=['u','g','r','i','z'], frames=['ima','weight']) # Write below the same lines for the test galaxies
#main(ra=[40.7434583], dec=[-0.9525556], name=['SPRC4'],  bands=['u','g','r','i','z'], frames=['ima','weight']) # Write below the same lines for the test galaxies
#main(ra=[312.0236], dec=[0.0688], name=['SPRC69'],  bands=['u','g','r','i','z'], frames=['ima','weight']) # Write below the same lines for the test galaxies
#main(ra=[8.0408333], dec=[1.1435], name=['SPRC73'],  bands=['u','g','r','i','z'], frames=['ima','weight']) # Write below the same lines for the test galaxies
#main(ra=[12.05075], dec=[-0.2154444], name=['SPRC74'],  bands=['u','g','r','i','z'], frames=['ima','weight']) # Write below the same lines for the test galaxies
'''
main(ra=[20.37216], dec=[0.6246389], name=['SPRC76'],  bands=['u','g','r','i','z'], frames=['ima','weight']) # Write below the same lines for the test galaxies
exit()
main(ra=[29.7432917], dec=[-0.4898055], name=['SPRC77'],  bands=['u','g','r','i','z'], frames=['ima','weight']) # Write below the same lines for the test galaxies
main(ra=[348.13375], dec=[-0.1102778], name=['SPRC185'],  bands=['u','g','r','i','z'], frames=['ima','weight']) # Write below the same lines for the test galaxies
main(ra=[12.04], dec=[-0.9122778], name=['SPRC186'],  bands=['u','g','r','i','z'], frames=['ima','weight']) # Write below the same lines for the test galaxies
main(ra=[53.5253333], dec=[1.0944444], name=['SPRC188'],  bands=['u','g','r','i','z'], frames=['ima','weight']) # Write below the same lines for the test galaxies
main(ra=[320.913125], dec=[-0.3764167], name=['SPRC234'],  bands=['u','g','r','i','z'], frames=['ima','weight']) # Write below the same lines for the test galaxies
main(ra=[354.9900417], dec=[0.1351944], name=['SPRC238'],  bands=['u','g','r','i','z'], frames=['ima','weight']) # Write below the same lines for the test galaxies
main(ra=[313.47275], dec=[-0.96817], name=['SPRC275'],  bands=['u','g','r','i','z'], frames=['ima','weight']) # Write below the same lines for the test galaxies
'''


''' 
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Download Stripe82") 
    parser.add_argument("ra", help="Right ascension [deg]", type=float)
    parser.add_argument("dec", help="Declination [deg]", type=float)
    parser.add_argument("--name", help="Optional: Name of the object", type=str, default=None)
    parser.add_argument("--bands", help="Optional: Bands to be downloaded, e.g. gri", type=str, default='ugriz') 
    parser.add_argument("--frames", help="Optional: Types of frames to be downloaded separated by comma [ima,weight,cat], e.g. rec,weight", type=str, default='rec')     

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


