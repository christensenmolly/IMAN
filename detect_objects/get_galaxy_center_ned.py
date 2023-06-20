import argparse
from astropy.io import fits as pyfits
from astropy import wcs
from astroquery import ned
import numpy as np


def get_coords_by_name(name):
    # Function to retrieve equatorial coordinates of an object by its NED name
    try:
        result_table = ned.Ned.query_object(name)
        if len(result_table)==1:
            try:
                RA,DEC = float(result_table['RA(deg)']),float(result_table['DEC(deg)']) # an astropy.table.Table
            except:
                RA,DEC = float(result_table['RA']),float(result_table['DEC']) # an astropy.table.Table
            return RA,DEC
        else:
            return float('nan'),float('nan')
    except:
         return float('nan'),float('nan')
     
def convert_radec_to_image(input_image, RA, DEC):
    # Function to convert wcs coordinates to image
    HDU = pyfits.open(input_image)
    header = HDU[0].header

    w = wcs.WCS(header)
    world = np.array([[RA,DEC]], np.float_)
    pixcrd = w.wcs_world2pix(world, 1)
    xc,yc = pixcrd[0][0],pixcrd[0][1]

    HDU.close()
    return xc,yc


def convert_image_to_radec(input_image, xc, yc):
    # Function to convert image to wcs coordinates
    HDU = pyfits.open(input_image)
    header = HDU[0].header

    w = wcs.WCS(header)
    pixcrd = np.array([[xc,yc]], np.float_)
    world = w.wcs_pix2world(pixcrd, 1)      
    RA,DEC = world[0][0],world[0][1]

    HDU.close()
    return RA,DEC    

def main(input_image, name=None, RA=None, DEC=None, verbosity=True):
    # Main function to retrieve image coordinates for the target object
    if (name is None) and (RA is None and DEC is None):
        if verbosity: print('\033[91m ERROR! No object is specified! Exiting... \033[0m')
        exit()
    
    if name is not None:
        RA_found,DEC_found = get_coords_by_name(name)
        
        if not np.isnan(RA_found):                
            if verbosity: print('\033[92m Working with Ra=%f and Dec=%f found in NED for the object %s \033[0m' % (float(RA_found), float(DEC_found), name))
            RA = RA_found 
            DEC = DEC_found
        else:
            if verbosity: print('\033[93m WARNING: The object name %s that you have submitted is not currently recognized by the NED name interpreter! \033[0m' % (name))
            if RA is None and DEC is None:
                if verbosity: print('\033[91m ERROR! Ra and Dec are None. Exiting... \033[0m')
                exit()
            else:
                if verbosity: print('\033[92m Working with Ra=%f and Dec=%f you specified as the galaxy center. \033[0m' % (float(RA), float(DEC)))
    
    xc,yc = convert_radec_to_image(input_image, RA, DEC)
    if verbosity: print(' The image coordinates of the target: %.2f, %.2f' % (float(xc), float(yc)))

    return xc,yc


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get center of galaxy target (in image coordinates) by name or ra,dec")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("--name", help="Galaxy name (NED)",type=str,default=None) 
    parser.add_argument("--ra", nargs='?', const=1, help="RA (deg) if no name is specified.",type=float, default=None) 
    parser.add_argument("--dec", nargs='?', const=1, help="DEC (deg) if no name is specified.",type=float, default=None) 
    
    args = parser.parse_args()

    input_image = args.inputImage
    name = args.name
    RA = args.ra
    DEC = args.dec

    main(input_image, name=name, RA=RA, DEC=DEC)
