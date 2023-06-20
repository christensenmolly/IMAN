import sys
from astropy.io import fits as pyfits
import numpy as np
import add_keyw_to_header

def write_keyw(header, keyword, value):
    value = float(value)
    if np.isnan(value)==False:
        if keyword in header.keys():
          header[keyword] = value
        else:
          header.append((keyword,value),end=True)
    return header

def main(input_image, EXPTIME, GAIN, NCOMBINE, m0=None):
    hdulist = pyfits.open(input_image, do_not_scale_image_data=True, mode='update')
    prihdr = hdulist[0].header
    
    add_keyw_to_header.write_keyw(prihdr, 'EXPTIME', EXPTIME)
    add_keyw_to_header.write_keyw(prihdr, 'GAIN', GAIN)
    add_keyw_to_header.write_keyw(prihdr, 'NCOMBINE', NCOMBINE)
    if m0 is not None:
        write_keyw(prihdr, 'M0', m0)

    hdulist.flush()



if __name__ == '__main__':
    input_image = str(sys.argv[1])
    EXPTIME = float(sys.argv[2])
    GAIN = float(sys.argv[3])
    NCOMBINE = float(sys.argv[4])
    main(input_image, EXPTIME, GAIN, NCOMBINE)