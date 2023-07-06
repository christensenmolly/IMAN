import sys
from astropy.io import fits as pyfits
import numpy as np




def main(input_image, header_keyword, keyword_value):
    hdulist = pyfits.open(input_image, do_not_scale_image_data=True, mode='update')
    prihdr = hdulist[0].header

    prihdr[header_keyword] = keyword_value

    hdulist.flush()



if __name__ == '__main__':
    input_image = str(sys.argv[1])
    header_keyword = str(sys.argv[2])
    keyword_value = float(sys.argv[3])
    main(input_image, header_keyword, keyword_value)
