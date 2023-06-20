import numpy as np
from os.path import splitext
from astropy.io import fits
from modopt.base.types import check_npndarray
import sys

def write_to_fits(file_name, data):
    """Write FITS file
    This method writes the output image array data to a FITS file.
    Parameters
    ----------
    file_name : str
        Name of file with path
    data : np.ndarray
        Image data array
    """

    fits.PrimaryHDU(data).writeto(file_name)
    
    
if __name__ == "__main__":
    input_image = sys.argv[1]
    output_image = sys.argv[2]

    
    data = np.load(input_image)
    
    write_to_fits(output_image, data)