#!/usr/bin/python
#!/usr/bin/python
# DESCRIPTION:
# Script to merge several masks into one.
# Both region and fits formats are supported. Input files should be separated by comma.
# MINIMAL USAGE: python merge_masks.py [input_masks_separated_by_comma] [output_mask]


# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import math
import argparse



def main(input_masks, output_mask, verbosity=True):
    # Find out about type:
    if '.fits' in input_masks[0]:
        mask_type = 'fits'
    elif '.reg' in input_masks[0]:
        mask_type = 'reg'
    else:
        if verbosity: print('ERROR. This type of mask file is not supported! Exiting.')

   
    if mask_type=='reg':
        outfile = open(output_mask, 'w')
        outfile.write('image\n')
        for input_mask in input_masks:
                input_mask = input_mask.strip()
                with open(input_mask) as infile:
                    Lines = infile.readlines()
                    for Line in Lines:
                      if "global" not in Line and 'image' not in Line and '#' not in Line[0]:
                        outfile.write(Line)   
        outfile.close()                        
    elif mask_type=='fits':
        for k in range(len(input_masks)):
            file = input_masks[k].strip()
            HDU = pyfits.open(file)
            Data = HDU[0].data
            if k==0:
                FinalData = Data
                header_data = HDU[0].header
            else:
                FinalData = FinalData + Data
        hdu = pyfits.PrimaryHDU(FinalData, header_data)
        hdu.writeto(output_mask, clobber=True)            
    if verbosity: print('Done!')
    
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Merge masks")
    parser.add_argument("i", help="Input masks (region or fits files). Should be separated by comma.")
    parser.add_argument("o", nargs='?', const=1, help="Optional: Output mask (fits or reg)",type=str, default=None) 
    args = parser.parse_args()

    input_masks = args.i
    output_mask = args.o
    
    input_masks = input_masks.split(',')
    
    
    main(input_masks, output_mask)

    #mask(input_masks, output_mask)
