
from astropy.io import fits
import argparse
from astropy.stats import sigma_clipped_stats, SigmaClip
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint
import numpy as np
import photutils

def write_keyw(header, keyword, value):
    value = float(value)
    if np.isnan(value)==False:
        if keyword in header.keys():
          header[keyword] = value
        else:
          header.append((keyword,value),end=True)
    return header

def main(input_image, output_image, sigma_threshold=2.0, npixels=5, footprint_radius=1):
    hdulist = fits.open(input_image, do_not_scale_image_data=True, mode='update')
    head = hdulist[0].header
    data = hdulist[0].data

    sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
    threshold = detect_threshold(data, nsigma=sigma_threshold, sigma_clip=sigma_clip)
    segment_img = detect_sources(data, threshold, npixels=npixels)
    footprint = circular_footprint(radius=footprint_radius)
    mask = segment_img.make_source_mask(footprint=None)
    
    outHDU = fits.PrimaryHDU(1.*mask, header=head)
    outHDU.writeto('mask.fits', overwrite=True)
    
    
    mean,median,std = sigma_clipped_stats(data, sigma=3.0, mask=mask)
    backsub_data = data - median
    print('Subtracted median sky is: %.8f' % (median))
    

    write_keyw(head, 'MED_SKY', median)

    hdulist.flush()
    
    fits.writeto(output_image, backsub_data, head, overwrite=True)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Estimate the flat background (median) and subtract it")
    parser.add_argument("input_image", help="Input fits image")
    parser.add_argument("output_image", help="Output image")
    parser.add_argument("--sigma", nargs='?', const=1, help="Optional: The number of standard deviations per pixel",type=float,default=2.0) 
    parser.add_argument("--npixels", nargs='?', const=1, help="Optional: The minimum number of connected pixels to be masked", type=float,default=3) 
    parser.add_argument("--fp", nargs='?', const=1, help="Optional: Footprint radius. The local footprint used for the source dilation. ", type=int, default=1) 

    

    
    
    args = parser.parse_args()

    input_image = args.input_image
    output_image = args.output_image
    sigma_threshold = args.sigma
    npixels = args.npixels
    footprint_radius = args.fp

    
    main(input_image, output_image, sigma_threshold=sigma_threshold, npixels=npixels, footprint_radius=footprint_radius)
