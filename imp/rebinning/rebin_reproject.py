from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from reproject import reproject_interp
import argparse

import rebin_image


def main(reference_image, input_image, output_image):

    pix2sec_ref,notes = rebin_image.resolution(reference_image)
    pix2sec_inp,notes = rebin_image.resolution(input_image)
    factor = float(pix2sec_inp) / float(pix2sec_ref)
    '''
    rebin_image.downsample(legacy_image, factor, output_image='legacy_rebin.fits', set_wcs=True, print_mes=True, norm=False)
    
    rebin_image.rebin('legacy_rebin.fits', sdss_image, output_image='sdss_rebin.fits', hdu_ref=0, hdu_inp=0, preserve_bad_pixels=True, no_interp=False)
    '''

    hdu1 = fits.open(reference_image)[0]
    hdu2 = fits.open(input_image)[0]
    array, footprint = reproject_interp(hdu2, hdu1.header)
    
    hdu_out = fits.PrimaryHDU(data=array/(factor**2), header=hdu1.header)
    hdu_out.writeto(output_image, overwrite=True)
    #'''

#main('PGC10108_sky_subtr_g.fits', 'PGC10108_g_trim.fits', 'sdss_rebin_new.fits')    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Rebinning (New version based on reprojected)")
    parser.add_argument("reference_image", help="Reference image")
    parser.add_argument("input_image", help="Input image")    
    parser.add_argument("output_image", help="Output image")
 
    args = parser.parse_args()

    reference_image = args.reference_image
    input_image = args.input_image
    output_image = args.output_image

    main(reference_image, input_image, output_image)
    
