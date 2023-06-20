
from astroscrappy import detect_cosmics
from astropy.io import fits
import argparse

# For ARCTIC, see https://www.apo.nmsu.edu/arc35m/Instruments/ARCTIC/#3p4


def main(input_image, output_image, ron, gain, satlevel, fwhm, psfsize):
    # See https://github.com/astropy/astroscrappy/blob/main/astroscrappy/astroscrappy.pyx
    data = fits.getdata(input_image, ext=0)
    head = fits.getheader(input_image, ext=0)
    
    mask,_clean = detect_cosmics(data, inmask=None, sigclip=4.5,
        sigfrac=0.3, objlim=5.0, gain=gain,
        readnoise=ron, satlevel=satlevel,
        niter=4, sepmed=True,
        cleantype='meanmask', fsmode='median',
        psfmodel='gauss', psffwhm=fwhm ,
        psfsize =psfsize,
        psfk=None, psfbeta=4.765, verbose=False)
    
    fits.writeto(output_image, _clean, head, overwrite=True)
    fits.writeto('cosmic_rays_mask.fits', mask * 1, head, overwrite=True)
    print('Done!')
    return output_image



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Remove cosmic rays")
    parser.add_argument("input_image", help="Input fits image")
    #parser.add_argument("model_file", help="Input model file")
    parser.add_argument("output_image", help="Output image")
    parser.add_argument("--ron", nargs='?', const=1, help="Optional: Read-out noise (e-)",type=float,default=6.5) 
    parser.add_argument("--gain", nargs='?', const=1, help="Optional: Gain (e-/ADU)",type=float,default=1.0)
    parser.add_argument("--satlevel", nargs='?', const=1, help="Optional: Saturation level (ADU)",type=float,default=65536.0)
    parser.add_argument("--fwhm", nargs='?', const=1, help="Optional: FWHM (pix)",type=float,default=2.5)
    

    
    
    args = parser.parse_args()

    input_image = args.input_image
    output_image = args.output_image
    ron = args.ron
    gain = args.gain
    satlevel = args.satlevel
    fwhm = args.fwhm
    
    psfsize = int(3.* fwhm)

    if (psfsize % 2) == 0:
        psfsize +=1

    
    main(input_image, output_image, ron, gain, satlevel, fwhm, psfsize)
    
    

