

import sys
import astropy.io.fits as fits
import os

LOCAL_DIR = "/imred"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]

sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/sky_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/masking'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))

import photometric_calibration
import remove_cosmic_rays
import rebin_image
import determine_sky
import auto_masking
import convert_reg_to_mask
import backg_subtr
import crop_image
import photutil_bkg
import glob
import numpy as np
import argparse

### CCD and IMAGE PARAMETERS:
fwhm = 1.0 # in arcsec
polynomial_degree = 1
x_l=3
y_l=3
x_r=2048
y_r=2048
satlevel = 65536.




def main(steps, cals_path=None, science_path=None, science_prefix=None, bias_prefix=None, dark_prefix=None, flat_prefix=None, input_images=None):
    if 1 in steps:
        if cals_path is None:
            cals_path = str(input('\n Provide the path to the directory with calibrations:'))

        if science_path is None:
            science_path = str(input('\n Provide the path to the directory with science images:'))
            
        if science_prefix is None:
            science_prefix = str(input('\n Provide the prefix for science images:'))

        if flat_prefix is None:
            flat_prefix = str(input('\n Provide the prefix for flat fields:'))            
            
       
        # Cals
        output_images = photometric_calibration.main(cals_path, science_path, science_prefix, bias_prefix=bias_prefix, dark_prefix=dark_prefix, flat_prefix=flat_prefix)
        
        return output_images

    if 2 in steps:
        input_images = glob.glob('*_reduced.fits')
        
        # Cropping
        crop_output_images = []
        for input_image in input_images:
            crop_image.main(input_image, x_l, y_l, x_r, y_r, output_image=None, hdu=0)
            crop_output_images.append(input_image.split('.fits')[0] + '_crop.fits')
        output_images = crop_output_images
        return output_images

    if 3 in steps:
        input_images = glob.glob('*_reduced_crop.fits')
        
        # Rays
        ray_output_images = []
        for input_image in input_images:
            hdulist = fits.open(input_image)
            img = hdulist[0].data
            header = hdulist[0].header
            gain = float(header['GTGAIN11'])
            ron = float(header['GTRON11'])
            
            
            pixelscale,note = rebin_image.resolution(input_image)
            
            output = remove_cosmic_rays.main(input_image, input_image.split('/')[-1].split('.fits')[0]+'_rayremoved.fits', ron, gain, satlevel, fwhm/float(pixelscale), 3.*2.*fwhm/float(pixelscale))
            
            ray_output_images.append(output)
        output_images = ray_output_images
        return output_images


        
    if 4 in steps:
        input_images = glob.glob('*_reduced_crop_rayremoved.fits')
        # Sky subtraction

        for input_image in input_images:
          do_image = str(input('\n Do you want to process %s? (yes):' % (input_image)) or 'yes')
          if do_image=='yes':
            #auto_masking.main(input_image, output_region_file='general_mask.reg', snr=2., min_pix=5, region_type='polygon', sextr_setup=None, galaxy_ellipse=None, offset_size=1.3, offset_pix=0., verbosity=True)
            #exit()
            
            ##backg_subtr.main(input_image, input_image.split('/')[-1].split('.fits')[0]+'_skysub.fits', sigma_threshold=2.0, npixels=3, footprint_radius=1)
            
            #convert_reg_to_mask.mask(input_image, 'general_mask.reg', output_image='mask.fits', output_mask=None, mask_value=1, mask_DN=None, verbosity=True)
            
            #determine_sky.sky_subtraction(input_image, 'segm.fits', polynomial_degree=1, output_image=input_image.split('/')[-1].split('.fits')[0]+'_skysub.fits', output_sky=None, hdu_inp=0, sampling_factor=1., sigma_smooth=None, verbosity=True, sky_value='mode')
            #exit()
            
            happy = 'no'
            while happy!='yes':
                polynomial_degree = int(input('\n Enter the polynomial degree for initial sky subtraction (1):') or 1)
                sigma = float(input('\n Enter the sigma threshold for sky subtraction (3.0):') or 3.0)
                box_size = int(input('\n Enter the box size for sky subtraction (50):') or 50)
                filter_size = int(input('\n Enter the filter size for sky subtraction (3):') or 3)
                
                photutil_bkg.main(input_image, output_image=input_image.split('/')[-1].split('.fits')[0]+'_skysub.fits', sigma=sigma, box_size=box_size, filter_size=filter_size, polynomial_degree=polynomial_degree)
                
                happy = str(input('\n Are you happy with this? (yes):') or 'yes')

            
    if 5 in steps:
        input_images = glob.glob('*_astro.fits')
        do_image = str(input('\n Do you want to process %s? (yes):' % (",".join(str(x) for x in input_images))) or 'yes')
        
        if do_image!='yes':
            input_images = str(input('\n Specify images to stack using , '))
            input_images = input_images.split(',')




        if True:
            for k in range(len(input_images)):
                print(k, input_images[k])
            
            reference_index = int(input('\nChoose the index for the reference image to which all other images will be alligned to. Default 0.: ') or 0)
            
            print('The reference image is %s' % (input_images[reference_index]))
            
            imgs = []
            for k in range(len(input_images)):
                rebin_image.rebin(input_images[reference_index], input_images[k], output_image=input_images[k].split('/')[-1].split('.fits')[0]+'_rebin.fits', hdu_ref=0, hdu_inp=0, preserve_bad_pixels=True, no_interp=False)

                hdulist = fits.open(input_images[k].split('/')[-1].split('.fits')[0]+'_rebin.fits')
                imgs.append(hdulist[0].data)
                if k==0:
                    header = hdulist[0].header
            
            
            
        do_median = str(input('\n Do you want to median these images? (yes):') or 'yes')
        
        if do_median=='yes':
            final_science = np.median(imgs, axis = 0)
        else:
            final_science = np.mean(imgs, axis = 0)
            
        outHDU = fits.PrimaryHDU(final_science, header=header)
            
        outHDU.writeto('final_image.fits', overwrite=True)        





if __name__ == '__main__':
    print("Photometric Data Reduction:")
    print("Step 1: Initial data reduction to correct for calibration frames.")
    print("Step 2: Automated cropping of the galaxy image [x_l, y_l: x_r, y_r] to remove a frame with bad pixels.")
    print("Step 3: Automated removement of cosmic rays.")
    print("Step 4: Sky subtraction - under human supervision.")
    print("Step 5: Median stacking.")
    parser = argparse.ArgumentParser(description="Photometric Data Reduction:")
    parser.add_argument("--steps", nargs='?', const=1, help="Specify steps separated by comma", type=str, default='1,2,3,4,5')
    parser.add_argument("--cals", nargs='?', const=1, help="Provide the path to calibrations", type=str, default=None)    
    parser.add_argument("--sci", nargs='?', const=1, help="Provide the path to science frames", type=str, default=None) 
    parser.add_argument("--sci_prefix", nargs='?', const=1, help="Provide the prefix for science frames", type=str, default=None) 
    parser.add_argument("--flats_prefix", nargs='?', const=1, help="Provide the prefix for flat fields", type=str, default='*flat.g*')
    parser.add_argument("--bias_prefix", nargs='?', const=1, help="Provide the prefix for bias frames", type=str, default='*bias*.fits')
    parser.add_argument("--dark_prefix", nargs='?', const=1, help="Provide the prefix for dark frames", type=str, default='*dark*.fits')
    
    args = parser.parse_args()
    
    
    steps = np.array(args.steps.split(','), int)
    cals_path = args.cals
    science_path = args.sci
    science_prefix = args.sci_prefix
    flat_prefix = args.flats_prefix
    bias_prefix = args.bias_prefix
    dark_prefix = args.dark_prefix
    
    main(steps=steps, cals_path=cals_path, science_path=science_path, science_prefix=science_prefix, bias_prefix=bias_prefix, dark_prefix=dark_prefix, flat_prefix=flat_prefix, input_images=None)





    
#output_images = main(steps=[1])
#output_images = main(steps=[2], input_images=output_images)
#output_images = main(steps=[3], input_images=output_images)

#input_images = glob.glob('*_rayremoved.fits')
#main(steps=[4], input_images=input_images)  #  sigma_clip = SigmaClip(sigma=3.0), Background2D(data, (50, 50), filter_size=(3, 3),



