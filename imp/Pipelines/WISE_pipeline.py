import sys
import astropy.io.fits as fits
import os
import subprocess
import astroquery

LOCAL_DIR = "/imp/Pipelines"
IMAN_DIR = os.path.dirname(__file__).rpartition(LOCAL_DIR)[0]


sys.path.append(os.path.join(IMAN_DIR, 'imp/rebinning'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/sky_fitting'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/masking'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/cropping'))
sys.path.append(os.path.join(IMAN_DIR, 'imp'))
sys.path.append(os.path.join(IMAN_DIR, 'misc_funcs'))
sys.path.append(os.path.join(IMAN_DIR, 'detect_objects'))
sys.path.append(os.path.join(IMAN_DIR, 'imp/psf'))
sys.path.append(os.path.join(IMAN_DIR, 'deca_lite'))

#import photometric_calibration
#import remove_cosmic_rays
#import rebin_image
#import determine_sky
import auto_masking
import convert_reg_to_mask
import determine_sky
import crop_image
#import photutil_bkg
import glob
import numpy as np
import argparse

#from IMP import IMP_main
import IMP
import interpolate_nans
import get_galaxy_center_ned
import fit_PSF_PA
import get_galaxy_ellipse
import rotate_psf
import deca_lite
import crop_galaxy_image

import platform
opsyst = platform.system()


# WISE W1 PSF FWHM:
FWHM = 6./1.375	# in pixels




def get_info_about_galaxy(galaxy_name):
            # Retrieve galaxy's coordinates'
            RA_object,DEC_object = get_galaxy_center_ned.get_coords_by_name(galaxy_name)



            # Get 2MASS galaxy parameters
            coords_cen = coordinates.SkyCoord(RA_object, DEC_object, unit=('deg','deg'), frame='fk5')
            table = astroquery.ipac.irsa.query_region(coords_cen, catalog='fp_xsc', spatial='Cone', radius= 10. * u.arcsec)


            dist = list(table['dist'])
            if dist!=[]:
                min_dist_index = dist.index(min(dist))
                sma_object = 2.5 * float(table['r_k20fe'][min_dist_index])	# in arcsec
                smb_object = 1.5 * float(sma_object * table['k_ba'][min_dist_index])	# in arcsec
                PA_object = float(table['k_phi'][min_dist_index]) + 90.
                q_ratio = smb_object / sma_object
                d25 = sma_object * 2. / 60. 	# in arcmin
            else:
                d25 = float('nan')

            if np.isnan(d25)==True:
                    if np.isnan(float(data['logd25'][k]))==True:
                        sma_object = 60.
                        smb_object = 60.
                        PA_object = 0.
                        q_ratio = 1.
                        d25 = 2.
                    else:
                        sma_object = 10**data['logd25'][k] / 10. * 60.	# in arcsec
                        smb_object = 1.5 * float(sma_object / 10**data['logr25'][k])
                        PA_object = float(data['pa'][k]) + 90.
                        q_ratio = smb_object / sma_object
                        d25 = sma_object * 2. / 60. 	# in arcmin

            return RA_object,DEC_object,sma_object,smb_object,PA_object,q_ratio,d25


def fix_images(galaxy_name, input_int_image, input_unc_image):
                    # New input images in the current directory:
                    input_int = galaxy_name + '_WISE_W1_int.fits'
                    input_unc = galaxy_name + '_WISE_W1_unc.fits'

                    # Fixing nans:
                    hdulist_inp = pyfits.open(input_int_image)
                    header = hdulist_inp[0].header
                    scidata = hdulist_inp[0].data
                    np.putmask(scidata, scidata==float('inf'), float('nan'))
                    new_data = interpolate_nans.main(scidata, 5, 0.5, kernel_size=int(ceil(FWHM)), method='localmean')

                    hdulist_unc = pyfits.open(input_unc_image)
                    header_unc = hdulist_unc[0].header
                    scidata_unc = hdulist_unc[0].data
                    np.putmask(scidata_unc, scidata==float('inf'), float('nan'))
                    new_data_unc = interpolate_nans.main(scidata_unc, 5, 0.5, kernel_size=int(ceil(FWHM)), method='localmean')

                    # Fixing header:
                    if 'COMMENT' in header:
                        del header['COMMENT']
                    if 'HISTORY' in header:
                        del header['HISTORY']
                    if '' in header:
                        del header['']

                    w = wcs.WCS(header)
                    world = np.array([[RA_object,DEC_object]], np.float_)
                    pixcrd = w.wcs_world2pix(world, 1)
                    xc = pixcrd[0][0]
                    yc = pixcrd[0][1]

                    # Save new images:
                    hdu = pyfits.PrimaryHDU(new_data, header)
                    hdu.writeto(input_int, clobber=True)

                    hdu_unc = pyfits.PrimaryHDU(new_data_unc, header_unc)
                    hdu_unc.writeto(input_unc, clobber=True)


def find_PA(file_Image, input_psf_file, RA, DEC, show_running=True):
    xc_psf = 120.
    yc_psf = 120.


    # Detect the brightest peaks:
    if show_running==True:
      print('\t Detecting peaks...')

    hdulist = pyfits.open(file_Image) # open FITS file
    Data = hdulist[0].data
    header = hdulist[0].header

    header = hdulist[0].header
    if 'COMMENT' in header:
		    del header['COMMENT']
    if 'HISTORY' in header:
		    del header['HISTORY']
    if '' in header:
        del header['']

    w = wcs.WCS(header)

    world = np.array([[RA,DEC]], np.float_)
    pixcrd = w.wcs_world2pix(world, 1)
    xc_gal = int(pixcrd[0][0])
    yc_gal = int(pixcrd[0][1])


    mean, median, std = sigma_clipped_stats(Data, sigma=3.0, iters=5)
    threshold = median + 3.*std

    tbl = find_peaks(Data, threshold, box_size=30, npeaks=50)
    tbl.sort('peak_value')
    tbl.reverse()

    if show_running==True:
      print('Done!')
    hdulist_segm = pyfits.open('segm.fits') # open FITS file
    Data_segm = hdulist_segm[0].data
    ny,nx = np.shape(Data_segm)

    sky_level = median
    for kk in range(len(tbl)):
      x = int(tbl['x_peak'][kk])
      y = int(tbl['y_peak'][kk])
      if Data_segm[y,x]!=Data_segm[yc_gal,xc_gal]:
        xc = x
        yc = y
        (yy,xx) = np.where(Data_segm==Data_segm[yc, xc])
        xmin = min(xx)
        ymin = min(yy)
        xmax = max(xx)
        ymax = max(yy)
        if xmax-xmin>300 or ymax-ymin>300:
            continue
        np.putmask(Data_segm, Data_segm==Data_segm[yc, xc], 0.)

        outHDU = pyfits.PrimaryHDU(Data_segm)
        outHDU.writeto('mask_psf.fits',clobber=True)
        #exit()
        if show_running==True:
            print('\tEstimating Angle...')
        Angle = fit_PSF_PA.main_fast(file_Image, input_psf_file, xc_psf, yc_psf, xc, yc, sky_level, xmin, ymin, xmax, ymax)
        break
    return Angle

def add_keyw(input_image, output_image, xc, yc, sky_level, sky_std, Angle, fit_sky=1):
    #http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4h.html
    hdulist = pyfits.open(input_image)
    primhdr = pyfits.getheader(input_image, 0)
    inframe = hdulist[0].data
    ySize, xSize = inframe.shape

    outframe = np.copy(inframe)

    primhdr['SCALE'] = 1.375
    primhdr['FWHM'] = 4.42	# sma = 6.08 arc (W1), http://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec4_4c.html
    primhdr['GAIN'] = 3.20
    primhdr['RDNOISE'] = 3.09
    primhdr['NCOMBINE'] = 1
    primhdr['EXPTIME'] = 7.7
    primhdr['M0'] = 20.5	# Was 20.752!!!!
    primhdr['SKY_LEVEL'] = sky_level
    primhdr['SKY_STD'] = sky_std
    primhdr['SKY_SUBTR'] = fit_sky
    primhdr['XC'] = xc
    primhdr['YC'] = yc
    primhdr['PSF_angle'] = Angle

    hdu = pyfits.PrimaryHDU(outframe, primhdr)
    hdu.writeto(output_image, clobber=True)


def main(galaxy_name, steps):

    if -2 in steps:
        # Get the coordinates general angular characteristics of the galaxy:
        os.chdir(galaxy_name)
        RA_object,DEC_object,sma_object,smb_object,PA_object,q_ratio,d25 = get_info_about_galaxy(galaxy_name)
        #TODO: write these to json
        os.chdir('..')


    if -1 in steps:
        # Automatically download the needed WISE frames:
        #TODO: Should we do this in the galaxy directory?
        os.chdir(galaxy_name)
        input_int_image,input_unc_image = download_wise_frames.main(data['pgc'][k], RA_object, DEC_object, 20., int_dir='Int', unc_dir = 'Unc', crop=True) #TODO: Find this script
        os.chdir('..')


    if 0 in steps:
        # Fix nans and header
        os.chdir(galaxy_name)
        wise_image = glob.glob('null*-int-*.fits')[0]
        wise_unc = glob.glob('null*-unc-*.fits')[0]
        fix_images(galaxy_name, wise_image, wise_unc)


    if 1 in steps:
        # General Masking
        os.chdir(galaxy_name)
        wise_image = galaxy_name + '_WISE_W1_int.fits'
        wise_unc = galaxy_name + '_WISE_W1_unc.fits'
        auto_masking.main(wise_image, output_region_file='general_mask.reg', snr=1.5, min_pix=2, region_type='polygon', sextr_setup='wise.sex', galaxy_ellipse=None, offset_size=1.2, offset_pix=0., verbosity=True) #TODO: Check this function - was wise_masking.sex. Compare with the commented line below!
        convert_reg_to_mask.mask(wise_image, 'general_mask.reg', output_image=None, output_mask='general_mask.fits', mask_value=1, verbosity=True, mask_DN=None)
        os.chdir('..')

        #imp_masking.main(input_int, output_region_file='mask_tmp.reg', snr=2., min_pix=5, region_type='polygon', sextr_setup='wise.sex', galaxy_ellipse=None, offset_size=1., offset_pix=0., verbosity=True)


    if 2 in steps:
        # Sky subtraction
        os.chdir(galaxy_name)
        wise_image = galaxy_name + '_WISE_W1_int.fits'
        wise_unc = galaxy_name + '_WISE_W1_unc.fits'
        #auto_masking.main(wise_image, output_region_file='general_mask.reg', snr=1.5, min_pix=2, region_type='polygon', sextr_setup='wise_masking.sex', galaxy_ellipse=None, offset_size=1.2, offset_pix=0., verbosity=True)
        determine_sky.sky_subtraction(wise_image, 'general_mask.fits', polynomial_degree=0, output_image='sky_subtr.fits', output_sky=None, hdu_inp=0, sampling_factor=1., sigma_smooth=None, verbosity=True, sky_value='mode')
        os.chdir('..')


    if 3 in steps:
        # PSF PA fiting
        os.chdir(galaxy_name)
        # Find the PSF angle of the PSF image and rotate the PSF image accordingly:
        #TODO: Read in RA_object and DEC_object from the json-file
        #TODO: path to input_psf_file
        Angle = find_PA('sky_subtr.fits', input_psf_file, RA_object, DEC_object, show_running=False)
        rotate_psf.rotate_psf(input_psf_file, 120., 120., Angle)
        os.chdir('..')


    if 4 in steps:
        # Get galaxy ellipse and final masking
        os.chdir(galaxy_name)
        xc,yc = get_galaxy_center_ned.main('sky_subtr.fits', name=None, RA=RA_object, DEC=DEC_object, verbosity=True)
        [xc,yc],sma,smb,PA = get_galaxy_ellipse.main('sky_subtr.fits', 'segm.fits', xc=xc, yc=yc)

        # Do final masking: inner and outer
        # TODO get sky_std from json
        IMP.IMP_main('sky_subtr.fits', 'final_masking', region_file='galaxy_ellipse.reg', psf_image='rot_psf.fits', output_image='galaxy_clean.fits', mask_image='mask.fits',sky_std=sky_std, user_interact=True, verbosity=True) # 'galaxy_mask.reg' will be created
        os.chdir('..')


    if 5 in steps:
        # Add keywords:
        #TODO: get nx,ny
        add_keyw('galaxy_clean.fits', 'galaxy_clean_galf.fits', xc, yc, 0., sky_std, Angle, 1)


    if 6 in steps:
        # Fitting:

        #TODO: get
        m0
        pixel_scale
        Distance
        Scale
        galaxy_name
        # Get borders for fitting:
        [[x_min, y_min], [x_max, y_max]] = crop_galaxy_image.ellipse_borders([xc,yc], sma, smb, PA, nx=nx, ny=ny, square=False)

        deca_lite.main('sky_subtr.fits', wise_unc, 'rot_psf.fits', 'mask.fits',
            m0, pixel_scale,
            xc=xc, yc=yc,
            comps = ['edgedisk','sersic'],
            comp_names = ['disk','bulge'],
            sampling=1,
            C0 = None,
            initial_pars=None,
            xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max,
            verbosity=True,
            plot_graphs=True,
            mask_cen_for_pure_disks=False,
            Distance=Distance,
            Scale=Scale,
            Aext=0.,
            Kcorr=0.,
            object_name=galaxy_name,
            run_sersic=True)



    # END #############














    if 3 in steps:
        # Determine the size of the galaxy
        print('Use a circle DS9 region to outline the edge of the galaxy. The center of the circe should coincide with the galaxy center (~the brightest central pixel). Save the region in a file called galaxy.reg')

        f=open('galaxy.reg', 'a').close()

        if opsyst=='Linux':
            ds9Proc = subprocess.Popen(["ds9", 'sky_subtr.fits',
                                                "-scale", "histequ"])
            ds9Proc.wait()
        elif opsyst=='Darwin':
            ds9Proc = subprocess.Popen(["/Applications/SAOImageDS9.app/Contents/MacOS/ds9", 'sky_subtr.fits',
                                                "-scale", "histequ"])
            ds9Proc.wait()



    if 5 in steps:
        # Do final masking: inner and outer
        os.chdir(galaxy_name)
        hdulist = pyfits.open('sky_subtr.fits')
        header = hdulist[0].header

        sky_std = header['Sky_std']

        IMP_main('sky_subtr.fits', 'final_masking', region_file='galaxy_ellipse.reg', psf_image='psf.fits', output_image='galaxy_clean.fits', mask_image='final_mask.fits',sky_std=sky_std, user_interact=True, verbosity=True)

        os.chdir('..')


















    if 22 in steps:
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
        input_images = glob.glob('*_reduced_crop_rayremoved_skysub.fits')
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
    #print("Photometric Data Reduction:")
    #print("Step 1: Initial data reduction to correct for calibration frames.")
    #print("Step 2: Automated cropping of the galaxy image [x_l, y_l: x_r, y_r] to remove a frame with bad pixels.")
    #print("Step 3: Automated removement of cosmic rays.")
    #print("Step 4: Sky subtraction - under human supervision.")
    #print("Step 5: Median stacking.")
    
    parser = argparse.ArgumentParser(description="Pipeline for WISE preparation and fitting:")
    parser.add_argument("--steps", nargs='?', const=1, help="Specify steps separated by comma", type=str, default='1,2,3,4,5')
    parser.add_argument("--galaxy", nargs='?', const=1, help="Provide a galaxy name", type=str, default=None)    
    #parser.add_argument("--sci", nargs='?', const=1, help="Provide the path to science frames", type=str, default=None) 
    #parser.add_argument("--sci_prefix", nargs='?', const=1, help="Provide the prefix for science frames", type=str, default=None) 
    #parser.add_argument("--flats_prefix", nargs='?', const=1, help="Provide the prefix for flat fields", type=str, default='*flat.g*')
    #parser.add_argument("--bias_prefix", nargs='?', const=1, help="Provide the prefix for bias frames", type=str, default='*bias*.fits')
    #parser.add_argument("--dark_prefix", nargs='?', const=1, help="Provide the prefix for dark frames", type=str, default='*dark*.fits')
    
    args = parser.parse_args()
    
    
    steps = np.array(args.steps.split(','), int)
    galaxy = args.galaxy
    
    main(galaxy, steps)




