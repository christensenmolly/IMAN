#!/usr/bin/python

# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import math
import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from photutils import detect_sources, segmentation
import argparse
import pyregion
import os
import astropy.units as u
import scipy.ndimage as ndi
from skimage import data, filters, measure, morphology
from skimage.measure import label, regionprops, regionprops_table
#import mask_objects


def unmask_galaxy(segm, xc=None, yc=None, min_radius=10.):
    # Find the masked area which corresponds to the galaxy and unmask it. All objects within at least 2*min_radius will be unmasked.

    segm_without_galaxy = np.copy(segm)
    
    
    ny,nx = np.shape(segm_without_galaxy)
    
    segm_only_galaxy = np.zeros(shape=(ny,nx))
    
    if xc is None or yc is None:
        xc = nx/2.
        yc = ny/2.

    labels = []
    for y in range(int(yc-2.*min_radius), int(yc+2.*min_radius)):
        for x in range(int(xc-2.*min_radius), int(xc+2.*min_radius)):
            try:
                lable = segm_without_galaxy[y,x]
                if lable not in labels and lable!=0:
                    labels.append(lable)
            except:
                z=1 # Beyond the image borders

    for label in labels:
        segm_without_galaxy[segm_without_galaxy == label] = 0.
        segm_only_galaxy[segm == label] = 1.

    hdu = pyfits.PrimaryHDU(segm_without_galaxy)
    hdu.writeto('no_galaxy.fits', clobber=True)


    hdu = pyfits.PrimaryHDU(segm_only_galaxy)
    hdu.writeto('with_galaxy.fits', clobber=True)

    return segm_without_galaxy, segm_only_galaxy  



def create_ellipse_region(input_image, xc, yc, sma, smb, PA, file, coord_format='image'):
    if coord_format=='image':
        f_galaxy = open(file, 'w')
        f_galaxy.write('%s\n' % ('image') )
        f_galaxy.write('ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # color=red\n' % (xc, yc, sma, smb, PA))
        f_galaxy.close()
    else:
        region = 'image;ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f)' % (xc, yc, sma, smb, PA)
        r = pyregion.parse(region)
        
        f = pyfits.open(input_image)
        r2 = pyregion.parse(r).as_imagecoord(f[0].header)
        
        f_galaxy = open(file, 'w')
        f_galaxy.write('%s\n' % ('fk5') )
        f_galaxy.write('ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # color=red\n' % (xc, yc, sma, smb, PA))
        f_galaxy.close()        

def convert_segm_to_boolean(mask):
    return np.ma.make_mask(mask)

def main(input_image, segm_image, I_lim=None, output_region='galaxy_ellipse.reg', xc=None, yc=None, min_radius=10., verbosity=True, method='iso', filter_size=10., use_input_coords=False, radius_factor=1.):
    '''
    Function to determine the ellipse which fits the outermost galaxy isophote at signal-to-noise ratio equal snr.
    The galaxy has a center of xc, yc. If they are None, then the center of the image is taken as the galaxy center.
    min_radius says that the galaxy semi-major axis should be larger than this when searching for the outermost ellipse.
    '''

    hdulist = pyfits.open(input_image)
    data = hdulist[0].data  

    if method=='SE':
            hdulist_segm = pyfits.open(segm_image)
            segm = hdulist_segm[0].data    
        
        
            ny,nx = np.shape(data)
            if xc is None or yc is None:
                xc = nx/2.
                yc = ny/2.

            # Unmask the target galaxy
            segm_without_galaxy, segm_only_galaxy  = unmask_galaxy(segm, xc=xc, yc=yc, min_radius=min_radius)
            
            
            segmap_float = ndi.uniform_filter(np.float64(segm_only_galaxy), size=filter_size)

            segm_only_galaxy = segmap_float > 0.5 #WARNING: NEW
            segm_only_galaxy = segm_only_galaxy.astype(int)

            # Convert to boolean mask
            mask = convert_segm_to_boolean(segm_without_galaxy)

            # Coarsely determine the galaxy ellipse.
            cat = segmentation.source_properties(data, segm_only_galaxy, mask=mask) # data should be sky subtracted!!!

            r = 3.    # approximate isophotal extent
            #apertures = []
            position = []; a = []; b = []; theta = []
            for obj in cat:
                #print(obj)
                position.append((obj.xcentroid.value, obj.ycentroid.value))
                a.append(obj.semimajor_axis_sigma.value * r)
                b.append(obj.semiminor_axis_sigma.value * r)
                theta.append(obj.orientation.to(u.rad).value)
                print(position,a,b,theta)
            try:
                # Choose the largest object in the table is the galaxy
                ind = np.argmax(a)

                PA = np.degrees(theta[ind]%(2*math.pi))
                xc_fit = position[ind][0]
                yc_fit = position[ind][1]
                sma = a[ind]
                smb = b[ind]
            except:
                PA = 0.
                xc_fit = nx/2.
                yc_fit = ny/2.
                sma = nx/4.+1.
                smb = ny/4.    

            os.remove('no_galaxy.fits')
            os.remove('with_galaxy.fits')

            
    if method=='iso':
        import create_cleaned_smooth_image
        cleaned_img, cleaned_mask = create_cleaned_smooth_image.main(input_image, I_lim, filter_size_img=None, filter_size_mask=filter_size, output_mask = 'mask_test.fits', output_image = 'img_test.fits')
    
        props = measure.regionprops(cleaned_mask)
        properties = ['centroid', 'orientation', 'axis_major_length', 'axis_minor_length']
        
        Largest_radius = 0.
        for index in range(1, cleaned_mask.max()):
            (xccc,yccc) = getattr(props[index], 'centroid')
            SSMA = getattr(props[index], 'axis_major_length')/2.
            SSMB = getattr(props[index], 'axis_minor_length')/2.
            #print(xccc, yccc, xc, yc)
            if getattr(props[index], 'axis_major_length')>Largest_radius and math.sqrt((xccc-xc)**2 + (yccc-yc)**2)<math.sqrt(SSMA*SSMB):
                Largest_radius = SSMA
                Largest_index = index
        
        (xc_fit,yc_fit) = getattr(props[Largest_index], 'centroid')
        sma = getattr(props[Largest_index], 'axis_major_length')/2.
        smb = getattr(props[Largest_index], 'axis_minor_length')/2.
        PA = 90.-np.degrees(getattr(props[Largest_index], 'orientation'))


    '''

    
    exit()

    from photutils.segmentation import SourceCatalog
    cat = SourceCatalog(data, segm)
    tbl = cat.to_table()
    tbl['xcentroid'].info.format = '.2f'  # optional format
    tbl['ycentroid'].info.format = '.2f'
    tbl['kron_flux'].info.format = '.2f'
    print(tbl)
    exit()
    

    import matplotlib.pyplot as plt
    from astropy.visualization import SqrtStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    norm = ImageNormalize(stretch=SqrtStretch())
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
    ax1.imshow(segm_only_galaxy, origin='lower', cmap='Greys_r', norm=norm)
    ax1.set_title('Data')
    cmap = segm.make_cmap(seed=123)
    ax2.imshow(segm, origin='lower', cmap=cmap, interpolation='nearest')
    ax2.set_title('Segmentation Image')
    plt.show()
    '''

    sma = sma * radius_factor
    smb = smb * radius_factor
    
    if verbosity: print('Galaxy ellipse:') 
    if verbosity: print('\txc,yc: %.1f, %.1f' % (xc_fit, yc_fit)) 
    if verbosity: print('\tsma [pix]: %.1f' % (sma))
    if verbosity: print('\tsmb [pix]: %.1f' % (smb))
    if verbosity: print('\tell: %.2f' % (1.-smb/sma))    
    if verbosity: print('\tPA [Degrees: Up=90, Right=0, counterclockwise]: %.1f' % (PA)) 
    
    if use_input_coords==False:
        create_ellipse_region(input_image, xc_fit, yc_fit, sma, smb, PA, output_region)
        
        return [xc_fit,yc_fit],sma,smb,PA
    else:
        create_ellipse_region(input_image, xc, yc, sma, smb, PA, output_region)
        
        return [xc,yc],sma,smb,PA        
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Determine galaxy ellipse")
    parser.add_argument("inputImage", help="Input image")
    parser.add_argument("segmImage", help="Segmentation image with marked objects")
    parser.add_argument("--xcen", nargs='?', const=1, help="Optional: Galaxy center x",type=float,default=None)  
    parser.add_argument("--ycen", nargs='?', const=1, help="Optional: Galaxy center y",type=float,default=None)  
    parser.add_argument("--bckg", nargs='?', const=1, help="Optional: Sky background",type=float,default=0.)  
    parser.add_argument("--min_radius", nargs='?', const=1, help="Optional: Minumum radius of the galaxy to be determined",type=float,default=10.)  
    args = parser.parse_args()

    input_image = args.inputImage
    segm_image = args.segmImage
    xcen = args.xcen
    ycen = args.ycen
    bckg = args.bckg
    min_radius = args.min_radius    
    
    main(input_image, segm_image, xc=xcen, yc=ycen, sky_background=bckg, min_radius=min_radius)
