from astropy.io import fits
import astroalign


### JUST STARTED!!!
# Ideas:
# https://pypi.org/project/astrometry/#get-started - only works with coordinates
# astroalign - must provide a reference reference image!











data = fits.getdata('skybackg_sibtr.fit', ext=0)
reference_img = fits.getdata(’path_to_ref_image’,ext=0)
reference_head = fits.getheader(’path_to_ref_image’,ext=0)
aligned_image = astroalign.register(data,reference_img)
fits.writeto(’path_to_aligned_image.fits’,aligned_image,reference_head,overwrite=True)
