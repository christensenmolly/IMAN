#! /usr/bin/env python

# Author: Victor Terron (c) 2015
# Email: `echo vt2rron1iaa32s | tr 132 @.e`
# License: GNU GPLv3

""" Python interface to IRAF's imcombine (via PyRAF) """

import sys
import tempfile
import os

try:
    import pyraf.iraf
except:
    z=1
#from pyraf.iraf import noao, artdata  # 'noao.artdat


def imcombine(images, output_path, flow, fhigh):
    """ Mim-max average combine all the images in the set.

    The method invokes the IRAF task 'imcombine' on the FITS images contained
    in 'images', an iterable. The type of combining operation performed on the
    pixels (the 'combine' parameter when 'imcombine' is invoked) is 'average',
    while the rejection operation (the 'reject' parameter) is 'minmax'. That is
    the reason why 'flow' and 'fhigh', which determine the fraction of low and
    high pixels, respectively, that are to be rejected, are required as
    arguments. If existing, the output image is silently overwritten.

    The WCS in the image is used to derive the offsets.

    """

    if not len(images):
        raise ValueError("no FITS images given")

    if not 0 <= flow <= 1 or not 0 <= fhigh <= 1:
        raise ValueError("both 'nlow' and 'fhigh' must be in the range [0,1]")

    # Number of low and high pixels are to reject; round to nearest integer
    nlow  = int(round(len(images) * flow))
    nhigh = int(round(len(images) * fhigh))

    # If only two images are combined, we always discard the highest pixel
    # (since some stars could be discernible in the case of flat-fields,
    # for example). When there are three or more images we also discard the
    # lowest and highest pixels, unless the specified fraction is zero: in
    # that case, it is understood as to mean "do not discard anything".
    if len(images) == 2 and not nhigh and fhigh:
        nhigh = 1
    elif len(images) >= 3:
        if not nlow and flow:
            nlow = 1
        if not nhigh and fhigh:
            nhigh = 1

    # There is what seems to be a bug in PyRAF that makes it impossible to
    # imcombine more than approximately forty images, so we save the list
    # of input images to a temporary file, use it as the input to IRAF's
    # imcombine and delete it before leaving the method.
    
    
    input_fd, input_path = tempfile.mkstemp(suffix = '.lst', text = True)
    #print(input_fd, input_path)
    #exit()
    try:
        for path in images:
            os.write(input_fd, str.encode('{0}\n'.format(path)))
        os.close(input_fd)

        if os.path.exists(output_path):
            os.unlink(output_path)

        pyraf.iraf.imcombine('@' + input_path,
                             output_path,
                             combine = 'average',   # average by default
                             reject = 'none',
                             nlow = nlow,
                             nhigh = nhigh,
                             offsets = 'wcs',
                             )
    finally:
        os.unlink(input_path)

if __name__ == "__main__":

    if len(sys.argv) < 3:
        msg = "usage: {0} INPUT_IMGS... OUTPUT_IMG".format(sys.argv[0])
        sys.exit(msg)

    input_paths  = sys.argv[1:-1]
    output_path = sys.argv[-1]
    imcombine(input_paths, output_path, 0.1, 0.1)