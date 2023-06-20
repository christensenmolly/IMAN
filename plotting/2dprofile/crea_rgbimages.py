from distutils import version
import os
import warnings

import tempfile
import shutil

import numpy as np
from astropy import log
from astropy.io import fits

from aplpy import image_util
from aplpy import math_util

def _data_stretch(image, vmin=None, vmax=None, pmin=0.25, pmax=99.75, \
                  stretch='linear', vmid=None, exponent=2):

    min_auto = not math_util.isnumeric(vmin)
    max_auto = not math_util.isnumeric(vmax)

    if min_auto or max_auto:
        auto_v = image_util.percentile_function(image)
        vmin_auto, vmax_auto = auto_v(pmin), auto_v(pmax)

    if min_auto:
        log.info("vmin = %10.3e (auto)" % vmin_auto)
        vmin = vmin_auto
    else:
        log.info("vmin = %10.3e" % vmin)

    if max_auto:
        log.info("vmax = %10.3e (auto)" % vmax_auto)
        vmax = vmax_auto
    else:
        log.info("vmax = %10.3e" % vmax)

    image = (image - vmin) / (vmax - vmin)

    data = image_util.stretch(image, stretch, exponent=exponent, midpoint=vmid)

    data = np.nan_to_num(data)
    data = np.clip(data * 255., 0., 255.)

    return data.astype(np.uint8)


def make_rgb_image(data, output, indices=(0, 1, 2), \
                   vmin_r=None, vmax_r=None, pmin_r=0.25, pmax_r=99.75, \
                   stretch_r='linear', vmid_r=None, exponent_r=2, \
                   vmin_g=None, vmax_g=None, pmin_g=0.25, pmax_g=99.75, \
                   stretch_g='linear', vmid_g=None, exponent_g=2, \
                   vmin_b=None, vmax_b=None, pmin_b=0.25, pmax_b=99.75, \
                   stretch_b='linear', vmid_b=None, exponent_b=2, \
                   embed_avm_tags=True):
    '''
    Make an RGB image from a FITS RGB cube or from three FITS files.

    Parameters
    ----------

    data : str or tuple or list
        If a string, this is the filename of an RGB FITS cube. If a tuple
        or list, this should give the filename of three files to use for
        the red, green, and blue channel.

    output : str
        The output filename. The image type (e.g. PNG, JPEG, TIFF, ...)
        will be determined from the extension. Any image type supported by
        the Python Imaging Library can be used.

    indices : tuple, optional
        If data is the filename of a FITS cube, these indices are the
        positions in the third dimension to use for red, green, and
        blue respectively. The default is to use the first three
        indices.

    vmin_r, vmin_g, vmin_b : float, optional
        Minimum pixel value to use for the red, green, and blue channels.
        If set to None for a given channel, the minimum pixel value for
        that channel is determined using the corresponding pmin_x argument
        (default).

    vmax_r, vmax_g, vmax_b : float, optional
        Maximum pixel value to use for the red, green, and blue channels.
        If set to None for a given channel, the maximum pixel value for
        that channel is determined using the corresponding pmax_x argument
        (default).

    pmin_r, pmin_r, pmin_g : float, optional
        Percentile values used to determine for a given channel the
        minimum pixel value to use for that channel if the corresponding
        vmin_x is set to None. The default is 0.25% for all channels.

    pmax_r, pmax_g, pmax_b : float, optional
        Percentile values used to determine for a given channel the
        maximum pixel value to use for that channel if the corresponding
        vmax_x is set to None. The default is 99.75% for all channels.

    stretch_r, stretch_g, stretch_b : { 'linear', 'log', 'sqrt', 'arcsinh', 'power' }
        The stretch function to use for the different channels.

    vmid_r, vmid_g, vmid_b : float, optional
        Baseline values used for the log and arcsinh stretches. If
        set to None, this is set to zero for log stretches and to
        vmin - (vmax - vmin) / 30. for arcsinh stretches

    exponent_r, exponent_g, exponent_b : float, optional
        If stretch_x is set to 'power', this is the exponent to use.

    embed_avm_tags : bool, optional
        Whether to embed AVM tags inside the image - this can only be done for
        JPEG and PNG files, and only if PyAVM is installed.
    '''

    try:
        from PIL import Image
    except ImportError:
        try:
            import Image
        except ImportError:
            raise ImportError("The Python Imaging Library (PIL) is required to make an RGB image")

    if isinstance(data, basestring):

        image = fits.getdata(data)
        image_r = image[indices[0], :, :]
        image_g = image[indices[1], :, :]
        image_b = image[indices[2], :, :]

        # Read in header
        header = fits.getheader(data)

        # Remove information about third dimension
        header['NAXIS'] = 2
        for key in ['NAXIS', 'CTYPE', 'CRPIX', 'CRVAL', 'CUNIT', 'CDELT', 'CROTA']:
            for coord in range(3, 6):
                name = key + str(coord)
                if name in header:
                    header.__delitem__(name)

    elif (type(data) == list or type(data) == tuple) and len(data) == 3:

        filename_r, filename_g, filename_b = data
        image_r = fits.getdata(filename_r)
        image_g = fits.getdata(filename_g)
        image_b = fits.getdata(filename_b)

        # Read in header
        header = fits.getheader(filename_r)

    else:
        raise Exception("data should either be the filename of a FITS cube or a list/tuple of three images")

    #log.info("Red:")
    image_r = Image.fromarray(_data_stretch(image_r, \
                                            vmin=vmin_r, vmax=vmax_r, \
                                            pmin=pmin_r, pmax=pmax_r, \
                                            stretch=stretch_r, \
                                            vmid=vmid_r, \
                                            exponent=exponent_r))

    #log.info("Green:")
    image_g = Image.fromarray(_data_stretch(image_g, \
                                            vmin=vmin_g, vmax=vmax_g, \
                                            pmin=pmin_g, pmax=pmax_g, \
                                            stretch=stretch_g, \
                                            vmid=vmid_g, \
                                            exponent=exponent_g))

    #log.info("Blue:")
    image_b = Image.fromarray(_data_stretch(image_b, \
                                            vmin=vmin_b, vmax=vmax_b, \
                                            pmin=pmin_b, pmax=pmax_b, \
                                            stretch=stretch_b, \
                                            vmid=vmid_b, \
                                            exponent=exponent_b))

    img = Image.merge("RGB", (image_r, image_g, image_b))
    img = img.transpose(Image.FLIP_TOP_BOTTOM)

    if output!=None:
        img.save(output)
        
        if embed_avm_tags:

            try:
                import pyavm
            except ImportError:
                warnings.warn("PyAVM 0.9.1 or later is not installed, so AVM tags will not be embedded in RGB image")
                return

            if version.LooseVersion(pyavm.__version__) < version.LooseVersion('0.9.1'):
                warnings.warn("PyAVM 0.9.1 or later is not installed, so AVM tags will not be embedded in RGB image")
                return

            from pyavm import AVM

            if output.lower().endswith(('.jpg', '.jpeg', '.png')):
                avm = AVM.from_header(header)
                avm.embed(output, output)
            else:
                warnings.warn("AVM tags will not be embedded in RGB image, as only JPEG and PNG files are supported")
    else:
        return img
