import argparse
import numpy as np
from ctypes import c_float, c_double
from mtolib.io_mto import read_fits_file
import mtolib.main as mto
from mtolib import _ctype_classes as ct
from mtolib.tree_filtering import init_double_filtering
from types import SimpleNamespace


def main(filename, out='segm_mto.fits', par_out='parameters.csv', soft_bias=None, gain=-1, bg_mean=None, bg_variance=-1, alpha=1e-6, move_factor=0.3, min_distance=0.0, verbosity=True):

    if verbosity: print('Creating segmentation map using MTO (Teeninga et al. 2011) ...')
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', action='store_true')
    parser.add_argument('-out', action='store_true')
    parser.add_argument('-par_out', action='store_true')
    parser.add_argument('-soft_bias', action='store_true')
    parser.add_argument('-gain', action='store_true')
    parser.add_argument('-bg_mean', action='store_true')
    parser.add_argument('-bg_variance', action='store_true')
    parser.add_argument('-alpha', action='store_true')
    parser.add_argument('-move_factor', action='store_true')
    parser.add_argument('-min_distance', action='store_true')
    parser.add_argument('-verbosity', action='store_true')
        
    params = parser.parse_args()
    '''
    #print(soft_bias,gain, bg_mean, bg_variance, alpha, move_factor, min_distance, verbosity)
    #exit()
    

    params = SimpleNamespace(filename = filename, out=out, par_out = par_out, soft_bias=soft_bias, gain=gain, bg_mean=bg_mean, bg_variance=bg_variance, alpha=alpha, move_factor=move_factor, min_distance=min_distance, verbosity=verbosity)
    
    '''
    params.filename = filename
    params.out = out
    params.par_out = par_out
    params.soft_bias = soft_bias
    params.gain = gain
    params.bg_mean = bg_mean
    params.bg_variance = bg_variance
    params.alpha = alpha
    params.move_factor = move_factor
    params.min_distance = min_distance
    params.verbosity = verbosity
    '''

    # Warn if using default soft bias
    if params.soft_bias is None:
        params.soft_bias = 0.0

    img = read_fits_file(params.filename)

    if verbosity:
        print("\n---Image dimensions---")
        print("Height = ", img.shape[0])
        print("Width = ", img.shape[1])
        print("Size = ", img.size)

    # Set the pixel type based on the type in the image
    params.d_type = c_float
    if np.issubdtype(img.dtype, np.float64):
        params.d_type = c_double
        init_double_filtering(params)

    # Initialise CTypes classes
    ct.init_classes(params.d_type)


    # Pre-process the image
    processed_image,bg_mean,bg_variance,gain = mto.preprocess_image(img, params, n=2)

    # Build a max tree
    mt = mto.build_max_tree(processed_image, params)

    # Filter the tree and find objects
    id_map, sig_ancs = mto.filter_tree(mt, processed_image, params)


    # Relabel objects for clearer visualisation
    id_map = mto.relabel_segments(id_map, shuffle_labels=False)

    # Generate output files
    mto.generate_image(img, id_map, params)
    mto.generate_parameters(img, id_map, sig_ancs, params)   
    if verbosity: print('Done!')
    return gain
    
    
#main('1_r_trim.fits')
