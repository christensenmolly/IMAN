#!/usr/bin/env python
# -*- coding: utf8 -*-


import numpy as np


            

def main(array, max_iter, tol, kernel_size=1, method='localmean'):
    # Initialize arrays
    filled = np.empty( [array.shape[0], array.shape[1]], dtype=np.float64)
    kernel = np.empty( (2*kernel_size+1, 2*kernel_size+1), dtype=np.float64)

    # Indices where array is NaN
    inans, jnans = np.nonzero( np.isnan(array) )
    
    # Number of NaN elements
    n_nans = len(inans)
    
    # Arrays which contain replaced values to check for convergence
    replaced_new = np.zeros( n_nans, dtype=np.float64)
    replaced_old = np.zeros( n_nans, dtype=np.float64)

    # Depending on kernel type, fill kernel array
    if method == 'localmean':
        for i in range(2*kernel_size+1):
            for j in range(2*kernel_size+1):
                kernel[i,j] = 1.

    elif method == 'idw':
        kernel = np.array([[0, 0.5, 0.5, 0.5,0],
                  [0.5,0.75,0.75,0.75,0.5],
                  [0.5,0.75,1,0.75,0.5],
                  [0.5,0.75,0.75,0.5,1],
                  [0, 0.5, 0.5 ,0.5 ,0]])
        #print kernel, 'kernel'

    else:
        raise ValueError( 'method not valid. Should be one of `localmean`.')
    
    # Fill new array with input elements
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            filled[i,j] = array[i,j]

    # Make several passes
    # until we reach convergence
    for it in range(max_iter):
        
        # for each NaN element
        for k in range(n_nans):
            i = inans[k]
            j = jnans[k]
            
            # Initialize to zero
            filled[i,j] = 0.0
            n = 0
            
            # Loop over the kernel
            for I in range(2*kernel_size+1):
                for J in range(2*kernel_size+1):
                   
                    # if we are not out of the boundaries
                    if i+I-kernel_size < array.shape[0] and i+I-kernel_size >= 0:
                        if j+J-kernel_size < array.shape[1] and j+J-kernel_size >= 0:
                                                
                            # if the neighbour element is not NaN itself.
                            if filled[i+I-kernel_size, j+J-kernel_size] == filled[i+I-kernel_size, j+J-kernel_size] :
                                
                                # do not sum itself
                                if I-kernel_size != 0 and J-kernel_size != 0:
                                    
                                    # convolve kernel with original array
                                    filled[i,j] = filled[i,j] + filled[i+I-kernel_size, j+J-kernel_size]*kernel[I, J]
                                    n = n + 1

            # Divide value by effective number of added elements
            if n != 0:
                filled[i,j] = filled[i,j] / n
                replaced_new[k] = filled[i,j]
            else:
                filled[i,j] = 0.0
                
        # Check if mean square difference between values of replaced
        # elements is below a certain tolerance
        if np.mean( (replaced_new-replaced_old)**2 ) < tol:
            break
        else:
            for l in range(n_nans):
                replaced_old[l] = replaced_new[l]
    
    return filled
