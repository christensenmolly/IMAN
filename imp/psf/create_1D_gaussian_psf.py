import numpy as np

def gaussian_profile(I0, FWHM, radius):
    r = np.arange(1, radius)
    sigma = FWHM/2.354
    return r, I0 * np.exp(-r**2/(2.*sigma**2))
    
    
    
def main(I0, FWHM, radius, file='gaussian_1D.txt'):
    r,I = gaussian_profile(I0, FWHM, radius)
    
    I = I/np.sum(I)
    
    f_res = open(file, 'w') #### Output model is saved in a text-file as well
    f_res.truncate(0)
    f_res.write("# sma[pix]\tflux[DN]\tflux_err[DN]\n")
    
    for k in range(len(r)):
        f_res.write("%7.2f\t%15.8e\t%15.8e\n" % (r[k], I[k], 0.))
    
    f_res.close()

main(1., 5., 16)
