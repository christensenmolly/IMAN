
# Import the necessary modules
from astropy.io import fits as pyfits
import numpy as np
import math
import argparse
#    noise^2 = object_flux + sky + rdnoise^2
#
# Since sigma_adu = sigma_e/gain, we can go from
#    noise(e-)^2 = object_flux(e-) + sky(e-) + rdnoise^2
# to
#    noise(adu)^2 = object_flux(adu)/gain + sky(adu)/gain + rdnoise^2/gain^2
# or just
#    noise(adu)^2 = (object_flux(adu) + sky(adu))/gain + rdnoise^2/gain^2
# (assuming that read noise is in units of e-, as is usual)
#
# Exposure time and number of averaged images can be accounted for by including them in 
# the effective gain:
#    gain_eff = (gain * t_exp * N_combined)
# HOWEVER, in this case we also have to account for the multiple readouts, which means
# that the read noise term is multiplied by N_combined, so that we end up with
#    noise(adu)^2 = (object_flux(adu) + sky(adu))/gain_eff + N_combined * rdnoise^2/gain_eff^2
# (where "adu" can be adu/sec if t_exp != 1)



def create_weight_image(input_image, output_image='weight.fits', input_sky=0., gain=1., t_exp=1., N_combined=1, rdnoise=10.):
        print('Creating weight image as 1/sigma...')
        HDU = pyfits.open(input_image)
        Data = HDU[0].data
        header = HDU[0].header
        ySize, xSize = Data.shape
        
        try:
            HDUsky = pyfits.open(input_sky)
            Sky = HDUHDUsky[0].data
        except:
            Sky = float(input_sky)

        gain_eff = (gain * t_exp * N_combined)
        totalFlux = Data + Sky
        totalFlux[totalFlux<0.] = 0.

        noise_squared = totalFlux/gain_eff + N_combined*(rdnoise**2/(gain_eff**2))
        weight = 1.0 / np.sqrt(noise_squared)
        
        outHDU = pyfits.PrimaryHDU(weight, header=header)
        outHDU.writeto(output_image, clobber=True)   
        print('Done!')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky background estimation")
    parser.add_argument("inputImage", help="Input fits image in ADU")
    parser.add_argument("--outputImage", help="Optional: Output weight image",type=str,default='weight.fits') 
    parser.add_argument("--input_sky", nargs='?', const=1, help="Optional: Input sky (subtracted): image of value", type=str, default='0.')
    parser.add_argument("--gain", nargs='?', const=1, help="Optional: gain (e-/ADU)", type=float, default=1.)
    parser.add_argument("--t_exp", nargs='?', const=1, help="Optional: one exposure time t_exp (sec)", type=float, default=1.)
    parser.add_argument("--n", nargs='?', const=1, help="Optional: number of combined images", type=int, default=1)
    parser.add_argument("--ron", nargs='?', const=1, help="Optional: read-out-noise (e-)", type=float, default=10.)
    
    
    args = parser.parse_args()

    input_image = args.inputImage
    output_image = args.outputImage
    input_sky = args.input_sky
    gain = args.gain
    t_exp = args.t_exp
    N_combined = args.n
    rdnoise = args.ron

    create_weight_image(input_image, output_image=output_image, input_sky=input_sky, gain=gain, t_exp=t_exp, N_combined=N_combined, rdnoise=rdnoise)
