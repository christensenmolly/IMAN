#!/usr/bin/python
# DESCRIPTION:
# Script to rebin image to a reference image.
# MINIMAL USAGE: python  schlafly_extinction.py [input_object] [wavelength_nm]
# EXAMPLE:  python3 schlafly_extinction.py NGC891 0.55
# EXAMPLE:  python3 schlafly_extinction.py 38.583792,32.505611 0.55

from astroquery.irsa_dust import IrsaDust
import astropy.coordinates as coord
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import sys
import argparse
# wavelengths are in um 

def main(pts, wavelengths=None, frame='fk5', units='deg', output_file=None):
  #input_coords = "22h57m57.5s +26d09m00.09s"
  Lambdas = []
  A_Lambdas = []
  for k in range(len(pts)):
    try:
        coords = str(pts[k][0]) + ' ' + str(pts[k][1])
        if units==None:
            C = coord.SkyCoord(coords, frame=frame)
        else:
            C = coord.SkyCoord(coords, frame=frame, unit='deg')
    except:
      C = pts[k][0]
    table = IrsaDust.get_extinction_table(C)
    Lambdas.append(table["LamEff"])
    A_Lambdas.append(table["A_SandF"])
    

  Lambdas = np.array(Lambdas,float)
  A_Lambdas = np.array(A_Lambdas,float)


  F = []
  for k in range(len(Lambdas)):
      lambdas = Lambdas[k]
      A_lambdas = A_Lambdas[k]

      xx = np.arange(min(lambdas),max(lambdas),0.001)
      f2 = interp1d(lambdas,A_lambdas)  

      #https://en.wikipedia.org/wiki/Photometric_system
      A_B = f2(0.445)	# B-band
      A_V = f2(0.551)	# V-band



      # For the UV fluxes:
      #http://arxiv.org/pdf/1301.1427v1.pdf
      #http://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html  - Convert GALEX counts to AB-mag
      E_BV = A_B - A_V
      R_FUV = 4.89 # Alternatively, 7.9
      R_NUV = 7.24 # Alternatively, 8.0
      A_NUV = R_NUV * E_BV
      A_FUV = R_FUV * E_BV
      lambda_FUV = 0.1528
      lambda_NUV = 0.2271

      lambdas = np.insert(lambdas, 0, lambda_NUV)
      lambdas = np.insert(lambdas, 0, lambda_FUV)


      
      A_lambdas = np.insert(A_lambdas, 0, A_NUV)
      A_lambdas = np.insert(A_lambdas, 0, A_FUV)
      
     
      xx = np.arange(min(lambdas),max(lambdas),0.001)
      f2 = interp1d(lambdas,A_lambdas) 

      F.append(f2)
      #plt.plot(xx,f2(xx),'o')
      #plt.show()
  
  if wavelengths==None:
    return F
  else:
    Extinctions = []
    if output_file!=None:
      f = open(output_file, 'w')


    for k in range(len(F)):
      EXT = F[k](wavelengths[k])
      Extinctions.append(EXT)


    if output_file!=None:
        for i in range(len(wavelengths[0])):
          for k in range(len(F)):
            if k!= len(F)-1:
                f.write(str(round(Extinctions[k][i],3))+'\t')
            else:
                f.write(str(round(Extinctions[k][i],3))+'\n')
      
    if output_file!=None:
      f.close()
    return np.array(Extinctions, float)





'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A module to retrieve extinction from Schlafly and Finkbeiner (2011) for a specified object")
    parser.add_argument("input_object", help="Input_object (NED name or coordinates separated by comma)")    
    parser.add_argument("wavelength", help="Wavelength in nm",type=float)
    args = parser.parse_args()

    input_object = args.input_object
    wavelength = args.wavelength

    if ',' in input_object:
        [ra,dec] = input_object.split(',')
        ext = main([[str(ra),str(dec)]], wavelengths=[[wavelength]], frame='fk5', units='deg', output_file=None)
    else:
        ext = main([[input_object]], wavelengths=[[wavelength]], frame='fk5', units='deg', output_file=None)
    print('A at %.3f nm: %.3f mag' % (wavelength, ext[0][0]))
'''





'''
res = main([['22h57m57.5s','+26d09m00.09s'], ['22h57m57.5s','+26d09m00.09s']])
print res
'''
#res = main([[38.583792,32.505611]],wavelengths=[0.655])
#print(res)


'''
Ext = main([['NGC973'],['NGC5907']],wavelengths=[[0.5517,0.54],[0.3,0.5]],output_file='ext.dat')
print(Ext)
'''





