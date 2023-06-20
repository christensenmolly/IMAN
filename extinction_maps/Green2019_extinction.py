from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarQuery
import astropy.units as units
import numpy as np
import os
from scipy.interpolate import interp1d
PATH_TO_MODULE = os.path.dirname(__file__)


def convert_l_b(RA=None, DEC=None):
    c = SkyCoord(ra=RA*units.degree, dec=DEC*units.degree, frame='fk5').galactic
    return c.l.degree,c.b.degree


def return_ext_curve(Rv=3.1):
    # Here Rv_* is Ab/E(B-V)SFD
    lambda_eff,Rv_21,Rv_31,Rv_41,Rv_51 = np.loadtxt(PATH_TO_MODULE+'/Table_6_SF2011.dat', usecols=[1,2,3,4,5], unpack=True, skiprows = 1, dtype=float, delimiter='\t')
    if Rv==2.1:
        Ab_EBV = Rv_21

    if Rv==3.1:
        Ab_EBV = Rv_31

    if Rv==4.1:
        Ab_EBV = Rv_41        

    if Rv==5.1:
        Ab_EBV = Rv_51        

    inds = lambda_eff.argsort()
    Ab_EBV = Ab_EBV[inds]
    lambda_eff = lambda_eff[inds]

    f = interp1d(lambda_eff, Ab_EBV) 
    return f
    


    
def main(pts, wavelength):
    # wavelength in nm
    bayestar = BayestarQuery(max_samples=1)
    res = []
    for k in range(len(pts)):
        #print k
        [ll,bb,distance] = pts[k]
        coords = SkyCoord(ll*units.deg, bb*units.deg, distance=distance*units.pc, frame='galactic')
        if True:
                egr = bayestar(coords, mode='best')
                f = return_ext_curve(Rv=3.1)
                res.append([egr, egr*f(wavelength*10000.)])
        else:
                res.append([float('nan'),float('nan')])
    return res

#l,b = convert_l_b(RA=0., DEC=60.)
#print(main([[l,b, 60000000000.]], 0.655))