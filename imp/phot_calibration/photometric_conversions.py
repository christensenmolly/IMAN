#!/usr/bin/python
# DESCRIPTION:
# Script to convert from one photometric system to the other



def transformation_from_gaia_to_sdss(G, GRP, GBP):
    # https://gea.esac.esa.int/archive/documentation/GDR2/Data_processing/chap_cu5pho/sec_cu5pho_calibr/ssec_cu5pho_PhotTransf.html
    x = GBP - GRP
    g = G - 0.13518 + 0.46245*x + 0.25171*x**2 - 0.021349*x**3
    r = G + 0.12879 - 0.24662*x + 0.027464*x**2 + 0.049465*x**3
    i = G + 0.29676 - 0.64728*x + 0.10141*x**2
    return g,r,i

def transformation_from_panstarrs_to_sdss(u, g, r, i, z, y):
    # http://iopscience.iop.org/article/10.3847/0004-637X/822/2/66/pdf
    x = g - i
    g = g + 0.01808 + 0.13595 * x - 0.01941 * x**2 + 0.00183 * x**3
    i = i - 0.01170 + 0.00400 * x - 0.00066 * x**2 + 0.00058 * x**3
    r = r + 0.01836 + 0.03577 * x - 0.02612 * x**2 + 0.00558 * x**3
    u = u - 0.04438 + 2.26095 * x + 0.13387 * x**2 - 0.27099 * x**3
    z = z + 0.01062 - 0.07529 * x + 0.03592 * x**2 - 0.00890 * x**3
    y = y - 0.08924 + 0.20878 * x - 0.10360 * x**2 + 0.02441 * x**3
    return u,g,r,i,z,y


def transformation_from_UBVRI_to_sdss(U, B, V, R, I):
    #https://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php
    # Jester et al. (2005)
    g = V + 0.60*(B-V) - 0.12
    r = V - 0.42*(B-V) + 0.11
    u = 1.28*(U-B) + 1.13 + g
    i = r - 0.91*(R-I) + 0.20
    z = r - 1.72*(R-I) + 0.41
    
    # Correct for shifts between SDSS and AB-system
    u = u - 0.04 # See http://www.sdss3.org/dr8/algorithms/fluxcal.php#SDSStoAB
    z = z + 0.02 # See http://www.sdss3.org/dr8/algorithms/fluxcal.php#SDSStoAB
    
    # Now all mags are in AB!
    return u,g,r,i,z


def transformation_from_sdss_to_UBVRI(u, g, r, i, z):
    #https://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php
    # Jester et al. (2005)
    B = g + 0.39*(g-r) + 0.21
    V = g - 0.59*(g-r) - 0.01
    R = V - 1.09*(r-i) - 0.22
    I = R - 1.00*(r-i) - 0.21
    U = B + 0.78*(u-g) - 0.88
    return U,B,V,R,I # Not AB!


def transformation_from_sdss_to_BVRI(u,g,r,i,z):
    #https://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php
    # Lupton (2005)
    B = u - 0.8116*(u - g) + 0.1313#;  sigma = 0.0095
    B = g + 0.3130*(g - r) + 0.2271#;  sigma = 0.0107

    V = g - 0.2906*(u - g) + 0.0885#;  sigma = 0.0129
    V = g - 0.5784*(g - r) - 0.0038#;  sigma = 0.0054

    R = r - 0.1837*(g - r) - 0.0971#;  sigma = 0.0106
    R = r - 0.2936*(r - i) - 0.1439#;  sigma = 0.0072

    I = r - 1.2444*(r - i) - 0.3820#;  sigma = 0.0078
    I = i - 0.3780*(i - z)  -0.3974#;  sigma = 0.0063
    return B,V,R,I



def get_mag(u, g, r, i, z, y, U, B, V, R, I, band):
  if band=='u':
    return u
  if band=='g':
    return g
  if band=='r':
    return r
  if band=='i':
    return i
  if band=='z':
    return z
  if band=='y':
    return y
  if band=='U':
    return U
  if band=='B':
    return B
  if band=='V':
    return V
  if band=='R':
    return R
  if band=='I':
    return I  

  
'''
G = 5.5821
GRP = 5.1393
GBP = 5.9337
x = transformation_from_gaia_to_sdss(G, GRP, GBP)
print(x)
'''
