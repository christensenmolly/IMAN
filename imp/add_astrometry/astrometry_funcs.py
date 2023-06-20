


def calc_from_cd(cd1_1, cd1_2, cd2_1, cd2_2):

        # TODO: Check if first coordinate in CTYPE is latitude
        # if (ctype EQ 'DEC-') or (strmid(ctype, 1) EQ 'LAT')  then $
        #    cd = reverse(cd,1)

        det = cd1_1*cd2_2 - cd1_2*cd2_1
        if det < 0:
            sgn = -1
        else:
            sgn = 1
        ## if det > 0:
        ##     raise ValueError("Astrometry is for a right-handed coordinate system")

        if (cd2_1 == 0.0) or (cd1_2 == 0.0):
            # Unrotated coordinates?
            xrot = 0.0
            yrot = 0.0
            cdelt1 = cd1_1
            cdelt2 = cd2_2
        else:
            xrot = math.atan2(sgn * cd1_2, sgn * cd1_1)
            yrot = math.atan2(-cd2_1, cd2_2)

            cdelt1 = sgn * math.sqrt(cd1_1**2 + cd1_2**2)
            cdelt2 = math.sqrt(cd1_1**2 + cd2_1**2)

        return xrot, yrot, cdelt1, cdelt2

def convert_pix_to_wcs(I,J,I0,J0,X0,Y0,cdelt1,cdelt2,crota2=0.0):
  if crota2!=0.0:
    crota2 = crota2 * math.pi / 180	# in radians

  sinrot = math.sin(crota2)
  cosrot = math.cos(crota2)
  pcmatrix = [cosrot, -sinrot, sinrot, cosrot]
  cd1_1 = cdelt1*pcmatrix[0]
  cd1_2 = cdelt2*pcmatrix[1]
  cd2_1 = cdelt1*pcmatrix[2]
  cd2_2 = cdelt2*pcmatrix[3]
  
  X = X0 + cd1_1* (I - I0) + cd1_2* (J - J0)
  Y = Y0 + cd2_1* (I - I0) + cd2_2* (J - J0)
  return X,Y

def define_wcs(input_image, pix_coords,cel_coords):
  from astropy.io import fits as pyfits
  shutil.copy(input_image, 'galaxy.fits')
  hdulist = pyfits.open('galaxy.fits')
  data = hdulist[0].data
  
  #http://www.atnf.csiro.au/people/mcalabre/WCS/Intro/WCS06.html
  I0 = pix_coords[0][0]
  J0 = pix_coords[0][1]

  I1 = pix_coords[1][0]
  J1 = pix_coords[1][1]

  I2 = pix_coords[2][0]
  J2 = pix_coords[2][1]

  X0 = cel_coords[0][0]
  Y0 = cel_coords[0][1]
  
  X1 = cel_coords[1][0]
  Y1 = cel_coords[1][1]

  X2 = cel_coords[2][0]
  Y2 = cel_coords[2][1]
  
  #(I1-I0) * C11 + (J1-J0) * C12 = X1-X0
  #(I2-I0) * C11 + (J2-J0) * C12 = X2-X0

  a = np.array([[(I1-I0),(J1-J0)], [(I2-I0),(J2-J0)]])
  b = np.array([X1-X0,X2-X0])
  
  [C11,C12] = np.linalg.solve(a, b)

  #(I1-I0) * C21 + (J1-J0) * C22 = Y1-Y0
  #(I2-I0) * C21 + (J2-J0) * C22 = Y2-Y0
  
  #a = np.array([[(I1-I0),(J1-J0)], [(I2-I0),(J2-J0)]])
  b = np.array([Y1-Y0,Y2-Y0])
  
  [C21,C22] = np.linalg.solve(a, b)


  xrot, yrot, cdelt1, cdelt2 = calc_from_cd(C11, C12, C21, C22)
  
  from astropy.io import fits as pyfits
  from astropy import wcs

  w = wcs.WCS(naxis=2)

  # what is the center pixel of the XY grid.
  w.wcs.crpix = [I0,J0]

  # what is the galactic coordinate of that pixel.
  w.wcs.crval = [X0, Y0]

  # what is the pixel scale in lon, lat.
  w.wcs.cdelt = np.array([cdelt1, cdelt2])

  # you would have to determine if this is in fact a tangential projection. 
  w.wcs.ctype = ["RA-TAN", "DEC-TAN"]

  pc1_1 = C11 / cdelt1
  pc1_2 = C12 / cdelt1
  pc2_1 = C21 / cdelt2
  pc2_2 = C22 / cdelt2



  # write the HDU object WITH THE HEADER
  header = w.to_header()

  header['PC1_1'] = pc1_1
  header['PC1_2'] = pc1_2
  header['PC2_1'] = pc2_1
  header['PC2_2'] = pc2_2

  hdu = pyfits.PrimaryHDU(data, header=header)
  hdu.writeto('galaxy.fits', clobber=True)


def add_wcs_rot(inputname,outputname,x_ref,y_ref,RA_ref,DEC_ref,cdelt1,cdelt2,angle):
    hdulist = pyfits.open(inputname)
    header = pyfits.getheader(inputname, 0)
    outframe = pyfits.getdata(inputname, 0)


    header['EQUINOX'] = 2.000000000000E+03
    header['RADECSYS'] = 'FK5'
    
    header['CTYPE1'] = 'RA---TAN'
    header['CUNIT1'] = 'deg'    
    header['CRVAL1'] = RA_ref
    header['CRPIX1'] = x_ref

    
    header['CTYPE2'] = 'DEC---TAN'
    header['CUNIT2'] = 'deg'  
    header['CRVAL2'] = DEC_ref
    header['CRPIX2'] = y_ref

    header['CDELT1'] = cdelt1
    header['CDELT2'] = cdelt2
    header['CROTA2'] = angle    


    hdu = pyfits.PrimaryHDU(data=outframe, header=header)
    hdu.writeto(outputname,clobber=True)   

def get_angle(p0, p1=np.array([0,0]), p2=None):
    # http://stackoverflow.com/questions/13226038/calculating-angle-between-two-lines-in-python
    ''' compute angle (in degrees) for p0p1p2 corner
    Inputs:
        p0,p1,p2 - points in the form of [x,y]
    '''
    if p2 is None:
        p2 = p1 + np.array([1, 0])
    v0 = np.array(p0) - np.array(p1)
    v1 = np.array(p2) - np.array(p1)

    angle = np.math.atan2(np.linalg.det([v0,v1]),np.dot(v0,v1))
    return np.degrees(angle)

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_wcs_angle(pts_wcs_ref, wcs_stars, pts_wcs_wrong):
  # pts_wcs_wrong - arbitrary converted wcs coordinates of the selected points
  # wcs_stars - selected stars from the catalogue which correspond to the given points

  Angles = []

  for k in range(len(pts_wcs_wrong)):
    for i in range(len(wcs_stars[k])):
      Angles.append(get_angle(np.array(wcs_stars[k][i]), p1=np.array(pts_wcs_ref), p2=np.array(pts_wcs_wrong[k])))

  
  
  n, bins,patches = plt.hist(Angles, bins=360., normed=True, fc='k', alpha=0.3)

  elem = np.argmax(n)
  PA =  bins[elem]
  PA_prec = find_nearest(np.array(Angles,float),bins[elem])

  return PA_prec
