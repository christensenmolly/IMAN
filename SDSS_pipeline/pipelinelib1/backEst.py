#! /usr/bin/env python
# Example: python ~/programs/Sergey/pipeline-master/pipelinelib/backEst.py frame-g-003813-1-0201.fits --method annulus --sky_subtr yes --x_obj 179.63042 --y_obj 43.947405 --manual yes
import sys
import subprocess
import pyfits
from pylab import *
import itertools
import os
from os.path import exists
from os import remove
from scipy.spatial import cKDTree
from scipy.optimize import fmin_tnc, fmin
import arithm_operations
import argparse

from datetime import datetime

# Disable astropy logging except for warnings and errors
from astropy import log
log.setLevel("WARNING")

import warnings
warnings.filterwarnings("ignore")

startTime = datetime.now()

#---------------------------
# Parameters for the method to estimate flat sky background
N_box = 10000	# Number of random placed boxes
N_tr = 5	# Number of trials
Box_size = 5	# Side length of the box
#---------------------------
Cov_coeff = 1.2	# Coefficient to enlarge the size of the masked obects detected bt Sextractor, i.e. R_new = R * Cov_coeff  

Annulus_coeff = 2	# Coefficient such that the annulus will have Sma = Annulus_coeff * Sma_sextr
R_find = 10.	# Search the object inside this radius

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


class ObjParams:
    def __init__(self, xCen, yCen, ellA, ellB, PA):
        self.cen = Point(xCen, yCen)
        self.ellA = ellA
        self.ellB = ellB
        self.PA = PA

# Colors to highlight the output text
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''
        
def numpy_weighted_median(data, weights=None):
    """Calculate the weighted median of an array/list using numpy."""
    # Taken from https://pypi.python.org/pypi/weightedstats/0.2
    import numpy as np
    if weights is None:
        return np.median(np.array(data).flatten())
    data, weights = np.array(data).flatten(), np.array(weights).flatten()
    if any(weights > 0):
        sorted_data, sorted_weights = map(np.array, zip(*sorted(zip(data, weights))))
        midpoint = 0.5 * sum(sorted_weights)
        if any(weights > midpoint):
            return (data[weights == np.max(weights)])[0]
        cumulative_weight = np.cumsum(sorted_weights)
        below_midpoint_index = np.where(cumulative_weight <= midpoint)[0][-1]
        if cumulative_weight[below_midpoint_index] == midpoint:
            return np.mean(sorted_data[below_midpoint_index:below_midpoint_index+2])
        return sorted_data[below_midpoint_index+1]








def polyfit2d(x, y, z, order, linear):

    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
        if linear & (i != 0.) & (j != 0.):
            G[:, k] = 0
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m
    
# Values to two-dimensional polynomial fit. 
# Based uppon code provided by Joe Kington
def polyval2d(x, y, m):

    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z




def ellipse_annulus(output_image,inner_ellipse, outer_ellipse, inframe, mask, xSize, ySize, header, sky_subtr):
	from astropy.stats import sigma_clipped_stats
	def define_ellipse(cen, ellA, ellB, ellPA):
	    cospa = cos(radians(ellPA))
	    sinpa = sin(radians(ellPA))
	    # Check if whole ellipse is inside of the image
	    # and obtain size of the ellipse in xy plane
	    xMax = 0.0
	    yMax = 0.0
	    xMin = 1e10
	    yMin = 1e10
	    for e in linspace(0, 4*pi, 1000):
		cose = cos(e)
		sine = sin(e)
		x = cen.x + ellA * cose * cospa - ellB * sine * sinpa
		y = cen.y + ellB * sine * cospa + ellA * cose * sinpa
		if x > xMax:
		    xMax = x
		if y > yMax:
		    yMax = y
		if x < xMin:
		    xMin = x
		if y < yMin:
		    yMin = y
	    xMin = max(0, int(round(xMin)))
	    xMax = min(xSize, int(round(xMax)))
	    yMin = max(0, int(round(yMin)))
	    yMax = min(ySize, int(round(yMax)))
	    
	    
	    focusR = (ellA ** 2.0 - ellB ** 2.0) ** 0.5
	    focus10 = Point(cen.x + focusR, cen.y)  # Unrotated
	    focus20 = Point(cen.x - focusR, cen.y)  #
	    focus1 = rot_point(focus10, cen, radians(ellPA))
	    focus2 = rot_point(focus20, cen, radians(ellPA))
	    return xMin,xMax,yMin,yMax,focus1,focus2

	cen_in, ellA_in, ellB_in, ellPA_in = inner_ellipse
	cen_out, ellA_out, ellB_out, ellPA_out = outer_ellipse
	
	## Inner ellipse
	xMin_in,xMax_in,yMin_in,yMax_in,focus1_in,focus2_in = define_ellipse(cen_in, ellA_in, ellB_in, ellPA_in)

	## Outer ellipse
	xMin_out,xMax_out,yMin_out,yMax_out,focus1_out,focus2_out = define_ellipse(cen_out, ellA_out, ellB_out, ellPA_out)
	
	xMin = min([xMin_in,xMin_out])
	yMin = min([yMin_in,yMin_out])

	xMax = max([xMax_in,xMax_out])
	yMax = max([yMax_in,yMax_out])

        # Find pixels inside of the ellipse annulus
        dEll_in = 2 * ellA_in
        dEll_out = 2 * ellA_out
        
        I_an = []
        I_new_mask = (mask>0.)
	I_dan = np.copy(inframe)
	X = []
	Y = []
	Z = []
        for x in xrange(xMin, xMax):
            for y in xrange(yMin, yMax):
                dFocus1_in = hypot(x-focus1_in.x, y-focus1_in.y)
                dFocus2_in = hypot(x-focus2_in.x, y-focus2_in.y)
                dPoint_in = dFocus1_in + dFocus2_in

                dFocus1_out = hypot(x-focus1_out.x, y-focus1_out.y)
                dFocus2_out = hypot(x-focus2_out.x, y-focus2_out.y)
                dPoint_out = dFocus1_out + dFocus2_out                
                if dPoint_out < dEll_out and dPoint_in > dEll_in and mask[y,x]==0.:
		  try:
			I_an.append(inframe[y,x])
			I_new_mask[y,x] = False
			I_dan[y,x] = 1000.
			X = np.append(X, float(x))
			Y = np.append(Y, float(y))
			Z = np.append(Z, float(inframe[y,x]))
		  except:
		    zz=1
	print 'Number of pixels in the annulus is', len(I_an)
	outHDU = pyfits.PrimaryHDU(I_dan)
	outHDU.writeto('annulus.fits',clobber=True)
	mean, median, std = sigma_clipped_stats(I_an, sigma=3.0, iters=5)
	#print mean, median, std
	sky = median

	'''
	#### TO DO. DON'T USE THIS OPTION BELOW
	degree = degree + 1
	if degree>0:
	  # Fit a 3rd degree, 2d polynomial
	  m = polyfit2d(X,Y,Z,degree,True)
	  
	  # Evaluate it on a grid...
	  nx, ny = 20, 20
	  xx, yy = np.meshgrid(np.linspace(0, xSize, nx),
			      np.linspace(0, ySize, ny))
	  zz = polyval2d(xx, yy, m)

	  # Make a sky array and subtract it from the object frame
	  sky = np.copy(inframe)
	  strip = np.arange(1.0,np.float(xSize+1.0))
	  for i in range(0,ySize):
	      sky[i,:] = polyval2d(strip,i,m)
	'''
	if sky_subtr==True:
	  sframe = inframe - sky
	else:
	  sframe = inframe

	header['Sky_med'] = sky
	header['Sky_std'] = std

	header['Ann_x0'] = cen_in.x
	header['Ann_y0'] = cen_in.y
	header['Ann_sma'] = ellA_in
	header['Ann_smb'] = ellB_in
	header['Ann_PA'] = ellPA_in
	outHDU = pyfits.PrimaryHDU(sframe,header=header)
	outHDU.writeto(output_image,clobber=True)  
	return sky,std



def sky_in_ell_annulus(input_image,output_image,mask_image,ell_pars,N,sky_subtr,manual):  
    hdulist_ima = pyfits.open(input_image)
    data_image = hdulist_ima[0].data
    header = hdulist_ima[0].header

    hdulist_mask = pyfits.open(mask_image)
    data_mask = hdulist_mask[0].data

    ySize,xSize = np.shape(data_image)
    
    '''
    mask_astropy = np.zeros_like(mask_infr, dtype=bool)
    nx, ny =mask_infr.shape[1], mask_infr.shape[0]
    for k in range(ny):
      for i in range(nx):
	if mask_infr[k,i]>0:
	  mask_astropy[k,i] = True
	else:
	  mask_astropy[k,i] = False
    '''

    if np.isnan(ell_pars[0])==False:
	# Inner ellipse (given)
	xc_in = ell_pars[0]
	yc_in = ell_pars[1]
	sma_in = ell_pars[2]
	smb_in = ell_pars[3] #*1.5 for NGC4013
	PA_in = ell_pars[4]
	inner_ellipse = [Point(xc_in,yc_in),sma_in,smb_in,PA_in]
    else:
	xc_in = xSize/2.
	yc_in = ySize/2.
	sma_in = 50
	smb_in = 50
	PA_in = 45
	inner_ellipse = [Point(xc_in,yc_in),sma_in,smb_in,PA_in]      



      
    # Find width of the annulus: the annulus should include N pixels
    ann_width = ceil(0.5*(-(sma_in+smb_in)+sqrt((sma_in+smb_in)**2+4.*N*2/math.pi)))
    print 'Width of the annulus is', ann_width
    
    # Outer ellipse (calculated)
    xc_out = xc_in
    yc_out = yc_in
    sma_out = sma_in + ann_width
    smb_out = smb_in + ann_width #*1.5 for NGC4013
    PA_out = PA_in
    outer_ellipse = [Point(xc_out,yc_out),sma_out,smb_out,PA_out]
    
    if manual==True:
	    crea_annulus(inner_ellipse,outer_ellipse)
	    inner_ellipse,outer_ellipse = read_annulus()
    
    sky,std = ellipse_annulus(output_image,inner_ellipse, outer_ellipse, data_image, data_mask, xSize, ySize, header, sky_subtr)
    return sky,std
    

def flat_sky_subtr(input_image):
    from astropy.stats import sigma_clipped_stats
    
    hdulist = pyfits.open('segm.fits') # open a FITS file
    mask_infr = hdulist[0].data
    '''
    ny,nx = np.shape(mask_infr)
    mask_astropy = np.zeros_like(mask_infr, dtype=bool)

    for k in range(ny):
      for i in range(nx):
	if mask_infr[k,i]>0.:
	  mask_astropy[k,i] = True
    '''
    mask_astropy = (mask_infr>0.)

    # Taken from Mratin-Navarro et al. (2012) and Fliri & Trujillo (2016)
    # This method is valid for frames with flat sky background (no gradient)
    
    ## Subtract two images produced by SExtractor => background
    #arithm_operations.main(input_image,'objects.fits','sub','backgr.fits')
    
    ## Read in the file with the background:
    #hdulist = pyfits.open('backgr.fits') # open a FITS file
    hdulist = pyfits.open(input_image) # open a FITS file
    scidata = hdulist[0].data	
    ny,nx = np.shape(scidata)
    
    ## Measure median values in N_box randomly placed boxes
    I_MED = []
    I_STD = []
    XC = []
    YC = []
    ii = 0
    while ii<N_tr:
	k = 0
	I_med = []
	I_std = []
	while k<N_box:
		I = []
		xc = randint(Box_size,nx-Box_size-1)
		yc = randint(Box_size,ny-Box_size-1)

		x = arange(xc-Box_size,xc+Box_size,1)
		y = arange(yc-Box_size,yc+Box_size,1)
		
		breaker = False
		for m in y:
		    for n in x:
			  if mask_infr[m,n]==0.:
				I.append(scidata[m,n])
			  else:
				breaker = True 
				break
		    if breaker:
		      break
		if breaker==False:
		  I_med.append(np.median(I))
		  I_std.append(np.std(I))
		  XC.append(xc)
		  YC.append(yc)
		  del I
		  k = k + 1
	I_mean, I_median, std = sigma_clipped_stats(I_med, sigma=3.0, iters=5)
	I_me, I_medi, std = sigma_clipped_stats(I_std, sigma=3.0, iters=5)
	I_MED.append(I_median)
	I_STD.append(I_medi)
	ii = ii + 1
    #import matplotlib.pylab as plt
    #plt.plot(XC,YC,'o')
    #plt.show()
    ## Remove the file with the background
    #os.remove("backgr.fits")
    
    #mean, median, std = sigma_clipped_stats(I_med, sigma=3.0, iters=5)
    I_median = np.mean(I_MED)
    I_std = np.mean(I_STD)
    
    '''
    mean, median, std = sigma_clipped_stats(scidata, sigma=3.0, iters=10, mask=mask_astropy)
    print median, std    
    
    from photutils.background import Background
    bkg = Background(scidata, (Box_size, Box_size), 1, method='sextractor', mask=mask_astropy)
    print bkg.background_median,bkg.background_rms_median

    import matplotlib.pylab as plt
    from astropy.visualization import SqrtStretch
    from astropy.visualization.mpl_normalize import ImageNormalize
    norm = ImageNormalize(stretch=SqrtStretch())
    bkg3 = bkg.background * ~mask_astropy
    plt.imshow(bkg3, origin='lower', cmap='Greys_r', norm=norm)    
    plt.show()
    '''
    return I_median,I_std   

    
    


def rot_point(p, orig, angle):
    """ performs rotation of a point (x,y) around point (x0,y0)
        by angle radians"""
    x1 = cos(angle)*(p.x-orig.x) - sin(angle)*(p.y-orig.y) + orig.x
    y1 = sin(angle)*(p.x-orig.x) + cos(angle)*(p.y-orig.y) + orig.y
    return Point(x1, y1)


def call_SE(fitsFile):
    # detect the name of SExtractor
    if subprocess.call("which sex >/dev/null", shell=True) == 0:
        callString = "sex %s " % fitsFile
        CallString = "sex %s " % fitsFile
    elif subprocess.call("which sextractor >/dev/null", shell=True) == 0:
        callString = "sextractor %s " % fitsFile
        CallString = "sextractor %s " % fitsFile
    else:
        print "SExtractor was not found. Exiting..."
        exit(1)
    
    # Cold catalogue
    callString += "-c %s/cold.sex" % (os.path.dirname(__file__))
    callString += " -VERBOSE_TYPE=QUIET"
    subprocess.call(callString, shell=True)
    os.rename('field.cat','field_cold.cat')
    os.rename('segm.fits','segm_cold.fits')

    # Hot catalogue
    CallString += "-c %s/hot.sex" % (os.path.dirname(__file__))
    CallString += " -VERBOSE_TYPE=QUIET"
    subprocess.call(CallString, shell=True)
    os.rename('field.cat','field_hot.cat')
    os.rename('segm.fits','segm_hot.fits')
    
    arithm_operations.main('segm_cold.fits','segm_hot.fits','add','segm.fits')
    filenames = ['field_hot.cat', 'field_cold.cat']
    with open('field.cat', 'w') as outfile:
	for fname in filenames:
	    with open(fname) as infile:
		for line in infile:
		    outfile.write(line)    
    remove('segm_hot.fits')
    remove('segm_cold.fits')
    remove('field_hot.cat')
    remove('field_cold.cat')   

def get_SE_results():
    objects = []
    for line in open("field.cat"):
        if line.startswith("#"):
            continue
        params = line.split()
        n = int(params[0])
        xCen = float(params[1])
        yCen = float(params[2])
        kron = float(params[8])
        ellA = kron * float(params[4]) * Cov_coeff
        ellB = kron * float(params[5]) * Cov_coeff
        PA = float(params[6])
        objects.append(ObjParams(xCen, yCen, ellA, ellB, PA))
    return objects

def find_object(x_obj,y_obj):
    for line in open("field.cat"):
        if line.startswith("#"):
            continue
        params = line.split()
        n = int(params[0])
        xCen = float(params[1])
        yCen = float(params[2])
        if sqrt((xCen-x_obj)**2+(yCen-y_obj)**2)<R_find and float(params[10])<0.1:
	  kron = float(params[8])
	  ellA = kron * float(params[4])
	  ellB = kron * float(params[5])
	  PA = float(params[6])
	  break
    try:
      return xCen,yCen,ellA,ellB,PA
    except:
      return float('nan'),float('nan'),float('nan'),float('nan'),float('nan')
  
def create_mask(objects, xSize, ySize):
    mask = zeros((ySize, xSize))
    for obj in objects:
        cospa = cos(radians(obj.PA))
        sinpa = sin(radians(obj.PA))
        # Check if whole ellipse is inside of the image
        # and obtain size of the ellipse in xy plane
        xMax = 0.0
        yMax = 0.0
        xMin = 1e10
        yMin = 1e10
        for e in linspace(0, 4*pi, 1000):
            cose = cos(e)
            sine = sin(e)
            x = obj.cen.x + obj.ellA * cose * cospa - obj.ellB * sine * sinpa
            y = obj.cen.y + obj.ellB * sine * cospa + obj.ellA * cose * sinpa
            if x > xMax:
                xMax = x
            if y > yMax:
                yMax = y
            if x < xMin:
                xMin = x
            if y < yMin:
                yMin = y
        xMin = max(0, int(round(xMin)))
        xMax = min(xSize, int(round(xMax)))
        yMin = max(0, int(round(yMin)))
        yMax = min(ySize, int(round(yMax)))
        focusR = (obj.ellA ** 2.0 - obj.ellB ** 2.0) ** 0.5
        focus10 = Point(obj.cen.x + focusR, obj.cen.y)  # Unrotated
        focus20 = Point(obj.cen.x - focusR, obj.cen.y)  #
        focus1 = rot_point(focus10, obj.cen, radians(obj.PA))
        focus2 = rot_point(focus20, obj.cen, radians(obj.PA))
        # Find pixels inside of the ellipse
        dEll = 2 * obj.ellA
        for x in xrange(xMin, xMax):
            for y in xrange(yMin, yMax):
                dFocus1 = hypot(x-focus1.x, y-focus1.y)
                dFocus2 = hypot(x-focus2.x, y-focus2.y)
                dPoint = dFocus1 + dFocus2
                if dPoint < dEll:
                    mask[y, x] = 1
    
    outHDU = pyfits.PrimaryHDU(mask)
    outHDU.writeto('mask.fits',clobber=True)
    return mask


def buildDiagram(mask):
    maskedPixels = where(mask > 0)
    maskKDTree = cKDTree(data=zip(maskedPixels[0], maskedPixels[1]),
                         leafsize=100)
    diagram = zeros_like(mask)
    listOfAllCoords = list(itertools.product(range(mask.shape[0]), range(mask.shape[1])))
    distances, nearestPoints = maskKDTree.query(listOfAllCoords)
    diagram = reshape(distances, (mask.shape[0], mask.shape[1]))
    diagram = diagram / sum(diagram)
    diagramHDU = pyfits.PrimaryHDU(diagram)
    if exists("diagram.fits"):
        remove("diagram.fits")
    diagramHDU.writeto("diagram.fits")
    return diagram


def create_2d_poly_of_n_degree(n):
    """Function returns function that represents 2d polynomial
    of n degree"""
    def pk(coeffs, xy):
        r = coeffs[0]
        for k in xrange(1, n+1):
            r += coeffs[k] * xy[0]**k + coeffs[k+n] * xy[1]**k
        return r
    return pk


def fitness_function(data, poly, coeffs, diagram):
    xSize, ySize = data.shape
    yCen = int(ySize/2.0)
    xCen = int(xSize/2.0)
    coords = meshgrid(xrange(ySize), xrange(xSize))
    coords[0] = coords[0] - yCen
    coords[1] = coords[1] - xCen
    difference = sum(diagram * (data - poly(coeffs, coords))**2) / (xSize*ySize)
    return difference


def fit_sky_by_nd_poly(data, diagram, degree):
    #meanSky = average(data, weights=diagram)
    meanSky = numpy_weighted_median(data, weights=diagram)
    if True:
      from astropy.stats import sigma_clipped_stats
      hdulist = pyfits.open('segm.fits') # open a FITS file
      mask_infr = hdulist[0].data
      mask_astropy = (mask_infr>0.)

      mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5, mask=mask_astropy)
    else:
      std = math.sqrt(np.average((data-meanSky)**2, weights=diagram))

    print "Free coeff:", meanSky
    # if the degree of the polynomial is equal to zero
    # than return just weighted average of the sky level
    if degree == 0:
        return meanSky
    coeffs = [meanSky]
    # if the degree is not equal to zero.
    # this cycle fitts the data with polynomials of
    # increasing degree (from 1 to degree=degree). At each
    # inerration the results of previous iteration (k-1 degree poly)
    # are being used as initial guess. Coefficients for higher degrees
    # of x and y are set to be zero as initial guesses
    for k in xrange(1, degree+1):
        print "Fitting by %i degree polynomial" % k
        kdPoly = create_2d_poly_of_n_degree(k)
        fitFunc = lambda c: fitness_function(data, kdPoly, c, diagram)
        # unpack results of previous iteration
        freeCoeff = coeffs[0]
        xCoeffs = list(coeffs[1:k])
        yCoeffs = list(coeffs[k:])
        coeffs = None
        # new highest degree coefficients are zeros
        xCoeffs.append(0.0)
        yCoeffs.append(0.0)
        # assemble all coefficients are new initial guesses
        coeffs = [freeCoeff]
        coeffs.extend(xCoeffs)
        coeffs.extend(yCoeffs)
        # run gradient descent
        print "coeffs before fmin", coeffs, xCoeffs, yCoeffs
        res = fmin_tnc(fitFunc, x0=coeffs, approx_grad=True, maxfun=1000)
        coeffs = res[0]
        print "res:", coeffs
    return coeffs,std


def get_mean_sky(data, diagram):
    #meanSky = average(data, weights=diagram)
    meanSky = numpy_weighted_median(data, weights=diagram)
    
    if True:
      from astropy.stats import sigma_clipped_stats
      hdulist = pyfits.open('segm.fits') # open a FITS file
      mask_infr = hdulist[0].data
      mask_astropy = (mask_infr>0.)

      mean, median, std = sigma_clipped_stats(data, sigma=3.0, iters=5, mask=mask_astropy)
    else:
      std = math.sqrt(np.average((data-meanSky)**2, weights=diagram))
    return meanSky,std


def subtract_back(inFitsName, outFitsName, degree, coeffs, std, sky_subtr):
    inHDU = pyfits.open(inFitsName)
    data = inHDU[0].data
    header = inHDU[0].header
    xSize, ySize = data.shape
    yCen = int(ySize/2.0)
    xCen = int(xSize/2.0)
    if degree == 0:       
        Sky = coeffs
        if sky_subtr=='True':
	  outData = data - coeffs
	else:
	  outData = data
    else:
	Sky = coeffs[0]
	if sky_subtr=='True':
	  poly = create_2d_poly_of_n_degree(degree)
	  coords = meshgrid(xrange(ySize), xrange(xSize))
	  coords[0] = coords[0] - yCen
	  coords[1] = coords[1] - xCen
	  sky_data = poly(coeffs, coords)       
	  outData = data - sky_data
	else:
	  outData = data
        
    header['Sky_med'] = Sky
    header['Sky_std'] = std
    outHDU = pyfits.PrimaryHDU(outData, header=header)
    outHDU.writeto(outFitsName,clobber=True)

def crea_annulus(inner_ellipse,outer_ellipse):
    [cen_in,sma_in,smb_in,PA_in] = inner_ellipse
    [cen_out,sma_out,smb_out,PA_out] = outer_ellipse

    fout = open("annulus.reg", "w")
    fout.truncate(0)
    fout.write("# Region file format: DS9 version 4.1\n")
    fout.write('global color=green dashlist=8 3 width=1 font="helvetica 10 ')
    fout.write('normal roman" select=1 highlite=1 dash=0 fixed=0 ')
    fout.write('edit=1 move=1 delete=1 include=1 source=1\n')
    fout.write('image\n')
    #fout.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # inner ellipse\n" % (cen_in.x, cen_in.y,
    #                                                             sma_in, smb_in,
    #                                                             PA_in))
    #fout.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # outer ellipse\n" % (cen_out.x, cen_out.y,
    #                                                             sma_out, smb_out,
    #                                                             PA_out))
    fout.write("ellipse(%1.1f,%1.1f,%1.1f,%1.1f,%1.1f,%1.1f,%1.1f) # ellipse annulus\n" % (cen_in.x, cen_in.y,
                                                                 sma_in, smb_in, sma_out, smb_out,
                                                                 PA_in))    
    
    
    
    fout.close()      
    ds9Proc = subprocess.Popen(["ds9", input_image,
                                "-regions", "annulus.reg",
                                "-scale", "log"])
    ds9Proc.wait()      

def read_annulus():
   ellipses = []
   for line in open("annulus.reg"):
        if "ellipse" in line:
	  params = line.split(",")
	  cen = Point(float(params[0].split('(')[1]),
		      float(params[1]))
	  sma_in = float(params[2])
	  smb_in = float(params[3])
	  sma_out = float(params[4])
	  smb_out = float(params[5])
	  PA = float(params[6].split(')')[0])
   inner_ellipse = [cen,sma_in,smb_in,PA]
   outer_ellipse = [cen,sma_out,smb_out,PA]
   return inner_ellipse,outer_ellipse
 
def estimate_sky(fitsFileName, cleanName = "back_clean.fits", method = 'kd-tree', degree=2, sky_subtr=True, x_obj = None, y_obj = None, coords = 'world', manual=False):
    # Function to find sky background
    # cleanName - the output sky-subtracted image (if sky_subtr=True, by default)
    #             or the same image as the input but with added values of the sky level (median and std dev) to the header (if sky_subtr=False)
    # degree - the polinomial degree of the fitted background
    # method - which method to use:
    #          'random' - randomly placed boxes with the parameters on the top of this script
    #          'kd-tree' - use kd-tree diagram
    #          'annulus' - use the annulus around the given object. In this case you should specify x_obj and y_obj
    # x_obj  - RA or x-coordinate in pixels of the object around which one wants to find sky background 
    # y_obj  - DEC or y-coordinate in pixels of the object around which one wants to find sky background 
    # coords - Specify units of x_obj and y_obj: 'world' or 'pixels'

    ####
    #xc = 1207.9678
    #yc = 1163.7023
    #sma = 108.357
    #smb = 23.793
    #theta = 68.63
    #kron_r = 3.50
    #ann_width = 30.
    #ell_pars = [xc,yc,sma*kron_r* Annulus_coeff,smb*kron_r* Annulus_coeff,theta]
    #sky_in_ell_annulus(sys.argv[1],'segm.fits',ell_pars,ann_width)
    #exit()
    ####

    print '\t\t\t',bcolors.OKBLUE+ "*****************************************" + bcolors.ENDC
    print '\t\t\t',bcolors.OKBLUE+ "      SKY BACKGROUND ESTIMATION 2016     " + bcolors.ENDC
    print '\t\t\t',bcolors.OKBLUE+ "*****************************************" + bcolors.ENDC
    
    print bcolors.OKGREEN+ 'The chosen method is ' + bcolors.ENDC, method 
    print bcolors.OKGREEN+'The polinomial degree is ' + bcolors.ENDC, degree 
    print bcolors.OKGREEN+'Sky background will be subtracted in the output image ' + bcolors.ENDC, str(sky_subtr)     
    
    ## Open the input file
    HDU = pyfits.open(fitsFileName)
    data = HDU[0].data
    header = HDU[0].header
    ySize, xSize = data.shape
    
    ## Run SExtractor 
    call_SE(fitsFileName)
    
    if method=='kd-tree':
	## Define objects from SExtractor output catalog
	objects = get_SE_results()
	
	## Create mask for detected objects
	mask = create_mask(objects, xSize, ySize)
	
	
	diagram = buildDiagram(mask)
	if degree == 0:
	    sky,std = get_mean_sky(data, diagram)
	else:
	    sky,std = fit_sky_by_nd_poly(data, diagram, degree)

	subtract_back(fitsFileName, cleanName, degree, sky, std, sky_subtr)
    elif method=='random':
	sky,std = flat_sky_subtr(fitsFileName)
	subtract_back(fitsFileName, cleanName, 0, sky, std, sky_subtr)
    elif method=='annulus':
	if coords=='world':
	  from astropy.wcs import WCS
	  w = WCS(fitsFileName)
	  [[x_obj,y_obj]] = w.wcs_world2pix([[x_obj,y_obj]], 1)
	xc,yc,sma,smb,theta = find_object(x_obj,y_obj)
	if np.isnan(xc)==True:
	  print bcolors.FAIL+'The object (%f,%f) is not found!' % (x_obj,y_obj) + bcolors.ENDC
	  if manual==False:
	    exit()
	else:
	  print bcolors.OKGREEN+'The object (%f,%f) is found!' % (xc,yc) + bcolors.ENDC    
	ell_pars = [xc, yc, sma* Annulus_coeff, smb* Annulus_coeff, theta]
	Number_of_an_pixels = N_box
	sky,std = sky_in_ell_annulus(fitsFileName,cleanName,'segm.fits',ell_pars,Number_of_an_pixels, sky_subtr, manual)
    else:
	print bcolors.FAIL+'The input method is not found!' + bcolors.ENDC
	exit()
    if isinstance(sky, (list, tuple, np.ndarray)):
      sky = sky[0]
    print bcolors.OKGREEN+ 'Median sky level is' + bcolors.ENDC, sky
    print bcolors.OKGREEN+ 'Background std deviation is' + bcolors.ENDC, std
    print 'Done!'
    print 'Ellapsed time', datetime.now() - startTime
    return sky


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sky background estimation")
    parser.add_argument("inputImage", help="Input fits image with the centered object")
    parser.add_argument("--outputImage", nargs='?', const=1, help="Optional: Input the name of the output file",type=str,default='back_clean.fits') 
    parser.add_argument("--method", nargs='?', const=1, help="Optional: Specify the method: random, kd-tree, annulus",type=str,default='random')
    parser.add_argument("--degree", nargs='?', const=1, help="Optional: Input the polinomial order to fit the sky",type=int,default=2)
    parser.add_argument("--sky_subtr", nargs='?', const=1, help="Optional: Do you want to subtract sky from the original image (yes,no)",type=str,default='yes')
    parser.add_argument("--x_obj", nargs='?', const=1, help="Optional: Input the x-coordinate (or RA) of the object where to estimate the sky",type=str,default='none') 
    parser.add_argument("--y_obj", nargs='?', const=1, help="Optional: Input the y-coordinate (or DEC) of the object where to estimate the sky",type=str,default='none') 
    parser.add_argument("--coords", nargs='?', const=1, help="Optional: Input the system of coordinates: world or pixels",type=str,default='world')
    parser.add_argument("--manual", nargs='?', const=1, help="Optional: Only in case of the annulus method (yes,no)",type=str,default='no')
    args = parser.parse_args()

    input_image = args.inputImage
    output_image = args.outputImage
    method = args.method
    degree = args.degree
    if method!='kd-tree':
      degree = 0
    
    if args.sky_subtr=='yes':
      sky_subtr = True
    else:
      sky_subtr = False
    if args.x_obj=='none':
      x_obj = None
    else:
      x_obj = float(args.x_obj)
    if args.y_obj=='none':
      y_obj = None
    else:
      y_obj = float(args.y_obj)
    if args.manual=='yes':
      manual = True
    else:
      manual = False

    coords = args.coords
    sky = estimate_sky(input_image, output_image, method, degree, sky_subtr, x_obj, y_obj, coords, manual)
    

