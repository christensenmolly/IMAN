
#************* INI CROP  *******************
Ini_crop_factor = 6.
#*******************************************

#************* SKY SUBTRACTION *************
method_to_use = 'kd-tree'	# 'annulus', 'random', 'kd-tree'
Cov_coeff = 1.2	# Coefficient to enlarge the size of the masked obects detected by Sextractor, i.e. R_new = R * Cov_coeff  
sky_subtr = True
degree = 2
# Parameters for the method to estimate flat sky background
N_box = 10000	# Number of randomly placed boxes
N_tr = 5	# Number of trials
Box_size = 5	# Side length of the box (pix)

# Elliptical annulus method
Annulus_coeff = 1.0	# Coefficient such that the annulus will have Sma = Annulus_coeff * Sma_sextr
Number_of_an_pixels = 100000	# Number of annulus pixels for estimating sky background
R_find = 10.	# Search the object inside this radius
#*******************************************


#*************** ROTATION ******************
manual_method = 'line'
inner_level = 10.	# in std of background
outer_level = 3.	# in std of background
#*******************************************


#*************** CENTERING *****************
exend_factor = 2.	#2.5 - for 25.5, 5 - optimal
#*******************************************

#*************** MASKING *******************
mask_comp = 'no'
norm = False
exend_factor_mask = 2.5
#*******************************************


#************ PSF SExtractor setup**********
box_psf = 51
beta = 3.4
window = 'moffat'

N_psf_stars = 10	# Maximal number of best PSF stars
star_cl = 0.98
Npix_around_star = 3	# Minimal width of annulus (in pixels) around PSF stars to estimate local background 

# READ MANUAL TO SEXTRACTOR
DETECT_MINAREA = 4
DETECT_THRESH = 3.0
FILTER_NAME = 'default.conv'
DEBLEND_NTHRESH = 5
DEBLEND_MINCONT = 0.00005
CLEAN_PARAM = 2.0
PHOT_AUTOPARAMS = (2.0, 4.0)
BACK_SIZE = 300
BACK_FILTERSIZE = 5
BACKPHOTO_TYPE = 'LOCAL'
MEMORY_PIXSTACK = 3000000
#*******************************************