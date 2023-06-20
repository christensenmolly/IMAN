# *************************** INTERNAL PARAMETERS FOR DECA ***************************

#====================================================================================
# DETAILS OF THE DECOMPOSITION
code = 'GALFIT'             # Decomposition code: 'GALFIT', 'GALFITM', or 'IMFIT'
minAlgorithm = 'LM'         # Mimimization algorithm to be used in case of IMFIT: 'LM' (Levenberg-Marquardt), 'NM' (Nelder-Mead), or 'GA' (Genetic algorithm)
                            # For GALFIT the only option is 'LM' so far
imfitPath = ''              # Path to imfit and makeimage binaries '/Users/mosenkov/Programs/'
galfitPath = ''             # Path to galfit binary
galfitmPath = ''            # Path to galfitm binary
addImfitStr = 'none'        # If you want to add some other text to imfit command
#------------------------------------------------------------------------------------
# GALFITM OPTIONS
param_method = '0'
fit_method = '0'

# CHEBYSHEV polynomials: orders
chb_coord_x = 1
chb_coord_y = 1

# BULGE (component 1)
chb_bulge_mag = 3
chb_bulge_reb = 2
chb_bulge_n = 2
chb_bulge_q = 1
chb_bulge_PA = 1

# DISC (component 2)
chb_disc_mag = 3
chb_disc_rs = 2
chb_disc_z0 = 2
chb_disc_q = 1
chb_disc_PA = 1
#------------------------------------------------------------------------------------
# GENETIC ALGORITHM DETAILS
numOfCores = 30		# Number of cores to scatter computations - the same as --n #TODO: REMOVE?
ZeroGenSize = 150	# Size of random zero generation
Popsize = 120		# Size of i-th generatio
selectNbest = 30	# Number of best organisms selected in the generation
addNew = 3		# Add this many random organisms in each generation
maxGenNumber = 100	# Maximum number of generation
fTol = 1e-5		# Relative fitness change for stop condition
fSpan = 5		# How many generations in a row have to have improvement less than fTol to raise the stop condition
saveGens = 'no'		# save image of the best model for each generation (yes/no)
runLM = 'yes'		# Run LM optimisation at the end. If no, just create a srcript to run it later (no/yes)
numOfLM = 5		# Number of LM model to run at the end
LMCores = 4		# Number of cores to use for each LM model
numOfFits = 5		# Number of repeated fits to estimate errors of parameters
#------------------------------------------------------------------------------------
# PSF 
conv_size = 51              # Convolution size (in pixels). If conv_size=None, then the size of the PSF kernel is used as the convolution box size. 
#====================================================================================


#====================================================================================
# SAVE RESULTS
res_directory = './results2'            # Specify path to the output directory where the results will be stored
path_to_save_out = 'deca_results.log'   # Specify a name for the log-file
del_files = False
#------------------------------------------------------------------------------------
# SHOW RESULTS
plot_hor = False            # Plot a horizontal photometric cut crossing the centre of the object: True or False
plot_azim_aver = True       # Plot an azimutahlly-averaged profile: True or False
plot_pix = False            # Plot a profile which represents all pixels in dependance on the distance between them and the object centre: True or False
plot_2D = True              # Plot reference image, model and residual image: True or False
plot_2d_allignment = 'line' # Specify allignment of the plot_2D plots: line or collumn
pic_extension = '.png'      # Choose extension (format) of the output images: .jpg, .png, .eps
#------------------------------------------------------------------------------------
# TEXT EDITOR
text_editor = 'gedit'       # Path to your favourite text editor you want to use during editing input files
#------------------------------------------------------------------------------------
# IMAGE VIEWER
image_viewer = 'okular'     # Path to yout favourite pdf viewer you want to use to check decomposition results
ds9_path = ''               # Path to SAO IMAGE DS9 program
#------------------------------------------------------------------------------------
# UNITS FOR OUTPUT RESULTS
mag_system = 'Vega'         # Magnitude system: Vega or AB
geom_units = 'arcsec'       # Geometrical units: pix, arcsec, or kpc
lum_units = 'mag'           # Luminosity units: ADU, mag, or Lsun
SB_units = 'mag/arcsec2'    # Surface brightness units: mag/arcsec2, ADU/pix2, or Lsun/pc^2
#====================================================================================



#====================================================================================
# DECA in process
ini_guess = 'Vika+2014'     # Choose method to find initial guess: 'Vika+2014' (see Vika et al. MNRAS 444,3603-3621(2014)) or 'DECA' (not recommended)
show_code_run = False       # Show run of the code: True or False
timeout = 600               # Timeout for Galfit and IMFIT: Decomposition will be interruped in [timeout] sec
timeoutM = 3600             # Timeout for GalfitM: Decomposition will be interruped in [timeoutM] sec
stop_deca = True            # Stop DECA if it faces some kind of error: True or False. Otherwise, an error code will be written in the log-file
#====================================================================================



#====================================================================================
# COSMOLOGY
H0 = 72.                    # THE HUBBLE CONSTANT IN km/s/Mpc
WM = 0.3                    # DARK MATTER 
WV = 0.7                    # VACUUM
#====================================================================================



#====================================================================================
# COMPONENTS
add_sky_comp = True	    # If True and sky background=0 (sky has been subtracted),
#			    then Sky component will be added (GALFIT) or --sky (IMFIT).
#			    If True and sky background!=0. (sky is still present), then
#			    Sky component will be created for both codes.
#			    If False and sky background=0, then no sky component will be added.
#			    If False and sky_background!=0, then Sky component will be created
#			    for both codes.
q_eon = 0.35	            # Maximal apparent flattening for edge-on galaxies
c0 = False		    # Fit boxyness: True or False
#coord_vary_perc = 5        # Let vary coordinates in %
imfit_disc_n = 1            # Disc vertical profile law: If n=1, then it is the sech^2 -law, if n->oo, then it is exp-law.


CENTRE_CONSTRS = {}
CENTRE_CONSTRS['COM_CENTRE'] = 'yes'    # set the same coordinates for all components (applicable only in GALFIT)
CENTRE_CONSTRS['MAX_OFFSET'] = '2pix'   # maximal offset from the input coordinates during the fitting in pix or arcsec, along the X and Y axes

# Constraints:
set_constraints = False # False - no constraints
                      # True - apply CONSTRS + CENTRE_CONSTRS
                      
                      
# Initial guess from scratch and constraints:
#COEFFS = collections.OrderedDict(); CONSTRS = collections.OrderedDict()
COEFFS = {}; CONSTRS = {}

# Recipes to get initial guess for a multicomponent model, based on the results of the single Sersic decomposition
# Constraints: CONSTRS[COMPONENT] = [...,[MIN_VALUE_k,MAX_VALUE_k],...] , where MIN_VALUE_k and MAX_VALUE_k are minimal and maximal values for the k-th parameter

#=========================== sersic ==========================
# For initial single Sersic decomposition (first guess) and elliptical galaxies
CONSTRS['sersic:sersic'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]
#=============================================================


#========================= exp_disc ==========================
# Model for late-type, bulgeless galaxies
COEFFS['exp_disc:exp_disc'] = ['+0.0','/1.68','+0.0','+0.0'] # [Mag, h, q, PA]
CONSTRS['exp_disc:exp_disc'] = [None,None,['0.05','1.0'],['0.0','360.0']] # [Mag, h, q, PA]
#=============================================================


#========================= eon_disc ==========================
# Model for edge-on bulgeless galaxies
COEFFS['eon_disc:eon_disc'] = ['+0.0','/1.68','+0.0'] # [Mag, h, PA]
CONSTRS['eon_disc:eon_disc'] = [None,None,['0.0','360.0']] # [Mag, h, PA]
#=============================================================


#=========================== ring ============================
# Model for pure-ring galaxies
COEFFS['ring:ring'] = ['+0.0','*0.2','*3.0','+0.0','+0.0'] # [m0, sigma_r, R_ring, q, PA]
CONSTRS['ring:ring'] = [None,None,None,['0.05','1.0'],['0.0','360.0']] # [m0, sigma_r, R_ring, q, PA]
#=============================================================


#====================== sersic+sersic ========================
# Model for galaxies which consist of two Sersic components (general model for spiral galaxies)
# Sersic1 (bulge)
COEFFS['sersic+sersic:sersic1'] = ['+0.75','*0.1','2.0','0.8','+0.0'] # [Mag, reb, n, q, PA]
CONSTRS['sersic+sersic:sersic1'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]

# Sersic2 (disc)
COEFFS['sersic+sersic:sersic2'] = ['+0.65','*1.0','1.0','+0.0','+0.0'] # [Mag, reb, n, q, PA]
CONSTRS['sersic+sersic:sersic2'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]
#==============================================================


#======================= sersic+exp_disc ======================
# Model for spiral galaxies
# Bulge
COEFFS['sersic+exp_disc:sersic'] = ['+0.75','*0.1','2.0','0.8','+0.0'] # [Mag, reb, n, q, PA]
CONSTRS['sersic+exp_disc:sersic'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]

# Disc
COEFFS['sersic+exp_disc:exp_disc'] = ['+0.65','/1.68','+0.0','+0.0'] # [Mag, h, q, PA]
CONSTRS['sersic+exp_disc:exp_disc'] = [None,None,['0.05','1.0'],['0.0','360.0']] # [Mag, h, q, PA]
#=============================================================


#======================= sersic+eon_disc =====================
# Model for edge-on spiral galaxies
# Bulge
COEFFS['sersic+eon_disc:sersic'] = ['+0.75','*0.1','2.0','0.8','+0.0'] # [Mag, reb, n, q, PA]
CONSTRS['sersic+eon_disc:sersic'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]

# Eon-disc
COEFFS['sersic+eon_disc:eon_disc'] = ['+0.65','/1.68','+0.0'] # [Mag, h, PA]
CONSTRS['sersic+eon_disc:eon_disc'] = [None,None,['0.0','360.0']] # [Mag, h, PA]
#=============================================================


#======================= agn+exp_disc ========================
# Model for late-type, bulgeless galaxies with AGN
# Agn
COEFFS['agn+exp_disc:agn'] = [None] # [Mag]
CONSTRS['agn+exp_disc:agn'] = [None] # [Mag]

# Disc
COEFFS['agn+exp_disc:exp_disc'] = ['+0.0','/1.68','+0.0','+0.0'] # [Mag, h, q, PA]
CONSTRS['agn+exp_disc:exp_disc'] = [None,None,['0.05','1.0'],['0.0','360.0']] # [Mag, h, q, PA]
#=============================================================



#======================= agn+eon_disc ========================
# Model for edge-on bulgeless galaxies with AGN
# Agn
COEFFS['agn+eon_disc:agn'] = [None] # [Mag]
CONSTRS['agn+eon_disc:agn'] = [None] # [Mag]

# Eon-disc
COEFFS['agn+eon_disc:eon_disc'] = ['+0.0','/1.68','+0.0'] # [Mag, h, PA]
CONSTRS['agn+eon_disc:eon_disc'] = [None,None,['0.0','360.0']] # [Mag, h, PA]
#=============================================================


#======================== agn+ring ===========================
# Model for ring galaxies
# Agn
COEFFS['agn+ring:agn'] = [None] # [Mag]
CONSTRS['agn+ring:agn'] = [None] # [Mag]

# Ring
COEFFS['agn+ring:ring'] = ['+0.0','*0.2','*3.0','+0.0','+0.0'] # [Mag, sigma_r, R_ring, q, PA]
CONSTRS['agn+ring:ring'] = [None,None,None,['0.05','1.0'],['0.0','360.0']] # [Mag, sigma_r, R_ring, q, PA]
#=============================================================


#======================= sersic+ring =========================
# Model for ring galaxies
# Host galaxy
COEFFS['sersic+ring:sersic'] = ['+0.1','*1.0','3.0','+0.0','+0.0'] # [Mag, reb, n, q, PA]
#CONSTRS['sersic+ring:sersic'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]
CONSTRS['sersic+ring:sersic'] = [None,None,None,None,['0.0','360.0']] # [Mag, reb, n, q, PA]
# Ring
COEFFS['sersic+ring:ring'] = ['+1.0','*0.2','*3.0','+0.0','+0.0'] # [Mag, sigma_r, R_ring, q, PA]
CONSTRS['sersic+ring:ring'] = [None,None,None,None,['0.0','360.0']] # [Mag, sigma_r, R_ring, q, PA]
#=============================================================



#====================== ferrer+sersic ========================
# Model for bulgeless bar galaxies
# Bar
COEFFS['ferrer+sersic:ferrer'] = ['+3.0','*0.2','4.0','2.0','+0.0','+0.0'] # [m0d, Rad, Alpha, Beta, q, PA]
CONSTRS['ferrer+sersic:ferrer'] = [None,None,None,None,['0.05','1.0'],['0.0','360.0']] # [m0d, Rad, Alpha, Beta, q, PA]

# Disc
COEFFS['ferrer+sersic:sersic'] = ['+0.65','*1.0','1.0','+0.0','+0.0'] # [Mag, reb, n, q, PA]
CONSTRS['ferrer+sersic:sersic'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]
#=============================================================
  


#=================== sersic+sersic+sersic ====================
# Model for galaxies which consist of three Sersic components (general model for spiral galaxies)
# Sersic1 (bulge)
COEFFS['sersic+sersic+sersic:sersic1'] = ['+1.0','*0.1','2.0','0.8','+0.0'] # [Mag, reb, n, q, PA]
CONSTRS['sersic+sersic+sersic:sersic1'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]

# Sersic2 (inner disc or bar)
COEFFS['sersic+sersic+sersic:sersic2'] = ['+0.75','*0.3','1.0','+0.0','+0.0'] # [Mag, reb, n, q, PA]
CONSTRS['sersic+sersic+sersic:sersic2'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]

# Sersic3 (outer disc)
COEFFS['sersic+sersic+sersic:sersic3'] = ['+1.0','*1.0','1.0','+0.0','+0.0'] # [Mag, reb, n, q, PA]
CONSTRS['sersic+sersic+sersic:sersic3'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]
#=============================================================



#==================== agn+sersic+exp_disc ====================
# Model for spiral galaxies with AGN
# Agn
COEFFS['agn+sersic+exp_disc:agn'] = [None] # [Mag]
CONSTRS['agn+sersic+exp_disc:agn'] = [None] # [Mag]

# Bulge
COEFFS['agn+sersic+exp_disc:sersic'] = ['+0.75','*0.1','2.0','0.8','+0.0'] # [Mag, reb, n, q, PA]
CONSTRS['agn+sersic+exp_disc:sersic'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]

# Disc
COEFFS['agn+sersic+exp_disc:exp_disc'] = ['+0.65','/1.68','+0.0','+0.0'] # [Mag, h, q, PA]
CONSTRS['agn+sersic+exp_disc:exp_disc'] = [None,None,['0.05','1.0'],['0.0','360.0']] # [Mag, h, q, PA]
#=============================================================



#==================== agn+sersic+eon_disc ====================
# Model for edge-on spiral galaxies with AGN
# Agn
COEFFS['agn+sersic+eon_disc:agn'] = [None] # [Mag]
CONSTRS['agn+sersic+eon_disc:agn'] = [None] # [Mag]

# Bulge
COEFFS['agn+sersic+eon_disc:sersic'] = ['+0.75','*0.1','2.0','0.8','+0.0'] # [Mag, reb, n, q, PA]
CONSTRS['agn+sersic+eon_disc:sersic'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]

# Eon-disc
COEFFS['agn+sersic+eon_disc:eon_disc'] = ['+0.65','/1.68','+0.0'] # [Mag, h, PA]
CONSTRS['agn+sersic+eon_disc:eon_disc'] = [None,None,['0.0','360.0']] # [Mag, h, PA]
#=============================================================


#================== sersic+ferrer+sersic =====================
# Model for spiral galaxies with a bar
# Bulge
COEFFS['sersic+ferrer+sersic:sersic1'] = ['+0.75','*0.3','2.0','+0.0','+0.0'] # [Mag, reb, n, q, PA]
CONSTRS['sersic+ferrer+sersic:sersic1'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]

# Bar
COEFFS['sersic+ferrer+sersic:ferrer'] = ['+3.0','*0.2','4.0','2.0','+0.0','+0.0'] # [m0d, Rad, Alpha, Beta, q, PA]
CONSTRS['sersic+ferrer+sersic:ferrer'] = [None,None,None,None,['0.05','1.0'],['0.0','360.0']] # [m0d, Rad, Alpha, Beta, q, PA]

# Disc
COEFFS['sersic+ferrer+sersic:sersic2'] = ['+0.65','*1.0','1.0','+0.0','+0.0'] # [Mag, reb, n, q, PA]
CONSTRS['sersic+ferrer+sersic:sersic2'] = [None,None,['0.3','15.0'],['0.05','1.0'],['0.0','360.0']] # [Mag, reb, n, q, PA]
#=============================================================
#====================================================================================
