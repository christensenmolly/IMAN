#!/usr/bin/env python
# Main script to run DECA
# *****************************************************************
# **       DECA -- DEComposition Analysis of galaxy images       **
# **          Astronomical Observatory, Ghent University         **
# *****************************************************************
# USAGE: python ~/MEGA/MyPrograms/IMAN/DECA/NEW/deca.py --i fit_sample.csv --g 11 all --n 1

# Import standard modules
import pylab
import sys
import os
import shutil
import math
import numpy as np
import scipy as sp
from numpy import *
from pylab import *
import subprocess
import argparse
import glob
import pyfits
from multiprocessing import Process
from joblib import Parallel, delayed
import filecmp
from astropy.stats import sigma_clipped_stats
import shelve
import urllib2
import datetime

import rebin_image




# Path to additional packages:
PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0].split('DECA')[0]

FNULL = open(os.devnull, 'w')

tmp_out = sys.stdout
warnings.filterwarnings("ignore")
from misc_functions import * 

# -----------------------------------------------------------------
# Import config file for DECA: default or specified by user
if os.path.exists('deca_setup.py')==False:
    OUT_INFO('\nThe default setup file is being used!')
    import deca_setup
else:
    current_dir = os.getcwd()
    sys.path.insert(0, current_dir)
    import deca_setup
    OUT_INFO('\nThe user-created setup file is being used!')
# -----------------------------------------------------------------


# Import DECA modules
import do_mode
import redo_mode
import do_mode_galaxies
import galfitm_input
import psf_fit
import deca_tk
import read_data
import imfit_input


# -----------------------------------------------------------------
# FUNCTION TO READ IN DECA SETUP FILE
def read_setup():    
    # MAIN_SETUP
    MAIN_SETUP = {}
    code = deca_setup.code
    minAlgorithm = deca_setup.minAlgorithm
    imfitPath = deca_setup.imfitPath
    galfitPath = deca_setup.galfitPath
    galfitmPath = deca_setup.galfitmPath
    addImfitStr = deca_setup.addImfitStr
    
    MAIN_SETUP['code'] = code
    MAIN_SETUP['minAlgorithm'] = minAlgorithm
    MAIN_SETUP['imfitPath'] = imfitPath
    MAIN_SETUP['galfitPath'] = galfitPath
    MAIN_SETUP['addImfitStr'] = addImfitStr
    
    # GENTIC ALGORITHM SETUP
    GA_SETUP = {}
    numOfCores = deca_setup.numOfCores
    ZeroGenSize = deca_setup.ZeroGenSize
    Popsize = deca_setup.Popsize
    selectNbest = deca_setup.selectNbest
    addNew = deca_setup.addNew
    maxGenNumber = deca_setup.maxGenNumber
    fTol = deca_setup.fTol
    fSpan = deca_setup.fSpan
    saveGens = deca_setup.saveGens
    runLM = deca_setup.runLM
    numOfLM = deca_setup.numOfLM
    LMCores = deca_setup.LMCores
    numOfFits = deca_setup.numOfFits
    
    GA_SETUP['numOfCores'] = numOfCores
    GA_SETUP['ZeroGenSize'] = ZeroGenSize
    GA_SETUP['Popsize'] = Popsize
    GA_SETUP['selectNbest'] = selectNbest
    GA_SETUP['addNew'] = addNew
    GA_SETUP['maxGenNumber'] = maxGenNumber
    GA_SETUP['fTol'] = fTol
    GA_SETUP['fSpan'] = fSpan
    GA_SETUP['saveGens'] = saveGens
    GA_SETUP['runLM'] = runLM
    GA_SETUP['numOfLM'] = numOfLM
    GA_SETUP['LMCores'] = LMCores
    GA_SETUP['numOfFits'] = numOfFits  
    
    return MAIN_SETUP,GA_SETUP
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO READ IN DECA INPUT FILE
def read_input_file(inputDECA, galaxy):
    NUMBER,FILTER,NAME,GALAXY_IMAGE,MASK_IMAGE,SIGMA_IMAGE,PSF_IMAGE,SAMPLING,DISTANCE,FITTING_PROC = \
    loadtxt(inputDECA, usecols=[0,1,2,3,4,5,6,7,8,9],dtype=str, unpack=True,skiprows=1,delimiter=',')

    idx_deca = numpy.where(NUMBER_d==galaxy)[0]

    IMAGES = [GALAXY_IMAGE[idx_deca],MASK_IMAGE[idx_deca],SIGMA_IMAGE[idx_deca],PSF_IMAGE[idx_deca]]
    INFO = [NUMBER[idx_deca],FILTER[idx_deca],NAME[idx_deca],SAMPLING[idx_deca],DISTANCE[idx_deca],FITTING_PROC[idx_deca]]
    return IMAGES,INFO
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO RUN DECA-TK
def run_deca_tk(observation_info, object_info, add_info, new_images):
	OUT_INFO('Tk mode will be launched')
	parameters_tk = shelve.open('parameters_tk.inp')

	[nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,Sky_level,SkySubtr,sampling,convolution_box,FILTER] = observation_info
	[xc,yc,name,distance,Aext,Kext] = object_info
	[number,fitting_proc,mode,model,region_to_fit] = add_info
	[input_image,sigma_image,psf_image,mask_image] = new_images


	inputParamsgeneral = {}
	inputParamsgeneral["image_fits"] = input_image
	inputParamsgeneral["weight_fits"] = sigma_image
	inputParamsgeneral["psf_fits"] = psf_image
	inputParamsgeneral["mask_fits"] = mask_image

	inputParamsgeneral["magZP"] = m0
	inputParamsgeneral["scale"] = scale
	inputParamsgeneral["GAIN"] = gain
	inputParamsgeneral["NCOMBINE"] = ncombine
	inputParamsgeneral["EXPTIME"] = exptime
	inputParamsgeneral["READNOISE"] = read_out_noise
	inputParamsgeneral["FWHM"] = fwhm
	inputParamsgeneral["Sky"] = Sky_level
	inputParamsgeneral["SkySubtr"] = SkySubtr
	inputParamsgeneral["Sampling"] = sampling

	inputParamsgeneral["RA"] = xc
	inputParamsgeneral["DEC"] = yc
	inputParamsgeneral["NAME"] = name
	inputParamsgeneral["Distance"] = distance
	inputParamsgeneral["Extinction"] = Aext
	inputParamsgeneral["Kcorr"] = Kext
	inputParamsgeneral["CoordsCheck"] = 1
	inputParamsgeneral["Filter_index"] = FILTERS_LIST('band').index(FILTER)

	parameters_tk["inputParamsgeneral"] = inputParamsgeneral
	parameters_tk.close()
	
	# Launch TKinter
	deca_tk.main(input_file='parameters_tk.inp')
# -----------------------------------------------------------------



# -----------------------------------------------------------------
# FUNCTION TO RUN GALFITM
def run_galfitm(control, del_files, FITTING_PROC, NUMBER, FILTER, GALAXY_IMAGE, SIGMA_IMAGE, PSF_IMAGE, MASK_IMAGE, REGION_TO_FIT, galaxy_name, k, current_dir, path_to_save_out, log_text):
	considered_galaxy = list(np.where(NUMBER == NUMBER[k])[0])    # get indices for this galaxy in the input table

	last_fits = []
	number_of_do = 0
	done_galaxies = []
	for ii in considered_galaxy:
            # Check if this filter should be taken into account
            if FITTING_PROC[ii]=='DO' or FITTING_PROC[ii]=='REDO':
                # Find out if this filter has been done and remember path to its output file
                done_galaxies.append(ii)
                number_of_do = number_of_do + 1
                arr = glob.glob(deca_setup.res_directory+'/%s_%s/%s/galfit_*' % (NUMBER[ii],galaxy_name,FILTER[ii]))
                fits = []
                try:
                            for iii in range(len(arr)):
                                fits.append(int(arr[iii].split('galfit_')[-1]))
                            last_fits.append(max(fits)) # Get last number of the fits
                except:
                            OUT_FAIL('The filter %s does not (yet?) have a decomposition!' % (FILTER[ii]))
                            return 1
        '''        
	if number_of_do==len(last_fits):
                # The number of 'DO' filters is equal to the done galaxies
                if len(set(last_fits))!=1:
                        # If all filters do not have the same number of the output directory 
                        #f_out_gm = open(current_dir+'/'+path_to_save_out,'a')
			#print >>f_out_gm, '%s\t%s\t%s\t%i' % (NUMBER[i],galaxy_name,FILTER[i],2) # status=2
			#f_out_gm.close()
			status = 8
			write_log_deca(current_dir+'/'+path_to_save_out, NUMBER[k], galaxy_name, '_'.join(filters), status, log_text)
			#return 1                    
	else:
                # Not all filters have been done
                #return 0
                z=1
        '''  
	if number_of_do!=len(last_fits):
            return 1
        
	#done_galaxies = considered_galaxy  ####TODO: This is not necessarilly so!
		
	try:
                filters = FILTER[done_galaxies]  
                if len(done_galaxies)>1: 
			gm_dir = deca_setup.res_directory+'/'+str(NUMBER[k])+'_'+galaxy_name+'/'+'_'.join(filters)

                        if FITTING_PROC[k]=='REDO' and os.path.exists(gm_dir):
                            # Remove the output galaxy directory which was created before
                            shutil.rmtree(gm_dir)

			if not os.path.exists(gm_dir):
			      os.makedirs(gm_dir)

			# 2. Copy the images
			COPIED_IMAGES = []
			FILTERS = []
			NX = []
			NY = []
			MAGZPS = []
			SCALES = []
			SAMPLINGS = []
			CONVOLUTION_BOXES = []
			new_images = []
			for i in done_galaxies:
			  IMAGES = [GALAXY_IMAGE[i],SIGMA_IMAGE[i],PSF_IMAGE[i],MASK_IMAGE[i]]
			  copied_images = copy_images(IMAGES, gm_dir+'/', FILTER[i])
			  COPIED_IMAGES.append(copied_images)
			  FILTERS.append(FILTER[i])
			  [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,Sky_level,SkySubtr,sampling,convolution_box],xc,yc,log_txt = read_header(GALAXY_IMAGE[i], MASK_IMAGE[i], PSF_IMAGE[i], None, log_txt)
			  MAGZPS.append(str(m0-2.5*log10(exptime)))
			  SCALES.append(scale)
			  SAMPLINGS.append(str(sampling))
			  CONVOLUTION_BOXES.append(str(convolution_box))
			  NX.append(nx)
			  NY.append(ny)
			  new_images.append(IMAGES)

			new_images = np.transpose(np.array(new_images))                           
			if len(set(NX))!=1 or len(set(NY))!=1 or len(set(SCALES))!=1 or len(set(SAMPLINGS))!=1 or len(set(CONVOLUTION_BOXES))!=1:
			    OUT_FAIL('Galaxy images have different dimensions! Exiting...')
			    exit()

			# 3. Copy the galfit files with initial guess for each band
			GALFIT_INPUT_FILES = []
			BANDS = []; COMP = []; PAR = []; CONSTR = []
			for i in done_galaxies:
				  arr = glob.glob(deca_setup.res_directory+'/%s_%s/%s/galfit_*' % (NUMBER[i],galaxy_name,FILTER[i]))
				  fits = []
				  for ii in range(len(arr)):
				    fits.append(int(arr[ii].split('_')[-1]))
				  last_fit = max(fits)
				  shutil.copy(arr[fits.index(last_fit)]+'/galfit.01',gm_dir+'/galfit_%s.01' % (FILTER[i]))
				  GALFIT_INPUT_FILES.append('galfit_%s.01' % (FILTER[i]))

				  if os.path.exists(arr[fits.index(last_fit)]+'/constraint.txt'):
					comp,par,constr = constraints.read_constraints(arr[fits.index(last_fit)]+'/constraint.txt')
					COMP.append(comp)
					PAR.append(par)
					CONSTR.append(constr)
					BANDS.append(FILTER[i]))

                        if COMP!=[]:
                            ff = open(gm_dir+'/constraint.txt', 'w')
                            ff.write('# Component/    parameter   constraint	Comment\n')
                            ff.write('# operation	(see below)   range\n\n')
                            for kk in range(len(COMP)):
                                for ii in range(len(COMP[kk])):
                                    ff.write('    %s              %s_%s        %s\n' % (COMP[kk][ii],PAR[kk][ii],BANDS[kk],CONSTR[kk][ii]))
                            ff.close()
                            constr_file = 'constraint.txt'
			else:
                            constr_file = None
		
			# 3. Change the directory
			os.chdir(gm_dir)

			# 4. Create galfitm input file
			f = open('galfitm.inp', "w") 
			sys.stdout = f  
			galfitm_input.header(COPIED_IMAGES, FILTERS, constr_file, 'model.fits', 1, nx, 1, ny, MAGZPS, SCALES, SAMPLINGS, CONVOLUTION_BOXES)
			galfitm_input.merge_components_data(GALFIT_INPUT_FILES, len(MAGZPS))
			sys.stdout = tmp_out
			f.close()

			os.chmod(r"galfitm.inp",0777)    
			#subprocess.call(deca_setup.galfitmPath+'galfitm' + " galfitm.inp", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)

                        if delay==False:
                            status,log_text = do_mode.launch_compl_decomposition(control, del_files, log_text, new_images, None, None, run_line=None)
                        else:
                            status,log_text = do_mode.build_initial_model(control, del_files, log_text, new_images, None, None, run_line=None)


			write_log_deca(current_dir+'/'+path_to_save_out, NUMBER[k], galaxy_name, '_'.join(filters), status, log_text)
			os.chdir(current_dir)
                else: 
		    z=1
	except:
		  os.chdir(current_dir)
		  #f_out_gm = open(path_to_save_out,'a')
		  #print >>f_out_gm, '%s\t%s\t%s\t%i' % (NUMBER[i],galaxy_name,FILTER[i],1)
		  #f_out_gm.close()
		  status = 9
		  
		  write_log_deca(current_dir+'/'+path_to_save_out, NUMBER[k], galaxy_name, '_'.join(filters), status, log_text)
# -----------------------------------------------------------------    


# -----------------------------------------------------------------  
# FUNCTION TO PRINT OUT RESULTS OF DECOMPOSITION
def write_log_deca(outFilePath, galaxy_number, galaxy_name, galaxy_filter, status, log_txt):
    f_out = open(outFilePath, 'a')
    header = '==========================================================\n'
    header = header + 'RESULTS for %s,%s,%s: %i\n' % (galaxy_number,galaxy_name,galaxy_filter,status)
    log_txt = header + log_txt
    log_txt = log_txt + '=========================================================='
    f_out.write('%s\n\n' % (log_txt))
    f_out.close()
# -----------------------------------------------------------------  
      

# -----------------------------------------------------------------
# MAIN FUNCTION TO DO EACH ITEM IN THE INPUT TABLE
def do_each_galaxy(k, data, last_model):
                galaxy_dir = deca_setup.res_directory+"/%s_%s/%s/%s_%s" % (str(int(data["NUMBER"][k])), data["NAME"][k], data["FILTER"][k], deca_setup.code.lower(),last_model[k])
                res_file = galaxy_dir + '/' + deca_setup.code.lower() + '.01'
                
                
                
                os.chdir(deca_setup.res_directory+"/%s_%s/%s" % (str(int(data["NUMBER"][k])), data["NAME"][k], data["FILTER"][k]))
                os.copy(res_file, deca_setup.code.lower() + '_err.inp')
                new_inp_file = deca_setup.code.lower() + '_err.inp'
                
                os.mkdir('./psf_unc')
                
                old_psf = data["PSF_IMAGE"][k].split('/')[-1].split('.fits')[0]+'_old.fits'
                new_psf = data["PSF_IMAGE"][k].split('/')[-1]
                shutil.move(new_psf, old_psf)
                remove_file('log.txt')
                remove_file('galfit.01')
                remove_file('imfit.01')
                
                err = ['plus', 'minus']
                
                for er in err:
                    # Rebin psf image
                    if er=='plus':
                        rebin_image.downsample(old_psf, 1.+data["PSF_STD_FWHM"][k], new_psf, set_wcs=True)
                    else:
                        rebin_image.downsample(old_psf, 1.-data["PSF_STD_FWHM"][k], new_psf, set_wcs=True)
                    
                    hdulist = pyfits.open(new_psf)
                    data = hdulist_in[0].data
                    ny,nx = np.shape(data)
                    
                    if deca_setup.code=='GALFIT':
                        f = open(res_file, 'r')
                        ff = open(new_inp_file, 'w')
                        lines = f.readlines()
                        for line in lines:
                            if '# Size of the convolution box (x y)' in line and 'I)' in line:
                                ff.write('I) %i     %i           # Size of the convolution box (x y)\n' % (nx,ny))
                            else:
                                ff.write(line)
                    elif deca_setup.code=='IMFIT':
                        f = open(res_file, 'r')
                        lines = f.readlines()
                        string = lines[0].split('#')[-1].strip()
                    f.close()
                    ff.close()
                    
                    if deca_setup.code=='GALFIT':
                        status,log_text = do_mode.run_galfit(new_inp_file, '')
                        if status==0:
                            shutil.move('galfit.01', './psf_unc/galfit_%s.01' % (er)) 
                    elif deca_setup.code=='IMFIT':
                        status,log_text = run_imfit(string, '')
                        if status==0:
                            shutil.move('imfit.01', './psf_unc/imfit_%s.01' % (er)) 
                    remove_file('log.txt')
                    remove_file('model.fits')
                shutil.move(old_psf, new_psf)                
# -----------------------------------------------------------------












if __name__ == '__main__':
    # Arguments:
    parser = argparse.ArgumentParser(description="DEComposition Analysis")   
    parser.add_argument("--i", default='sample.csv',
                    help="input file for DECA", required=True)
    
    parser.add_argument("--m", help="output file with built models")
    
    parser.add_argument("--g", default=None, nargs='+',
                    help="select galaxies by their number in the input file, along with their filters separated by comma (or type all for all filters given for the galaxies). Example: --g 1 all 2 u,g,r ")

    parser.add_argument("--n", default=1,
                    help="number of Python processes (1, by default)", required=False)

    # Parse arguments
    args = parser.parse_args()

    inputDECA = args.i
    
    models = args.m
    
    if args.g!=None:
      galaxies = args.g

      galaxy_numbers = []; galaxy_filters = []
      for k in range(len(galaxies)):
	if k%2==0:
	  galaxy_numbers.append(galaxies[k])
	else:
	  galaxy_filters.append(galaxies[k])
    else:
      galaxy_numbers = None
      galaxy_filters = None
   
    n_jobs = int(args.n)


    # Salutation
    OUT_TITLE("\n\t\t\t*****************************************")
    OUT_TITLE("\t\t\t                PSF ERRORS                ")
    OUT_TITLE("\t\t\t        A. Mosenkov & S. Savchenko          ")
    OUT_TITLE("\t\t\t*****************************************")

    # Read the setup file for DECA
    # Notice that if you run DECA from the directory in which deca_stup.py file is located,
    # it will be read by DECA instead of the default deca_setup.py file located in the DECA directory!
    MAIN_SETUP,GA_SETUP = read_setup()


    data,units = read_data.main(inputDECA, delimiter=',', header_line=0, units_line=None, skip_lines = [])
    
    try:
        PSF_STD_FWHM = np.array(data["PSF_STD_FWHM"],float)
    except:
        OUT_FAIL('The column %s is not found in the input file! Exiting' % ("PSF_STD_FWHM"))
        exit()


    last_model,model,numbers = np.loadtxt(models, usecols=[3,4,5],dtype=str, unpack=True,skiprows=1,delimiter='\t')
    Parallel(n_jobs=n_jobs)(delayed(do_each_galaxy)(i, data, last_model) for i in range(len(PSF_STD_FWHM)))
