#!/usr/bin/env python
# Module to do photometric decomposition
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **        © Astronomical Observatory, Ghent University         **
# *****************************************************************

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
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import glob
from PIL import Image

# Import astro modules
import pyfits
from photutils import datasets
from astropy.stats import sigma_clipped_stats
from photutils import aperture_photometry
from photutils import EllipticalAperture
from astropy.table import hstack
 
# Import DECA modules
import deca_setup
import models
import galfit_input
import imfit_input
import initial_guess
import save_results
import threading
import signal
import constraints

# Import IMP modules
sys.path.append(PATH_TO_PACKAGE+'IMP_NEW')
import convert_segm_to_region
import mask_indiv

tmp_out = sys.stdout
FNULL = open(os.devnull, 'w')

# Path to DECA package 
PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0].split('DECA')[0]

# deca_setup variables
galfitPath = deca_setup.galfitPath
imfitPath = deca_setup.imfitPath
show_code_run = deca_setup.show_code_run

# -----------------------------------------------------------------
# CLASS TO EXECUTE A COMMAND WITH TIMEOUT
class RunCmd(threading.Thread):
    # Taken from http://stackoverflow.com/questions/4789837/how-to-terminate-a-python-subprocess-launched-with-shell-true
    def __init__(self, cmd, timeout):
        threading.Thread.__init__(self)
        self.cmd = cmd
        self.timeout = timeout

    def run(self):
	if show_code_run==True:
	  self.p = subprocess.Popen(self.cmd, shell=True,preexec_fn=os.setsid)#, stdout=FNULL, stderr=subprocess.STDOUT)
	else:
	  self.p = subprocess.Popen(self.cmd, shell=True, stdout=FNULL, stderr=subprocess.STDOUT,preexec_fn=os.setsid)
        self.p.wait()

    def Run(self):
        self.start()
        self.join(self.timeout)

        if self.is_alive():
            #self.p.kill()#p.terminate()      #use self.p.kill() if process needs a kill -9
            os.killpg(os.getpgid(self.p.pid), signal.SIGTERM)
            self.join()
            return 1    # If timeout happenes, the process will be killed. This error code will be returned.
        else:
            return 0
# -----------------------------------------------------------------

'''
def RunCmd(cmd):
  from __future__ import print_function
  from threading import Timer
  #from subprocess import Popen, PIPE, STDOUT
  #from subprocess import call

  def terminate(process):
      if process.poll() is None:
	  subprocess.call('taskkill /F /T /PID ' + str(process.pid))
  # start process
  process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
		  bufsize=1, universal_newlines=True)

  # terminate process in 15 seconds
  timer = Timer(timeout, terminate, args=[process])
  timer.start()

  # print output
  for line in iter(process.stdout.readline, ''):
      print(line, end='')
  process.stdout.close()
  process.wait() # wait for the child process to finish
  timer.cancel()
'''

# -----------------------------------------------------------------
# FUNCTION TO CREATE A PDF-FILE WITH MAIN RESULTS PRESENTED
def crea_report(files, output_pdf):
    from reportlab.pdfgen import canvas
    from reportlab.lib.units import inch, cm
    from reportlab.lib.pagesizes import letter, A4
    imgDoc = canvas.Canvas(output_pdf)
    
    imgDoc.setPageSize(A4) # This is actually the default page size
    document_width, document_height = A4
    '''
    for file in files:
        imgDoc.drawImage(file, 10, 10, width=10*cm, preserveAspectRatio=True)
        imgDoc.showPage()
    '''
    for image in files:
        # Open the image file to get image dimensions
        Image_file = Image.open(image)
        image_width, image_height = Image_file.size
        image_aspect = image_height / float(image_width)

        # Determine the dimensions of the image in the overview
        print_width = document_width
        print_height = document_width * image_aspect

        # Draw the image on the current page
        # Note: As reportlab uses bottom left as (0,0) we need to determine the start position by subtracting the
        #       dimensions of the image from those of the document
        imgDoc.drawImage(image, document_width - print_width, document_height - print_height, width=print_width,
                         height=print_height)
        imgDoc.showPage()        
    imgDoc.save()    
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO LAUNCH GALFIT WITH TIMEOUT
def run_galfit(file, log_text):
  if show_code_run==True:
    #subprocess.call(galfitPath+'galfit ' + file, shell=True)
    status = RunCmd([galfitPath+'galfit '+file], deca_setup.timeout).Run()
  else:
    #subprocess.call(galfitPath+'galfit ' + file, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    status = RunCmd([galfitPath+'galfit '+file], deca_setup.timeout).Run()
  
  if status==1:
      misc_functions.OUT_FAIL('Timeout happened! Decomposition was terminated!')
      log_text = log_text + 'Timeout happened! Decomposition was terminated!\n'
  if not os.path.exists('model.fits'):
      misc_functions.OUT_FAIL('Galfit crashed.')
      log_text = log_text + 'Galfit crashed.\n'
  return status, log_text
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO LAUNCH IMFIT WITH TIMEOUT
def run_imfit(string, log_text):
  if show_code_run==True:
    #subprocess.call(string, shell=True)
    status = RunCmd([string], deca_setup.timeout).Run()
  else:
    #subprocess.call(string, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    status = RunCmd([string], deca_setup.timeout).Run()
  
  if status==1:
      misc_functions.OUT_FAIL('Timeout happened! Decomposition was terminated!')
      log_text = log_text + 'Timeout happened! Decomposition was terminated!\n'
  if not os.path.exists('model.fits'):
      misc_functions.OUT_FAIL('Imfit crashed.')
      log_text = log_text + 'Imfit crashed.\n'
  return status, log_text
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO LAUNCH GALFITM WITH TIMEOUT
def run_galfitm(file, log_text):
  if show_code_run==True:
    #subprocess.call(galfitmPath+'galfitm ' + file, shell=True)
    status = RunCmd([galfitmPath+'galfitm '+file], deca_setup.timeoutM).Run()
  else:
    #subprocess.call(galfitmPath+'galfitm ' + file, shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
    status = RunCmd([galfitmPath+'galfitm '+file], deca_setup.timeoutM).Run()
  
  if status==1:
      misc_functions.OUT_FAIL('Timeout happened! Decomposition was terminated!')
      log_text = log_text + 'Timeout happened! Decomposition was terminated!\n'
  if not os.path.exists('model.fits'):
      misc_functions.OUT_FAIL('GalfitM crashed.')
      log_text = log_text + 'GalfitM crashed.\n'
  return status, log_text
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO FIND THE NEAREST ELEMENT IN THE ARRAY
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO DETERMINE THE CODE OF THE INPUT MODEL FILE
def determine_code(model):
    f = open(model, 'r')
    
    galfit_model = False
    imfit_model = False
    for line in f:
        if 'A)' in line and '# Input data image (FITS file)' in line:
            galfit_model = True
            dec_code = 'galfit'
            break
        if 'FUNCTION' in line:
            imfit_model = True
            dec_code = 'imfit'
            break
    if galfit_model == False and imfit_model == False:
        print 'Cannot recognize the input file %s!' % (model)
        exit()
    return dec_code 
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO FIND FLUXES WITHIN AN ELLIPSE TO ESTIMATE TOTAL LUMINOSITY AND EFFECTIVE RADIUS
def find_eff_radius(data, mask_data, xc, yc, sma, smb, theta):
  Sma = [max([sma/30.,2.]),sma/20.,sma/10.,sma/2.,sma]
  Sma = np.array(Sma)
  Smb = Sma * smb/sma
  fluxes = []

  for k in range(len(Sma)):
    apertures = EllipticalAperture([(xc,yc)], Sma[k],Smb[k],radians(theta))
    fluxes.append(aperture_photometry(data, apertures, mask=mask_data))
  phot_table = hstack(fluxes)

  sum1 = float(phot_table['aperture_sum_1'])
  sum2 = float(phot_table['aperture_sum_2'])
  sum3 = float(phot_table['aperture_sum_3'])
  sum4 = float(phot_table['aperture_sum_4'])
  sum5 = float(phot_table['aperture_sum_5'])
  Sum = np.array([sum1,sum2,sum3,sum4,sum5])

  f2 = interp1d(Sma, Sum)
  x = np.linspace(max([sma/30.,2.]), sma, num=1000)
  re = x[find_nearest(f2(x),sum5/2.)]
  return re,sum5,data[int(yc),int(xc)]
# -----------------------------------------------------------------  
    
  
# -----------------------------------------------------------------
# FUNCTION TO REPEAT SERSIC DECOMPOSITION IF NEEDED
def sersic_crash(control, log_text, run_line=None):
    if deca_setup.code=='GALFIT':
                while not os.path.exists('model.fits'):
                    if control==True:
                            misc_functions.OUT_FAIL('Galfit crashed. Please change the input file!')
                            try:
                                shutil.move('galfit.01','galfit.inp')
                            except:
                                z=1
                            subprocess.call(deca_setup.text_editor + ' galfit.inp' , shell=True)
                            raw_input("Press Enter to continue...")
                            status,log_text = run_galfit('galfit.inp', log_text)
                            
                    else:
                        #log_text = log_text + 'Sersic decomposition crashed!\n'
                        return 5,log_text
                return 0,log_text                    
    
    if deca_setup.code=='IMFIT':                    
                while not os.path.exists('model.fits'):
                    if control==True:
                            misc_functions.OUT_FAIL('Imfit crashed. Please change the input file!')
                            subprocess.call(deca_setup.text_editor + ' imfit.inp' , shell=True)
                            raw_input("Press Enter to continue...")
                            status,log_text = run_imfit(run_line, log_text)
                    else:
                        #log_text = log_text + 'Sersic decomposition crashed!\n'
                        return 5,log_text                    
                return 0,log_text                      
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO CREATE NEW DIRECTORY WHERE NEW RESULTS WILL BE OUTPUT
def crea_new_dir():
            arr = glob.glob('%s_*' % (deca_setup.code.lower()))

            dir_numbers = []
            for ar in arr:
                dir_numbers.append(int(ar.split('_')[-1]))
            if dir_numbers==[]:
                new_dir = './%s_2' % (deca_setup.code.lower())
            else:
                new_dir = './%s_' % (deca_setup.code.lower()) + str(max(dir_numbers)+1) 

            return new_dir
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO LAUNCH A COMPLICATED DECOMPOSITION
def launch_compl_decomposition(control, del_files, log_text, new_images, scale, m0, run_line=None):                
                if deca_setup.code=='GALFIT':
                    input_model_file = 'galfit.inp'
                    output_model_file = 'galfit.01'
                elif deca_setup.code=='IMFIT':
                    input_model_file = 'imfit.inp'
                    output_model_file = 'imfit.01'
                elif deca_setup.code=='GALFITM':
                    input_model_file = 'galfitm.inp'
                    output_model_file = 'model.galfit.01.band'
                    
                repeat = True                
                while repeat==True:
                        if os.path.exists('model.fits'):
                            os.remove('model.fits')
                        
                        if deca_setup.code=='GALFIT':                            
                            status,log_text = run_galfit(input_model_file, log_text)
                        elif deca_setup.code=='IMFIT':
                            status,log_text = run_imfit(run_line, log_text)
                        elif deca_setup.code=='GALFITM':
                            status,log_text = run_galfitm(input_model_file, log_text)
                            
                        while not os.path.exists('model.fits'):
                            if control==True:
                                    print 'Decomposition crashed. Please change the input file!'
                                    try:
                                        shutil.move(output_model_file, input_model_file)
                                    except:
                                        z=1 
                                    subprocess.call(deca_setup.text_editor + ' %s' % (input_model_file) , shell=True)
                                    raw_input("Press Enter to continue...")
                                    if deca_setup.code=='GALFIT':                            
                                        status,log_text = run_galfit(input_model_file, log_text)
                                    elif deca_setup.code=='IMFIT':
                                        status,log_text = run_imfit(run_line, log_text)
                                    elif deca_setup.code=='GALFITM':
                                        status,log_text = run_galfitm(input_model_file, log_text)                            
                            else:
                                    if deca_setup.code!='GALFITM':
                                        #log_text = log_text + 'More complicated model decomposition crashed!\n'
                                        return 6, log_text
                                    else:
                                        log_text = log_text + 'Multi-band decomposition crashed!\n'
                                        return 7, log_text
                        
                        # For the first fit, when a user model is given, the output decomposition will be stored in the galfit_2 or imfit_2 directory!
                        '''
                        if os.path.exists('./galfit_1')==False and deca_setup.code=='GALFIT':
                            NEW_DIR = './galfit_2'
                        elif os.path.exists('./imfit_1')==False and deca_setup.code=='IMFIT':
                            NEW_DIR = './imfit_2'  
                        else:
                            NEW_DIR = None
                        '''
                        NEW_DIR = crea_new_dir()
                        new_dir,status = save_results.main(new_images,scale,m0,new_dir=NEW_DIR) ####TODO!!!
                        
                        if control==True:
                                        pictures = glob.glob(new_dir + '/*.eps') + glob.glob(new_dir + '/*.png')
                                        '''
                                        for picture in pictures:                                                                                            
                                            img = Image.open(picture)
                                        img.show()
                                        plt.clf()
                                        plt.close()
                                        plt.close('all')
                                        '''
                                        
                                        output_pdf = 'report.pdf'
                                        crea_report(pictures, output_pdf)
                                        if sys.platform=="darwin":
                                            os.system("open %s" % (output_pdf))
                                        else:
                                            subprocess.call(deca_setup.image_viewer + ' %s' % (output_pdf), shell=True)
                                        
                                        
                                        '''
                                        imgs = [Image.open(i) for i in pictures ]
                                        imgs_comb = np.vstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
                                        imgs_comb = Image.fromarray( imgs_comb)
                                        imgs_comb.show()
                                        '''
                                        try:
                                            with open('fit.log', 'r') as fin:
                                                print fin.read()
                                        except:
                                            if deca_setup.code=='GALFIT': 
                                                with open('galfit.log', 'r') as fin:
                                                    print fin.read()
                                            elif deca_setup.code=='IMFIT': 
                                                with open('imfit.log', 'r') as fin:
                                                    print fin.read()  
                                            elif deca_setup.code=='GALFITM': 
                                                with open('galfitm.log', 'r') as fin: ####??????????
                                                    print fin.read()
                                                    
                                        Happy = 'yes'
                                        Happy = str(raw_input('Are you happy with the results?  (%s): ' % (str(Happy)))) or str(Happy)
                                        if Happy == 'no':
                                            Change_Mask = 'no'
                                            Change_Mask = str(raw_input('Do you want to change the mask (not available with GALFITM) [no]?  (%s): ' % (str(Change_Mask)))) or str(Change_Mask)
                                            if Change_Mask=='yes' and deca_setup.code!='GALFITM':
                                                new_mask_file = new_images[3].split('.fits')[0] + '_' + str(int(crea_new_dir().split('_')[-1]))+'.fits'
                                                
                                                convert_segm_to_region.main(new_images[3], 'mask_' + str(int(crea_new_dir().split('_')[-1])) +'.reg', fits_slice = 0, scale=None, offset=None, ignore_value=None)
                                                print 'Please modify the mask file mask_%s.reg!' % (str(int(crea_new_dir().split('_')[-1])))
                                                ds9Proc = subprocess.Popen([deca_setup.ds9_path + "ds9", os.path.abspath(new_images[0]),"-regions", os.path.abspath('mask_' + str(int(crea_new_dir().split('_')[-1])) +'.reg'),"-scale", "log"])
                                                ds9Proc.wait()                                                 
                                                mask_indiv.mask(new_images[0], 'mask_' + str(int(crea_new_dir().split('_')[-1])) +'.reg', factor=1, output_image=None, output_mask=new_mask_file, obj_ellipse=None, mask_value=1, show_running=True, mask_DN=None)
                                                #os.remove('tmp_mask.reg')
                                                #exit()
                                                if deca_setup.code=='GALFIT':                            
                                                    fff = open(input_model_file, 'r')
                                                    ffff = open('tmp.inp', 'w')
                                                    lines = fff.readlines()
                                                    for k in range(len(lines)):
                                                        if 'F)' in lines[k] and '# Bad pixel mask (FITS image or ASCII coord list)' in lines[k]:
                                                            ffff.write( 'F) %s           # Bad pixel mask (FITS image or ASCII coord list)' % (new_mask_file))
                                                        else:
                                                            ffff.write(lines[k])
                                                    ffff.close()
                                                    fff.close()
                                                    shutil.move('tmp.inp', input_model_file)
                                                    
                                                elif deca_setup.code=='IMFIT':
                                                    run_line = run_line.replace(new_images[3], new_mask_file)
                                              
                                            Change_Model = '1'
                                            Change_Model = str(raw_input('Do you want to change the input [1] or output model [2]?  (%s): ' % (str(Change_Model)))) or str(Change_Model)
                                            if Change_Model=='1':
                                                if os.path.exists(input_model_file):
                                                    subprocess.call(deca_setup.text_editor + ' %s' % (input_model_file) , shell=True)
                                                    raw_input("Press Enter to continue...")
                                            else:
                                                if os.path.exists(output_model_file):
                                                    try:
                                                        shutil.move(output_model_file,input_model_file)
                                                    except:
                                                        z=1 
                                                    subprocess.call(deca_setup.text_editor + ' %s' % (input_model_file) , shell=True)
                                                    raw_input("Press Enter to continue...")                                            
                                            
                                            
                            
                                            repeat = True
                                            #new_dir,status = save_results.main(new_images,scale,m0)
                                            if del_files==True:
                                                save_results.remove_files_func(new_images, del_all=False)
                                            for file in ['fit.log','galfit.01','imfit.01','composed_model.fits','model.galfit.01','model.galfit.01.band']:
                                                if os.path.exists(file): 
                                                    os.remove(file)
                                                    
                                            arr = glob.glob('galfit_*.01')
                                            for file in arr:
                                                os.remove(file)
                                        else:
                                            repeat = False
                                            #new_dir,status = save_results.main(new_images,scale,m0)
                                            if del_files==True:
                                                os.remove(output_pdf)
                                                save_results.remove_files_func(new_images, del_all=True)
                                            for file in ['fit.log','galfit.01','imfit.01','composed_model.fits','model.galfit.01','model.galfit.01.band']:
                                                if os.path.exists(file): 
                                                    os.remove(file)
                                            
                                            arr = glob.glob('galfit_*.01')
                                            for file in arr:
                                                os.remove(file)
                                            
                                            if deca_setup.code!='GALFITM':  
                                                log_text = log_text + 'More complicated model decomposition: OK\n'
                                            else:  
                                                log_text = log_text + 'Multi-band decomposition: OK\n'
                                                
                                            return status,log_text
                                        os.remove(output_pdf)
                        else:
                                        repeat = False
                                        if del_files==True:
                                            save_results.remove_files_func(new_images, del_all=True)

                                        for file in ['fit.log','galfit.01','imfit.01','composed_model.fits','model.galfit.01','model.galfit.01.band']:
                                            if os.path.exists(file): 
                                                os.remove(file)

                                        arr = glob.glob('galfit_*.01')
                                        for file in arr:
                                                os.remove(file)

                                        if deca_setup.code!='GALFITM':  
                                                log_text = log_text + 'More complicated model decomposition: OK\n'
                                        else:  
                                                log_text = log_text + 'Multi-band decomposition: OK\n'
                                        return status,log_text
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO BUILD AN INITIAL MODEL
def build_initial_model(control, del_files, log_text, new_images, scale, m0, run_line=None):                
                if deca_setup.code=='GALFIT':
                    input_model_file = 'galfit.inp'
                    output_model_file = 'galfit.01'
                elif deca_setup.code=='IMFIT':
                    input_model_file = 'imfit.inp'
                    output_model_file = 'imfit.01'
                elif deca_setup.code=='GALFITM':
                    input_model_file = 'galfitm.inp'
                    output_model_file = 'model.galfit.01.band'

                repeat = True                
                while repeat==True:
                        if os.path.exists('model.fits'):
                            os.remove('model.fits')
                        
                        if deca_setup.code=='GALFIT':                            
                            subprocess.call("%sgalfit -o1 %s" % (galfitPath, input_model_file), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                        elif deca_setup.code=='IMFIT':
                            subprocess.call("%smakeimage %s -o model.fits --refimage %s --psf %s" % (imfitPath, input_model_file,new_images[0], new_images[2]), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                        elif deca_setup.code=='GALFITM':                            
                            subprocess.call("%sgalfitm -o1 %s" % (galfitmPath, input_model_file), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                        
                        while not os.path.exists('model.fits'):
                            if control==True:
                                    print 'Decomposition crashed. Please change the input file!'
                                    try:
                                        shutil.move(output_model_file,input_model_file)
                                    except:
                                        z=1 
                                    subprocess.call(deca_setup.text_editor + ' %s' % (input_model_file) , shell=True)
                                    raw_input("Press Enter to continue...")
                                    if deca_setup.code=='GALFIT':                            
                                        subprocess.call("%sgalfit -o1 %s" % (galfitPath, input_model_file), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                                    elif deca_setup.code=='IMFIT':
                                        subprocess.call("%smakeimage %s -o model.fits --refimage %s --psf %s" % (imfitPath, input_model_file,new_images[0], new_images[2]), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                                    elif deca_setup.code=='GALFITM':                            
                                        subprocess.call("%sgalfitm -o1 %s" % (galfitmPath, input_model_file), shell=True, stdout=FNULL, stderr=subprocess.STDOUT)                            
                            else:
                                    if deca_setup.code!='GALFITM':
                                        log_text = log_text + 'More complicated model decomposition crashed!\n'
                                        return 6,log_text
                                    else:
                                        log_text = log_text + 'Multi-band decomposition crashed!\n'
                                        return 9,log_text
                                    
                        # For the first fit, when a user model is given, the output decomposition will be stored in the galfit_2 or imfit_2 directory!
                        #NEW_DIR = crea_new_dir()
                        new_dir,status = save_results.main(new_images,scale,m0,new_dir='') #### TODO!!!!!!!!!!!!

                        if control==True:
                                        pictures = glob.glob(new_dir + '*.eps') + glob.glob(new_dir + '*.png')
                                        '''
                                        for picture in pictures:                                                                                            
                                            img = Image.open(picture)
                                        img.show()
                                        plt.clf()
                                        plt.close()
                                        plt.close('all')
                                        '''
                                        
                                        output_pdf = 'report.pdf'
                                        crea_report(pictures, output_pdf)
                                        if sys.platform=="darwin":
                                            os.system("open %s" % (output_pdf))
                                        else:
                                            subprocess.call(deca_setup.image_viewer + ' %s' % (output_pdf), shell=True)
                                        
                                        
                                        '''
                                        imgs = [Image.open(i) for i in pictures ]
                                        imgs_comb = np.vstack( (np.asarray( i.resize(min_shape) ) for i in imgs ) )
                                        imgs_comb = Image.fromarray( imgs_comb)
                                        imgs_comb.show()
                                        '''
 
                                        Happy = 'yes'
                                        Happy = str(raw_input('Are you happy with the results?  (%s): ' % (str(Happy)))) or str(Happy)
                                        arr = []
                                        if Happy == 'no':
                                            Change_Mask = 'no'
                                            Change_Mask = str(raw_input('Do you want to change the mask (not available with GALFITM) [no]?  (%s): ' % (str(Change_Mask)))) or str(Change_Mask)
                                            if Change_Mask=='yes' and deca_setup.code!='GALFITM':
                                                new_mask_file = new_images[3].split('.fits')[0] + '_' + str(int(crea_new_dir().split('_')[-1]))+'.fits'
                                                
                                                convert_segm_to_region.main(new_images[3], 'mask_' + str(int(crea_new_dir().split('_')[-1])) +'.reg', fits_slice = 0, scale=None, offset=None, ignore_value=None)
                                                print 'Please modify the mask file mask_%s.reg!' % (str(int(crea_new_dir().split('_')[-1])))
                                                ds9Proc = subprocess.Popen([deca_setup.ds9_path + "ds9", os.path.abspath(new_images[0]),"-regions", os.path.abspath('mask_' + str(int(crea_new_dir().split('_')[-1])) +'.reg'),"-scale", "log"])
                                                ds9Proc.wait()                                                 
                                                mask_indiv.mask(new_images[0], 'mask_' + str(int(crea_new_dir().split('_')[-1])) +'.reg', factor=1, output_image=None, output_mask=new_mask_file, obj_ellipse=None, mask_value=1, show_running=True, mask_DN=None)
                                                #os.remove('tmp_mask.reg')
                                                #exit()
                                                if deca_setup.code=='GALFIT':                            
                                                    fff = open(input_model_file, 'r')
                                                    ffff = open('tmp.inp', 'w')
                                                    lines = fff.readlines()
                                                    for k in range(len(lines)):
                                                        if 'F)' in lines[k] and '# Bad pixel mask (FITS image or ASCII coord list)' in lines[k]:
                                                            ffff.write( 'F) %s           # Bad pixel mask (FITS image or ASCII coord list)' % (new_mask_file))
                                                        else:
                                                            ffff.write(lines[k])
                                                    ffff.close()
                                                    fff.close()
                                                    shutil.move('tmp.inp', input_model_file)
                                                    
                                                elif deca_setup.code=='IMFIT':
                                                    run_line = run_line.replace(new_images[3], new_mask_file)
                                              
                                            Change_Model = '1'
                                            Change_Model = str(raw_input('Do you want to change the input [1] or output model [2]?  (%s): ' % (str(Change_Model)))) or str(Change_Model)
                                            if Change_Model=='1':
                                                if os.path.exists(input_model_file):
                                                    subprocess.call(deca_setup.text_editor + ' %s' % (input_model_file) , shell=True)
                                                    raw_input("Press Enter to continue...")
                                            else:
                                                if os.path.exists(output_model_file):
                                                    try:
                                                        shutil.move(output_model_file,input_model_file)
                                                    except:
                                                        z=1 
                                                    subprocess.call(deca_setup.text_editor + ' %s' % (input_model_file) , shell=True)
                                                    raw_input("Press Enter to continue...")                                            
                                            
                                            
                                            
                                            
                                            arr = glob.glob('%s_*.inp' % (deca_setup.code.lower()))

                                            file_numbers = []
                                            for ar in arr:
                                                file_numbers.append(int(ar.split('_')[-1]))
                
                                            
                                            shutil.copy(input_model_file, '%s_%s.inp' % (deca_setup.code.lower(), str(max(file_numbers)+1)))
                                            repeat = True

                                            for file in ['fit.log','galfit.01','imfit.01','composed_model.fits','model.galfit.01','model.galfit.01.band']+pictures:
                                                if os.path.exists(file): 
                                                    os.remove(file)
                                            arr = glob.glob('galfit_*.01')
                                            for file in arr:
                                                    os.remove(file)
                                                
                                        else:
                                            repeat = False
                                            for file in ['fit.log','galfit.01','imfit.01','composed_model.fits','model.galfit.01','model.galfit.01.band']+pictures+arr:
                                                if os.path.exists(file): 
                                                    os.remove(file)
                                            arr = glob.glob('galfit_*.01')
                                            for file in arr:
                                                    os.remove(file)
                                                    
                                            log_text = log_text + 'Input model: OK\n'
                                            return status,log_text
                                        os.remove(output_pdf)
                        else:
                                        repeat = False
                                        for file in ['fit.log','galfit.01','imfit.01','composed_model.fits','model.galfit.01','model.galfit.01.band']:
                                            if os.path.exists(file): 
                                                os.remove(file)
                                        arr = glob.glob('galfit_*.01')
                                        for file in arr:
                                                    os.remove(file)
                                                    
                                        log_text = log_text + 'Input model: OK\n'
                                        return status,log_text
# -----------------------------------------------------------------

                    
# -----------------------------------------------------------------
# MAIN FUNCTION
def main(new_images,observation_info,object_info,add_info,log_text,keys,file_with_galaxies,del_files=True):
  #*******************************************************************        
  [nx,ny,m0,scale,gain,ncombine,exptime,read_out_noise,fwhm,Sky_level,SkySubtr,sampling,convolution_box,FILTER] = observation_info
  [xc,yc,name] = object_info
  [number,fitting_proc,mode,model,region_to_fit] = add_info
  [input_image,sigma_image,psf_image,mask_image] = new_images
  [control,delay] = keys
  #*******************************************************************


  ## Read in the input galaxy image
  hdulist = pyfits.open(input_image)
  data = hdulist[0].data
  ny,nx = np.shape(data)

  ## Specify the region where to fit the object (cut-out)
  if region_to_fit!=None:                
                #  x1:x2;y1:y2 format
                xmin = int(region_to_fit.split(':')[0])  
                xmax = int(region_to_fit.split(':')[1].split(';')[0])  
                ymin = int(region_to_fit.split(';')[1].split(':')[0])  
                ymax = int(region_to_fit.split(':')[2])
  else:
                xmin = 0
                xmax = nx
                ymin = 0 
                ymax = ny
  
  ## Cut out the galaxy image and apply mask according to the specified region
  if mask_image!='none' and os.path.exists(mask_image):
      ## Read the input mask and convert it to astropy mask
      hdulist_mask = pyfits.open(mask_image)
      data_mask = hdulist_mask[0].data

      mask_astropy = np.zeros_like(np.array(data,dtype=float), dtype=bool)

      for k in range(ny):
	for i in range(nx):
	  if data_mask[k,i]!=0.:
	    mask_astropy[k,i] = True
  else:
      mask_astropy = np.zeros_like(np.array(data, dtype=float), dtype=bool) 

  data = data[ymin:ymax,xmin:xmax]
  mask_astropy = mask_astropy[ymin:ymax,xmin:xmax]



  # -----------------------------------------------------------------
  ## TWO WAYS: THE INPUT MODEL IS UNKNOWN ('none') OR SPECIFIED BY USER
  if '.' not in model:
    # Build model
    if True:    # TODO: Replace if True with try!
            if (os.path.exists('./galfit_1')==False and deca_setup.code=='GALFIT') or (os.path.exists('./imfit_1')==False and deca_setup.code=='IMFIT'):
                # -----------------------------------------------------------------
                #  ***************FIT WITH A SINGLE SERSIC FUNCTION****************
                ## Coarsely estimate the single sersic parameters
                # The galfit_01 and imfit_01 directories are reserved for single Sersic decomposition!                                              
                
                if file_with_galaxies==None:
                    ## Define the single Sersic model for one object:
                    model = [model]
                    PARS = models.define_model_galaxies(number_of_objects=1)
                    
                    ## Find coarse parameters for the initial guess
                    re,lum,I0 = find_eff_radius(data=data, mask_data=mask_astropy, xc=xc-xmin, yc=yc-ymin, sma=xc-xmin, smb=yc-ymin, theta=0.)

                    # Default parameters:
                    PARS["First_Ser_Mb"][0] = m0 - 2.5*log10(lum)
                    PARS["First_Ser_reb"][0] = re
                    PARS["First_Ser_n"][0] = 2.5
                    PARS["First_Ser_PA"][0] = 90. 
                    PARS["First_Ser_q"][0] = 1.
                    PARS["First_Ser_meb"][0] = initial_guess.meb_bulge_f(PARS["First_Ser_reb"][0]*scale,PARS["First_Ser_Mb"][0],PARS["First_Ser_n"][0],PARS["First_Ser_q"][0])
                else:
                    # Read in input objects from the region file
                    [xc,yc,ellA,ellB,ellPA,Name,Model,Redshift,AExt,KCorr],[xc_star,yc_star] = misc_functions.READ_REGION_OBJECTS(file_with_galaxies, input_image)
                    model = Model
                    xc = np.array(xc) - xmin
                    yc = np.array(yc) - ymin
                    xc_star = np.array(xc_star) - xmin
                    yc_star = np.array(yc_star) - ymin
                    
                    PARS = models.define_model_galaxies(len(xc))
                    PARS_stars = models.define_model_stars(len(xc_star))
                
                
                run_line = None

                ## Create input model and launch decomposition
                if deca_setup.code=='GALFIT':
                    ## Create constraints galfit file                    
                    if file_with_galaxies==None:
                        constraint_file, status, log_text = constraints.create_constraint_file('sersic', scale, log_text)
                        if status!=0:
                                return status,log_text
                    else:
                        constraint_file, status, log_text = constraints.create_constraint_file(len(model)*['sersic'], scale, log_text)
                        if status!=0:
                                return status, log_text


                    ## Create input galfit file
                    f = open('galfit.inp', "w") 
                    sys.stdout = f  
                    galfit_input.header(input_image,sigma_image,psf_image,mask_image,constraint_file,'model.fits',xmin,xmax,ymin,ymax,m0-2.5*log10(exptime),scale,int(sampling),convolution_box)
                    
                    if file_with_galaxies==None:
                        galfit_input.Sersic(component=1,xc=xc,yc=yc,meb=PARS["First_Ser_meb"][0],
                        Mb=PARS["First_Ser_Mb"][0],reb=PARS["First_Ser_reb"][0],
                        n=PARS["First_Ser_n"][0],q=PARS["First_Ser_q"][0],PA=PARS["First_Ser_PA"][0])
                        Ncomp = 1
                    else:
                        Ncomp = 0
                        for k in range(len(xc)):
                            ## Coarsely estimate the single sersic parameters
                            aperture = EllipticalAperture([(xc[k],yc[k])], ellA[k],ellB[k],radians(ellPA[k]))
                            phot_table = aperture_photometry(data, aperture, mask=mask_astropy)
                            if SkySubtr==1:
                                lum = phot_table['aperture_sum'][0]
                            else:
                                bkg_sum = Sky_level * aperture.area()
                                lum = phot_table['aperture_sum'][0] - bkg_sum

                            ## Define the single Sersic model:
                            galfit_input.Sersic(component=k+1,xc=xc[k],yc=yc[k],meb=float('nan'),
                                Mb = m0 - 2.5*log10(lum), reb=ellA[k]/5.,
                                n=2.5, q=ellB[k]/ellA[k], PA=ellPA[k])
                            Ncomp = Ncomp + 1

                        for k in range(len(xc_stars)):    
                            ## Coarsely estimate the single sersic parameters
                            aperture = EllipticalAperture([(xc_stars[k],xc_stars[k])], 2.*fwhm,2.*fwhm,0.)
                            phot_table = aperture_photometry(data, aperture, mask=mask_astropy)
                            if SkySubtr==1:
                                lum = phot_table['aperture_sum'][0]
                            else:
                                bkg_sum = Sky_level * aperture.area()
                                lum = phot_table['aperture_sum'][0] - bkg_sum

                            ## Define the single Sersic model:
                            galfit_input.Agn(component=len(xc)+k+1,xc=xc_stars[k],yc=yc_stars[k],mag=m0 - 2.5*log10(lum))
                            Ncomp = Ncomp + 1
		                        
                    
                    
                    galfit_input.Sky(component=Ncomp+1,sky_level=Sky_level,fit_sky=SkySubtr)
                    sys.stdout = tmp_out
                    f.close()

                    os.chmod(r"galfit.inp",0777)
                                
                    ## Fit data with the single Sersic model
                    if os.path.exists('model.fits'):
                        os.remove('model.fits')
                    status, log_text = run_galfit('galfit.inp', log_text)
                    
                elif deca_setup.code=='IMFIT':
                    ## Create input imfit file
                    run_line = imfit_input.imfit_exec_line(new_images,observation_info,object_info,add_info,'imfit.inp')
                    print run_line

                    f = open('imfit.inp', "w") 
                    sys.stdout = f
                    
                    if file_with_galaxies==None:
                        imfit_input.header(run_line,xc,yc)
                        imfit_input.Sersic(component=1,xc=xc,yc=yc,PA=PARS["First_Ser_PA"][0], ell=1.-PARS["First_Ser_q"][0],n=PARS["First_Ser_n"][0],I_e=10**(0.4*(m0-PARS["First_Ser_meb"][0]))*(scale**2),r_e=PARS["First_Ser_reb"][0],m0=m0,pix2sec=scale,model='sersic',write_coords=False)

                    
                        imfit_input.Sky(component=2,sky_level=Sky_level,fit_sky=SkySubtr)
                        sys.stdout = tmp_out
                        f.close()
                    else:
                        imfit_input.header(run_line,0.,0.)

                        for k in range(len(xc)):
                            ## Coarsely estimate the single sersic parameters
                            aperture = EllipticalAperture([(xc[k],yc[k])], ellA[k],ellB[k],radians(ellPA[k]))
                            phot_table = aperture_photometry(data, aperture, mask=mask_astropy)
                            if SkySubtr==1:
                                lum = phot_table['aperture_sum'][0]
                            else:
                                bkg_sum = Sky_level * aperture.area()
                                lum = phot_table['aperture_sum'][0] - bkg_sum

                            ## Define the single Sersic model:
                            imfit_input.Sersic(component=k+1,xc=xc[k],yc=yc[k],meb=float('nan'),
                                Mb = m0 - 2.5*log10(lum),reb=ellA[k]/5.,
                                n=2.,q=ellB[k]/ellA[k],PA=ellPA[k],m0=m0,pix2sec=scale,write_coords=False)

                        
                        imfit_input.Sky(component=len(xc)+1,sky_level=Sky_level)
                        sys.stdout = tmp_out
                        f.close()


                    os.chmod(r"imfit.inp",0777)
                    
                    ## Fit data with the single Sersic model
                    if os.path.exists('model.fits'):
                        os.remove('model.fits')
                    #print run_line
                    log_text = log_text + run_line + '\n'
                    status, log_text = run_imfit(run_line, log_text)
                
                    
                ## Check if decomposition failed
                status,log_text = sersic_crash(control, log_text, run_line=run_line)
                if status==5:
                        return 5,log_text


                if file_with_galaxies==None:
                    ## Read and save results
                    if deca_setup.code=='GALFIT':                
                        hdulist_model = pyfits.open('model.fits')
                        model_header = hdulist_model[2].header                
                        PARS = models.read_results_First_Ser(PARS,'galfit.01',scale)
                        new_dir,status = save_results.main(new_images,scale,m0,new_dir='./galfit_1')
                    elif deca_setup.code=='IMFIT':              
                        PARS = models.read_results_First_Ser_IMFIT(PARS,'imfit.01',m0,scale)
                        new_dir,status = save_results.main(new_images,scale,m0,new_dir='./imfit_1')
                else:
                    ## Read and save results
                    if deca_setup.code=='GALFIT':
                        PARS = models.read_results_First_Ser(PARS,'galfit.01',scale)
                        PARS_stars = models.read_results_stars(PARS_stars,'galfit.01')
                        new_dir,status = save_results.main(new_images,scale,m0,new_dir='./galfit_1')
                    elif deca_setup.code=='IMFIT':
                        PARS = models.read_results_First_Ser_IMFIT(PARS,'imfit.01',m0,scale)
                        PARS_stars = models.read_results_stars_IMFIT(PARS_stars,'imfit.01')
                        new_dir,status = save_results.main(new_images,scale,m0,new_dir='./imfit_1')
                
                
                
                
                if del_files==True:
                    save_results.remove_files_func(new_images, del_all=False)
                
                for file in ['fit.log','galfit.01','imfit.01','composed_model.fits']:
                    if os.path.exists(file): 
                        os.remove(file)
                        
                if (model[0]=='sersic' and len(model)==1) or (model[0]=='sersic' and len(set(the_list))==1):
                    log_text = log_text + 'Sersic decomposition: OK\n'
                    return status,log_text
                log_text = log_text + 'Sersic decomposition: OK\n'
            '''
            if file_with_galaxies==None:                
              if (os.path.exists('./galfit_1')==True and deca_setup.code=='GALFIT') or (os.path.exists('./imfit_1')==True and deca_setup.code=='IMFIT'):
                ## Define the single Sersic model:
                PARS = models.define_model()
                
                if deca_setup.code=='GALFIT':                            
                    PARS = models.read_results_First_Ser(1,PARS,'./galfit_1/galfit.01',scale)
                elif deca_setup.code=='IMFIT':              
                    PARS = models.read_results_First_Ser_IMFIT(1,PARS,'./imfit_1/imfit.01',m0,scale)
            '''
            # -----------------------------------------------------------------





            # -----------------------------------------------------------------
            #  ***************FIT TO A MORE COMPLICATED MODEL************
            for kk in range(len(model)):
                if model[kk] == 'none':
                    if PARS["First_Ser_q"][kk]<=0.3:
                        if PARS["First_Ser_n"][kk]<=1.2:
                            model[kk] = 'eon_disc'
                        else:
                            model[kk] = 'sersic+eon_disc'
                    else:
                        if PARS["First_Ser_n"][kk]<=1.2:
                            model[kk] = 'exp_disc'
                        else:
                            model[kk] = 'sersic+sersic'
        
            
            # First guess
            if deca_setup.ini_guess=='Vika+2014':
                PARS = initial_guess.initial_guess_from_scratch(PARS,scale,model)  
            else:
                PARS = initial_guess.initial_guess_deca(new_images,observation_info,object_info,add_info,PARS,model)

            run_line = None
            if deca_setup.code=='GALFIT':            
                ## Create constraints galfit file
                constraint_file, status, log_text = constraints.create_constraint_file(model, scale, log_text)
                if status!=0:
                        return status,log_text
                
                # Create input galfit file
                f = open('galfit.inp', "w") 
                sys.stdout = f  
                galfit_input.header(input_image,sigma_image,psf_image,mask_image,constraint_file,'model.fits',xmin,xmax,ymin,ymax,m0-2.5*log10(exptime),scale,int(sampling),convolution_box)
                
                
                
                comp_numb = 1
                for k in range(len(PARS["First_Ser_meb"])):
                        # Define each component in the given model
                        if np.isnan(PARS["PSF_M"][k])==False:
                                galfit_input.Agn(component=comp_numb,xc=xc[k],yc=yc[k],mag=PARS["PSF_M"][k])
                                comp_numb = comp_numb + 1
                        if np.isnan(PARS["First_Ser_reb"][k])==False:
                                galfit_input.Sersic(component=comp_numb,xc=xc[k],yc=yc[k],meb=PARS["First_Ser_meb"][k],
                                    Mb=PARS["First_Ser_Mb"][k],reb=PARS["First_Ser_reb"][k],
                                    n=PARS["First_Ser_n"][k],q=PARS["First_Ser_q"][k],PA=PARS["First_Ser_PA"][k])
                                comp_numb = comp_numb + 1
                        if np.isnan(PARS["First_Disk_z0"][k])==True and np.isnan(PARS["First_Disk_hd"][k])==False:
                                galfit_input.ExpDisc(component=comp_numb,xc=xc[k],yc=yc[k],m0d=PARS["First_Disk_m0d"][k],
                                                Md=PARS["First_Disk_Md"][k],h=PARS["First_Disk_hd"][k],q=PARS["First_Disk_q"][k],PA=PARS["First_Disk_PA"][k])
                                comp_numb = comp_numb + 1
                        if np.isnan(PARS["First_Disk_z0"][k])==False:
                                galfit_input.EdgeDisc(component=comp_numb,xc=xc[k],yc=yc[k],m0d=PARS["First_Disk_m0d"][k],
                                                z0=PARS["First_Disk_z0"][k],h=PARS["First_Disk_hd"][k],PA=PARS["First_Disk_PA"][k])
                                comp_numb = comp_numb + 1

                        if np.isnan(PARS["FERRER_m0"][k])==False:
                                galfit_input.Ferrer(component=comp_numb,xc=xc[k],yc=yc[k],m0d=PARS["FERRER_m0"][k],rtr=PARS["FERRER_Rad"][k],alpha=PARS["FERRER_Alfa"],beta=PARS["FERRER_Beta"][k],q=PARS["FERRER_q"][k],PA=PARS["FERRER_PA"][k])
                                comp_numb = comp_numb + 1 

                        if np.isnan(PARS["Second_Ser_reb"][k])==False:
                                galfit_input.Sersic(component=comp_numb,xc=xc[k],yc=yc[k],meb=PARS["Second_Ser_meb"][k],
                                    Mb=PARS["Second_Ser_Mb"][k],reb=PARS["Second_Ser_reb"][k],
                                    n=PARS["Second_Ser_n"][k],q=PARS["Second_Ser_q"][k],PA=PARS["Second_Ser_PA"][k])
                                comp_numb = comp_numb + 1
                        if np.isnan(PARS["Third_Ser_reb"][k])==False:
                                galfit_input.Sersic(component=comp_numb,xc=xc[k],yc=yc[k],meb=PARS["Third_Ser_meb"][k],
                                    Mb=PARS["Third_Ser_Mb"][k],reb=PARS["Third_Ser_reb"][k],
                                    n=PARS["Third_Ser_n"][k],q=PARS["Third_Ser_q"][k],PA=PARS["Third_Ser_PA"][k])
                                comp_numb = comp_numb + 1
                        

                for k in range(len(xc_stars)):    
                            galfit_input.Agn(component=len(xc)+k+1,xc=xc_stars[k],yc=yc_stars[k],mag=PARS_stars["PSF_M"][k])
                            comp_numb = comp_numb + 1             
                
                
                # Add sky if needed           
                if deca_setup.add_sky_comp==True:
                    galfit_input.Sky(component=comp_numb,sky_level=Sky_level,fit_sky=SkySubtr)
                
                sys.stdout = tmp_out
                f.close()

                os.chmod(r"galfit.inp",0777)

            elif deca_setup.code=='IMFIT':
                # Create input imfit file
                run_line = imfit_input.imfit_exec_line(new_images,observation_info,object_info,add_info,'imfit.inp')
                log_text = log_text + run_line + '\n'

             
                f = open('imfit.inp', "w") 
                sys.stdout = f

                imfit_input.header(run_line,xc,yc)

                # Define each component in the given model                    
                comp_numb = 1
                for k in range(len(PARS["First_Ser_meb"])):
                    if np.isnan(PARS["PSF_M"][k])==False:
                        imfit_input.Agn(component=comp_numb,xc=xc,yc=yc,mag=PARS["PSF_M"][k],m0=m0,pix2sec=scale,model=model,write_coords=False)
                        comp_numb = comp_numb + 1
                    if np.isnan(PARS["First_Ser_reb"][k])==False:
                        imfit_input.Sersic(component=comp_numb,xc=xc,yc=yc,PA=PARS["First_Ser_PA"][k], ell=1.-PARS["First_Ser_q"][k],n=PARS["First_Ser_n"][k],I_e=10**(0.4*(m0-PARS["First_Ser_meb"][k]))*(scale**2),r_e=PARS["First_Ser_reb"][k],m0=m0,pix2sec=scale,model=model,write_coords=False)
                        comp_numb = comp_numb + 1
                    if np.isnan(PARS["First_Disk_z0"][k])==True and np.isnan(PARS["First_Disk_hd"][k])==False:
                        #imfit_input.BrokenExpDisc(component=comp_numb,xc=xc,yc=yc,PA=PARS["First_Disk_PA"],ell=1.-PARS["First_Disk_q"],I_0=10**(0.4*(m0-PARS["First_Disk_m0d"]))*(scale**2),h=PARS["First_Disk_hd"],m0=m0,pix2sec=scale,model=model,write_coords=False)
                        
                        imfit_input.ExpDisc(component=comp_numb,xc=xc,yc=yc,PA=PARS["First_Disk_PA"][k],ell=1.-PARS["First_Disk_q"][k],I_0=10**(0.4*(m0-PARS["First_Disk_m0d"][k]))*(scale**2),h=PARS["First_Disk_hd"][k],m0=m0,pix2sec=scale,model=model,write_coords=False)
                        
                        comp_numb = comp_numb + 1
                    if np.isnan(PARS["First_Disk_z0"][k])==False:
                        imfit_input.EdgeDisc(component=comp_numb,xc=xc,yc=yc,PA=PARS["First_Disk_PA"][k],L_0=(10**(0.4*(m0-PARS["First_Disk_m0d"][k]))) * (scale**2) / (2.*PARS["First_Disk_hd"][k]),h=PARS["First_Disk_hd"][k],z_0=PARS["First_Disk_z0"][k],m0=m0,pix2sec=scale,model=model,write_coords=False)      
                        comp_numb = comp_numb + 1
                    if np.isnan(PARS["Second_Ser_reb"][k])==False:
                        imfit_input.Sersic(component=comp_numb,xc=xc,yc=yc,PA=PARS["Second_Ser_PA"][k], ell=1.-PARS["Second_Ser_q"][k],n=PARS["Second_Ser_n"][k],I_e=10**(0.4*(m0-PARS["Second_Ser_meb"][k]))*(scale**2),r_e=PARS["Second_Ser_reb"][k],m0=m0,pix2sec=scale,model=model,write_coords=False,comp_number='2')
                        comp_numb = comp_numb + 1
                    if np.isnan(PARS["Third_Ser_reb"][k])==False:
                        imfit_input.Sersic(component=comp_numb,xc=xc,yc=yc,PA=PARS["Third_Ser_PA"][k], ell=1.-PARS["Third_Ser_q"][k],n=PARS["Third_Ser_n"][k],I_e=10**(0.4*(m0-PARS["Third_Ser_meb"][k]))*(scale**2),r_e=PARS["Third_Ser_reb"][k],m0=m0,pix2sec=scale,model=model,write_coords=False,comp_number='3')
                        comp_numb = comp_numb + 1
                    if np.isnan(PARS["Ring_m0"][k])==False:                   
                        imfit_input.GaussRing(component=comp_numb,xc=xc,yc=yc,PA=PARS["Ring_PA"][k],ell=1.-PARS["Ring_q"][k],A=10**(0.4*(m0-PARS["Ring_m0"][k]))*(scale**2),R_ring=PARS["Ring_R_r"][k],sigma_r=PARS["Ring_sigma_r"][k],m0=m0,pix2sec=scale,model=model,write_coords=False)
                        comp_numb = comp_numb + 1


                # Add sky if needed        
                imfit_input.Sky(component=2,sky_level=Sky_level,fit_sky=SkySubtr)
                
                sys.stdout = tmp_out
                f.close()

                os.chmod(r"imfit.inp",0777)
            if delay==False:
                status,log_text = launch_compl_decomposition(control,del_files,log_text,new_images,scale,m0, run_line=run_line)
            else:
                status,log_text = build_initial_model(control, del_files, log_text, new_images, scale, m0, run_line=None)
            return status,log_text
            # -----------------------------------------------------------------
            
    else:
        log_text = log_text + 'More complicated model decomposition crashed!\n'
        return 6,log_text
  # -----------------------------------------------------------------
  
  else:
            # -----------------------------------------------------------------
            # If user provides an input galfit/imfit file
            
            # Find out to which code this input file refers
            dec_code = determine_code(model)
            print 'The specified input file %s refers to the %s code' % (model, dec_code)
            run_line = None
            if deca_setup.code=='GALFIT':
                ff = open(model, 'r')
                fff = open("galfit.inp", 'w')
                for line in ff:
                    if '# Output data image block' not in line:
                        fff.write(line)
                    else:
                        fff.write('B) model.fits      # Output data image block\n')
                ff.close()
                fff.close()
            elif deca_setup.code=='IMFIT':
                shutil.move(model,'imfit.inp')
                # Create input imfit file
                run_line = imfit_input.imfit_exec_line(new_images,observation_info,object_info,add_info,'imfit.inp')
                print run_line
                log_text = log_text + run_line + '\n'
            if delay==False:                
                status,log_text = launch_compl_decomposition(control,del_files,log_text,new_images,scale,m0,run_line=run_line)
            else:
                status,log_text = build_initial_model(control, del_files, log_text, new_images, scale, m0, run_line=None)
            return status,log_text
            # -----------------------------------------------------------------
      
