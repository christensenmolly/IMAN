#!/usr/bin/env python
# Module to create final pdf file 
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **        © Astronomical Observatory, Ghent University         **
# *****************************************************************

# Import standard modules
import os
import numpy as np
import pyfits
import subprocess
import shutil
import collections
import sys
import glob
import argparse
from multiprocessing import Process
from joblib import Parallel, delayed
FNULL = open(os.devnull, 'w')

# Import DECA modules
DECA_PATH = os.path.dirname(__file__)
sys.path.append(DECA_PATH+'/deca_tk_lib')
import read_data
import galfit_parser
import galfit_input
import imfit_input
import imfit_parser
import misc
import tex_creator
import deca
import misc_functions
import create_final_table

PATH_TO_SCRIPT = os.path.dirname(os.path.realpath(__file__))
PATH_TO_PACKAGE = os.path.split(PATH_TO_SCRIPT)[0].split('DECA')[0]

sys.path.append(PATH_TO_PACKAGE+'Cosmo_funcs')
import cosmo_calc_NED_func

# -----------------------------------------------------------------
# Import config file for DECA: default or specified by user
if os.path.exists('deca_setup.py')==False:
    print('\nThe default setup file is being used!')
    import deca_setup
else:
    current_dir = os.getcwd()
    sys.path.insert(0, current_dir)
    import deca_setup
    print('\nThe user-created setup file is being used!')
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO CALCULATE LUMINOSITIES AND LUMINOSITY RATIOS
def find_luminosities(composed_image, components=None, galfitm=False):
  if galfitm==False:
   if components==None:
      hdu = pyfits.open(composed_image)
      number_of_layers = len(hdu)
      fluxes = []
      for k in range(number_of_layers):
	  if k==1:
	    Flux = np.sum(hdu[k].data)
	  if k>=4:
	    fluxes.append(np.sum(hdu[k].data))
      luminosities = []
      for k in range(len(fluxes)):
	luminosities.append(fluxes[k]/Flux)
      luminosities.append(Flux)
      return luminosities
   else:
      hdu = pyfits.open(composed_image)

      fluxes = []
      for k in components:
	    fluxes.append(np.sum(hdu[k].data))
      Flux = np.sum(fluxes)
      luminosities = []
      for k in components:
	luminosities.append(fluxes[k]/Flux)
      luminosities.append(Flux)
      return luminosities
  else:
   if components==None:
      hdu = pyfits.open(composed_image)
      header = hdu[0].header
      number_of_bands = header['number_of_bands']
      number_of_comps = header['number_of_comps']
      

      fluxes = []
      Flux = []
      for k in range(number_of_bands):            
	    Flux.append(np.sum(hdu[k*(4+number_of_comps)+1].data))
            fl = []
            for i in range(number_of_comps):
                fl.append(np.sum(hdu[k*(4+number_of_comps)+4+i].data))
            fluxes.append(fl)
      
      luminosities = []
      for k in range(number_of_bands):
          lum = []
          for i in range(number_of_comps):
            lum.append(fluxes[i]/Flux[k])
            lum.append(Flux[k])
          luminosities.append(lum)
      return luminosities
   else:
      hdu = pyfits.open(composed_image)
      header = hdu[0].header
      number_of_bands = header['number_of_bands']
      number_of_comps = header['number_of_comps']
      

      fluxes = []
      Flux = []
      for k in range(number_of_bands):            
	    Flux.append(np.sum(hdu[k*(4+number_of_comps)+1].data))
            fl = []
            for i in components:
                fl.append(np.sum(hdu[k*(4+number_of_comps)+4+i].data))
            fluxes.append(fl)
      
      luminosities = []
      for k in range(number_of_bands):
          lum = []
          for i in components:
            lum.append(fluxes[i]/Flux[k])
            lum.append(Flux[k])
          luminosities.append(lum)
      return luminosities
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO COMPILE TEX FILE
def compile_tex():
    if deca_setup.show_code_run==False:
                if deca_setup.pic_extension == '.eps':
                    subprocess.call('latex results.tex', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('latex results.tex', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('latex results.tex', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('dvips -Ppdf results.dvi -o', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('ps2pdf results.ps results.pdf', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('rm -f *.aux *.log ', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                else:
                    subprocess.call('pdflatex results.tex', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('pdflatex results.tex', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('pdflatex results.tex', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('rm -f *.aux *.log ', shell=True, stdout=FNULL, stderr=subprocess.STDOUT)   
    else:
                if deca_setup.pic_extension == '.eps':
                    subprocess.call('latex results.tex', shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('latex results.tex', shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('latex results.tex', shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('dvips -Ppdf results.dvi -o', shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('ps2pdf results.ps results.pdf', shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('rm -f *.aux *.log ', shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)
                else:
                    subprocess.call('pdflatex results.tex', shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('pdflatex results.tex', shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('pdflatex results.tex', shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)
                    subprocess.call('rm -f *.aux *.log ', shell=True)#, stdout=FNULL, stderr=subprocess.STDOUT)                     
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO DO EACH GALAXY
def do_each_galaxy(kkk, data, fil, last_model, model, numbers, geom_units, lum_units, SB_units, galfitm):
            Bands = misc_functions.FILTERS_LIST(out='band')
            Waves = misc_functions.FILTERS_LIST(out='wave')
            
            current_dir = os.getcwd()    
            
            print 'Object', ' ,'.join(kkk)
            galaxy_dir = os.path.abspath(deca_setup.res_directory+"/%s_%s/%s/%s_%s" % (str(int(data["NUMBER"][kkk[0]])), data["NAME"][kkk[0]], fil[kkk[0]], deca_setup.code.lower(), last_model[kkk[0]]))
            os.chdir(galaxy_dir)
                    
            if not os.path.exists('./tex'):
                        os.mkdir('./tex')
            os.chdir('./tex')

            object_name = '%s,%s,%s' % (str(int(data["NUMBER"][kkk[0]])), data["NAME"][kkk[0]], data["FILTER"][kkk[0]])
            output_file = galaxy_dir + '/'+ deca_setup.code.lower() + '.01'
            pictures = glob.glob(galaxy_dir + '/*%s' % (deca_setup.pic_extension)) 


            if '.reg' not in data["NAME"][k]:
                luminosities = find_luminosities(galaxy_dir + '/composed_model.fits', galfitm=galfitm)
                        
                m0s = []; scales = []; Filters = []; Scales = []; DL_Mpcs = []; AEXTs = []; KCORRs = [] 
                for k in kkk:
                    m0,scale = create_final_table.read_head(data["GALAXY_IMAGE"][k])
                    m0s.append(m0)
                    scales.append(scale)

                    try:
                        Filter = Waves[Bands.index(data["FILTER"][k])]
                    except:
                        misc_functions.OUT_WARNING('Filter %s is not found! Set to None!' % (data["FILTER"][k]))
                        Filter = None
                    Filters.append(Filter)                
                

                    try:
                            Scale = data["DISTANCE"][k] * 1000. / 206265.
                            DL_Mpc =  data["DISTANCE"][k]
                    except:
                            try:
                                DA_Mpc,Scale,DL_Mpc = cosmo_calc_NED_func.main(data["REDSHIFT"][k], H0=deca_setup.H0, WM=deca_setup.WM, WV=deca_setup.WV)
                            except:
                                if geom_units=='kpc' or lum_units=='mag' or lum_units=='Lsun' or SB_units=='Lsun/pc2':
                                    misc_functions.OUT_FAIL('Nor DISTANCE nor REDSHIFT columns are found in the input file %s! Exiting...' % (deca_input_file))
                                else:
                                    Scale = None
                                    DL_Mpc = None
                    Scales.append(Scale)
                    DL_Mpcs.append(DL_Mpc)
                    AEXTs.append(data["AEXT"][k])
                    KCORRs.append(data["KCORR"][k])

                    tex_creator.main(object_name,output_file,deca_setup.code,luminosities,pictures,m0s,scales,DL_Mpcs,Scales,Filters,AEXTs,KCORRs,geom_units,lum_units,SB_units)
                    
                    compile_tex()
                    os.chdir(current_dir)
            else:

                    [xc,yc,ellA,ellB,ellPA,Name,Model,Redshift,AExt,KCorr],[xc_star,yc_star] = misc_functions.READ_REGION_OBJECTS(data["NAME"][kkk[0]], data["GALAXY_IMAGE"][kkk[0]])
                    numb_of_objects = len(xc)
                    
                    Model_list = []
                    comp_number = 0
                    for kk in range(len(Model)):
                        Model_tmp = Model[kk].split('+')
                        aa = []
                        for ii in range(len(Model_tmp)):
                            aa.append(comp_number + ii)
                        Model_list.append(aa)
                        comp_number = comp_number + len(Model_tmp)
                    
                     
                    
                    
                    
                    for i in range(numb_of_objects):
                        SScale = []; DDL_Mpc = []; luminosities = []
                        m0s = []; scales = []; Filters = []; Scales = []; DL_Mpcs = []; AEXTs = []; KCORRs = [] 
                        for k in kkk:
                            [xc,yc,ellA,ellB,ellPA,Name,Model,Redshift,AExt,KCorr],[xc_star,yc_star] = misc_functions.READ_REGION_OBJECTS(data["NAME"][k], data["GALAXY_IMAGE"][k])
                            try:
                                    DA_Mpc,Scale,DL_Mpc = cosmo_calc_NED_func.main(Redshift[i], H0=deca_setup.H0, WM=deca_setup.WM, WV=deca_setup.WV)
                            except:
                                    if geom_units=='kpc' or lum_units=='mag' or lum_units=='Lsun' or SB_units=='Lsun/pc2':
                                        misc_functions.OUT_FAIL('Nor DISTANCE nor REDSHIFT columns are found in the input file %s! Exiting...' % (deca_input_file))
                                    else:
                                        Scale = None
                                        DL_Mpc = None

                            SScale.append(Scale)
                            DDL_Mpc.append(DL_Mpc)

                            try:
                                Filter = Waves[Bands.index(data["FILTER"][k])]
                            except:
                                misc_functions.OUT_WARNING('Filter %s is not found! Set to None!' % (data["FILTER"][k]))
                                Filter = None
                            Filters.append(Filter)                               
                            
                            luminosities.append(find_luminosities(galaxy_dir + '/composed_model.fits', components=Model_list[i]), galfitm=galfitm)
                            m0,scale = create_final_table.read_head(data["GALAXY_IMAGE"][k])
                            m0s.append(m0)
                            scales.append(scale)

                            Scales.append(Scale)
                            DL_Mpcs.append(DL_Mpc)
                            AEXTs.append(AExt)
                            KCORRs.append(KCorr)
                                              
                        if i==0:
                            tex_creator.main(Name[i],output_file,deca_setup.code,luminosities,pictures,m0s,scales,DL_Mpcs,Scales,Filters,AEXTs,KCORRs,geom_units,lum_units,SB_units, continue_file=False, last_object == False)
                        else
                            if i==numb_of_objects-1:
                                tex_creator.main(Name[i],output_file,deca_setup.code,luminosities,pictures,m0s,scales,DL_Mpcs,Scales,Filters,AEXTs,KCORRs,geom_units,lum_units,SB_units, continue_file=True, last_object == False)                                 
                            else:
                                tex_creator.main(Name[i],output_file,deca_setup.code,luminosities,pictures,m0s,scales,DL_Mpcs,Scales,Filters,AEXTs,KCORRs,geom_units,lum_units,SB_units, continue_file=True, last_object == True)                    


            compile_tex()
            os.chdir(current_dir)
# -----------------------------------------------------------------
        
if __name__ == '__main__':
    # Arguments:
    parser = argparse.ArgumentParser(description="Script to read DECA results and put them into a pdf file")   
    parser.add_argument("input", help="input file for DECA")
    
    parser.add_argument("models", help="output file with built models")

    parser.add_argument("--geom_units", default='pix',
                    help="units for output geometrical parameters (pix, by default). Can be pix, arcsec, kpc", required=False)   

    parser.add_argument("--lum_units", default='mag',
                    help="units for output luminosities (mag, by default). Can be mag, ADU, Lsun", required=False)   

    parser.add_argument("--SB_units", default='mag/arcsec2',
                    help="units for output Surface Brightnesses (mag/arcsec2, by default). Can be mag/arcsec2, ADU/pix2, Lsun/pc2", required=False)   

    parser.add_argument("--n", default=1,
                    help="number of Python processes (1, by default)", required=False)

    parser.add_argument("--merge", action="store_true", default=False,
                    help="Merge all pdf files for individual table items into one pdf which will be stored in ./tex directory")

    parser.add_argument("--galfitm", action="store_true", default=False,
                        help="if you want to mark items for which GALFITM decomposition has been done")
    
    # Parse arguments
    args = parser.parse_args()
    
    deca_input_file = args.input
    models = args.models
    geom_units = args.geom_units
    lum_units = args.lum_units
    SB_units = args.SB_units
    n_jobs = int(args.n)
    merge = args.merge
    galfitm = args.galfitm
    
    data,units = read_data.main(deca_input_file, delimiter=',', header_line=0, units_line=None, skip_lines = [])
    fil,last_model,model,numbers = np.loadtxt(models, usecols=[2,3,4,5],dtype=str, unpack=True,skiprows=1,delimiter='\t')    

    # Find items to work with
    filter_indices = []
    for k in range(len(last_model)):
        if last_model[k]!='None':
            filter_ind = []
            for i in range(len(last_model)):
                if NUMBER[k]==NUMBER[i] and FILTER[i] in fil[k]:
                    filter_ind.append(i)
            filter_indices.append(filter_ind)    
    
    Parallel(n_jobs=n_jobs)(delayed(do_each_galaxy)(i,data,fil,last_model,model,numbers,geom_units,lum_units,SB_units,galfitm) for i in filter_indices)
    
    if merge:
        current_dir = os.getcwd()
        res_directory = os.path.abspath(deca_setup.res_directory) 
        if not os.path.exists('./tex'):
                    os.mkdir('./tex')
        os.chdir('./tex')
        ff = open('results.tex', 'w')
        print >>ff, tex_creator.header()
        
        lines_table = []
        for k in range(len(data["NUMBER"])):
            if last_model[k]!='None':
                #data["FILTER"][k] = 'gri' ####TODO:REMOVE
                tex_dir = os.path.abspath(res_directory+"/%s_%s/%s/%s_%s/tex" % (str(int(data["NUMBER"][k])), data["NAME"][k], data["FILTER"][k], deca_setup.code.lower(),last_model[k]))
                fff = open(tex_dir + '/results.tex', 'r')
                lines = fff.readlines()
                for i in range(len(lines)):
                    if '\\section{' in lines[i]:
                        start_line = i
                    if '\\end{center}' in lines[i]:
                        end_line = i
                        break
                for i in range(start_line, end_line, 1):
                    lines_table.append(lines[i])
                
                lines_table.append('\n')
                
                fff.close()
                
        for k in range(len(lines_table)):
            print >>ff, lines_table[k]
    
        end_line = tex_creator.end()
        print >>ff, end_line                
        ff.close()        
        compile_tex()
        os.chdir(current_dir)
        print 'Done!'
        
