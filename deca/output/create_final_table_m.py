#!/usr/bin/env python
# Module to create final table based on GalfitM decompositions
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **        © Astronomical Observatory, Ghent University         **
# *****************************************************************

# Import standard modules
import os
import numpy as np
import glob
import argparse
import sys
import collections
import pyfits

# Import DECA modules
DECA_PATH = os.path.dirname(__file__)
sys.path.append(DECA_PATH+'/deca_tk_lib')
import deca
import read_data
import galfit_parser
import galfit_input
import imfit_input
import imfit_parser
import misc
import misc_functions

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
# FUNCTION TO FIND WHERE A FUNCTION ENDS
def find_end_function(lines, k_start):
    k_end = None
    for k in range(k_start,len(lines),1):
        if 'FUNCTION' in lines[k] and deca_setup.code=='IMFIT':
            k_end = k -1
            break
        if '#  Component type' in line and '0)' in line and deca_setup.code=='GALFIT':
            k_end = k - 2
            
    if k_end==None and deca_setup.code=='IMFIT':
        k_end = len(lines)

    if k_end==None and deca_setup.code=='GALFIT':
        k_end = len(lines)-4
        
    return k_end
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO READ SOME KEYWORDS IN THE FITS-HEADER
def read_head(galaxy_image):
    hdulist = pyfits.open(galaxy_image)
    prihdr = hdulist[0].header
    
    m0 = misc_functions.READ_KEYWORD(prihdr, 'M0', default_value=20.5, log_txt='')[0]
    
    scale,note = misc_functions.GET_RESOLUTION(galaxy_image)        
    scale = float(scale)

    if 'SCALE' in prihdr:
        scale,log_txt,status_tmp = READ_KEYWORD(prihdr, 'SCALE', default_value=scale, log_txt=log_txt)
    
    return m0,scale    
# -----------------------------------------------------------------    


# -----------------------------------------------------------------
# MAIN FUNCTION
def main(deca_input_file, models, components, geom_units, lum_units, SB_units):
    # Read in input DECA file
    data,units = read_data.main(deca_input_file, delimiter=',', header_line=0, units_line=None, skip_lines = [])
    fil,last_model,model,numbers = np.loadtxt(models, usecols=[2,3,4,5],dtype=str, unpack=True,skiprows=1,delimiter='\t')

    
    
    FITTING_PROC = data["FITTING_PROC"]
    MODEL = data["MODEL"]

    
    ff = open('%s_table_m.dat' % (deca_setup.res_directory), 'w')
    
    components = components.split('+')

    GALFIT_FUNCTIONS = []
    for i in range(len(components)):
                functions,Galfit_functions = galfit_input.read_list_of_funcs(functions=[components[i]])
                if i==0:
                    for k in range(len(Galfit_functions)):
                        Galfit_functions[k] = Galfit_functions[k]+'/1'
                    GALFIT_FUNCTIONS = Galfit_functions    
                else:
                    c = collections.Counter(components[0:i])
                    number = int(c[components[i]])
                    for k in range(len(Galfit_functions)):
                        Galfit_functions[k] = Galfit_functions[k]+'/'+str(number+1)    
                    GALFIT_FUNCTIONS = GALFIT_FUNCTIONS + Galfit_functions
    PARS = collections.OrderedDict()

    for i in range(len(GALFIT_FUNCTIONS)):
        PARS[GALFIT_FUNCTIONS[i]] = []




    CompNames, ComponentFunctions = galfit_parser.read_list_of_funcs()
    geom_pars = galfit_input.geom_pars()
    lum_pars = galfit_input.lum_pars()
    SB_pars = galfit_input.SB_pars()


    # Find items to work with
    filter_indices = []
    for k in range(len(last_model)):
        if last_model[k]!='None':
            filter_ind = []
            for i in range(len(last_model)):
                if NUMBER[k]==NUMBER[i] and FILTER[i] in fil[k]:
                    filter_ind.append(i)
            filter_indices.append(filter_ind)

    Bands = misc_functions.FILTERS_LIST(out='band')
    Waves = misc_functions.FILTERS_LIST(out='wave')
                                        
    for k in range(len(filter_indices)):
            galaxy_dir = deca_setup.res_directory+"/%s_%s/%s/galfitm_%s" % (str(int(data["NUMBER"][filter_indices[k][0]])), data["NAME"][filter_indices[k][0]], fil[filter_indices[k][0]],last_model[filter_indices[k][0]])
            file = galaxy_dir + '/model.galfit.01.band'

            # Read in m0 and scale from the galfitm file:
	    with open(file) as f:
		  lines = f.readlines()
            for line in lines:
                if 'J) ' in line and '# Magnitude photometric zeropoint' in line:
                    m0s = list(np.array(line.split('J) ')[-1].split('#')[0].split(','),float))
                if 'K) ' in line and '# Plate scale (dx dy)   [arcsec per pixel]' in line:
                    scale = float(line.split()[1])                    
            
            for i in filter_indices[k]:
                m0 = m0s[i]
                try:
                    Filter = Waves[Bands.index(data["FILTER"][i])]
                except:
                    misc_functions.OUT_WARNING('Filter %s is not found! Set to None!' % (data["FILTER"][i]))
                    Filter = None

                
                try:
                    Scale = data["DISTANCE"][i] * 1000. / 206265.
                    DL_Mpc =  data["DISTANCE"][i]
                except:
                    try:
                        DA_Mpc,Scale,DL_Mpc = cosmo_calc_NED_func.main(data["REDSHIFT"][i], H0=deca_setup.H0, WM=deca_setup.WM, WV=deca_setup.WV)
                    except:
                        if geom_units=='kpc' or lum_units=='mag' or lum_units=='Lsun' or SB_units=='Lsun/pc2':
                            misc_functions.OUT_FAIL('Nor DISTANCE nor REDSHIFT columns are found in the input file %s! Exiting...' % (deca_input_file))
                        else:
                            Scale = None
                            DL_Mpc = None


                
                
            TMP_PARS, func_names = galfitm_parser.to_deca_format(file, CompNames, ComponentFunctions)


            TMP_PARS_funcs = []
            TMP_PARS_vals = []
            TMP_PARS_units = []
            for i in range(len(TMP_PARS)):
                        TMP_PARS_funcs.append(TMP_PARS[i].split('/')[0]+'/'+TMP_PARS[i].split('/')[1].split(':')[0])
                        TMP_PARS_vals.append(TMP_PARS[i].split(':')[3].split(',')[0])
                        if str(TMP_PARS[i].split(':')[2].split('/')[0]) in geom_pars:
                            TMP_PARS_units.append('geom')                    
                        elif str(TMP_PARS[i].split(':')[2].split('/')[0]) in lum_pars:
                            TMP_PARS_units.append('lum')                      
                        elif str(TMP_PARS[i].split(':')[2].split('/')[0]) in SB_pars:
                            TMP_PARS_units.append('SB')
                        else:
                            TMP_PARS_units.append('')

            for i in range(len(GALFIT_FUNCTIONS)):
                    try:
                        val = TMP_PARS_vals[TMP_PARS_funcs.index(GALFIT_FUNCTIONS[i])]
                        unit_cat = TMP_PARS_units[TMP_PARS_funcs.index(GALFIT_FUNCTIONS[i])]
                        val = val.split()
                        Value = []
                        for v in val:
                            value = misc.units_converter(float(v),unit_cat,'pix','mag/arcsec2','mag',geom_units,SB_units,lum_units,m0,scale,DL_Mpc,Scale,Filter,data["AEXT"][k],data["KCORR"][k])
                            Value.append(value)
                        PARS[GALFIT_FUNCTIONS[i]].append(Value) 
                    except:
                        PARS[GALFIT_FUNCTIONS[i]].append([float('nan')]) 
            
            #else:
            #    for i in range(len(GALFIT_FUNCTIONS)):                        
            #            PARS[GALFIT_FUNCTIONS[i]].append(float('nan')) 




    
    # Header
    ff.write('# Units of the geometrical parameters %s are %s\n' % (geom_pars, geom_units))
    ff.write('# Units of the luminosities %s are %s\n' % (lum_pars, lum_units))
    ff.write('# Units of the surface brightnesses %s are %s\n' % (SB_pars, SB_units))
    
    columns = ["NUMBER"] + ["NAME"] + ["FILTER"] + list(PARS.keys())
    for i in range(len(columns)):
        if i!=len(columns)-1:
            ff.write(columns[i]+'\t')
        else:
            ff.write(columns[i]+'\n')

    
    for k in in range(len(filter_indices)):
        NN = 0
        for i in filter_indices[k]:
            ff.write('%s\t%s\t%s\t' % (str(int(data["NUMBER"][filter_indices[k][0]])), data["NAME"][filter_indices[k][0]],  data["FILTER"][i]) )
            keys = list(PARS.keys())
            for l in range(len(keys)):
                if l!=len(keys)-1:
                        ff.write(str(PARS[keys[l]][k][NN])+'\t')
                else:
                        ff.write(str(PARS[keys[l]][k][NN])+'\n')
            NN = NN + 1
    ff.close()
    print 'Done!'
# -----------------------------------------------------------------        
        
if __name__ == '__main__':
    # Arguments:
    parser = argparse.ArgumentParser(description="Script to read DECA results and put them into a table")   
    parser.add_argument("i",help="input file for DECA")
    
    parser.add_argument("m", help="output file with built models")
    
    parser.add_argument("--c", default='sersic',
                    help="unified model which should be created for each galaxy (sersic, by default)", required=False)

    parser.add_argument("--geom_units", default='pix',
                    help="units for output geometrical parameters (pix, by default). Can be pix, arcsec, kpc", required=False)   

    parser.add_argument("--lum_units", default='mag',
                    help="units for output luminosities (mag, by default). Can be mag, ADU, Lsun", required=False)   

    parser.add_argument("--SB_units", default='mag/arcsec2',
                    help="units for output Surface Brightnesses (mag/arcsec2, by default). Can be mag/arcsec2, ADU/pix2, Lsun/pc2", required=False)   

    # Parse arguments
    args = parser.parse_args()
    
    deca_input_file = args.i
    models = args.m
    components = args.c
    geom_units = args.geom_units
    lum_units = args.lum_units
    SB_units = args.SB_units
    
    main(deca_input_file, models, components=components, geom_units=geom_units, lum_units=lum_units, SB_units=SB_units)
