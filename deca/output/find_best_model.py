#!/usr/bin/env python
# Module to find the best model (with lowest chi2) from the recieved decomposition models
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **        © Astronomical Observatory, Ghent University         **
# *****************************************************************

# Import standard modules
import os
import sys
import os
import glob
from math import *

# Import DECA modules
import read_data

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
# MAIN FUNCTION
def read_chi2(file):
    f = open(file, 'r')
    lines = f.readlines()
    for line in lines:
        if 'Chi^2/nu = ' in line:
            chi2 = float(line.split('Chi^2/nu = ')[-1].split(',')[0])
            break
        if 'Reduced Chi^2 = 'in line:
            chi2 = float(line.split('Reduced Chi^2 = ')[-1].split('\n')[0])
            break            
    f.close()
    return chi2

def find_nearest(array, value):
    n = [abs(i-value) for i in array]
    idx = n.index(min(n))
    return idx

def main(deca_input_file, galfitm=False):
    # Read in input DECA file
    data,units = read_data.main(deca_input_file, delimiter=',', header_line=0, units_line=None, skip_lines = [])
    FITTING_PROC = data["FITTING_PROC"]
    MODEL = data["MODEL"]
    
    ff = open('model_numbers.txt', 'w')
    

    
    # Write first line with column names
    ff.write('NUMBER\tNAME\tFILTER\tBEST_MODEL_NUMBER\n')
    
    # Check if a MODEL has been built for each galaxy
    
    if galfitm==False:
        for k in range(len(MODEL)):
                galaxy_dir = deca_setup.res_directory+"/%s_%s/%s" % (str(int(data["NUMBER"][k])), data["NAME"][k], data["FILTER"][k])

                arr = glob.glob(galaxy_dir + '/%s_*' % (deca_setup.code.lower()))
                dir_numbers = []
                Chi2 = []
                for ar in arr:
                    dir_numbers.append(int(ar.split('_')[-1]))
                    chi2 = read_chi2( ar + '/' + deca_setup.code.lower() + '.01' )
                    Chi2.append(chi2)
                
                #best_model = dir_numbers[Chi2.index(min(Chi2))]
                best_model = dir_numbers[find_nearest(Chi2, 1.)]


                ff.write('%s\t%s\t%s\t%s\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], data["FILTER"][k],best_model))
    else:
        for k in range(len(MODEL)):
                considered_galaxy = list(np.where(data["NUMBER"] == data["NUMBER"][k] && (data["FITTING_PROC"][k]=='DO' || data["FITTING_PROC"][k]=='REDO'))[0])
                done_galaxies = considered_galaxy  ####TODO: This is not necessarilly so!
                filters = data["FILTER"][done_galaxies]
                galaxy_dir = deca_setup.res_directory+'/'+str(data["NUMBER"][k])+'_'+galaxy_name+'/'+'_'.join(filters)
                arr = glob.glob(galaxy_dir + '/galfitm_*')            
            
                dir_numbers = []
                Chi2 = []
                for ar in arr:
                    dir_numbers.append(int(ar.split('_')[-1]))
                    chi2 = read_chi2( ar + '/' + deca_setup.code.lower() + '.01' )
                    Chi2.append(chi2)
                
                #best_model = dir_numbers[Chi2.index(min(Chi2))]
                best_model = dir_numbers[find_nearest(Chi2, 1.)]

                if data["FILTER"][k] in filters:
                    ff.write('%s\t%s\t%s\t%s\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], galaxy_dir.split('/')[-1], best_model))
                else:
                    ff.write('%s\t%s\t%s\t%s\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], 'None', 'None'))

    ff.close()
# -----------------------------------------------------------------
        
        
if __name__ == '__main__':
    deca_input_file = str(sys.argv[1]) 
    galfitm = str(sys.argv[2])

    main(deca_input_file, galfitm=galfitm)
