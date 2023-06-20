#!/usr/bin/env python
# Module to get information on the last models
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **        © Astronomical Observatory, Ghent University         **
# *****************************************************************

# Import standard modules
import os
import numpy as np
import sys

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
def main(deca_input_file, last_models, galfitm=False):
    # Read in input DECA file
    data,units = read_data.main(deca_input_file, delimiter=',', header_line=0, units_line=None, skip_lines = [])
    fil,last_model = np.loadtxt(last_models, usecols=[2,3],dtype=str, unpack=True,skiprows=1,delimiter='\t')

    
    
    FITTING_PROC = data["FITTING_PROC"]
    MODEL = data["MODEL"]
    
    ff = open('final_models.txt', 'w')
    

    
    # Write first line with column names
    ff.write('NUMBER\tNAME\tFILTER\tMODEL_NUMBER\tMODEL\tOUTPUT_COMPONENTS\n')
    
    # Check if a MODEL has been built for each galaxy
    
    for k in range(len(MODEL)):
        if last_model[k]!='None':
            if galfitm==False:
                galaxy_dir = deca_setup.res_directory+"/%s_%s/%s/%s_%s" % (str(int(data["NUMBER"][k])), data["NAME"][k], data["FILTER"][k], deca_setup.code.lower(),last_model[k])
                file = galaxy_dir + '/' + deca_setup.code.lower() + '.01'
            else:
                galaxy_dir = deca_setup.res_directory+"/%s_%s/%s/galfitm_%s" % (str(int(data["NUMBER"][k])), data["NAME"][k], fil[k], last_model[k])
                file = galaxy_dir + '/galfitm.01'
            
            
            
            f = open(file, 'r')
            lines = f.readlines()
            components = []
            if deca_setup.code == 'IMFIT':
                for line in lines:
                    if 'FUNCTION' in line:
                        components.append(line.split()[1])
            elif deca_setup.code == 'GALFIT' or deca_setup.code == 'GALFITM':
                for line in lines:
                    if '#  Component type' in line and '0)' in line:
                        components.append(line.split()[1])
                        
            model = ''
            numbers = ''
            for i in range(len(components)):
                if i == 0:
                    model = components[i] 
                    numbers = str(i+1)
                else:
                    model = model + '+' + components[i] 
                    numbers = numbers + ',' + str(i+1)
            if galfitm==False:        
                ff.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], data["FILTER"][k], last_model[k], model, numbers))
            else:        
                ff.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], fil[k], last_model[k], model, numbers))
        else:
            if galfitm==False:   
                ff.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], data["FILTER"][k], 'None', 'None', 'None'))
            else:   
                ff.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], fil[k], 'None', 'None', 'None'))
    ff.close()
# -----------------------------------------------------------------
        
        
if __name__ == '__main__':
    deca_input_file = str(sys.argv[1]) 
    last_models = str(sys.argv[2])
    galfitm = str(sys.argv[3])
    
    main(deca_input_file, last_models, galfitm)
