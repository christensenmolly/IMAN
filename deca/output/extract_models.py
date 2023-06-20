#!/usr/bin/env python
# Module to extract models specified by user 
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **        © Astronomical Observatory, Ghent University         **
# *****************************************************************

# Import standard modules
import os
import numpy as np

# Import DECA modules
import read_data
import collections
import sys

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
        if '#  Component type' in lines[k] and '0)' in lines[k] and deca_setup.code=='GALFIT':
            k_end = k - 2
            
    if k_end==None and deca_setup.code=='IMFIT':
        k_end = len(lines)

    if k_end==None and deca_setup.code=='GALFIT':
        k_end = len(lines)-4
        
    return k_end
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# MAIN FUNCTION
def main(deca_input_file, models, components='sersic_1+sersic_2', galfitm=False):
    # Read in input DECA file
    data,units = read_data.main(deca_input_file, delimiter=',', header_line=0, units_line=None, skip_lines = [])
    fil,last_model,model,numbers = np.loadtxt(models, usecols=[2,3,4,5],dtype=str, unpack=True,skiprows=1,delimiter='\t')

    
    
    FITTING_PROC = data["FITTING_PROC"]
    MODEL = data["MODEL"]
    
    ff = open('extracted_models.txt', 'w')
    

    components = components.split('+')
    

    
    COMPONENTS = collections.OrderedDict()
    for k in range(len(components)):
        COMPONENTS[components[k]] = collections.OrderedDict()
    
    if galfitm==False:
        for k in range(len(MODEL)):
                if k!=0:
                    ff.write('-------------------------------------------------\n')            
                
                if last_model[k]!='None':
                    ff.write('GALAXY: %s\t%s\t%s\t0\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], data["FILTER"][k]))
                    galaxy_dir = deca_setup.res_directory+"/%s_%s/%s/%s_%s" % (str(int(data["NUMBER"][k])), data["NAME"][k], data["FILTER"][k], deca_setup.code.lower(),last_model[k])
                    file = galaxy_dir + '/' + deca_setup.code.lower() + '.01'
                    f = open(file, 'r')
                    lines = f.readlines()
                    TEXT = []
                    if deca_setup.code == 'IMFIT':
                        for i in range(len(lines)):
                            if 'FUNCTION' in lines[i]:
                                i_end = find_end_function(lines, i)
                                text = lines[i:i_end+1]
                                text = map(''.join, text) + '\n'
                                TEXT.append(text)

                                
                                
                    elif deca_setup.code == 'GALFIT':
                        for i in range(len(lines)):
                            if '#  Component type' in lines[i] and '0)' in lines[i]:
                                i_end = find_end_function(lines, i)
                                text = lines[i:i_end+1]
                                text = ''.join(text)
                                TEXT.append(text)                                     
                    
                    Text = ''
                    for number in numbers[k].split(","):
                        Text = Text + TEXT[int(number)-1]
                        
                            
                    ff.write('%s' % (Text))
                else:
                    ff.write('GALAXY: %s\t%s\t%s\t1\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], data["FILTER"][k]))
                ff.write('-------------------------------------------------\n\n')
    else:
        for k in range(len(MODEL)):
                if k!=0:
                    ff.write('-------------------------------------------------\n')            
                
                if last_model[k]!='None':
                    ff.write('GALAXY: %s\t%s\t%s\t0\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], fil[k]))
                    galaxy_dir = deca_setup.res_directory+"/%s_%s/%s/galfitm_%s" % (str(int(data["NUMBER"][k])), data["NAME"][k], fil[k], last_model[k])
                    file = galaxy_dir + '/galfitm.01'
                    f = open(file, 'r')
                    lines = f.readlines()
                    TEXT = []

                    for i in range(len(lines)):
                            if '#  Component type' in lines[i] and '0)' in lines[i]:
                                i_end = find_end_function(lines, i)
                                text = lines[i:i_end+1]
                                text = ''.join(text)
                                TEXT.append(text)                                     
                    
                    Text = ''
                    for number in numbers[k].split(","):
                        Text = Text + TEXT[int(number)-1]
                        
                            
                    ff.write('%s' % (Text))
                else:
                    ff.write('GALAXY: %s\t%s\t%s\t1\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], fil[k]))
                ff.write('-------------------------------------------------\n\n')    
    
    
    ff.close()
# -----------------------------------------------------------------
        
        
if __name__ == '__main__':
    deca_input_file = str(sys.argv[1]) 
    models = str(sys.argv[2])
    galfitm = str(sys.argv[3])
    main(deca_input_file, models, galfitm=galfitm)
