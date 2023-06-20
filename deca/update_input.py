#!/usr/bin/env python
# -*- coding: utf8 -*-
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **        © Astronomical Observatory, Ghent University         **
# *****************************************************************

# Import standard modules
import os

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


def main(deca_input_file, deca_input_file_new):
    # Read in input DECA file
    data,units = read_data.main(deca_input_file, delimiter=',', header_line=0, units_line=None, skip_lines = [])
    FITTING_PROC = data["FITTING_PROC"]
    MODEL = data["MODEL"]
    
    f = open(deca_input_file, 'r')
    ff = open(deca_input_file_new, 'w')
    
    lines = f.readlines()
    
    # Write first line with column names
    ff.write(lines[0])
    
    # Check if a MODEL has been built for each galaxy
    
    for k in range(len(MODEL)):
        galaxy_dir = deca_setup.res_directory+"/%s_%s/%s" % (data["NUMBER"][k], data["NAME"][k], data["FILTER"][k])
        if MODEL[k]=='sersic':
            if not os.path.exists(galaxy_dir + '/' + deca_setup.code.lower() + '_1/' + deca_setup.code.lower() + '.01'):
                # Add this to the new file
                ff.write(lines[kk+1])

        else:
            if not os.path.exists(galaxy_dir + '/' + deca_setup.code.lower() + '_2/' + deca_setup.code.lower() + '.01'):
                # Add this to the new file
                ff.write(lines[kk+1])

    f.close()
    ff.close()
    

    ff = open(deca_input_file_new, 'r')
    lines = ff.readlines()
    if len(lines)<=1:
        print 'All items in the input table have been decomposed! No new table was created!'
        ff.close()
        os.remove(deca_input_file_new)
    else:
        ff.close()
        
if __name__ == '__main__':
    deca_input_file = str(sys.argv[1]) 
    deca_input_file_new = str(sys.argv[2])

    main(deca_input_file, deca_input_file_new)
