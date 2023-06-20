#!/usr/bin/env python
# Module to read in last model in the output directory
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **        © Astronomical Observatory, Ghent University         **
# *****************************************************************

# Import standard modules
import os
import sys
import os
import glob
import argparse

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
def main(deca_input_file, galfitm=False):
    # Read in input DECA file
    data,units = read_data.main(deca_input_file, delimiter=',', header_line=0, units_line=None, skip_lines = [])
    FITTING_PROC = data["FITTING_PROC"]
    MODEL = data["MODEL"]
    
    ff = open('model_numbers.txt', 'w')
    

    
    # Write first line with column names
    ff.write('NUMBER\tNAME\tFILTER\tLAST_MODEL_NUMBER\n')
    
    # Check if a MODEL has been built for each galaxy
    
    for k in range(len(MODEL)):
        if galfitm==False:
            galaxy_dir = deca_setup.res_directory+"/%s_%s/%s" % (str(int(data["NUMBER"][k])), data["NAME"][k], data["FILTER"][k])
            arr = glob.glob(galaxy_dir + '/%s_*' % (deca_setup.code.lower()))

            dir_numbers = []
            for ar in arr:
                dir_numbers.append(int(ar.split('_')[-1]))
            if dir_numbers==[]:
                last_model = 'None'
            else:
                last_model = str(max(dir_numbers))

            ff.write('%s\t%s\t%s\t%s\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], data["FILTER"][k],last_model))
        else:
            considered_galaxy = list(np.where(data["NUMBER"] == data["NUMBER"][k] && (data["FITTING_PROC"][k]=='DO' || data["FITTING_PROC"][k]=='REDO'))[0])
            done_galaxies = considered_galaxy  ####TODO: This is not necessarilly so!
            filters = data["FILTER"][done_galaxies]
            galaxy_dir = deca_setup.res_directory+'/'+str(data["NUMBER"][k])+'_'+galaxy_name+'/'+'_'.join(filters)
            arr = glob.glob(galaxy_dir + '/galfitm_*')
            
            dir_numbers = []
            for ar in arr:
                dir_numbers.append(int(ar.split('_')[-1]))
            if dir_numbers==[]:
                last_model = 'None'
            else:
                last_model = str(max(dir_numbers))
                
            if data["FILTER"][k] in filters:
                ff.write('%s\t%s\t%s\t%s\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], galaxy_dir.split('/')[-1], last_model))
            else:
                ff.write('%s\t%s\t%s\t%s\n' % (str(int(data["NUMBER"][k])), data["NAME"][k], 'None', 'None'))
            #k = max(done_galaxies)

    ff.close()
# -----------------------------------------------------------------
        
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Function to read in last model in the output directory.")     
    parser.add_argument("deca_input_file", help="input file for DECA", required=True)
    parser.add_argument("--galfitm", action="store_true", default=False,
                        help="if you want to mark items for which GALFITM decomposition has been done")
    args = parser.parse_args()
    deca_input_file = args.deca_input_file
    galfitm = args.galfitm
    
    main(deca_input_file, galfitm=galfitm)
