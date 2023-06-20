#!/usr/bin/python
# DESCRIPTION:
# Script to parse SExtractor output catalog by specifying its columns
# EXAMPLE: X_IMAGE,Y_IMAGE,A_IMAGE,B_IMAGE,KRON_RADIUS,THETA_IMAGE,CLASS_STAR = find_sex_column('field.cat', ['X_IMAGE','Y_IMAGE','A_IMAGE','B_IMAGE','KRON_RADIUS','THETA_IMAGE','CLASS_STAR'], float)


import numpy as np



def find_sex_column(sex_catalog, col_names, dtypes=str):
  print('Reading SExtractor catalogue ...')
  f = open(sex_catalog, 'r')

  column_numbers = []
  lines = f.readlines()
  header_lines = [name for name in lines if (name[0] in '#')]

  for col_name in col_names:
    found_column = False
    for line in header_lines:
      if '#' in line and col_name == line.split()[2]:	
        cols = line.split()
        column_numbers.append(int(cols[1])-1)
        found_column = True
    if found_column == False:
      print('Column with the name %s was not found. Exiting ...' % (col_name))
      exit()
  

  
  if isinstance(dtypes, list):
    COLUMNS = []
    columns = np.loadtxt(sex_catalog, usecols=column_numbers, unpack=True, dtype=str)
    
    for k in range(len(columns)):
      COLUMNS.append(columns[k])
      COLUMNS[k] = np.array(COLUMNS[k], dtype=dtypes[k])
    
  else:
    columns = np.loadtxt(sex_catalog, usecols=column_numbers, unpack=True, dtype=dtypes)
    COLUMNS = columns
  f.close()

  return COLUMNS

# Read Sextractor catalogue  
#X_IMAGE,Y_IMAGE,A_IMAGE,B_IMAGE,KRON_RADIUS,THETA_IMAGE,CLASS_STAR = find_sex_column('field.cat', #['X_IMAGE','Y_IMAGE','A_IMAGE','B_IMAGE','KRON_RADIUS','THETA_IMAGE','CLASS_STAR'], float)
#print(X_IMAGE)
