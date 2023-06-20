#!/usr/bin/python
# Module to read in a data file
import sys
import collections

# -----------------------------------------------------------------
# MAIN FUNCTION
def main(data_file, delimiter=',', header_line=0, units_line=None, skip_lines = []):
    # Open the file
    f = open(data_file, 'r')
    
    # Read in the file
    lines = f.readlines()
    
    # Determine the header which will be used to give names to the columns
    header = lines[header_line].split(delimiter)
    header[-1] = header[-1].split('\n')[0]
    if units_line!=None:
      unit = lines[units_line].split(delimiter)
    else:
      unit = None

    # Create the dictionary with the names from the header
    data = collections.OrderedDict()
    units = collections.OrderedDict()
    for i in range(len(header)):
      data[header[i]] = []
      if unit!=None:
	units[header[i]] = unit[i]
      else:
	units[header[i]] = None
    
      
    
    # Read in each line except for the header and the skip_lines
    if unit!=None:
      start = 2
    else:
      start = 1
    #print header

    for k in range(header_line+start,len(lines)):
     if k not in skip_lines:
      line = lines[k].split(delimiter)
      
      for i in range(len(line)):
	if line[i]=='':
	  line[i]= float('nan')
	else:
	  try:
	    line[i] = float(line[i])
	  except:
	    line[i] = str(line[i].split('\n')[0])
	data[header[i]].append(line[i])
    f.close()
    return data,units
# -----------------------------------------------------------------
