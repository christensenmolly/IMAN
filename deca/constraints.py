#!/usr/bin/env python
# Module to put constraints on decomposition parameters
# *****************************************************************
# **         DECA -- DECompoistion Analysis of galaxies          **
# **        © Astronomical Observatory, Ghent University         **
# *****************************************************************
import sys
import deca_setup
import math
#tmp_out = sys.stdout


# -----------------------------------------------------------------
# FUNCTION TO CREATE CONSTRAINTS FILE FOR GALFIT
def create_constraint_file(model, scale, log_text):
 constraints = str()
 if deca_setup.set_constraints==True:
  for kk in range(len(model)):  
        if 'ring' in model[kk]:
            log_text = log_text + 'ERROR: Ring model is not found in GALFIT!\n'
            return None, 3, log_text
    
	CONSTRS = deca_setup.CONSTRS
	
	model_split = model[kk].split('+')
	for k in range(len(model_split)):
            # Offset
            if 'pix' in deca_setup.CENTRE_CONSTRS['MAX_OFFSET']:
                offset = int(math.ceil(float(deca_setup.CENTRE_CONSTRS['MAX_OFFSET'].split('pix')[0])))
            elif 'arcsec' in deca_setup.CENTRE_CONSTRS['MAX_OFFSET']:
                offset = int(math.ceil(float(deca_setup.CENTRE_CONSTRS['MAX_OFFSET'].split('arcsec')[0]) / scale))
            else:
                log_text = log_text + 'ERROR: Problem with constraints!\n'
                return None, 4, log_text
            
            constraints = constraints + '    %i	  	   x	     %i  %i\n' % (k+1,-offset, offset)
            constraints = constraints + '    %i	  	   y	     %i  %i\n' % (k+1,-offset, offset)
            
            if 'sersic' in model_split[k]:
                pars = ['mag','re','n','q','pa']
                for i in range(len(pars)):
                    try:
                        par_constr = CONSTRS[model[kk]+':'+model_split[k]][i]
                    except:
                        par_constr = CONSTRS[model[kk]+':'+model_split[k]+str(k+1)][i]
                    if par_constr!=None:
                        min_value = float(par_constr[0])
                        max_value = float(par_constr[1])
                        constraints = constraints + '    %i	  	   %s	     %f  to  %f\n' % (k+1, pars[i], min_value, max_value)

            if model_split[k]=='exp_disc':
                pars = ['mag','h','q','pa']
                for i in range(len(pars)):
                    par_constr = CONSTRS[model[kk]+':'+model_split[k]][i]
                    if par_constr!=None:
                        min_value = float(par_constr[0])
                        max_value = float(par_constr[1])
                        constraints = constraints + '    %i	  	   %s	     %f  to  %f\n' % (k+1, pars[i], min_value, max_value)

            if model_split[k]=='eon_disc':
                pars = ['mag','h','pa']
                for i in range(len(pars)):
                    par_constr = CONSTRS[model[kk]+':'+model_split[k]][i]
                    if par_constr!=None:
                        min_value = float(par_constr[0])
                        max_value = float(par_constr[1])
                        constraints = constraints + '    %i	  	   %s	     %f  to  %f\n' % (k+1, pars[i], min_value, max_value)


            if model_split[k]=='agn':
                pars = ['mag']
                for i in range(len(pars)):
                    par_constr = CONSTRS[model[kk]+':'+model_split[k]][i]
                    if par_constr!=None:
                        min_value = float(par_constr[0])
                        max_value = float(par_constr[1])
                        constraints = constraints + '    %i	  	   %s	     %f  to  %f\n' % (k+1, pars[i], min_value, max_value)


            if model_split[k]=='ring':
                z=1
                '''
                pars = ['mag','sigma_r','R_ring','q','pa']
                for i in range(len(pars)):
                    par_constr = CONSTRS[model[kk]+':'+model_split[k]][i]
                    if par_constr!=None:
                        min_value = float(par_constr[0])
                        max_value = float(par_constr[1])
                        constraints = constraints + '    %i	  	   %s	     %f  to  %f\n' % (k+1, pars[i], min_value, max_value)
                '''

            if model_split[k]=='ferrer':
                pars = ['mag','rs','alpha','beta','q','pa']
                for i in range(len(pars)):
                    par_constr = CONSTRS[model[kk]+':'+model_split[k]][i]
                    if par_constr!=None:
                        min_value = float(par_constr[0])
                        max_value = float(par_constr[1])
                        constraints = constraints + '    %i	  	   %s	     %f  to  %f\n' % (k+1, pars[i], min_value, max_value)


	if deca_setup.CENTRE_CONSTRS['COM_CENTRE'] == 'yes' and len(model_split)>=2:               
            for k in range(len(model_split)):
                constraints = constraints + '    %i-%i	  	   x	     0.0  to  0.0\n' % (k+1,k+2)
                constraints = constraints + '    %i-%i	  	   y	     0.0  to  0.0\n' % (k+1,k+2) 
	
	
  f = open('constraint.txt', 'w')
  print >>f, '%s'% (constraints)
  f.close()
   
  log_text = log_text + 'Constraints file has been successfully created\n'
  return 'constraint.txt', 0, log_text
 else:
  return 'none', 0, log_text
# -----------------------------------------------------------------


# -----------------------------------------------------------------
# FUNCTION TO READ CONSTRAINTS FILE
def read_constraints(constraints_file):
    f = open('constraint.txt', 'r')
    lines = f.readlines()
    comp = []; par = []; constr = []
    for line in lines:
        if line.strip()[0] != '#' and line!='\n':
            s = line.strip().split('#')[0].strip().split()
            comp.append(s[0])
            par.append(s[1])
            constr.append(" ".join(s[2:]))
    return comp,par,constr            
# -----------------------------------------------------------------            
