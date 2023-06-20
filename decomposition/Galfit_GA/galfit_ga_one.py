#! /usr/bin/env python
"""
Genetic algorithms applied to GALFIT optimization problem
"""

from random import random
from math import sqrt

from libs.pygene.gene import FloatGene, FloatGeneMax, FloatGeneRandom
from libs.pygene.organism import Organism, MendelOrganism
from libs.pygene.population import Population

import time
import sys
import os
import pyfits
import subprocess

import read_galfit
import shutil

FNULL = open(os.devnull, 'w')

#*** Colour fonts ***
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''


# MAIN PARAMETERS FOR GENETIC ALGORITHMS:
#*****************************************
#geneRandMin = 0.0
#geneRandMax = 10.0
geneMutProb = 0.1
geneMutAmt = .2         # only if not using FloatGeneRandom

'''
popInitSize = 1000
popChildCull = 250
popChildCount = 250
popIncest = 700           # number of best parents to add to children
popNumMutants = 0.1     # proportion of mutants each generation
popNumRandomOrganisms = 50  # number of random organisms per generation
'''

'''
Generations = 100
popInitSize = 200
popChildCull = 50
popChildCount = 190
popIncest = 0           # number of best parents to add to children
popNumMutants = 0.1     # proportion of mutants each generation
popNumRandomOrganisms = 10  # number of random organisms per generation
'''
#test
NumberOfEqualChi2 = 10
Generations = 100
popInitSize = 200
popChildCull = 50
popChildCount = 160
popIncest = 20           # number of best parents to add to children
popNumMutants = 0.1     # proportion of mutants each generation
popNumRandomOrganisms = 20  # number of random organisms per generation


mutateOneOnly = False

BaseGeneClass = FloatGene
BaseGeneClass = FloatGeneMax
#BaseGeneClass = FloatGeneRandom

OrganismClass = MendelOrganism
#OrganismClass = Organism

mutateAfterMating = True

crossoverRate = 0.65
#*****************************************

def move_file(file,out_dir):
	try:
		shutil.move(file,'./'+out_dir+'/'+file)
	except:
		z = 1

def copy_file(file,out_dir):
	try:
		shutil.copy(file,'./'+out_dir+'/'+file)
	except:
		z = 1
		
def remove_file(file):
	try:
		os.remove(file)
	except:
		z = 1
		
		
# Timer
start = time.time()
print bcolors.OKBLUE+'\n\n************ GENETIC ALGORITHMS WITH GALFIT ************' + bcolors.ENDC
inputGalfFile = str(sys.argv[1])

o_flag = False
out_flag = False
out_dir = 'none'
k = -1
for a in sys.argv:
    k = k + 1
    if a == '--galf':
        o_flag = True
    if a == '-o':	
        out_flag = True
	out_dir = sys.argv[k+1]
	if not os.path.exists(out_dir):
	    os.makedirs(out_dir)

'''
# List all files in the directory:
files = []
for file in os.listdir("./"):
	files.append(file)
'''

# Create output 
f_gen = open('generations_res.txt', "w")
f_gen.close()

if out_flag == True:
    move_file('generations_res.txt',out_dir)
    
# Read in the input files from the galfit file
image_file,model_file,sigma_file,psf_file,mask_file = read_galfit.read_input_files(inputGalfFile)
print '\nTHE INPUT IMAGES ARE:'
print '\tGalaxy image: %s\n\tSigma image: %s\n\tPsf image: %s\n\tMask image: %s\n' % (image_file,sigma_file,psf_file,mask_file)

# Copy files to the out_dir if specified:
if out_flag == True:
    for file in [image_file,model_file,sigma_file,psf_file,mask_file,inputGalfFile]:
        copy_file(file,out_dir)


# Change the directory to out_dir:
os.chdir(out_dir)

    
# Just to find free parameters and remove appeared file model_file:
galfPars,number_of_pars = read_galfit.main('no',inputGalfFile,model_file)
#print galfPars
os.remove('model.txt')

genome = {}
for key in sorted(galfPars):
	minim = galfPars[key][1]
	maxim = galfPars[key][2]
	class GalfGene(BaseGeneClass):     
	  randMin = minim
	  randMax = maxim

	  mutProb = geneMutProb
	  mutAmt = geneMutAmt
    
	genome[key] = GalfGene
	

# Create header for the output file
f_gen = open('generations_res.txt', "a")
list_of_pars = ''
for key in sorted(galfPars):
  list_of_pars = list_of_pars + '\t' + key

print >>f_gen, "%s%s\t%s" % ('#',list_of_pars,'Chi2')  
f_gen.close()



class TSPSolution(OrganismClass):
	genome = genome

	mutateOneOnly = mutateOneOnly

	crossoverRate = crossoverRate

	numMutants = popNumMutants

	def fitness(self):
	    new_par = {}
	    for key in genome:
	      new_par[key] = self[key]
	    read_galfit.input_new('ini',inputGalfFile,'model.txt',new_par,'non-free')
	    subprocess.call("galfit model.txt", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
	    hdulist = pyfits.open(model_file)
	    prihdr = hdulist[2].header
	    Chi2 = float(prihdr['CHI2NU'])
	    return Chi2



class GalfitPopulation(Population):
        #print 'here'
	initPopulation = popInitSize
	species = TSPSolution

	# cull to this many children after each generation
	childCull = popChildCull

	# number of children to create after each generation
	childCount = popChildCount

	# number of best parents to add in with next gen
	incest = popIncest

	mutants = popNumMutants

	numNewOrganisms = popNumRandomOrganisms

	mutateAfterMating = mutateAfterMating
	
	#par,number_of_pars = read_galf.main('ini',GalfitFile,model_file)
	#subprocess.call("galfit " + str(model_file), shell=True)
        #print 'here'    

    
# create initial population
pop = GalfitPopulation()
#pop.gen()

CHI2 = []

try:
        generations = 0
        while generations<Generations:
            # execute a generation
            pop.gen()
            generations += 1
            best = pop.organisms[0]
            print "Generation %s: Best Chi2=%s" % (generations, pop.best().get_fitness())
            f_gen = open('generations_res.txt', "a")

	    best_par = {}
	    pars_out = ''
	    for key in sorted(genome):
	      best_par[key] = best[key]
	      if pars_out =='':
		pars_out = str(best[key])
	      else:
		pars_out = pars_out + '\t'+ str(best[key])
	    print >>f_gen, "%s\t%s\t%s" % (generations,pars_out,pop.best().get_fitness())    
	    read_galfit.input_new('ini',inputGalfFile,'model.txt',best_par,'non-free')
	    subprocess.call("galfit model.txt", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
	    shutil.move(model_file,'best_model.fits')
            f_gen.close()
            last_chi2 = pop.best().get_fitness()
            CHI2.append(last_chi2)
            
	    indices = [i for i, x in enumerate(CHI2) if x == last_chi2]
	    if len(indices)>=NumberOfEqualChi2:
		print 'All last %i generations had Chi2=%f. Terminated according to the convergence condition!' % (NumberOfEqualChi2,last_chi2)
                break


except KeyboardInterrupt:
        pass
print "Executed", generations, "generations."


# get the best solution
solution = pop.best()  
    
end = time.time()
# Output the total time of the fitting:
duration = end - start # in seconds
hms = time.strftime('%H:%M:%S', time.gmtime(duration))
time_print = "\nFinished fit in %.1f s (%s)\n" % (float(end - start),hms)
print bcolors.OKGREEN+time_print+bcolors.ENDC

if o_flag==True:
  print "Galfit running..."
  read_galfit.input_new('ini',inputGalfFile,'best_model.txt',best_par,'free')
  subprocess.call("galfit best_model.txt", shell=True, stdout=FNULL, stderr=subprocess.STDOUT)
  try:
    hdulist = pyfits.open(model_file)
    prihdr = hdulist[2].header
    Chi2 = float(prihdr['CHI2NU'])
    print "Galfit output file is: %s" % (model_file)
    print "Galfit Chi2 is: %f\n" % (Chi2)
    print "Genetic Algorithms output file is: %s" % ('best_model.fits')
    print "Genetic Algorithms Chi2 is: %f" % (last_chi2)
  except:
    print bcolors.FAIL+'GALFIT failed'+bcolors.ENDC
    print "Genetic Algorithms output file is: %s" % ('best_model.fits')
    print "Genetic Algorithms Chi2 is: %f" % (last_chi2)
    remove_file('best_model.txt')

# Remove unnecessary files:
remove_file('model.txt')
remove_file('fit.log')
if out_flag == True:
    for file in [image_file,model_file,sigma_file,psf_file,mask_file,inputGalfFile]:
        remove_file(file)



 
# Back to the initial directory: 
os.chdir("..")

'''
# Copying new created files to the directory if specified:
if out_flag==True and out_dir!='none':
	files_new = []
	for file in os.listdir("./"):
		files_new.append(file)

	crea_files = list(set(files_new) - set(files))

	for file in crea_files:
		shutil.move(file,'./'+out_dir+'/'+file)
	if (crea_files==[]) or ('galfit' in crea_files):
		move_file('best_model.fits',out_dir)
		move_file('generations_res.txt',out_dir)
		move_file(model_file,out_dir)

'''
    
