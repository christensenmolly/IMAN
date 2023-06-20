#! /usr/bin/env python

from os import getcwd
from os.path import exists
import sys


def parse_imfit_line(line):
    """ Function parses line of imfit data file and
    returns parameters"""
    params = line.split()
    name = params[0]
    value = float(params[1]) # Value of the parameter is the second entry of the line
    isFixed = False
    # (after the name)
    # Lets now find range of values. Its have to contain the coma and must
    # be the second enrty (otherwise imfin wont work).
    # Some parameters can be fixed, so we have to check this possibility at first
    if (len(params) == 2):
        # No bounds specified at all
        lowerLim = upperLim = None
    elif ("fixed" in params[2]):
        lowerLim = upperLim = None
        isFixed = True
    else:
        rangeParams = params[2].split(",")
        lowerLim = float(rangeParams[0])
        upperLim = float(rangeParams[1])
	    
    return ImfitParameter(name, value, lowerLim, upperLim, isFixed)


class ImfitParameter(object):
    """ Just a container of parameter instance:
    parameter name, value and its range"""
    def __init__(self, name, value, lowerLim, upperLim, fixed):
        self.name = name
        self.value = value
        self.lowerLim = lowerLim
        self.upperLim = upperLim
        self.fixed = fixed

    def tostring(self, fixAll):
        if fixAll or self.fixed:
            return "{:12}{:6.2f}     fixed\n".format(self.name, self.value)
        elif (self.lowerLim is None) and (not self.fixed):
            return "{:12}{:6.2f}\n".format(self.name, self.value)
        else:
            return "{:12}{:6.2f}{:12.2f},{:1.2f}\n".format(self.name, self.value, self.lowerLim, self.upperLim)
    def change_value(self, newValue):
        self.value = newValue
        if not self.fixed:
            if self.lowerLim > newValue:
                self.lowerLim = newValue - newValue*0.001
            if self.upperLim < newValue:
                self.upperLim = newValue + newValue*0.001


class ImfitFunction(object):
    """ Class represents imfit function with
    all its parameters their ranges"""
    def __init__(self, funcName, ident):
        # ident is a unical number of the function
        # funcName is just a type of the function and it can be 
        # the same for different galaxy components
        self.name = funcName
        self.ident = ident
        self.sameNameIdx = None
        self.uname = "%s.%i" % (funcName, ident) # func unique name
        self.params = []
    def add_parameter(self, newParameter):
        self.params.append(newParameter)
    def num_of_params(self):
        return len(self.params)
    def get_par_by_name(self, name):
        for par in self.params:
            if par.name == name:
                return par


class ImfitModel(object):
    """Imfit functions and their parameters"""
    def __init__(self, modelFileName):
        if modelFileName is None:
            return
        print("Reading '%s':" % (modelFileName))
        # Read imfit input file
        self.listOfFunctions = []
        self.allFuncNames = []
        self.numberOfParams = 0
        funcName = None
        ident = -1
        for line in open(modelFileName):
            sLine = line.strip()
            if sLine.startswith("#"):
                # It is a comment line, just skip it
                continue
            if len(sLine) == 0:
                "Empty line"
                continue
            if "#" in sLine:
                # Drop the comment part of the line if exists
                sLine = sLine[:sLine.index("#")].strip()
            if sLine.startswith("X0"):
                x0 = parse_imfit_line(sLine)
            elif sLine.startswith("Y0"):
                y0 = parse_imfit_line(sLine)
            elif sLine.startswith("FUNCTION"):
                # New function is found.
                ident += 1
                # If we are working already with some function, then
                # the list of parameters for this function is over and we can
                # add it to the function list
                if funcName is not None:
                    self.listOfFunctions.append(currentFunction)
                    self.allFuncNames.append(funcName)
                funcName = sLine.split()[1]
                currentFunction = ImfitFunction(funcName, ident)
                currentFunction.add_parameter(x0)
                currentFunction.add_parameter(y0)
                currentFunction.sameNameIdx = self.allFuncNames.count(funcName) + 1
                self.numberOfParams += 2
            else:
                # If line does not contain nor coordinates nor function name
                # then in has to be a parameter line
                param = parse_imfit_line(sLine)
                currentFunction.add_parameter(param)
                self.numberOfParams += 1
        # append the last function
        self.listOfFunctions.append(currentFunction)
        self.allFuncNames.append(funcName)
        # Print some statistics
        print("  %i functions found (%i parameters)\n" % (len(self.listOfFunctions), self.numberOfParams))        

    def get_func_by_uname(self, uname):
        for func in self.listOfFunctions:
            if uname == func.uname:
                return func

    def create_input_file(self, fileName="model.imfit", fixAll=False):
        fout = open(fileName, "w")
        fout.truncate(0)
        for func in self.listOfFunctions:
            fout.write(func.get_par_by_name("X0").tostring(fixAll))
            fout.write(func.get_par_by_name("Y0").tostring(fixAll))
            fout.write("FUNCTION " + func.name+"\n")
            for par in func.params[2:]:
                fout.write(par.tostring(fixAll))
        fout.close()
        # print("Model was saved to '%s'\n" % (fileName))
        return fileName

    def check_boundaries(self, resModel):
        """ Method takes other model object and checks if its parameter
        values are close to the parameter limiths of the self model"""
        badParams = []
        for selfFunc, resFunc in zip(self.listOfFunctions, resModel.listOfFunctions):
            for selfParam, resParam in zip(selfFunc.params, resFunc.params):
                if selfParam.fixed:
                    # Fixed parameter, so it does not have boundaries. Nothing to check.
                    continue
                parRange = selfParam.upperLim - selfParam.lowerLim
                eps = parRange / 1000.0
                if (abs(resParam.value-selfParam.upperLim)<eps) or (abs(resParam.value-selfParam.lowerLim)<eps):
                    badParams.append("%s(%i): %s" % (selfFunc.name, selfFunc.ident, selfParam.name))
        return badParams

    def model_to_text(self, genNumber, ftns, textFile):
        """ Method saves current values of model parameters to a text file """
        if not exists(textFile):
            fout = open(textFile, "w", buffering=1)
            # Create a header as a first line of a file
            fout.write("# genNumber   fintess")
            for func in self.listOfFunctions:
                for param in func.params:
                    fout.write("  %s.%s" % (func.uname, param.name))
            fout.write("\n")
        else:
            fout = open(textFile, "a", buffering=1)
        fout.write("%i   %1.3f" % (genNumber, ftns))
        for func in self.listOfFunctions:
            for param in func.params:
                fout.write("  %9.3f" % param.value)
        fout.write("\n")

    def to_deca_format(self):
        outList = []
        for func in self.listOfFunctions:
            for param in func.params:
                line = "%s:%s:%s/%i:%1.4f," % (func.name, param.name, param.name,
                                               func.sameNameIdx, param.value)
                if (param.lowerLim is None) and (not param.fixed):
                    # No boundaries given
                    line += ","
                elif param.fixed:
                    line += "*,*"
                else:
                    line += "%1.4f,%1.4f" % (param.lowerLim, param.upperLim)
                outList.append(line)
        return outList,self.listOfFunctions
    '''
    def return_comps(self):
        outList = []
        for func in self.listOfFunctions:
	    Component = OrderedDict()
            for param in func.params:
	      Component.update({func.name+':'+param.name+'/'+func.sameNameIdx:[]})
    '''

def parse_deca_line(line):
    isFixed = False
    name = line.split(":")[1]
    valueString = line.split(":")[-1]
    value = float(valueString.split(",")[0])
    lowerLimString = valueString.split(",")[1]
    upperLimString = valueString.split(",")[2]
    if lowerLimString == "":
        # no bounds given
        lowerLim = upperLim = None
    elif "*" in lowerLimString:
        lowerLim = upperLim = None
        isFixed = True
    else:
        lowerLim = float(lowerLimString)
        upperLim = float(upperLimString)
    return ImfitParameter(name, value, lowerLim, upperLim, isFixed)
    
    
def from_deca_format(decaList, fileName):
    """ Function creates imfit file from the given list
    of deca-formatted lines"""

    fout = open(fileName, "w")
    fout.truncate(0)
    
    currentName = None
    allFuncs = {}
    for line in decaList:
        funcName = line.split(":")[0]
        idx = line.split(":")[2].split("/")[1]
        token = funcName + idx
        if token not in allFuncs:
            allFuncs[token] = []
        allFuncs[token].append(line)

    for funcList in allFuncs.values():
        for line in funcList:
            if "X0" in line:
                par = parse_deca_line(line)
                fout.write(par.tostring(fixAll=False))
            if "Y0" in line:
                par = parse_deca_line(line)
                fout.write(par.tostring(fixAll=False))
        fout.write("FUNCTION %s\n" % funcList[0].split(":")[0])
        for line in funcList:
            if ("X0" not in line) and ("Y0" not in line):
                par = parse_deca_line(line)
                fout.write(par.tostring(fixAll=False))


#model = ImfitModel(sys.argv[1])
#decaList =  model.to_deca_format()

#from_deca_format(decaList, "tmp.imfit")
