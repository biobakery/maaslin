#!/usr/bin/env python
#####################################################################################
#Copyright (C) <2012>
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in the
#Software without restriction, including without limitation the rights to use, copy,
#modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
#and to permit persons to whom the Software is furnished to do so, subject to
#the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies
#or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# This file is a component of the MaAsLin (Multivariate Associations Using Linear Models), 
# authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Timothy Tickle, ttickle@hsph.harvard.edu).
#####################################################################################

__author__ = "Timothy Tickle"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle"]
__license__ = ""
__version__ = ""
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@sph.harvard.edu"
__status__ = "Development"

import argparse
import csv
from operator import itemgetter
import re
import sys

#Helper function which returns a boolean indicator of an input string being parsable as an int
def funcIsInt(strInt):
  try:
    int(strInt)
    return True
  except:
    return False

#Helper function that gets the index of the name and gives the last value of the list for - or the first value depending on the position
# This supports the ranging in the read.config files
#If no range is given then the result is just one index of the given name
def funcGetIndices(lsFeature, lsFunctionNames):
  if(len(lsFeature)) == 1:
    return lsFeatureNames.index(lsFeature)
  if(len(lsFeature)) == 2:
    iIndices = []
    iPosition = 1
    for sFeature in lsFeature:
      if(sFeature==""):
        if(iPosition==1):
          iIndices.append(2)
        elif(iPosition==2):
          iIndices.append(len(lsFunctionNames)-1)
      else:
        iIndices.append(lsFeatureNames.index(sFeature))
      iPosition = iPosition + 1
    return iIndices

#Constants
#The line indicating the rows to read
c_MatrixName = "Matrix:"
c_DataMatrix = "Abundance"
c_strRows = "Read_PCL_Rows:"

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "PCLToGraphlanCoreGene.py",
    description = """Converts PCL files to Graphlan core gene files.""" )

#Arguments
argp.add_argument("strInputPCL", metavar = "PCLFile", type = argparse.FileType("r"), help ="Input PCl file used in maaslin")
argp.add_argument("strInputRC", metavar = "RCFile", type = argparse.FileType("r"), help ="Input read config file used in maaslin")
argp.add_argument("strOutputCoreGene", metavar = "CoreGeneFile", type = argparse.FileType("w"), help ="Output core gene file for graphlan")

args = argp.parse_args( )

#Read in read config table and get the rows/columns to use
#Indicates if we are reading a data matrix
fIsData = False
#Holds the indices ranges
#List of lists,each internal list hold 1 or 2 indices, if two it indicates a range from the first to the second
llsIndices = []
csvRC = open(args.strInputRC,'r') if isinstance(args.strInputRC, str) else args.strInputRC
fRC = csv.reader(csvRC, delimiter=" ")
for sLine in fRC:
  #Get the row indexs or names
  if len(sLine):
    if sLine[0] == c_MatrixName:
      fIsData = sLine[1] == c_DataMatrix
    if sLine[0] == c_strRows:
      if fIsData:
        llsIndices = [sIndexRange.split("-") for sIndexRange in sLine[1].split(",")]
        break
csvRC.close()

#Read in the PCL file and parse the file names to core genes format
csvPCL = open(args.strInputPCL,'r') if isinstance(args.strInputPCL, str) else args.strInputPCL
fPCL = csv.reader(csvPCL,delimiter="\t")
#The first column of the csv file
lsFeatureNames = [sLine[0] for sLine in fPCL]
csvPCL.close()

#If the indices are names switch with numbers otherwise subtract 1 because they are ment for R
liConvertedRangedIndices = [sIndex-1 if funcIsInt(sIndex) else funcGetIndices(sIndex,lsFeatureNames) for sIndex in llsIndices]
llsIndices = None

#If there are any ranges, reduce to lists of indices
liConvertedIndices = []
for lsIndices in liConvertedRangedIndices:
  lsIndices.sort()
  iLenIndices = len(lsIndices)
  if iLenIndices > 2:
    print "Error, received more than 2 indices in a range. Stopped."
    exit()
  liConvertedIndices.extend(lsIndices if iLenIndices == 1 else range(lsIndices[0],lsIndices[1]+1))
liConvertedRangedIndices = None

#Collapse all indices to a set which is then sorted
liConvertedIndices = sorted(list(set(liConvertedIndices)))

#Reduce name of features to just bugs indicated by indices
lsFeatureNames = itemgetter(*liConvertedIndices)(lsFeatureNames)
liConvertedIndices = None

#Change the bug names to the correct formatting (clades seperated by .)
lsFeatureNames = sorted(lsFeatureNames)
lsFeatureNames = [re.sub("^[A-Za-z]__","",sBug) for sBug in lsFeatureNames]
lsFeatureNames = [[re.sub("\|*[A-Za-z]__|\|",".",sBug)] for sBug in lsFeatureNames]

#If this is an OTU, append the number and the genus level together for a more descriptive termal name
lsFeatureNamesModForOTU = []
for sBug in lsFeatureNames:
  lsBug = sBug[0].split(".")
  if(len(lsBug))> 1:
    if(lsBug[-1].isdigit()):
      lsBug[-2]=lsBug[-2]+"_"+lsBug[-1]
      lsBug = lsBug[0:-1]
    lsFeatureNamesModForOTU.append([".".join(lsBug)])
  else:
    lsFeatureNamesModForOTU.append([lsBug[0]])

#Output core gene file
csvCG = open(args.strOutputCoreGene,'w') if isinstance(args.strOutputCoreGene, str) else args.strOutputCoreGene
fCG = csv.writer(csvCG)
fCG.writerows(lsFeatureNamesModForOTU)
csvCG.close()
