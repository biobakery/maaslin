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

#Constants
sClass = "class"

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MaaslinToGraphlanAnnotation.py",
    description = """Converts summary files to graphlan annotation files.""" )

#Arguments
argp.add_argument("strInputSummary", metavar = "SummaryFile", type = argparse.FileType("r"), help ="Input summary file produced by maaslin")
argp.add_argument("strInputHeader", metavar = "HeaderFile", type = argparse.FileType("r"), help ="Input header file to append to the generated annotation file.")
argp.add_argument("strOutputAnnotation", metavar = "AnnotationFile", type = argparse.FileType("w"), help ="Output annotation file for graphlan")

args = argp.parse_args( )

#Read in the summary file and transform to class based descriptions
csvSum = open(args.strInputSummary,'r') if isinstance(args.strInputSummary, str) else args.strInputSummary
fSum = csv.reader(csvSum, delimiter="\t")
#Skip header (until i do this a better way)
fSum.next()

#Extract just class info
lsClassData = [[sLine[2],sClass,sLine[1]] for sLine in fSum]

csvSum.close()

#Manipulate names to graphlan complient names (clades seperated by .)
lsClassData = sorted(lsClassData, key=itemgetter(0))
lsClassData = [[re.sub("^[A-Za-z]__","",sBug[0])]+sBug[1:] for sBug in lsClassData]
lsClassData = [[re.sub("\|*[A-Za-z]__|\|",".",sBug[0])]+sBug[1:] for sBug in lsClassData]

#Read in the annotation header file
csvHdr = open(args.strInputHeader,'r') if isinstance(args.strInputHeader, str) else args.strInputHeader
fHdr = csv.reader(csvHdr, delimiter="\t")

#Output annotation file
csvAnn = open(args.strOutputAnnotation,'w') if isinstance(args.strOutputAnnotation, str) else args.strOutputAnnotation
fAnn = csv.writer(csvAnn, delimiter="\t")
fAnn.writerows(fHdr)
fAnn.writerows(lsClassData)
csvAnn.close()
csvHdr.close()
