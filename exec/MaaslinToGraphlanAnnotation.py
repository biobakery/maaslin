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
import math
from operator import itemgetter
import re
import string
import sys

#def funcGetColor(fNumeric,fMax):
#  if fNumeric>0:
#    return("#"+str(int(99*fNumeric/fMax)).zfill(2)+"0000")
#  if fNumeric<0:
#    return("#00"+str(int(99*abs(fNumeric/fMax))).zfill(2)+"00")
#  return("#000000")

def funcGetColor(fNumeric):
  if fNumeric>0:
    return sRingPositiveColor
  else:
    return sRingNegativeColor

def funcGetAlpha(fNumeric,fMax):
  return max(abs(fNumeric/fMax),dMinAlpha)

#Constants
sAnnotation = "annotation"
sAnnotationColor = "annotation_background_color"
sClass = "class"
sRingAlpha = "ring_alpha"
dMinAlpha = .075
sRingColor = "ring_color"
sRingHeight = "ring_height"
#sRingHeightMin = 0.5
sStandardizedRingHeight = "1.01"
sRingLabel = "ring_label"
sRingLabelSizeWord = "ring_label_font_size"
sRingLabelSize = 10
sRingLineColor = "#999999"
sRingPositiveWord = "Positive"
sRingPositiveColor = "#990000"
sRingNegativeWord = "Negative"
sRingNegativeColor = "#009900"
sRingLineColorWord = "ring_separator_color"
sRingLineThickness = "0.5"
sRingLineThicknessWord = "ring_internal_separator_thickness"
sCladeMarkerColor = "clade_marker_color"
sCladeMarkerSize = "clade_marker_size"
sHighlightedMarkerSize = "10"
c_dMinDoubleValue = 0.00000000001

#Set up arguments reader
argp = argparse.ArgumentParser( prog = "MaaslinToGraphlanAnnotation.py",
    description = """Converts summary files to graphlan annotation files.""" )

#### Read in information
#Arguments
argp.add_argument("strInputSummary", metavar = "SummaryFile", type = argparse.FileType("r"), help ="Input summary file produced by maaslin")
argp.add_argument("strInputCore", metavar = "CoreFile", type = argparse.FileType("r"), help ="Core file produced by Graphlan from the maaslin pcl")
argp.add_argument("strInputHeader", metavar = "HeaderFile", type = argparse.FileType("r"), help ="Input header file to append to the generated annotation file.")
argp.add_argument("strOutputAnnotation", metavar = "AnnotationFile", type = argparse.FileType("w"), help ="Output annotation file for graphlan")

args = argp.parse_args( )

#Read in the summary file and transform to class based descriptions
csvSum = open(args.strInputSummary,'r') if isinstance(args.strInputSummary, str) else args.strInputSummary
fSum = csv.reader(csvSum, delimiter="\t")
#Skip header (until i do this a better way)
fSum.next()

#Extract associations (Metadata,taxon,coef,qvalue)
lsAssociations = [[sLine[1],sLine[2],sLine[4],sLine[7]] for sLine in fSum]
csvSum.close()

#### Read in default graphlan settings provided by maaslin
#Read in the annotation header file
csvHdr = open(args.strInputHeader,'r') if isinstance(args.strInputHeader, str) else args.strInputHeader
fHdr = csv.reader(csvHdr, delimiter="\t")

#Begin writting the output
#Output annotation file
csvAnn = open(args.strOutputAnnotation,'w') if isinstance(args.strOutputAnnotation, str) else args.strOutputAnnotation
fAnn = csv.writer(csvAnn, delimiter="\t")
fAnn.writerows(fHdr)
csvHdr.close()

#If no associatiosn were found
if(len(lsAssociations)==0):
  csvAnn.close()

else:
  #### Fix name formats
  #Manipulate names to graphlan complient names (clades seperated by .)
  lsAssociations = sorted(lsAssociations, key=itemgetter(1))
  lsAssociations = [[sBug[0]]+[re.sub("^[A-Za-z]__","",sBug[1])]+sBug[2:] for sBug in lsAssociations]
  lsAssociations = [[sBug[0]]+[re.sub("\|*[A-Za-z]__|\|",".",sBug[1])]+sBug[2:] for sBug in lsAssociations]

  #If this is an OTU, append the number and the genus level together for a more descriptive termal name
  lsAssociationsModForOTU = []
  for sBug in lsAssociations:
    lsBug = sBug[1].split(".")
    if(len(lsBug))> 1:
      if(lsBug[-1].isdigit()):
        lsBug[-2]=lsBug[-2]+"_"+lsBug[-1]
        lsBug = lsBug[0:-1]
      lsAssociationsModForOTU.append([sBug[0]]+[".".join(lsBug)]+sBug[2:])
    else:
      lsAssociationsModForOTU.append([sBug[0]]+[lsBug[0]]+sBug[2:])

  #Extract just class info
  #lsClassData = [[sLine[2],sClass,sLine[1]] for sLine in fSum]

  ### Make rings
  #Setup rings
  dictRings = dict([[enumData[1],enumData[0]] for enumData in enumerate(set([lsData[0] for lsData in lsAssociationsModForOTU]))])

  #Ring graphlan setting: rings represent a metadata that associates with a feature
  #Rings have a line to help differetiate them
  lsRingSettings = [[sRingLabel,lsPair[1],lsPair[0]] for lsPair in dictRings.items()]
  lsRingLineColors = [[sRingLineColorWord,lsPair[1],sRingLineColor] for lsPair in dictRings.items()]
  lsRingLineThick = [[sRingLineThicknessWord,lsPair[1],sRingLineThickness] for lsPair in dictRings.items()]
  lsRingLineLabelSize = [[sRingLabelSizeWord,lsPair[1], sRingLabelSize] for lsPair in dictRings.items()]

  #Create coloring for rings color represents the directionality of the relationship
  dMaxCoef = max([abs(float(sAssociation[2])) for sAssociation in lsAssociationsModForOTU])
  lsRingColors = [[lsAssociation[1], sRingColor, dictRings[lsAssociation[0]], funcGetColor(float(lsAssociation[2]))] for lsAssociation in lsAssociationsModForOTU]
  lsRingAlpha = [[lsAssociation[1], sRingAlpha, dictRings[lsAssociation[0]], funcGetAlpha(float(lsAssociation[2]), dMaxCoef)] for lsAssociation in lsAssociationsModForOTU]

  #Create height for rings representing the log tranformed q-value?
  dMaxQValue = max([-1*math.log(max(float(sAssociation[3]), c_dMinDoubleValue)) for sAssociation in lsAssociationsModForOTU])
  #lsRingHeights = [[lsAssociation[1], sRingHeight, dictRings[lsAssociation[0]], ((-1*math.log(max(float(lsAssociation[3]), c_dMinDoubleValue)))/dMaxQValue)+sRingHeightMin] for lsAssociation in lsAssociationsModForOTU]
  lsRingHeights = [[lsAssociation[1], sRingHeight, dictRings[lsAssociation[0]], sStandardizedRingHeight] for lsAssociation in lsAssociationsModForOTU]

  #### Marker
  # Marker colors (mainly to make legend
  lsMarkerColors = [[lsAssociation[1], sCladeMarkerColor, funcGetColor(float(lsAssociation[2]))] for lsAssociation in lsAssociationsModForOTU]
  lsMarkerSizes = [[lsAssociation[1], sCladeMarkerSize, sHighlightedMarkerSize] for lsAssociation in lsAssociationsModForOTU]

  #### Make internal highlights
  #Highlight the associated clades
  lsUniqueAssociatedTaxa = sorted(list(set([lsAssociation[1] for lsAssociation in lsAssociationsModForOTU])))

  lsHighlights = []
  sABCPrefix = ""
  sListABC = string.ascii_lowercase
  iListABCIndex = 0
  for lsHighlight in lsUniqueAssociatedTaxa:
    lsTaxa = lsHighlight.split(".")
    sLabel = sABCPrefix+sListABC[iListABCIndex]+":"+lsTaxa[-1] if len(lsTaxa) > 2 else lsTaxa[-1]
    lsHighlights.append([lsHighlight, sAnnotation, sLabel])
    iListABCIndex = iListABCIndex + 1
    if iListABCIndex > 25:
      iListABCIndex = 0
      sABCPrefix = sABCPrefix + sListABC[len(sABCPrefix)]

  #Read in the core file
  csvCore = open(args.strInputCore,'r') if isinstance(args.strInputCore, str) else args.strInputCore
  fSum = csv.reader(csvCore, delimiter="\t")

  #Add in all phylum just incase they were not already included here
  lsAddSecondLevel = list(set([sUnique[0].split(".")[1] for sUnique in fSum if len(sUnique[0].split(".")) > 1]))
  lsHighlights.extend([[sSecondLevel, sAnnotation, sSecondLevel] for sSecondLevel in lsAddSecondLevel])
  lsHighlightColor = [[lsHighlight[0], sAnnotationColor,"b"] for lsHighlight in lsHighlights]

  #### Write the remaining output annotation file
  fAnn.writerows(lsRingSettings)
  fAnn.writerows(lsRingLineColors)
  fAnn.writerows(lsRingColors)
  fAnn.writerows(lsRingAlpha)
  fAnn.writerows(lsRingLineThick)
  fAnn.writerows(lsRingLineLabelSize)
  fAnn.writerows(lsRingHeights)
  fAnn.writerows(lsMarkerColors)
  fAnn.writerows(lsMarkerSizes)
  fAnn.writerows([[sRingPositiveWord, sCladeMarkerColor, sRingPositiveColor]])
  fAnn.writerows([[sRingNegativeWord, sCladeMarkerColor, sRingNegativeColor]])
  fAnn.writerows(lsHighlights)
  fAnn.writerows(lsHighlightColor)
  csvAnn.close()
