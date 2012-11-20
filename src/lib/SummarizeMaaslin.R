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

inlinedocs <- function(
##author<< Curtis Huttenhower <chuttenh@hsph.harvard.edu> and Timothy Tickle <ttickle@hsph.harvard.edu>
##description<< Creates a summary of association detail files.
) { return( pArgs ) }

#Logging class
suppressMessages(library(logging, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))

# Get logger
c_logrMaaslin	<- getLogger( "maaslin" )

funcSummarizeDirectory = function(
### Summarizes the massline detail files into one file based on significance.
astrOutputDirectory,
### The output directory to find the MaAsLin results.
strBaseName,
### The prefix string used in maaslin to start the detail files.
astrSummaryFileName,
### The summary file's name, should be a path not a file name
astrKeyword,
### The column name of the data to check significance before adding a detail to the summary
afSignificanceLevel
### The value of significance the data must be at or below to be included in the summary (0.0 is most significant; like p-values)
){
  #Store significant data elements
  dfSignificantData = NULL

  #Get detail files in output directory
  astrlsDetailFiles = list.files(astrOutputDirectory, pattern=paste(strBaseName,"-","[[:print:]]*",c_sDetailFileSuffix,sep=""), full.names=TRUE)
  logdebug(format(astrlsDetailFiles),c_logrMaaslin)

  #For each file after the first file
  for(astrFile in astrlsDetailFiles)
  {
    #Read in data and reduce to significance
    dfDetails = read.table(astrFile, header=TRUE, sep=c_cTableDelimiter)
    dfDetails = dfDetails[which(dfDetails[astrKeyword] <= afSignificanceLevel),]
    #Combine with other data if it exists
    if(is.null(dfSignificantData))
    {
      dfSignificantData = dfDetails
    } else {
      dfSignificantData = rbind(dfSignificantData,dfDetails)
    }
  }
  #Write data to file
  unlink(astrSummaryFileName)
  if(is.null(dfSignificantData))
  {
    funcWrite("No data files found to combine.",astrSummaryFileName)
  } else {
    funcWriteTable( dfSignificantData, astrSummaryFileName, fAppend = FALSE )
  }
} 
