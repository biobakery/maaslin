####################################
# Summary: Takes the individual detail data from different factors
# And condenses them to one file for reading
# Author: Timothy Tickle
# Start Date: 12-23-2011 current
####################################

inlinedocs <- function(
##author<< Timothy Tickle <ttickle@hsph.harvard.edu>
) { return( pArgs ) }

#Logging class
suppressMessages(library(logging, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
#source(file.path("input","maaslin","src","Utility.R"))
#source("Utility.R")

# Get logger
c_logrMaaslin	<- getLogger( "maaslin" )

### Summarizes the massline detail files into one file based on significance.
### astrOutputDirectory The output directory to find the MaAsLin results.
### strBasename The prefix string used in maaslin to strart the detail files.
### astrSummaryFileName The summary file's name, should be a path not a file name
### astrKeyword The column name of the data to check significance before adding a detail to the summary
### afSignificanceLevel The value of significance the data must be at or below to be included in the summary (0.0 is most significant; like p-values)
funcSummarizeDirectory = function(astrOutputDirectory, strBaseName, astrSummaryFileName, astrKeyword, afSignificanceLevel)
{
  #Store significant data elements
  dfSignificantData = NULL

  #Get detail files in output directory
  astrlsDetailFiles = list.files(astrOutputDirectory, pattern=paste(strBaseName,"-","[0-9A-Za-z]*",c_DETAIL_FILE_SUFFIX,sep=""), full.names=TRUE)
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
