####################################
# Summary: Takes the individual detail data from different factors
# And condenses them to one file for reading
# Author: Timothy Tickle
# Start Date: 12-23-2011
####################################

library(logging)

c_logrMaaslin	<- getLogger( "maaslin" )

#Transposes a data file
#Assumes the file data is rectangular/square
funcSummarizeDirectory = function(astrOutputDirectory, strBaseName, astrSummaryFileName, astrKeyword, afSignificanceLevel, acharDelimiter="\t")
{
  c_logrMaaslin$debug("funcSummarizeDirectory")
  c_logrMaaslin$debug(format(astrOutputDirectory))
  c_logrMaaslin$debug("astrSummaryFileName")
  c_logrMaaslin$debug(format(astrSummaryFileName))
  c_logrMaaslin$debug("afSignificanceLevel")
  c_logrMaaslin$debug(format(afSignificanceLevel))

  filePrefix = paste(strBaseName,"-",sep="")
  aiPrefixLength = nchar(filePrefix)
  c_DETAIL_FILE_SUFFIX = ".txt"
  aiSuffixLength = nchar(c_DETAIL_FILE_SUFFIX)

  #Store significant data elements
  vData = c()

  #Store detail file paths
  astrlsDetailFiles = c()

  #Get files in output directory
  astrlsDirectoryFiles = list.files(astrOutputDirectory, full.names=FALSE)
  c_logrMaaslin$debug(format(astrlsDirectoryFiles))

  #Reduce to detail files
  for(astrLowestPathLevel in astrlsDirectoryFiles)
  {
#    astrLowestPathLevel = strsplit(astrRandomFilePath,"/")[[1]]
#    astrLowestPathLevel = astrLowestPathLevel[length(astrLowestPathLevel)]
    aiLowestPathLevelLength = nchar(astrLowestPathLevel)
    if(aiLowestPathLevelLength>aiPrefixLength)
    {
      if((substr(astrLowestPathLevel,1,aiPrefixLength)==filePrefix)&&
         (substr(astrLowestPathLevel,aiLowestPathLevelLength-aiSuffixLength+1,aiLowestPathLevelLength)==c_DETAIL_FILE_SUFFIX))
      {
        astrlsDetailFiles = c(astrlsDetailFiles,paste(astrOutputDirectory,"/",astrLowestPathLevel,sep=""))
      }
    }
  }

  #Evaluate header on first file
  #Will hold the position of the keyword to use to evaluate significance
  aiKeywordPosition = -1

  #Read in file line by line into a list
  astrHeader = "No data files found to combine."
  if(!is.null(astrlsDetailFiles[1]))
  {
    astrlsFileContents = readLines(astrlsDetailFiles[1])
    
    #Holds the initial header to make sure all data file are consistent in format
    astrHeader = astrlsFileContents[1]

    #Elements of the header to search through
    astrlsHeaderElements = strsplit(astrHeader,acharDelimiter)[[1]]
    #Search through header for keyword element position
    aiCurrentHeaderElementPosition = 0
    for(astrHeaderElement in astrlsHeaderElements)
    {
      aiCurrentHeaderElementPosition = aiCurrentHeaderElementPosition + 1
      if(astrHeaderElement == astrKeyword)
      {
        aiKeywordPosition = aiCurrentHeaderElementPosition
        break
      }
    }

    #Check to see if the keyword was found, if not, skip file
    if(aiKeywordPosition == -1)
    {
      print(paste("The following file did not have the keyword in the header, this file was skipped in the summary file. File:",astrFile," Keyword:",astrKeyword," Header:",astrHeader,sep=""))
      return(-1)
    }

    #Read through file and store lines that pass significance
    for(astrLine in astrlsFileContents[2:length(astrlsFileContents)])
    {
      astrlsHeaderElements = strsplit(astrLine,acharDelimiter)[[1]]
      if(astrlsHeaderElements[aiCurrentHeaderElementPosition]<=afSignificanceLevel)
      {
        vData = c(vData,astrLine)
      }
    }

    #For each file after the first file
    for(astrFile in astrlsDetailFiles[2:length(astrlsDetailFiles)])
    {
      #Read in file line by line into a list
      astrlsFileContents = readLines(astrFile)
      #If file has no data entries, pass
      if(length(astrlsFileContents)<2)
      {
        print(paste("The following file did not have data lines, this file was skipped in the summary file. File:",astrFile,sep=""))
        pass
      }
      #If file has an inconsistent header, pass
      if(astrlsFileContents[1] != astrHeader)
      {
        print(paste("The following file had an inconsistent header, this file was skipped in the summary file. File:",astrFile,sep=""))
        pass
      }

      #Read through file and store lines that pass significance
      for(astrLine in astrlsFileContents[2:length(astrlsFileContents)])
      {
        astrlsHeaderElements = strsplit(astrLine,acharDelimiter)[[1]]
        if(as.numeric(astrlsHeaderElements[aiCurrentHeaderElementPosition])<=afSignificanceLevel)
        {
          vData = c(vData,astrLine)
        }
      }
    }
  }

  #Sort data and add header for readability
  vData = sort(vData)
  vData = c(astrHeader,vData)

  #Write data vector
  aiDataLength = length(vData)
  write(x=paste(vData,collapse="\n"), file=astrSummaryFileName, append=FALSE, sep="")
}
