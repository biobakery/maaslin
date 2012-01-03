#!/usr/bin/env Rscript

####################################
# Summary: Core Code for MaAsLin
# Author: Timothy Tickle
# Start Date: 10-26-2011
####################################

library( optparse )

pArgs <- OptionParser( usage = "%prog [options] <output.txt> <data.tsv> <data.read.config> <data.R> [source.R]*" )
lsArgs <- parse_args( pArgs, positional_arguments = TRUE )
if( length( lsArgs$args ) < 4 ) {
	stop( print_help( pArgs ) ) }
strOutputTXT	<- lsArgs$args[1]
strInputTSV		<- lsArgs$args[2]
strInputRC		<- lsArgs$args[3]
strInputR		<- lsArgs$args[4]
astrSourceR		<- lsArgs$args[5:length( lsArgs$args )]

#Remove any previously stored variables
#rm(list=ls(all=TRUE))

#Constants
c_dMinSamp = 0.1
c_dFreq = 0.01
#Input
c_strMatrixData		<- "Abundance"
c_strMatrixMetadata	<- "Metadata"
# Output
c_iMFA = 30
c_dHeight = 9
# Summary
c_strKeywordEvaluatedForInclusion = "Q-value"
c_strProcessFunction = "processFunction"

#Get command line arguments
inputFile = strInputRC
outputDirectory = paste(dirname(strOutputTXT), "/", sep = "" )
customDataProcessFunction = strInputR

strBase <- sub(".read.config", "", basename(strInputRC))
strSummaryFileName = paste(outputDirectory,strBase,"_Summary.txt",sep="")

print(astrSourceR)
print(inputFile)
print(outputDirectory)
print(customDataProcessFunction)

#Libraries
for( strR in astrSourceR ) {
	source( strR ) }

#Indictate start
print("Start MaAsLin")
print(lsArgs)

funcSourceScript = function(astrFunctionPath)
{
  #If is specified, set up the custom func clean variable
  #Check to see if is NA (or not given)
  if(!is.na(astrFunctionPath))
  {
    #Check to make sure the file exists
    if(file.exists(astrFunctionPath))
    {
      noImpute = NULL
      invert = FALSE
      #Read in the file
      source(astrFunctionPath)
      #If the script exists in the file, run.
      if(exists(c_strProcessFunction,mode="function"))
      {
        print(paste("Preprocessing script is loaded an available. Script Name:",c_strProcessFunction," Script File:",astrFunctionPath,sep=""))
        return( list(
          processFunction   = get( c_strProcessFunction ),
          significanceLevel = significanceLevel,
          invert            = invert,
          noImpute          = noImpute) )
      } else {
        print(paste("MaAsLin Error: Attempted to read in function but was not successful. Function Name: ",c_strProcessFunction,sep=""))
        return(-1)
      }
    #Handle when the file does not exist
    } else {
      print(paste("MaAsLin Error: A custom data manipulation script was indicated but was not found at the file path: ",astrFunctionPath,sep=""))
      return(-1)
    }
  }
}

if(!file.exists(inputFile))
{
  print("The following file does not exist please place in the appropriate directory or supply an appriate .pcl file for conversion. File:",inputFile,sep="")
  return(-1)
}

#Read file
inputFileData = funcReadMatrices(inputFile, strInputTSV, log=TRUE)

#Dimensions of the datasets
metaData = dim(inputFileData[[c_strMatrixMetadata]])
data = dim(inputFileData[[c_strMatrixData]])

#Check data before the merge
dataToWrite = list(Metadata = inputFileData[[c_strMatrixMetadata]])
saveToFile = c(paste(outputDirectory,"metadata.tsv",sep=""))
configFileName = c(paste(outputDirectory,"metadata.read.config",sep=""))
funcWriteMatrices(dataFrameList=dataToWrite, saveFileList=saveToFile, configureFileName=configFileName, acharDelimiter="\t")

dataToWrite = list(Data = inputFileData[[c_strMatrixData]])
saveToFile = c(paste(outputDirectory,"data.tsv",sep=""))
configFileName = c(paste(outputDirectory,"data.read.config",sep=""))
funcWriteMatrices(dataFrameList=dataToWrite, saveFileList=saveToFile, configureFileName=configFileName, acharDelimiter="\t")

#Merge data files together
frmeData = merge(inputFileData[[c_strMatrixMetadata]],inputFileData[[c_strMatrixData]],by.x=0,by.y=0)

#Reset rownames
row.names(frmeData) = frmeData[[1]]
frmeData = frmeData[-1]

#Record the data as it has been read
dataToWrite = list(Merged = frmeData)
saveToFile = c(paste(outputDirectory,"read-Merged.tsv",sep=""))
configFileName = c(paste(outputDirectory,"read-Merged.read.config",sep=""))
funcWriteMatrices(dataFrameList=dataToWrite, saveFileList=saveToFile, configureFileName=configFileName, acharDelimiter="\t")

#Data needed for the MaAsLin environment
#List of lists (one entry per file)
#Is contained by a container of itself
#lslsData = list()
#List
lsData = c()

#List of metadata indicies
aiMetadata = c(1:metaData[2])
lsData$aiMetadata = aiMetadata
#List of genetics indicies
aiGenetics = c()
lsData$aiGenetics = aiGenetics
#List of data indicies
aiData = c(1:data[2])+metaData[2]
lsData$aiData = aiData
#Add a list to hold qc metrics and counts
lsData$lsQCCounts$aiDataInitial = aiData
lsData$lsQCCounts$aiMetadataInitial = aiMetadata

#Raw data
lsData$frmeRaw = frmeData

lsScript = funcSourceScript(customDataProcessFunction)

#Clean the data and update the current data list to the cleaned data list
lsRet = funcClean( frmeData=frmeData, funcDataProcess=lsScript$processFunction, aiMetadata=aiMetadata, aiGenetics=aiGenetics, aiData=aiData, lsQCCounts=lsData$lsQCCounts, astrNoImpute=lsScript$noImpute )
print("lsRet")
print(lsRet)
#Update the variables after cleaning
lsRet$frmeRaw = frmeData
lsData = lsRet
frmeData = lsData$frmeData
aiMetadata = lsData$aiMetadata
aiGenetics = lsData$aiGenetics
aiData = lsData$aiData
lsData$lsQCCounts$aiDataCleaned = aiData
lsData$lsQCCounts$aiMetadataCleaned = aiMetadata

#Add List of metadata string names
astrMetadata = colnames(frmeData)[aiMetadata]
lsData$astrMetadata = astrMetadata

#Record the data after cleaning
dataToWrite = list(Cleaned = frmeData)
saveToFile = c(paste(outputDirectory,"read_cleaned.tsv",sep=""))
configFileName = c(paste(outputDirectory,"read_cleaned.read.config",sep=""))
funcWriteMatrices(dataFrameList=dataToWrite, saveFileList=saveToFile, configureFileName=configFileName, acharDelimiter="\t")

#Log file
strData = paste(outputDirectory,strBase,".txt",sep="")

#These variables will be used to count how many features get analysed
lsData$lsQCCounts$iBoosts = 0
lsData$lsQCCounts$iBoostErrors = 0
lsData$lsQCCounts$iNoTerms = 0
lsData$lsQCCounts$iLms = 0

#Run analysis
alsRetBugs = funcBugs( frmeData, funcBugHybrid, lsData, aiMetadata, aiGenetics, aiData, strData,
  lsScript$significanceLevel, lsScript$invert, dirname(strOutputTXT), astrScreen = c() )
aiBugs = alsRetBugs$aiReturnBugs
lsQCCounts = alsRetBugs$lsQCCounts

#Numeric vector of Metadata indexes
aiUMD <- intersect( aiMetadata, which( colnames( frmeData ) %in% lsData$astrMetadata ) )

#Output a summary file of analysis process
strProcessFileName = c(paste(outputDirectory,"ProcessQC.txt",sep=""))
unlink(strProcessFileName)
funcWrite( paste("Initial Metadata Matrix Size: Rows ",metaData[1],"  Columns ",metaData[2],sep=""), strProcessFileName )
funcWrite( paste("Initial Data Matrix Size: Rows ",data[1],"  Columns ",data[2],sep=""), strProcessFileName )
funcWrite( paste("\nInitial Data Count: ",length(lsQCCounts$aiDataInitial),sep=""), strProcessFileName )
funcWrite( paste("Initial Metadata Count: ",length(lsQCCounts$aiMetadataInitial),sep=""), strProcessFileName )
funcWrite( paste("Data Count after preprocess: ",length(lsQCCounts$aiAfterPreprocess),sep=""), strProcessFileName )
funcWrite( paste("Removed for missing metadata: ",length(lsQCCounts$iMissingMetadata),sep=""), strProcessFileName )
funcWrite( paste("Removed for missing data: ",length(lsQCCounts$iMissingData),sep=""), strProcessFileName )
funcWrite( paste("Data with outliers: ",length(lsQCCounts$aiSumOutlierPerDatum[lsQCCounts$aiSumOutlierPerDatum>0]),sep=""), strProcessFileName )
funcWrite( paste("Metadata count which survived clean: ",length(lsQCCounts$aiMetadataCleaned),sep=""), strProcessFileName )
funcWrite( paste("Data count which survived clean: ",length(lsQCCounts$aiDataCleaned),sep=""), strProcessFileName )
funcWrite( paste("\nBoostings: ",lsQCCounts$iBoosts,sep=""), strProcessFileName )
funcWrite( paste("Boosting Errors: ",lsQCCounts$iBoostErrors,sep=""), strProcessFileName )
funcWrite( paste("LMs with no terms suriving boosting: ",lsQCCounts$iNoTerms,sep=""), strProcessFileName )
funcWrite( paste("LMs performed: ",lsQCCounts$iLms,sep=""), strProcessFileName )
if(!is.null(lsQCCounts$lsQCCustom))
{
  funcWrite("Custom preprocess QC data: ", strProcessFileName )
  funcWrite(lsQCCounts$lsQCCustom, strProcessFileName )
} else {
  funcWrite("No custom preprocess QC data.", strProcessFileName )
}
funcWrite( "\n#Details###########################", strProcessFileName )
funcWrite("\nInitial Data Count: ", strProcessFileName )
funcWrite(lsQCCounts$aiDataInitial, strProcessFileName )
funcWrite("\nInitial Metadata Count: ", strProcessFileName )
funcWrite(lsQCCounts$aiMetadataInitial, strProcessFileName )
funcWrite("\nData Count after preprocess: ", strProcessFileName )
funcWrite(lsQCCounts$aiAfterPreprocess, strProcessFileName )
funcWrite("\nRemoved for missing metadata: ", strProcessFileName )
funcWrite(lsQCCounts$iMissingMetadata, strProcessFileName )
funcWrite("\nRemoved for missing data: ", strProcessFileName )
funcWrite(lsQCCounts$iMissingData, strProcessFileName )
funcWrite("\nOutlier Count per Datum: ", strProcessFileName )
funcWrite(lsQCCounts$aiSumOutlierPerDatum, strProcessFileName )
funcWrite("\nMetadata which survived clean: ", strProcessFileName )
funcWrite(lsQCCounts$aiMetadataCleaned, strProcessFileName )
funcWrite("\nData which survived clean: ", strProcessFileName )
funcWrite(lsQCCounts$aiDataCleaned, strProcessFileName )

#Run MFA and plot covariance of factors
if( length( aiBugs ) )
{
    print("MFA:in")
    lsMFA <- funcMFA( frmeData, aiUMD, aiBugs )
    print("MFA:out")
    if( class( lsMFA ) != "try-error" )
    {
        print("PlotMFA:in")
        funcPlotMFA( lsMFA, paste(outputDirectory,strBase,sep="") )
        print("PlotMFA:out")
    }
}

#Summarize output files based on a keyword and a significance threshold
#Look for less than or equal to the threshold (approapriate for p-value and q-value type measurements)
funcSummarizeDirectory(astrOutputDirectory=outputDirectory,
                       strBaseName = strBase,
                       astrSummaryFileName=strSummaryFileName, 
                       astrKeyword=c_strKeywordEvaluatedForInclusion, 
                       afSignificanceLevel=lsScript$significanceLevel)
