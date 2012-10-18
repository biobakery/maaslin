#!/usr/bin/env Rscript

inlinedocs <- function(
##author<< Curtis Huttenhower <chuttenh@hsph.harvard.edu> and Timothy Tickle <ttickle@hsph.harvard.edu>
) { return( pArgs ) }

### Logging class
suppressMessages(library( logging, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
### Class for commandline argument processing
suppressMessages(library( optparse, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))

### Create command line argument parser
pArgs <- OptionParser( usage = "%prog [options] <output.txt> <data.tsv> <data.read.config> <data.R> [source.R]*" )

# Settings for MaAsLin
## Maximum false discovery rate
pArgs <- add_option( pArgs, c("-d", "--fdr"), type="double", action="store", dest="dSignificanceLevel", default=0.25, metavar="significance", help="The threshold to use for significance for the generated q-values (BH FDR). Anything equal to or lower than this is significant.")
## Minimum feature relative abundance filtering
pArgs <- add_option( pArgs, c("-r", "--minRelativeAbundance"), type="double", action="store", dest="dMinAbd", default=0.0001, metavar="minRelativeAbundance", help="The minimum relative abundance allowed in the data. Values below this are removed and imputed as the median of the sample data.")
## Minimum feature prevalence filtering
pArgs <- add_option( pArgs, c("-p", "--minPrevalence"), type="double", action="store", dest="dMinSamp", default=0.1, metavar="minPrevalence", help="The minimum percentage of samples a feature can have abudance in before being removed.")

# Arguments used in validation of MaAsLin
## Model selection (enumerate) c("none","boost","forward","backward")
pArgs <- add_option( pArgs, c("-s", "--selection"), type="character", action="store", dest="strModelSelection", default="boost", metavar="model selection", help="Indicates which of the model selection techniques to use.")
## Argument indicating which method should be ran (enumerate) c("lefse","wilcoxon","spearman","gunifrac","lm","lasso")
pArgs <- add_option( pArgs, c("-m", "--method"), type="character", action="store", dest="strMethod", default="lm", metavar="method", help="Indicates which of the statistical analysis methods to run.")
## Argument indicating which link function is used c("none","neg_binomial","quasi","asinsqrt")
pArgs <- add_option( pArgs, c("-l", "--link"), type="character", action="store", dest="strTransform", default="asinsqrt", metavar="method", help="Indicates which link or transformation to use with a glm, if glm is not selected this argument will be set to none.")

#Miscellaneouse arguments
## Argument to control logging (enumerate)
strDefaultLogging = "INFO"
pArgs <- add_option( pArgs, c("-v", "--verbosity"), type="character", action="store", dest="fVerbosity", default=strDefaultLogging, metavar="verbosity", help="Logging verbosity")
### Argument for inverting background to black
pArgs <- add_option( pArgs, c("-i", "--invert"), type="logical", action="store_true", dest="fInvert", default="FALSE", metavar="invert", help="When given, flag indicates to invert the background of figures to black.")
### Selection Frequency
pArgs <- add_option( pArgs, c("-f","--selectionFrequency"), type="double", action="store", dest="dSelectionFrequency", default= NA, metavar="selectionFrequency", help="Selection Frequency")

### Parse arguments
lsArgs <- parse_args( pArgs, positional_arguments = TRUE )
logdebug("lsArgs", c_logrMaaslin)
logdebug(paste(lsArgs,sep=" "), c_logrMaaslin)

###Default configurations
#c_astrConfigurationValues <- c("noImpute", "invertPlots", "significanceLevel", "selectionFrequency", "processFunction")
#c_lsConfigurationDefaults <- list(NULL, lsArgs$options$fInvert, lsArgs$options$dSignificanceLevel, NA, NULL)

# Constants
### Minimum Relative Abundance (###TODO fix shadowing)
c_dMinSamp <- lsArgs$options$dMinSamp

#TODO
xNoImpute = NULL

### Allowable values for logging
c_lsLoggingValues = c("INFO","WARNING","ERROR","CRITICAL","LOG","EXCEPTION")
#If logging is not an allowable value, inform user and set to INFO
if(length(intersect(c_lsLoggingValues, c(lsArgs$options$fVerbosity))) == 0)
{
  print(paste("Maaslin::Error. Did not understand the value given for logging, please use any of the following: ",c_lsLoggingValues,"."))
  print(paste("Maaslin::Warning. Setting logging value to \"",strDefaultLogging,"\"."))
}

### Create logger
c_logrMaaslin <- getLogger( "maaslin" )
addHandler( writeToConsole, c_logrMaaslin )
setLevel( lsArgs$options$fVerbosity, c_logrMaaslin )

#Get positional arguments
if( length( lsArgs$args ) < 4 ) { stop( print_help( pArgs ) ) }
### Output file name
strOutputTXT <- lsArgs$args[1]
### Input TSV data file
strInputTSV <- lsArgs$args[2]
### Input Read config file
strInputRC <- lsArgs$args[3]
### Input optional R Script
strInputR <- lsArgs$args[4]
### External libraries to source
astrSourceR <- lsArgs$args[5:length( lsArgs$args )]
### Source all libraries
for( strR in astrSourceR ){source( strR )}

# Get analysis method options
# includes data transformations, model selection/regularization, regression models/links
lsArgs$options$strModelSelection = tolower(lsArgs$options$strModelSelection)
if(!lsArgs$options$strModelSelection %in% c("none","boost","forward","backward"))
{
  logerror(paste("Received an invalid value for the selection argument, received '",lsArgs$options$strModelSelection,"'"), c_logrMaaslin)
  stop( print_help( pArgs ) )
}
lsArgs$options$strMethod = tolower(lsArgs$options$strMethod)
if(!lsArgs$options$strMethod %in% c("wilcoxon","spearman","lm","lasso"))
{
  logerror(paste("Received an invalid value for the method argument, received '",lsArgs$options$strMethod,"'"), c_logrMaaslin)
  stop( print_help( pArgs ) )
}
lsArgs$options$strTransform = tolower(lsArgs$options$strTransform)
if(!lsArgs$options$strTransform %in% c("none","neg_binomial","quasi","asinsqrt"))
{
  logerror(paste("Received an invalid value for the transform/link argument, received '",lsArgs$options$strTransform,"'"), c_logrMaaslin)
  stop( print_help( pArgs ) )
}
# Make sure that the link happens for lm only
if(!lsArgs$options$strMethod == "lm" && lsArgs$options$strTransform %in% c("neg_binomial","quasi"))
{
  logdebug(paste("Analysis method selected was not lm, the link option '",lsArgs$options$strTransform,"' was not used.", sep=""), c_logrMaaslin)
  lsArgs$options$strTransform = "none"
}
# If lasso is selected, do not use a regularization technique. This will happen in the lasso call
if(lsArgs$options$strMethod == "lasso")
{
  logdebug(paste("Lasso was selected so no model selection ocurred outside the lasso call."), c_logrMaaslin)
  llsArgs$options$strModelSelection = "none"
}

# Get analysis modules
afuncVariableAnalysis = funcGetAnalysisMethods(lsArgs$options$strModelSelection,lsArgs$options$strTransform,lsArgs$options$strMethod)

#Indicate start
logdebug("Start MaAsLin", c_logrMaaslin)
#Log commandline arguments
logdebug("Commandline Arguments", c_logrMaaslin)
logdebug(lsArgs, c_logrMaaslin)

### Output directory for the study based on the requested output file
outputDirectory = dirname(strOutputTXT)
### Base name for the project based on the read.config name
strBase <- sub(".read.config", "", basename(strInputRC))

### Sources in the custom script
### If the custom script is not there then
### defaults are used and no custom scripts are ran
funcSourceScript <- function(strFunctionPath)
{
  #If is specified, set up the custom func clean variable
  #If the custom script is null then return 
  if(is.null(strFunctionPath)){return(NULL)}

  #Check to make sure the file exists
  if(file.exists(strFunctionPath))
  {
    #Read in the file
    source(strFunctionPath)
  } else {
    #Handle when the file does not exist
    stop(paste("MaAsLin Error: A custom data manipulation script was indicated but was not found at the file path: ",strFunctionPath,sep=""))
  }
}

#Read file
inputFileData = funcReadMatrices(strInputRC, strInputTSV, log=TRUE)

#Dimensions of the datasets
liMetaData = dim(inputFileData[[c_strMatrixMetadata]])
liData = dim(inputFileData[[c_strMatrixData]])

# If the QC internal file does not exist, make
strQCDir = file.path(outputDirectory,"QC")
dir.create(strQCDir, showWarnings = FALSE)

# Write metadata matrix before merge
funcWriteMatrices(dataFrameList=list(Metadata = inputFileData[[c_strMatrixMetadata]]), saveFileList=c(file.path(strQCDir,"metadata.tsv")), configureFileName=c(file.path(strQCDir,"metadata.read.config")), acharDelimiter="\t")

# Write data matrix before merge
funcWriteMatrices(dataFrameList=list(Data = inputFileData[[c_strMatrixData]]), saveFileList=c(file.path(strQCDir,"data.tsv")), configureFileName=c(file.path(strQCDir,"data.read.config")), acharDelimiter="\t")

#Merge data files together
frmeData = merge(inputFileData[[c_strMatrixMetadata]],inputFileData[[c_strMatrixData]],by.x=0,by.y=0)

#Reset rownames
row.names(frmeData) = frmeData[[1]]
frmeData = frmeData[-1]

#Record the data as it has been read
funcWriteMatrices(dataFrameList=list(Merged = frmeData), saveFileList=c(file.path(strQCDir,"read-Merged.tsv")), configureFileName=c(file.path(strQCDir,"read-Merged.read.config")), acharDelimiter="\t")

#Data needed for the MaAsLin environment
#List of lists (one entry per file)
#Is contained by a container of itself
#lslsData = list()
#List
lsData = c()

#List of metadata indicies
aiMetadata = c(1:liMetaData[2])
lsData$aiMetadata = aiMetadata
#List of genetics indicies
aiGenetics = c()
lsData$aiGenetics = aiGenetics
#List of data indicies
aiData = c(1:liData[2])+liMetaData[2]
lsData$aiData = aiData
#Add a list to hold qc metrics and counts
lsData$lsQCCounts$aiDataInitial = aiData
lsData$lsQCCounts$aiMetadataInitial = aiMetadata

#Raw data
lsData$frmeRaw = frmeData

#Load script if it exists, stop on error
funcProcess <- NULL
if(!is.null(funcSourceScript(strInputR))){funcProcess <- get(c_strCustomProcessFunction)}

#Clean the data and update the current data list to the cleaned data list
lsRet = funcClean( frmeData=frmeData, funcDataProcess=funcProcess, aiMetadata=aiMetadata, aiGenetics=aiGenetics, aiData=aiData, lsQCCounts=lsData$lsQCCounts, astrNoImpute=xNoImpute )
logdebug("lsRet", c_logrMaaslin)
logdebug(format(lsRet), c_logrMaaslin)
#Update the variables after cleaning
lsRet$frmeRaw = frmeData
lsRet$lsQCCounts$aiDataCleaned = lsRet$aiData
lsRet$lsQCCounts$aiMetadataCleaned = lsRet$aiMetadata

#Add List of metadata string names
astrMetadata = colnames(lsRet$frmeData)[lsRet$aiMetadata]
lsRet$astrMetadata = astrMetadata

#Record the data after cleaning
funcWriteMatrices(dataFrameList=list(Cleaned = lsRet$frmeData), saveFileList=c(file.path(strQCDir,"read_cleaned.tsv")), configureFileName=c(file.path(strQCDir,"read_cleaned.read.config")), acharDelimiter="\t")

#Log file
strData = file.path(outputDirectory,paste(strBase,".txt",sep=""))

#These variables will be used to count how many features get analysed
lsRet$lsQCCounts$iBoosts = 0
lsRet$lsQCCounts$iBoostErrors = 0
lsRet$lsQCCounts$iNoTerms = 0
lsRet$lsQCCounts$iLms = 0

#Run analysis
alsRetBugs = funcBugs( lsRet$frmeData, lsRet, lsRet$aiMetadata, lsRet$aiGenetics, lsRet$aiData, strData,
	lsArgs$options$dSelectionFrequency, lsArgs$options$dSignificanceLevel, lsArgs$options$fInvert,
        outputDirectory, astrScreen = c(), funcReg=afuncVariableAnalysis[[c_iSelection]],
        funcAnalysis=afuncVariableAnalysis[[c_iAnalysis]] )
aiBugs = alsRetBugs$aiReturnBugs

#Output a summary file of analysis process
#funcWriteQCReport(strProcessFileName=file.path(strQCDir,"ProcessQC.txt"), lsQCData=alsRetBugs$lsQCCounts, liDataDim=liData, liMetadataDim=liMetaData)

#Numeric vector of Metadata indexes or MFA
aiUMD <- intersect( lsRet$aiMetadata, which( colnames( lsRet$frmeData ) %in% lsRet$astrMetadata ) )

#Run MFA and plot covariance of factors
if( !length( aiBugs ) ) { aiBugs <- lsRet$aiData }
if( length( aiBugs ) )
{
  logdebug("MFA:in", c_logrMaaslin)
  lsMFA <- funcMFA( lsRet$frmeData, aiUMD, aiBugs )
  logdebug("MFA:out", c_logrMaaslin)
  if( class( lsMFA ) != "try-error" )
  {
    logdebug("PlotMFA:in", c_logrMaaslin)
    funcPlotMFA( fInvert=lsArgs$options$fInvert, lsMFA=lsMFA, tempSaveFileName=file.path(outputDirectory,strBase), funcPlotColors=lsRet$funcPlotColors, funcPlotPoints=lsRet$funcPlotPoints, funcPlotLegend=lsRet$funcPlotLegend )
    logdebug("PlotMFA:out", c_logrMaaslin)
  }
}

#Summarize output files based on a keyword and a significance threshold
#Look for less than or equal to the threshold (approapriate for p-value and q-value type measurements)
funcSummarizeDirectory(astrOutputDirectory=outputDirectory,
                       strBaseName=strBase,
                       astrSummaryFileName=file.path(outputDirectory,paste(strBase,c_sSummaryFileSuffix, sep="")), 
                       astrKeyword=c_strKeywordEvaluatedForInclusion, 
                       afSignificanceLevel=lsArgs$options$dSignificanceLevel)
