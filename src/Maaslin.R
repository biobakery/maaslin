#!/usr/bin/env Rscript
#######################################################################################
# This file is provided under the Creative Commons Attribution 3.0 license.
#
# You are free to share, copy, distribute, transmit, or adapt this work
# PROVIDED THAT you attribute the work to the authors listed below.
# For more information, please see the following web page:
# http://creativecommons.org/licenses/by/3.0/
#
# This file is a component of the MaAsLin (Multivariate Associations Using Linear Models), 
# authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Timothy Tickle, ttickle@hsph.harvard.edu).
#######################################################################################

funcGetAnalysisMethods <- function(
### Returns the appropriate functions for regularization, analysis, data transformation, and analysis object inspection.
### This allows modular customization per analysis step.
### To add a new method insert an entry in the switch for either the selection, transform, or method
### Insert them by using the pattern optparse_keyword_without_quotes = function_in_AnalysisModules
### Order in the return listy is curretly set and expected to be selection, transforms/links, analsis method
### none returns null
sModelSelectionKey,
### Keyword defining the method of model selection
sTransformKey,
### Keyword defining the method of data transformation
sMethodKey
### Keyword defining the method of analysis
){
  lRetMethods = list()
  #Insert selection methods here
  lRetMethods[[c_iSelection]] = switch(sModelSelectionKey,
    boost			= funcBoostModel,
    forward			= funcForwardModel,
    backward		= funcBackwardsModel,
    none			= NULL)

  #Insert transforms
  lRetMethods[[c_iTransform]] = switch(sTransformKey,
    asinsqrt		= funcArcsinSqrt,
    none			= funcNoTransform)

  #Insert analysis
  lRetMethods[[c_iAnalysis]] = switch(sMethodKey,
    neg_binomial	= funcBinomialMult,
    quasi			= funcQuasiMult,
    spearman		= funcSpearman,
    wilcoxon		= funcWilcoxon,
    lasso			= funcLasso,
    lm				= funcLM,
    none			= NULL)

  #Insert method to get results
  lRetMethods[[c_iResults]] = switch(sMethodKey,
    neg_binomial	= funcGetLMResults,
    quasi			= funcGetLMResults,
    spearman		= funcGetUnivariateResults,
    wilcoxon		= funcGetUnivariateResults,
    lasso			= funcGetLassoResults,
    lm				= funcGetLMResults,
    none			= NULL)

  return(lRetMethods)
  ### Returns a list of functions to be passed for regularization, data transformation, analysis,
  ### and custom analysis results introspection functions to pull from return objects data of interest
}

inlinedocs <- function(
##author<< Curtis Huttenhower <chuttenh@hsph.harvard.edu> and Timothy Tickle <ttickle@hsph.harvard.edu>
##description<< Main driver script. Should be called to perform MaAsLin Analysis.
) { return( pArgs ) }

### Logging class
suppressMessages(library( logging, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
### Class for commandline argument processing
suppressMessages(library( optparse, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))

### Necessary local import files
initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
strDir = file.path( dirname( script.name ), "lib" )
strSelf = basename( script.name )
for( strR in dir( strDir, pattern = "*.R$" ) )
{
  if( strR == strSelf ) {next}
  source( file.path( strDir, strR ) )
}

### Create command line argument parser
pArgs <- OptionParser( usage = "%prog [options] <output.txt> <data.tsv>" )

# Settings for MaAsLin
## Data configuration file
pArgs <- add_option( pArgs, c("-i", "--input_config"), type="character", action="store", dest="strInputConfig", metavar="data.read.config", help="Optional configuration file describing data input format.")
## Data manipulation/normalization file
pArgs <- add_option( pArgs, c("-I", "--input_process"), type="character", action="store", dest="strInputR", metavar="data.R", help="Optional configuration script normalizing or processing data.")

## Maximum false discovery rate
pArgs <- add_option( pArgs, c("-d", "--fdr"), type="double", action="store", dest="dSignificanceLevel", default=0.25, metavar="significance", help="The threshold to use for significance for the generated q-values (BH FDR). Anything equal to or lower than this is significant.")
## Minimum feature relative abundance filtering
pArgs <- add_option( pArgs, c("-r", "--minRelativeAbundance"), type="double", action="store", dest="dMinAbd", default=0.0001, metavar="minRelativeAbundance", help="The minimum relative abundance allowed in the data. Values below this are removed and imputed as the median of the sample data.")
## Minimum feature prevalence filtering
pArgs <- add_option( pArgs, c("-p", "--minPrevalence"), type="double", action="store", dest="dMinSamp", default=0.1, metavar="minPrevalence", help="The minimum percentage of samples a feature can have abudance in before being removed.")
## Fence for outlier, if not set Grubbs test is used
pArgs <- add_option( pArgs, c("-o", "--outlierFence"), type="double", action="store", dest="dOutlierFence", default=3.0, metavar="outlierFence", help="Outliers are defined as this number times the interquartile range added/subtracted from the 3rd/1st quartiles respectively. If set to 0, outliers are defined by the Grubbs test.")

# Arguments used in validation of MaAsLin
## Model selection (enumerate) c("none","boost","forward","backward")
pArgs <- add_option( pArgs, c("-s", "--selection"), type="character", action="store", dest="strModelSelection", default="boost", metavar="model_selection", help="Indicates which of the model selection techniques to use.")
## Argument indicating which method should be ran (enumerate) c("wilcoxon","spearman","lm","lasso","neg_binomial","quasi")
pArgs <- add_option( pArgs, c("-m", "--method"), type="character", action="store", dest="strMethod", default="lm", metavar="analysis_method", help="Indicates which of the statistical analysis methods to run.")
## Argument indicating which link function is used c("none","asinsqrt")
pArgs <- add_option( pArgs, c("-l", "--link"), type="character", action="store", dest="strTransform", default="asinsqrt", metavar="transform_method", help="Indicates which link or transformation to use with a glm, if glm is not selected this argument will be set to none.")

# Arguments to supress MaAsLin actions on certain data
## Do not perform model selection on the following data
pArgs <- add_option( pArgs, c("-F","--forced"), type="character", action="store", dest="strForcedPredictors", default=NULL, metavar="forced_predictors", help="Metadata features that will be forced into the model. Example 'Metadata2|Metadata6|Metadata10'")
## Do not impute the following
pArgs <- add_option( pArgs, c("-n","--noImpute"), type="character", action="store", dest="strNoImpute", default=NULL, metavar="no_impute", help="These data will not be imputed. Pipe delimited data feature names. Example 'Feature1|Feature4|Feature6'")

#Miscellaneouse arguments
### Argument to control logging (enumerate)
strDefaultLogging = "INFO"
pArgs <- add_option( pArgs, c("-v", "--verbosity"), type="character", action="store", dest="strVerbosity", default=strDefaultLogging, metavar="verbosity", help="Logging verbosity")
### Argument for inverting background to black
pArgs <- add_option( pArgs, c("-t", "--invert"), type="logical", action="store_true", dest="fInvert", default="FALSE", metavar="invert", help="When given, flag indicates to invert the background of figures to black.")
### Selection Frequency
pArgs <- add_option( pArgs, c("-f","--selectionFrequency"), type="double", action="store", dest="dSelectionFrequency", default= NA, metavar="selectionFrequency", help="Selection Frequency")

main <- function(
### The main function manages the following:
### 1. Optparse arguments are checked
### 2. A logger is created if requested in the optional arguments
### 3. The custom R script is sourced. This is the input *.R script named
### the same as the input *.pcl file. This script contains custom formating
### of data and function calls to the MFA visualization.
### 4. Matrices are written to the project folder as they are read in seperately as metadata and data and merged together.
### 5. Data is cleaned with custom filtering if supplied in the *.R script.
### 6. Transformations occur if indicated by the optional arguments
### 7. Standard quality control is performed on data
### 8. Cleaned metadata and data are written to output project for documentation.
### 9. A regularization method is ran (boosting by default).
### 10. An analysis method is performed on the model (optionally boostd model).
### 11. Data is summarized and PDFs are created for significant associations
### (those whose q-values {BH FDR correction} are <= the threshold given in the optional arguments.
pArgs
### Parsed commandline arguments
) {
lsArgs <- parse_args( pArgs, positional_arguments = TRUE )
#logdebug("lsArgs", c_logrMaaslin)
#logdebug(paste(lsArgs,sep=" "), c_logrMaaslin)

###Default configurations
#c_astrConfigurationValues <- c("noImpute", "invertPlots", "significanceLevel", "selectionFrequency", "processFunction")
#c_lsConfigurationDefaults <- list(NULL, lsArgs$options$fInvert, lsArgs$options$dSignificanceLevel, NA, NULL)

# Parse Piped parameters
lsForcedParameters = if(!is.null(lsArgs$strForcedPredictors)){unlist(strsplit(lsArgs$strForcedPredictors,"[|]"))}else{NULL}
xNoImpute = if(!is.null(lsArgs$strNoImpute)){unlist(strsplit(lsArgs$strNoImpute,"[|]"))}else{NULL}

#If logging is not an allowable value, inform user and set to INFO
if(length(intersect(names(loglevels), c(lsArgs$options$strVerbosity))) == 0)
{
  print(paste("Maaslin::Error. Did not understand the value given for logging, please use any of the following: DEBUG,INFO,WARN,ERROR."))
  print(paste("Maaslin::Warning. Setting logging value to \"",strDefaultLogging,"\"."))
}

### Create logger
c_logrMaaslin <- getLogger( "maaslin" )
addHandler( writeToConsole, c_logrMaaslin )
setLevel( lsArgs$options$strVerbosity, c_logrMaaslin )

#Get positional arguments
if( length( lsArgs$args ) != 2 ) { stop( print_help( pArgs ) ) }
### Output file name
strOutputTXT <- lsArgs$args[1]
### Input TSV data file
strInputTSV <- lsArgs$args[2]

# Get analysis method options
# includes data transformations, model selection/regularization, regression models/links
lsArgs$options$strModelSelection = tolower(lsArgs$options$strModelSelection)
if(!lsArgs$options$strModelSelection %in% c("none","boost","forward","backward"))
{
  logerror(paste("Received an invalid value for the selection argument, received '",lsArgs$options$strModelSelection,"'"), c_logrMaaslin)
  stop( print_help( pArgs ) )
}
lsArgs$options$strMethod = tolower(lsArgs$options$strMethod)
if(!lsArgs$options$strMethod %in% c("wilcoxon","spearman","lm","lasso","neg_binomial","quasi"))
{
  logerror(paste("Received an invalid value for the method argument, received '",lsArgs$options$strMethod,"'"), c_logrMaaslin)
  stop( print_help( pArgs ) )
}
lsArgs$options$strTransform = tolower(lsArgs$options$strTransform)
if(!lsArgs$options$strTransform %in% c("none","asinsqrt"))
{
  logerror(paste("Received an invalid value for the transform/link argument, received '",lsArgs$options$strTransform,"'"), c_logrMaaslin)
  stop( print_help( pArgs ) )
}
# If lasso is selected, do not use a regularization technique. This will happen in the lasso call
if(lsArgs$options$strMethod == "lasso")
{
  logdebug(paste("Lasso was selected so no model selection ocurred outside the lasso call."), c_logrMaaslin)
  lsArgs$options$strModelSelection = "none"
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
strBase <- sub("\\.[^.]*$", "", basename(strInputTSV))

### Sources in the custom script
### If the custom script is not there then
### defaults are used and no custom scripts are ran
funcSourceScript <- function(strFunctionPath)
{
  print("File:::")
  print(strFunctionPath)
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
inputFileData = funcReadMatrices(lsArgs$options$strInputConfig, strInputTSV, log=TRUE)
if(is.null(inputFileData[[c_strMatrixMetadata]])) { names(inputFileData)[1] <- c_strMatrixMetadata }
if(is.null(inputFileData[[c_strMatrixData]])) { names(inputFileData)[2] <- c_strMatrixData }

#Dimensions of the datasets
liMetaData = dim(inputFileData[[c_strMatrixMetadata]])
liData = dim(inputFileData[[c_strMatrixData]])

#Merge data files together
frmeData = merge(inputFileData[[c_strMatrixMetadata]],inputFileData[[c_strMatrixData]],by.x=0,by.y=0)
#Reset rownames
row.names(frmeData) = frmeData[[1]]
frmeData = frmeData[-1]

#Write QC files only in certain modes of verbosity
if( c_logrMaaslin$level <= loglevels["DEBUG"] ) {
  # If the QC internal file does not exist, make
  strQCDir = file.path(outputDirectory,"QC")
  dir.create(strQCDir, showWarnings = FALSE)
  # Write metadata matrix before merge
  funcWriteMatrices(dataFrameList=list(Metadata = inputFileData[[c_strMatrixMetadata]]), saveFileList=c(file.path(strQCDir,"metadata.tsv")), configureFileName=c(file.path(strQCDir,"metadata.read.config")), acharDelimiter="\t")
  # Write data matrix before merge
  funcWriteMatrices(dataFrameList=list(Data = inputFileData[[c_strMatrixData]]), saveFileList=c(file.path(strQCDir,"data.tsv")), configureFileName=c(file.path(strQCDir,"data.read.config")), acharDelimiter="\t")
  #Record the data as it has been read
  funcWriteMatrices(dataFrameList=list(Merged = frmeData), saveFileList=c(file.path(strQCDir,"read-Merged.tsv")), configureFileName=c(file.path(strQCDir,"read-Merged.read.config")), acharDelimiter="\t")
}

#Data needed for the MaAsLin environment
#List of lists (one entry per file)
#Is contained by a container of itself
#lslsData = list()
#List
lsData = c()

#List of metadata indicies
aiMetadata = c(1:liMetaData[2])
lsData$aiMetadata = aiMetadata
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
print("Inputs:::")
print(lsArgs$options$strInputR)

if(!is.null(funcSourceScript(lsArgs$options$strInputR))){funcProcess <- get(c_strCustomProcessFunction)}

#Clean the data and update the current data list to the cleaned data list
lsRet = funcClean( frmeData=frmeData, funcDataProcess=funcProcess, aiMetadata=aiMetadata, aiData=aiData, lsQCCounts=lsData$lsQCCounts, astrNoImpute=xNoImpute, dMinSamp = lsArgs$options$dMinSamp, dFence=lsArgs$options$dOutlierFence, funcTransform=afuncVariableAnalysis[[c_iTransform]])
logdebug("lsRet", c_logrMaaslin)
logdebug(format(lsRet), c_logrMaaslin)
#Update the variables after cleaning
lsRet$frmeRaw = frmeData
lsRet$lsQCCounts$aiDataCleaned = lsRet$aiData
lsRet$lsQCCounts$aiMetadataCleaned = lsRet$aiMetadata

#Add List of metadata string names
astrMetadata = colnames(lsRet$frmeData)[lsRet$aiMetadata]
lsRet$astrMetadata = astrMetadata

#Write QC files only in certain modes of verbosity
if( c_logrMaaslin$level <= loglevels["DEBUG"] ) {
	#Record the data after cleaning
	funcWriteMatrices(dataFrameList=list(Cleaned = lsRet$frmeData), saveFileList=c(file.path(strQCDir,"read_cleaned.tsv")), configureFileName=c(file.path(strQCDir,"read_cleaned.read.config")), acharDelimiter="\t") }

#These variables will be used to count how many features get analysed
lsRet$lsQCCounts$iBoosts = 0
lsRet$lsQCCounts$iBoostErrors = 0
lsRet$lsQCCounts$iNoTerms = 0
lsRet$lsQCCounts$iLms = 0

#Run analysis
alsRetBugs = funcBugs( lsRet$frmeData, lsRet, lsRet$aiMetadata, lsRet$aiData, strBase,
	lsArgs$options$dSelectionFrequency, lsArgs$options$dSignificanceLevel, lsArgs$options$dMinSamp, lsArgs$options$fInvert,
        outputDirectory, astrScreen = c(), funcReg=afuncVariableAnalysis[[c_iSelection]],
        funcAnalysis=afuncVariableAnalysis[[c_iAnalysis]], funcGetResults=afuncVariableAnalysis[[c_iResults]] )
aiBugs = alsRetBugs$aiReturnBugs

#Write QC files only in certain modes of verbosity
if( c_logrMaaslin$level <= loglevels["DEBUG"] ) {
	funcWriteQCReport(strProcessFileName=file.path(strQCDir,"ProcessQC.txt"), lsQCData=alsRetBugs$lsQCCounts, liDataDim=liData, liMetadataDim=liMetaData) }

#Numeric vector of Metadata indexes or MFA
aiUMD <- intersect( lsRet$aiMetadata, which( colnames( lsRet$frmeData ) %in% lsRet$astrMetadata ) )

#Run MFA and plot covariance of factors
if( !length( aiBugs ) ) { aiBugs <- lsRet$aiData }
if( length( aiBugs ) )
{
  logdebug("MFA:in", c_logrMaaslin)
  lsMFA <- funcMFA( lsRet$frmeData, lsArgs$options$dMinSamp, aiUMD, aiBugs)
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
}

# This is the equivalent of __name__ == "__main__" in Python.
# That is, if it's true we're being called as a command line script;
# if it's false, we're being sourced or otherwise included, such as for
# library or inlinedocs.
if( identical( environment( ), globalenv( ) ) &&
	!length( grep( "^source\\(", sys.calls( ) ) ) ) {
	main( pArgs ) }
