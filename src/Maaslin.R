#!/usr/bin/env Rscript
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
##description<< Main driver script. Should be called to perform MaAsLin Analysis.
) { return( pArgs ) }

funcGetAnalysisMethods <- function(
### Returns the appropriate functions for regularization, analysis, data transformation, and analysis object inspection.
### This allows modular customization per analysis step.
### To add a new method insert an entry in the switch for either the selection, transform, or method
### Insert them by using the pattern optparse_keyword_without_quotes = function_in_AnalysisModules
### Order in the return listy is curretly set and expected to be selection, transforms/links, analysis method
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
    boost = funcBoostModel,
    penalized = funcPenalizedModel,
    forward = funcForwardModel,
    backward = funcBackwardsModel,
    none = NA)

  #Insert transforms
  lRetMethods[[c_iTransform]] = switch(sTransformKey,
    asinsqrt = funcArcsinSqrt,
    none = funcNoTransform)

  #Insert untransform
  lRetMethods[[c_iUnTransform]] = switch(sTransformKey,
    asinsqrt = funcNoTransform,
    none = funcNoTransform)

  #Insert analysis
  lRetMethods[[c_iAnalysis]] = switch(sMethodKey,
    neg_binomial = funcBinomialMult,
    quasi = funcQuasiMult,
    univariate = funcDoUnivariate,
    lm = funcLM,
    none = NA)

  #Insert method to get results
  lRetMethods[[c_iResults]] = switch(sMethodKey,
    neg_binomial = funcGetLMResults,
    quasi = funcGetLMResults,
    univariate = funcGetUnivariateResults,
    lm = funcGetLMResults,
    none = NA)

  return(lRetMethods)
  ### Returns a list of functions to be passed for regularization, data transformation, analysis,
  ### and custom analysis results introspection functions to pull from return objects data of interest
}

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

# Input files for MaAsLin
## Data configuration file
pArgs <- add_option( pArgs, c("-i", "--input_config"), type="character", action="store", dest="strInputConfig", metavar="data.read.config", help="Optional configuration file describing data input format.")
## Data manipulation/normalization file
pArgs <- add_option( pArgs, c("-I", "--input_process"), type="character", action="store", dest="strInputR", metavar="data.R", help="Optional configuration script normalizing or processing data.")

# Settings for MaAsLin
## Maximum false discovery rate
pArgs <- add_option( pArgs, c("-d", "--fdr"), type="double", action="store", dest="dSignificanceLevel", default=0.25, metavar="significance", help="The threshold to use for significance for the generated q-values (BH FDR). Anything equal to or lower than this is significant.")
## Minimum feature relative abundance filtering
pArgs <- add_option( pArgs, c("-r", "--minRelativeAbundance"), type="double", action="store", dest="dMinAbd", default=0.0001, metavar="minRelativeAbundance", help="The minimum relative abundance allowed in the data. Values below this are removed and imputed as the median of the sample data.")
## Minimum feature prevalence filtering
pArgs <- add_option( pArgs, c("-p", "--minPrevalence"), type="double", action="store", dest="dMinSamp", default=0.1, metavar="minPrevalence", help="The minimum percentage of samples a feature can have abundance in before being removed. Also is the minimum percentage of samples a metadata can have NAs before being removed.")
## Fence for outlier, if not set Grubbs test is used
pArgs <- add_option( pArgs, c("-o", "--outlierFence"), type="double", action="store", dest="dOutlierFence", default=0, metavar="outlierFence", help="Outliers are defined as this number times the interquartile range added/subtracted from the 3rd/1st quartiles respectively. If set to 0 (default), outliers are defined by the Grubbs test.")
## Fixed (not random) covariates
pArgs <- add_option( pArgs, c("-R","--random"), type="character", action="store", dest="strRandomCovariates", default=NULL, metavar="fixed", help="These metadata will be treated as random covariates. Comma delimited data feature names. These features must be listed in the read.config file. Example '-R RandomMetadata1,RandomMetadata2'")
## Change the type of correction fo rmultiple corrections
pArgs <- add_option( pArgs, c("-T","--testingCorrection"), type="character", action="store", dest="strMultTestCorrection", default="BH", metavar="multipleTestingCorrection", help="This indicates which multiple hypothesis testing method will be used, available are holm, hochberg, hommel, bonferroni, BH, BY")

# Arguments used in validation of MaAsLin
## Model selection (enumerate) c("none","boost","penalized","forward","backward")
pArgs <- add_option( pArgs, c("-s", "--selection"), type="character", action="store", dest="strModelSelection", default="boost", metavar="model_selection", help="Indicates which of the variable selection techniques to use.")
## Argument indicating which method should be ran (enumerate) c("univariate","lm","neg_binomial","quasi")
pArgs <- add_option( pArgs, c("-m", "--method"), type="character", action="store", dest="strMethod", default="lm", metavar="analysis_method", help="Indicates which of the statistical inference methods to run.")
## Argument indicating which link function is used c("none","asinsqrt")
pArgs <- add_option( pArgs, c("-l", "--link"), type="character", action="store", dest="strTransform", default="asinsqrt", metavar="transform_method", help="Indicates which link or transformation to use with a glm, if glm is not selected this argument will be set to none.")

# Arguments to suppress MaAsLin actions on certain data
## Do not perform model selection on the following data
pArgs <- add_option( pArgs, c("-F","--forced"), type="character", action="store", dest="strForcedPredictors", default=NULL, metavar="forced_predictors", help="Metadata features that will be forced into the model seperated by commas. These features must be listed in the read.config file. Example '-F Metadata2,Metadata6,Metadata10'")
## Do not impute the following
pArgs <- add_option( pArgs, c("-n","--noImpute"), type="character", action="store", dest="strNoImpute", default=NULL, metavar="no_impute", help="These data will not be imputed. Comma delimited data feature names. Example '-n Feature1,Feature4,Feature6'")

#Miscellaneouse arguments
### Argument to control logging (enumerate)
strDefaultLogging = "INFO"
pArgs <- add_option( pArgs, c("-v", "--verbosity"), type="character", action="store", dest="strVerbosity", default=strDefaultLogging, metavar="verbosity", help="Logging verbosity")
### Argument for inverting background to black
pArgs <- add_option( pArgs, c("-O","--omitLogFile"), type="logical", action="store_true", default=FALSE, dest="fOmitLogFile", metavar="omitlogfile",help="Including this flag will stop the creation of the output log file.")
### Run maaslin without creating a log file
pArgs <- add_option( pArgs, c("-t", "--invert"), type="logical", action="store_true", dest="fInvert", default=FALSE, metavar="invert", help="When given, flag indicates to invert the background of figures to black.")
### Selection Frequency
pArgs <- add_option( pArgs, c("-f","--selectionFrequency"), type="double", action="store", dest="dSelectionFrequency", default=NA, metavar="selectionFrequency", help="Selection Frequency for boosting (max 1 will remove almost everything). Interpreted as requiring boosting to select metadata 100% percent of the time (or less if given a number that is less). Value should be between 1 (100&) and 0 (0%), NA (default is determined by data size).")
### All v All
pArgs <- add_option( pArgs, c("-a","--allvall"), type="logical", action="store_true", dest="fAllvAll", default=FALSE, metavar="compare_all", help="When given, the flag indicates that each fixed covariate that is not indicated as Forced is compared once at a time per data feature (bug). Made to be used with the -F option to specify one part of the model while allowing the other to cycle through a group of covariates. Does not affect Random covariates, which are always included when specified.")
pArgs <- add_option( pArgs, c("-N","--PlotNA"), type="logical", action="store_true", default=FALSE, dest="fPlotNA", metavar="plotNAs",help="Plot data that was originally NA, by default they are not plotted.")
### Alternative methodology settings
pArgs <- add_option( pArgs, c("-A","--pAlpha"), type="double", action="store", dest="dPenalizedAlpha", default=0.95, metavar="PenalizedAlpha",help="The alpha for penalization (1.0=L1 regularization, LASSO; 0.0=L2 regularization, ridge regression")

### Misc MFA plot arguments
pArgs <- add_option( pArgs, c("-c","--MFAFeatureCount"), type="integer", action="store", dest="iMFAMaxFeatures", default=3, metavar="maxMFAFeature", help="Number of features or number of bugs to plot (default=3; 3 metadata and 3 data).")
pArgs <- add_option( pArgs, c("-M","--MFAMetadataScale"), type="double", action="store", dest="dMFAMetadataScale", default=NULL, metavar="scaleForMetadata", help="A real number used to scale the metadata labels on the MFA plot (otherwise a default will be selected from the data).")
pArgs <- add_option( pArgs, c("-D","--MFADataScale"), type="double", action="store", dest="dMFADataScale", default=NULL, metavar="scaleForMetadata", help="A real number used to scale the metadata labels on the MFA plot (otherwise a default will be selected from the data).")
pArgs <- add_option( pArgs, c("-C", "--MFAColor"), type="character", action="store", dest="strMFAColor", default=NULL, metavar="MFAColorCovariate", help="A continuous metadata that will be used to color samples in the MFA ordination plot (otherwise a default will be selected from the data).")
pArgs <- add_option( pArgs, c("-S", "--MFAShape"), type="character", action="store", dest="strMFAShape", default=NULL, metavar="MFAShapeCovariate", help="A discontinuous metadata that will be used to indicate shapes of samples in the MFA ordination plot (otherwise a default will be selected from the data).")
pArgs <- add_option( pArgs, c("-P", "--MFAPlotFeatures"), type="character", action="store", dest="strMFAPlotFeatures", default=NULL, metavar="MFAFeaturesToPlot", help="Metadata and data features to plot (otherwise a default will be selected from the data). Comma Delimited.")

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

# Parse Piped parameters
lsForcedParameters = NULL
if(!is.null(lsArgs$options$strForcedPredictors))
{
  lsForcedParameters  = unlist(strsplit(lsArgs$options$strForcedPredictors,","))
}
xNoImpute = NULL
if(!is.null(lsArgs$options$strNoImpute))
{
  xNoImpute = unlist(strsplit(lsArgs$options$strNoImpute,"[,]"))
}
lsRandomCovariates = NULL
if(!is.null(lsArgs$options$strRandomCovariates))
{
  lsRandomCovariates = unlist(strsplit(lsArgs$options$strRandomCovariates,"[,]"))
}
lsFeaturesToPlot = NULL
if(!is.null(lsArgs$options$strMFAPlotFeatures))
{
  lsFeaturesToPlot = unlist(strsplit(lsArgs$options$strMFAPlotFeatures,"[,]"))
}

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
if(!lsArgs$options$strModelSelection %in% c("none","boost","penalized","forward","backward"))
{
  logerror(paste("Received an invalid value for the selection argument, received '",lsArgs$options$strModelSelection,"'"), c_logrMaaslin)
  stop( print_help( pArgs ) )
}
lsArgs$options$strMethod = tolower(lsArgs$options$strMethod)
if(!lsArgs$options$strMethod %in% c("univariate","lm","neg_binomial","quasi"))
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

if(!lsArgs$options$strMultTestCorrection %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))
{
  logerror(paste("Received an invalid value for the multiple testing correction argument, received '",lsArgs$options$strMultTestCorrection,"'"), c_logrMaaslin)
  stop( print_help( pArgs ) )
}

# Get analysis modules
afuncVariableAnalysis = funcGetAnalysisMethods(lsArgs$options$strModelSelection,lsArgs$options$strTransform,lsArgs$options$strMethod)

# Set up parameters for variable selection
lxParameters = list(dFreq=lsArgs$options$dSelectionFrequency, dPAlpha=lsArgs$options$dPenalizedAlpha)
if((lsArgs$options$strMethod == "lm")||(lsArgs$options$strMethod == "univariate")){ lxParameters$sFamily = "gaussian"
} else if(lsArgs$options$strMethod == "neg_binomial"){ lxParameters$sFamily = "binomial"
} else if(lsArgs$options$strMethod == "quasi"){ lxParameters$sFamily = "poisson"}

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

#Metadata and bug names
lsOriginalMetadataNames = names(inputFileData[[c_strMatrixMetadata]])
lsOriginalFeatureNames = names(inputFileData[[c_strMatrixData]])

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

if(!is.null(funcSourceScript(lsArgs$options$strInputR))){funcProcess <- get(c_strCustomProcessFunction)}

#Clean the data and update the current data list to the cleaned data list
lsRet = funcClean( frmeData=frmeData, funcDataProcess=funcProcess, aiMetadata=aiMetadata, aiData=aiData, lsQCCounts=lsData$lsQCCounts, astrNoImpute=xNoImpute, dMinSamp = lsArgs$options$dMinSamp, dMinAbd = lsArgs$options$dMinAbd, dFence=lsArgs$options$dOutlierFence, funcTransform=afuncVariableAnalysis[[c_iTransform]])
logdebug("lsRet", c_logrMaaslin)
logdebug(format(lsRet), c_logrMaaslin)
#Update the variables after cleaning
lsRet$frmeRaw = frmeData
lsRet$lsQCCounts$aiDataCleaned = lsRet$aiData
lsRet$lsQCCounts$aiMetadataCleaned = lsRet$aiMetadata

#Add List of metadata string names
astrMetadata = colnames(lsRet$frmeData)[lsRet$aiMetadata]
lsRet$astrMetadata = astrMetadata

# If plotting NA data reset the NA metadata indices to empty so they will not be excluded
if(lsArgs$options$fPlotNA)
{
  lsRet$liNaIndices = list()
}

#Write QC files only in certain modes of verbosity
if( c_logrMaaslin$level <= loglevels["DEBUG"] ) {
	#Record the data after cleaning
	funcWriteMatrices(dataFrameList=list(Cleaned = lsRet$frmeData[union(lsRet$aiMetadata,lsRet$aiData)]), saveFileList=c(file.path(strQCDir,"read_cleaned.tsv")), configureFileName=c(file.path(strQCDir,"read_cleaned.read.config")), acharDelimiter="\t") }

#These variables will be used to count how many features get analysed
lsRet$lsQCCounts$iBoosts = 0
lsRet$lsQCCounts$iBoostErrors = 0
lsRet$lsQCCounts$iNoTerms = 0
lsRet$lsQCCounts$iLms = 0

#Indicate if the residuals plots should occur
fDoRPlot=TRUE
#Should not occur for univariates
if(lsArgs$options$strMethod %in% c("univariate")){ fDoRPlot=FALSE }

#Run analysis
alsRetBugs = funcBugs( lsRet$frmeData, lsRet, lsRet$aiMetadata, lsRet$aiData, strBase, lsArgs$options$dSignificanceLevel, lsArgs$options$dMinSamp, lsArgs$options$fInvert,
        outputDirectory, astrScreen = c(), funcReg=afuncVariableAnalysis[[c_iSelection]], funcUnTransform=afuncVariableAnalysis[[c_iUnTransform]], lsForcedParameters,
        funcAnalysis=afuncVariableAnalysis[[c_iAnalysis]], lsRandomCovariates, funcGetResults=afuncVariableAnalysis[[c_iResults]], fDoRPlot=fDoRPlot, fOmitLogFile=lsArgs$options$fOmitLogFile, fAllvAll=lsArgs$options$fAllvAll, liNaIndices=lsRet$liNaIndices, lxParameters=lxParameters, strTestingCorrection=lsArgs$options$strMultTestCorrection )
aiBugs = alsRetBugs$aiReturnBugs

#Write QC files only in certain modes of verbosity
if( c_logrMaaslin$level <= loglevels["DEBUG"] ) {
	funcWriteQCReport(strProcessFileName=file.path(strQCDir,"ProcessQC.txt"), lsQCData=alsRetBugs$lsQCCounts, liDataDim=liData, liMetadataDim=liMetaData) }

#Numeric vector of Metadata indexes or MFA
aiUMD <- intersect( lsRet$aiMetadata, which( colnames( lsRet$frmeData ) %in% lsRet$astrMetadata ) )

#Run MFA and plot covariance of factors
if( !length( aiBugs ) ) { aiBugs <- lsRet$aiData }

#If a user defines a feature, make sure it is in the bugs/data indices
if(!is.null(lsFeaturesToPlot) || !is.null(lsArgs$options$strMFAColor) || !is.null(lsArgs$options$strMFAShape))
{
  lsCombinedFeaturesToPlot = unique(c(lsFeaturesToPlot,lsArgs$options$strMFAColor,lsArgs$options$strMFAShape))
  lsCombinedFeaturesToPlot = lsCombinedFeaturesToPlot[!is.null(lsCombinedFeaturesToPlot)]

  aiUMD = unique(c(aiUMD,which( colnames( lsRet$frmeData ) %in% intersect(lsCombinedFeaturesToPlot, lsOriginalMetadataNames))))
  aiBugs = unique(c(aiBugs,which( colnames( lsRet$frmeData ) %in% intersect(lsCombinedFeaturesToPlot, lsOriginalFeatureNames))))
}

#try(
  if( length( aiBugs ) )
  {
    logdebug("MFA:in", c_logrMaaslin)
    lsMFA <- funcMFA( lsRet$frmeData, lsArgs$options$dMinSamp, aiUMD, aiBugs)
    logdebug("MFA:out", c_logrMaaslin)

    if( class( lsMFA ) != "try-error" )
    {
      logdebug("PlotMFA:in", c_logrMaaslin)
      funcPlotMFA( lsMFA=lsMFA, frmeData=lsRet$frmeData, lsMetadata=lsOriginalMetadataNames, lsFeatures=lsOriginalFeatureNames, iMaxFeatures=lsArgs$options$iMFAMaxFeatures, strMFAColorCovariate=lsArgs$options$strMFAColor, strMFAShapeCovariate=lsArgs$options$strMFAShape, dMFAMetadataScale=lsArgs$options$dMFAMetadataScale, dMFADataScale=lsArgs$options$dMFADataScale, lsPlotFeatures=lsFeaturesToPlot, fPlotNA=lsArgs$options$fPlotNA, fInvert=lsArgs$options$fInvert, tempSaveFileName=file.path(outputDirectory,strBase), funcPlotColors=lsRet$funcPlotColors, funcPlotPoints=lsRet$funcPlotPoints, funcPlotLegend=lsRet$funcPlotLegend )
      logdebug("PlotMFA:out", c_logrMaaslin)
    }
  }
#)

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
