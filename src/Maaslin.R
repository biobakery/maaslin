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


### Install packages if not already installed
vDepLibrary = c("agricolae", "gam", "gamlss", "gbm", "glmnet", "inlinedocs", "logging", "MASS", "nlme", "optparse", "outliers", "penalized", "pscl", "robustbase", "testthat")
for(sDepLibrary in vDepLibrary)
{
  if(! require(sDepLibrary, character.only=TRUE) )
  {
    install.packages(pkgs=sDepLibrary, repos="http://cran.us.r-project.org")
  }
}

### Logging class
suppressMessages(library( logging, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
### Class for commandline argument processing
suppressMessages(library( optparse, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))


### Create command line argument parser
pArgs <- OptionParser( usage = "%prog [options] <output.txt> <data.tsv>" )

# Input files for MaAsLin
## Data configuration file
pArgs <- add_option( pArgs, c("-i", "--input_config"), type="character", action="store", dest="strInputConfig", metavar="data.read.config", help="Optional configuration file describing data input format.")
## Data manipulation/normalization file
pArgs <- add_option( pArgs, c("-I", "--input_process"), type="character", action="store", dest="strInputR", metavar="data.R", help="Optional configuration script normalizing or processing data.")

# Settings for MaAsLin
## Maximum false discovery rate
pArgs <- add_option( pArgs, c("-d", "--fdr"), type="double", action="store", dest="dSignificanceLevel", default=0.25, metavar="significance", help="The threshold to use for significance for the generated q-values (BH FDR). Anything equal to or lower than this is significant.  [Default %default]")
## Minimum feature relative abundance filtering
pArgs <- add_option( pArgs, c("-r", "--minRelativeAbundance"), type="double", action="store", dest="dMinAbd", default=0.0001, metavar="minRelativeAbundance", help="The minimum relative abundance allowed in the data. Values below this are removed and imputed as the median of the sample data.  [Default %default]")
## Minimum feature prevalence filtering
pArgs <- add_option( pArgs, c("-p", "--minPrevalence"), type="double", action="store", dest="dMinSamp", default=0.1, metavar="minPrevalence", help="The minimum percentage of samples a feature can have abundance in before being removed. Also is the minimum percentage of samples a metadata can have that are not NA before being removed.  [Default %default]")
## Fence for outlier, if not set Grubbs test is used
pArgs <- add_option( pArgs, c("-o", "--outlierFence"), type="double", action="store", dest="dOutlierFence", default=0, metavar="outlierFence", help="Outliers are defined as this number times the interquartile range added/subtracted from the 3rd/1st quartiles respectively. If set to 0 (default), outliers are defined by the Grubbs test.  [Default %default]")
## Significance for Grubbs test
pArgs <- add_option(pArgs, c("-G","--grubbsSig"), type="double", action="store", dest="dPOutlier", default=0.05, metavar="grubbsAlpha", help="This is the significance cuttoff used to indicate an outlier or not. The closer to zero, the more significant an outlier must be to be removed.  [Default %default]")
## Fixed (not random) covariates
pArgs <- add_option( pArgs, c("-R","--random"), type="character", action="store", dest="strRandomCovariates", default=NULL, metavar="fixed", help="These metadata will be treated as random covariates. Comma delimited data feature names. These features must be listed in the read.config file. Example '-R RandomMetadata1,RandomMetadata2'.  [Default %default]")
## Change the type of correction fo rmultiple corrections
pArgs <- add_option( pArgs, c("-T","--testingCorrection"), type="character", action="store", dest="strMultTestCorrection", default="BH", metavar="multipleTestingCorrection", help="This indicates which multiple hypothesis testing method will be used, available are holm, hochberg, hommel, bonferroni, BH, BY.  [Default %default]")
## Use a zero inflated model of the inference method indicate in -m
pArgs <- add_option( pArgs, c("-z","--doZeroInfated"), type="logical", action="store_true", default = FALSE, dest="fZeroInflated", metavar="fZeroInflated", help="If true, the zero inflated version of the inference model indicated in -m is used. For instance if using lm, zero-inflated regression on a gaussian distribution is used.  [Default %default].")

# Arguments used in validation of MaAsLin
## Model selection (enumerate) c("none","boost","penalized","forward","backward")
pArgs <- add_option( pArgs, c("-s", "--selection"), type="character", action="store", dest="strModelSelection", default="boost", metavar="model_selection", help="Indicates which of the variable selection techniques to use.  [Default %default]")
## Argument indicating which method should be ran (enumerate) c("univariate","lm","neg_binomial","quasi")
pArgs <- add_option( pArgs, c("-m", "--method"), type="character", action="store", dest="strMethod", default="lm", metavar="analysis_method", help="Indicates which of the statistical inference methods to run.  [Default %default]")
## Argument indicating which link function is used c("none","asinsqrt")
pArgs <- add_option( pArgs, c("-l", "--link"), type="character", action="store", dest="strTransform", default="asinsqrt", metavar="transform_method", help="Indicates which link or transformation to use with a glm, if glm is not selected this argument will be set to none.  [Default %default]")
pArgs <- add_option( pArgs, c("-Q","--NoQC"), type="logical", action="store_true", default=FALSE, dest="fNoQC", metavar="Do_Not_Run_QC", help="Indicates if the quality control will be ran on the metadata/data. Default is true.  [Default %default]")

# Arguments to suppress MaAsLin actions on certain data
## Do not perform model selection on the following data
pArgs <- add_option( pArgs, c("-F","--forced"), type="character", action="store", dest="strForcedPredictors", default=NULL, metavar="forced_predictors", help="Metadata features that will be forced into the model seperated by commas. These features must be listed in the read.config file. Example '-F Metadata2,Metadata6,Metadata10'.  [Default %default]")
## Do not impute the following
pArgs <- add_option( pArgs, c("-n","--noImpute"), type="character", action="store", dest="strNoImpute", default=NULL, metavar="no_impute", help="These data will not be imputed. Comma delimited data feature names. Example '-n Feature1,Feature4,Feature6'.  [Default %default]")

#Miscellaneouse arguments
### Argument to control logging (enumerate)
strDefaultLogging = "DEBUG"
pArgs <- add_option( pArgs, c("-v", "--verbosity"), type="character", action="store", dest="strVerbosity", default=strDefaultLogging, metavar="verbosity", help="Logging verbosity  [Default %default]")
### Run maaslin without creating a log file
pArgs <- add_option( pArgs, c("-O","--omitLogFile"), type="logical", action="store_true", default=FALSE, dest="fOmitLogFile", metavar="omitlogfile",help="Including this flag will stop the creation of the output log file.  [Default %default]")
### Argument for inverting background to black
pArgs <- add_option( pArgs, c("-t", "--invert"), type="logical", action="store_true", dest="fInvert", default=FALSE, metavar="invert", help="When given, flag indicates to invert the background of figures to black.  [Default %default]")
### Selection Frequency
pArgs <- add_option( pArgs, c("-f","--selectionFrequency"), type="double", action="store", dest="dSelectionFrequency", default=NA, metavar="selectionFrequency", help="Selection Frequency for boosting (max 1 will remove almost everything). Interpreted as requiring boosting to select metadata 100% percent of the time (or less if given a number that is less). Value should be between 1 (100%) and 0 (0%), NA (default is determined by data size).")
### All v All
pArgs <- add_option( pArgs, c("-a","--allvall"), type="logical", action="store_true", dest="fAllvAll", default=FALSE, metavar="compare_all", help="When given, the flag indicates that each fixed covariate that is not indicated as Forced is compared once at a time per data feature (bug). Made to be used with the -F option to specify one part of the model while allowing the other to cycle through a group of covariates. Does not affect Random covariates, which are always included when specified.  [Default %default]")
pArgs <- add_option( pArgs, c("-N","--PlotNA"), type="logical", action="store_true", default=FALSE, dest="fPlotNA", metavar="plotNAs",help="Plot data that was originally NA, by default they are not plotted.  [Default %default]")
### Alternative methodology settings
pArgs <- add_option( pArgs, c("-A","--pAlpha"), type="double", action="store", dest="dPenalizedAlpha", default=0.95, metavar="PenalizedAlpha",help="The alpha for penalization (1.0=L1 regularization, LASSO; 0.0=L2 regularization, ridge regression.  [Default %default]")
### Pass an alternative library dir
pArgs <- add_option( pArgs, c("-L", "--libdir"), action="store", dest="sAlternativeLibraryLocation", default=file.path( "","usr","share","biobakery" ), metavar="AlternativeLibraryDirectory", help="An alternative location to find the lib directory. This dir and children will be searched for the first maaslin/src/lib dir.")

#pArgs <- add_option( pArgs, c("-c","--MFAFeatureCount"), type="integer", action="store", dest="iMFAMaxFeatures", default=3, metavar="maxMFAFeature", help="Number of features or number of bugs to plot (default=3; 3 metadata and 3 data).")

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
### 10. An analysis method is performed on the model (optionally boosted model).
### 11. Data is summarized and PDFs are created for significant associations
### (those whose q-values {BH FDR correction} are <= the threshold given in the optional arguments.
pArgs
### Parsed commandline arguments
){
  lsArgs <- parse_args( pArgs, positional_arguments = TRUE )
  #logdebug("lsArgs", c_logrMaaslin)
  #logdebug(paste(lsArgs,sep=" "), c_logrMaaslin)

  # Parse parameters
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

  #If logging is not an allowable value, inform user and set to INFO
  if(length(intersect(names(loglevels), c(lsArgs$options$strVerbosity))) == 0)
  {
    print(paste("Maaslin::Error. Did not understand the value given for logging, please use any of the following: DEBUG,INFO,WARN,ERROR."))
    print(paste("Maaslin::Warning. Setting logging value to \"",strDefaultLogging,"\"."))
  }

  # Do not allow  mixed effect models and zero inflated models, don't have implemented
  if(lsArgs$options$fZeroInflated && !is.null(lsArgs$options$strRandomCovariates))
  {
    stop("MaAsLin Error:: The combination of zero inflated models and mixed effects models are not supported.")
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

  ### Necessary local import files
  ### Check to make sure the lib is in the expected place (where the script is)
  ### if not, then try the alternative lib location
  ### This will happen if, for instance the script is linked or
  ### on the path.
  # Get the first choice relative path
  initial.options <- commandArgs(trailingOnly = FALSE)
  script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
  strDir = file.path( dirname( script.name ), "lib" )
  # If this does not have the lib file then go for the alt lib
  if( !file.exists(strDir) )
  {
    lsPotentialListLocations = dir( path = lsArgs$options$sAlternativeLibraryLocation, pattern = "lib", recursive = TRUE, include.dirs = TRUE)
    if( length( lsPotentialListLocations ) > 0 )
    {
      sLibraryPath = file.path( "maaslin","src","lib" )
      iLibraryPathLength = nchar( sLibraryPath )
      for( strSearchDir in lsPotentialListLocations )
      {
        # Looking for the path where the end of the path is equal to the library path given earlier
        # Also checks before hand to make sure the path is atleast as long as the library path so no errors occur
        if ( substring( strSearchDir, 1 + nchar( strSearchDir ) - iLibraryPathLength ) == sLibraryPath )
        {
          strDir = file.path( lsArgs$options$sAlternativeLibraryLocation, strSearchDir )
          break
        }
      }
    }
  }
  
  strSelf = basename( script.name )
  for( strR in dir( strDir, pattern = "*.R$" ) )
  {
    if( strR == strSelf ) {next}
    source( file.path( strDir, strR ) )
  }

  # Get analysis modules
  afuncVariableAnalysis = funcGetAnalysisMethods(lsArgs$options$strModelSelection,lsArgs$options$strTransform,lsArgs$options$strMethod,lsArgs$options$fZeroInflated)

  # Set up parameters for variable selection
  lxParameters = list(dFreq=lsArgs$options$dSelectionFrequency, dPAlpha=lsArgs$options$dPenalizedAlpha)
  if((lsArgs$options$strMethod == "lm")||(lsArgs$options$strMethod == "univariate"))
  { lxParameters$sFamily = "gaussian"
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
  # Read in and merge files
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
  funcTransformData = afuncVariableAnalysis[[c_iTransform]]
  lsQCCounts = list(aiDataCleaned = c(), aiMetadataCleaned = c())
  lsRet = list(frmeData=frmeData, aiData=aiData, aiMetadata=aiMetadata, lsQCCounts=lsQCCounts, liNaIndices=c())

  viNotTransformedDataIndices = c()
  if(!lsArgs$options$fNoQC)
  {
    c_logrMaaslin$info( "Running quality control." )
    lsRet = funcClean( frmeData=frmeData, funcDataProcess=funcProcess, aiMetadata=aiMetadata, aiData=aiData, lsQCCounts=lsData$lsQCCounts, astrNoImpute=xNoImpute, dMinSamp = lsArgs$options$dMinSamp, dMinAbd = lsArgs$options$dMinAbd, dFence=lsArgs$options$dOutlierFence, funcTransform=funcTransformData, dPOutlier=lsArgs$options$dPOutlier)

    viNotTransformedDataIndices = lsRet$viNotTransformedData

    #If using a count based model make sure all are integer (QCing can add in numeric values during interpolation for example)
    if(lsArgs$options$strMethod %in% c_vCountBasedModels)
    {
      c_logrMaaslin$info( "Assuring the data matrix is integer." )
      for(iDataIndex in aiData)
      {
        lsRet$frmeData[ iDataIndex ] = round( lsRet$frmeData[ iDataIndex ] )
      }
    }
  } else {
    c_logrMaaslin$info( "Not running quality control, attempting transform." )
    ### Need to do transform if the QC is not performed
    iTransformed = 0
    for(iDataIndex in aiData)
    {
      if( ! funcTransformIncreasesOutliers( lsRet$frmeData[iDataIndex], funcTransformData ) )
      {
        lsRet$frmeData[iDataIndex]=funcTransformData(lsRet$frmeData[iDataIndex])
        iTransformed = iTransformed + 1
      } else {
        viNotTransformedDataIndices = c(viNotTransformedDataIndices, iDataIndex)
      }
    }
    c_logrMaaslin$info(paste("Number of features transformed = ", iTransformed))
  }

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
  alsRetBugs = funcBugs( frmeData=lsRet$frmeData, lsData=lsRet, aiMetadata=lsRet$aiMetadata, aiData=lsRet$aiData, aiNotTransformedData=viNotTransformedDataIndices, strData=strBase, dSig=lsArgs$options$dSignificanceLevel, fInvert=lsArgs$options$fInvert,
        strDirOut=outputDirectory, funcReg=afuncVariableAnalysis[[c_iSelection]], funcTransform=funcTransformData, funcUnTransform=afuncVariableAnalysis[[c_iUnTransform]], lsNonPenalizedPredictors=lsForcedParameters,
        funcAnalysis=afuncVariableAnalysis[[c_iAnalysis]], lsRandomCovariates=lsRandomCovariates, funcGetResults=afuncVariableAnalysis[[c_iResults]], fDoRPlot=fDoRPlot, fOmitLogFile=lsArgs$options$fOmitLogFile,
        fAllvAll=lsArgs$options$fAllvAll, liNaIndices=lsRet$liNaIndices, lxParameters=lxParameters, strTestingCorrection=lsArgs$options$strMultTestCorrection, 
        fIsUnivariate=afuncVariableAnalysis[[c_iIsUnivariate]], fZeroInflated=lsArgs$options$fZeroInflated )

  #Write QC files only in certain modes of verbosity
  if( c_logrMaaslin$level <= loglevels["DEBUG"] ) {
	funcWriteQCReport(strProcessFileName=file.path(strQCDir,"ProcessQC.txt"), lsQCData=alsRetBugs$lsQCCounts, liDataDim=liData, liMetadataDim=liMetaData)

    ### Write out the parameters used in the run
    unlink(file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite("Parameters used in the MaAsLin run", file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Optional input read.config file=",lsArgs$options$strInputConfig), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Optional R file=",lsArgs$options$strInputR), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("FDR threshold for pdf generation=",lsArgs$options$dSignificanceLevel), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Minimum relative abundance=",lsArgs$options$dMinAbd), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Minimum percentage of samples with measurements=",lsArgs$options$dMinSamp), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("The fence used to define outliers with a quantile based analysis. If set to 0, the Grubbs test was used=",lsArgs$options$dOutlierFence), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Ignore if the Grubbs test was not used. The significance level used as a cut-off to define outliers=",lsArgs$options$dPOutlier), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("These covariates are treated as random covariates and not fixed covariates=",lsArgs$options$strRandomCovariates), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("The type of multiple testing correction used=",lsArgs$options$strMultTestCorrection), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Zero inflated inference models were turned on=",lsArgs$options$fZeroInflated), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Feature selection step=",lsArgs$options$strModelSelection), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Statistical inference step=",lsArgs$options$strMethod), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Numeric transform used=",lsArgs$options$strTransform), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Quality control was run=",!lsArgs$options$fNoQC), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("These covariates were forced into each model=",lsArgs$options$strForcedPredictors), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("These features' data were not changed by QC processes=",lsArgs$options$strNoImpute), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Output verbosity=",lsArgs$options$strVerbosity), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Log file was generated=",!lsArgs$options$fOmitLogFile), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Data plots were inverted=",lsArgs$options$fInvert), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Ignore unless boosting was used. The threshold for the rel.inf used to select features=",lsArgs$options$dSelectionFrequency), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("All verses all inference method was used=",lsArgs$options$fAllvAll), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Ignore unless penalized feature selection was used. Alpha to determine the type of penalty=",lsArgs$options$dPenalizedAlpha), file.path(strQCDir,"Run_Parameters.txt"))
  }

  ### Write summary table
  # Summarize output files based on a keyword and a significance threshold
  # Look for less than or equal to the threshold (appropriate for p-value and q-value type measurements)
  # DfSummary is sorted by the q.value when it is returned
  dfSummary = funcSummarizeDirectory(astrOutputDirectory=outputDirectory,
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
