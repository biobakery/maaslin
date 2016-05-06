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


  #######################################################################################################
  #  The following changes were introduced  
  #   Added  the 'tools' library     - needed to find the suffix of a dataset name        
  #  George Weingart   george.weingart@gmail.com   4/28/15                  
  #######################################################################################################
 
### Load packages
vDepLibrary = c("agricolae", "gam", "gamlss", "gbm", "glmnet", "inlinedocs", "logging", "MASS", "nlme", "optparse", "outliers", "penalized", "pscl", "robustbase", "tools")
for(sDepLibrary in vDepLibrary)
{
  if(! require(sDepLibrary, character.only=TRUE) )
  {
    stop(paste("Please install the required package:",sDepLibrary,sep=" "))
  }
}

### Class for commandline argument processing
suppressMessages(library( optparse, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))

### Set the default run options
lsArgs = list()					 
lsArgs$strInputConfig = NULL
lsArgs$strInputR = NULL
lsArgs$dSignificanceLevel =  0.25	 
lsArgs$dMinAbd = 0.0001
lsArgs$dMinSamp = 0.1
lsArgs$dOutlierFence = 0
lsArgs$dPOutlier = 0.05	 
lsArgs$strRandomCovariates = NULL
lsArgs$strMultTestCorrection = "BH" 
lsArgs$fZeroInflated = FALSE
lsArgs$strModelSelection = "boost" 
lsArgs$strMethod = "lm"
lsArgs$strTransform = "asinsqrt"
lsArgs$fNoQC = FALSE
lsArgs$strForcedPredictors = NULL
lsArgs$strNoImpute = NULL
lsArgs$strVerbosity = "DEBUG"
lsArgs$fOmitLogFile = FALSE
lsArgs$fInvert = FALSE
lsArgs$dSelectionFrequency = NA
lsArgs$fAllvAll = FALSE
lsArgs$fPlotNA = FALSE
lsArgs$dPenalizedAlpha = 0.95
lsArgs$sAlternativeLibraryLocation = NULL
lsArgs$iLastMetadata = NULL  # Added by GW on 2015/04/28 to support lastmetadata row for pcl files
lsArgs$iFirstMetadata = NULL  # Added by GW on 2015/04/28 to support lastmetadata row for pcl files

### Create command line argument parser
pArgs <- OptionParser( usage = "%prog [options] <data.tsv> <outputdir>" )

# Input files for MaAsLin
## Data configuration file
pArgs <- add_option( pArgs, c("-i", "--input_config"), type="character", action="store", dest="strInputConfig", default=lsArgs$strInputConfig, metavar="data.read.config", help="Optional configuration file describing data input format.")
## Data manipulation/normalization file
pArgs <- add_option( pArgs, c("-I", "--input_process"), type="character", action="store", dest="strInputR", default=lsArgs$strInputR, metavar="data.R", help="Optional configuration script normalizing or processing data.")

# Settings for MaAsLin
## Maximum false discovery rate
pArgs <- add_option( pArgs, c("-d", "--fdr"), type="double", action="store", dest="dSignificanceLevel", default=lsArgs$dSignificanceLevel, metavar="significance", help="The threshold to use for significance for the generated q-values (BH FDR). Anything equal to or lower than this is significant.  [Default %default]")
## Minimum feature relative abundance filtering
pArgs <- add_option( pArgs, c("-r", "--minRelativeAbundance"), type="double", action="store", dest="dMinAbd", default=lsArgs$dMinAbd, metavar="minRelativeAbundance", help="The minimum relative abundance allowed in the data. Values below this are removed and imputed as the median of the sample data.  [Default %default]")
## Minimum feature prevalence filtering
pArgs <- add_option( pArgs, c("-p", "--minPrevalence"), type="double", action="store", dest="dMinSamp", default=lsArgs$dMinSamp, metavar="minPrevalence", help="The minimum percentage of samples in which a feature must have the minimum relative abundance in order not to be removed. Also this is the maximum percentage of samples for which a metadata can have NAs before being removed.  [Default %default]")
## Fence for outlier, if not set Grubbs test is used
pArgs <- add_option( pArgs, c("-o", "--outlierFence"), type="double", action="store", dest="dOutlierFence", default=lsArgs$dOutlierFence, metavar="outlierFence", help="Outliers are defined as this number times the interquartile range added/subtracted from the 3rd/1st quartiles respectively. If set to 0 (default), outliers are defined by the Grubbs test.  [Default %default]")
## Significance for Grubbs test
pArgs <- add_option(pArgs, c("-G","--grubbsSig"), type="double", action="store", dest="dPOutlier", default=lsArgs$dPOutlier, metavar="grubbsAlpha", help="This is the significance cuttoff used to indicate an outlier or not. The closer to zero, the more significant an outlier must be to be removed.  [Default %default]")
## Fixed (not random) covariates
pArgs <- add_option( pArgs, c("-R","--random"), type="character", action="store", dest="strRandomCovariates", default=lsArgs$strRandomCovariates, metavar="fixed", help="These metadata will be treated as random covariates. Comma delimited data feature names. These features must be listed in the read.config file. Example '-R RandomMetadata1,RandomMetadata2'.  [Default %default]")
## Change the type of correction fo rmultiple corrections
pArgs <- add_option( pArgs, c("-T","--testingCorrection"), type="character", action="store", dest="strMultTestCorrection", default=lsArgs$strMultTestCorrection, metavar="multipleTestingCorrection", help="This indicates which multiple hypothesis testing method will be used, available are holm, hochberg, hommel, bonferroni, BH, BY.  [Default %default]")
## Use a zero inflated model of the inference method indicate in -m
pArgs <- add_option( pArgs, c("-z","--doZeroInflated"), type="logical", action="store_true", default = lsArgs$fZeroInflated, dest="fZeroInflated", metavar="fZeroInflated", help="If true, the zero inflated version of the inference model indicated in -m is used. For instance if using lm, zero-inflated regression on a gaussian distribution is used.  [Default %default].")

# Arguments used in validation of MaAsLin
## Model selection (enumerate) c("none","boost","penalized","forward","backward")
pArgs <- add_option( pArgs, c("-s", "--selection"), type="character", action="store", dest="strModelSelection", default=lsArgs$strModelSelection, metavar="model_selection", help="Indicates which of the variable selection techniques to use.  [Default %default]")
## Argument indicating which method should be ran (enumerate) c("univariate","lm","neg_binomial","quasi")
pArgs <- add_option( pArgs, c("-m", "--method"), type="character", action="store", dest="strMethod", default=lsArgs$strMethod, metavar="analysis_method", help="Indicates which of the statistical inference methods to run.  [Default %default]")
## Argument indicating which link function is used c("none","asinsqrt")
pArgs <- add_option( pArgs, c("-l", "--link"), type="character", action="store", dest="strTransform", default=lsArgs$strTransform, metavar="transform_method", help="Indicates which link or transformation to use with a glm, if glm is not selected this argument will be set to none.  [Default %default]")
pArgs <- add_option( pArgs, c("-Q","--NoQC"), type="logical", action="store_true", default=lsArgs$fNoQC, dest="fNoQC", metavar="Do_Not_Run_QC", help="Indicates if the quality control will be ran on the metadata/data. Default is true.  [Default %default]")

# Arguments to suppress MaAsLin actions on certain data
## Do not perform model selection on the following data
pArgs <- add_option( pArgs, c("-F","--forced"), type="character", action="store", dest="strForcedPredictors", default=lsArgs$strForcedPredictors, metavar="forced_predictors", help="Metadata features that will be forced into the model seperated by commas. These features must be listed in the read.config file. Example '-F Metadata2,Metadata6,Metadata10'.  [Default %default]")
## Do not impute the following
pArgs <- add_option( pArgs, c("-n","--noImpute"), type="character", action="store", dest="strNoImpute", default=lsArgs$strNoImpute, metavar="no_impute", help="These data will not be imputed. Comma delimited data feature names. Example '-n Feature1,Feature4,Feature6'.  [Default %default]")

#Miscellaneouse arguments
### Argument to control logging (enumerate)
pArgs <- add_option( pArgs, c("-v", "--verbosity"), type="character", action="store", dest="strVerbosity", default=lsArgs$strVerbosity, metavar="verbosity", help="Logging verbosity  [Default %default]")
### Run maaslin without creating a log file
pArgs <- add_option( pArgs, c("-O","--omitLogFile"), type="logical", action="store_true", default=lsArgs$fOmitLogFile, dest="fOmitLogFile", metavar="omitlogfile",help="Including this flag will stop the creation of the output log file.  [Default %default]")
### Argument for inverting background to black
pArgs <- add_option( pArgs, c("-t", "--invert"), type="logical", action="store_true", dest="fInvert", default=lsArgs$fInvert, metavar="invert", help="When given, flag indicates to invert the background of figures to black.  [Default %default]")
### Selection Frequency
pArgs <- add_option( pArgs, c("-f","--selectionFrequency"), type="double", action="store", dest="dSelectionFrequency", default=lsArgs$dSelectionFrequency, metavar="selectionFrequency", help="Selection Frequency for boosting (max 1 will remove almost everything). Interpreted as requiring boosting to select metadata 100% percent of the time (or less if given a number that is less). Value should be between 1 (100%) and 0 (0%), NA (default is determined by data size).")
### All v All
pArgs <- add_option( pArgs, c("-a","--allvall"), type="logical", action="store_true", dest="fAllvAll", default=lsArgs$fAllvAll, metavar="compare_all", help="When given, the flag indicates that each fixed covariate that is not indicated as Forced is compared once at a time per data feature (bug). Made to be used with the -F option to specify one part of the model while allowing the other to cycle through a group of covariates. Does not affect Random covariates, which are always included when specified.  [Default %default]")
pArgs <- add_option( pArgs, c("-N","--PlotNA"), type="logical", action="store_true", default=lsArgs$fPlotNA, dest="fPlotNA", metavar="plotNAs",help="Plot data that was originally NA, by default they are not plotted.  [Default %default]")
### Alternative methodology settings
pArgs <- add_option( pArgs, c("-A","--pAlpha"), type="double", action="store", dest="dPenalizedAlpha", default=lsArgs$dPenalizedAlpha, metavar="PenalizedAlpha",help="The alpha for penalization (1.0=L1 regularization, LASSO; 0.0=L2 regularization, ridge regression.  [Default %default]")
### Pass an alternative library dir
pArgs <- add_option( pArgs, c("-L", "--libdir"), action="store", dest="sAlternativeLibraryLocation", default=lsArgs$sAlternativeLibraryLocation, metavar="AlternativeLibraryDirectory", help="An alternative location to find the lib directory. This dir and children will be searched for the first maaslin/src/lib dir.")


## Last Metadata Row    #Added by GW on 2015/04/28  to support last metadata row if the input was pcl
pArgs <- add_option( pArgs, c("--lastMetadata" ), type="integer", action="store", dest="iLastMetadata", default=lsArgs$iLastMetadata, metavar="LastMetadata", help="Last metadata row in the pcl input file (Number of the row) or Last metadata column in a tsv file.  [Default %default]")
## First Metadata Row   #Added by GW on 2015/04/28  to support last metadata row if the input was pcl
pArgs <- add_option( pArgs, c("--firstMetadata" ), type="integer", action="store", dest="iFirstMetadata", default=lsArgs$iFirstMetadata, metavar="FirstMetadata", help="First metadata row in the pcl input file (Number of the row) or first metadata column in a tsv file.  [Default %default]")



Maaslin <- function(
### The main function manages the following:
### 1. Arguments are checked
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
strInputTSV,
strOutputDIR,
strInputConfig=NULL,
strInputR = NULL,
dSignificanceLevel = 0.25,	
dMinAbd = 0.0001,
dMinSamp = 0.1,
dOutlierFence = 0,
dPOutlier = 0.05,
strRandomCovariates = NULL,
strMultTestCorrection = "BH",
fZeroInflated = FALSE,
strModelSelection = "boost",
strMethod = "lm",
strTransform = "asinsqrt",
fNoQC = FALSE,
strForcedPredictors = NULL,
strNoImpute = NULL,
strVerbosity = "DEBUG",
fOmitLogFile = FALSE,
fInvert = FALSE,
dSelectionFrequency = NA,
fAllvAll = FALSE,
fPlotNA = FALSE,
dPenalizedAlpha = 0.95,
sAlternativeLibraryLocation = NULL,
iLastMetadata = NULL,    #Added by GW on 2015/04/28  to support last metadata row if the input was pcl
iFirstMetadata = NULL)   #Added by GW on 2015/04/28  to support first metadata row if the input was pcl
{


  ### Logging class
  suppressMessages(library( logging, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))

  # Parse parameters
  lsForcedParameters = NULL
  if(!is.null(strForcedPredictors))
  {
    lsForcedParameters  = unlist(strsplit(strForcedPredictors,","))
  }
  xNoImpute = NULL
  if(!is.null(strNoImpute))
  {
    xNoImpute = unlist(strsplit(strNoImpute,"[,]"))
  }
  lsRandomCovariates = NULL
  if(!is.null(strRandomCovariates))
  {
    lsRandomCovariates = unlist(strsplit(strRandomCovariates,"[,]"))
  }
  
  #If logging is not an allowable value, inform user and set to INFO
  if(length(intersect(names(loglevels), c(strVerbosity))) == 0)
  {
    print(paste("Maaslin::Error. Did not understand the value given for logging, please use any of the following: DEBUG,INFO,WARN,ERROR."))
    print(paste("Maaslin::Warning. Setting logging value to \"",strVerbosity,"\"."))
  }

  
  # Do not allow  mixed effect models and zero inflated models, don't have implemented
  if(fZeroInflated && !is.null(strRandomCovariates))
  {
    stop("MaAsLin Error:: The combination of zero inflated models and mixed effects models are not supported.")
  }

  ### Create logger
  c_logrMaaslin <- getLogger( "maaslin" )
  addHandler( writeToConsole, c_logrMaaslin )
  setLevel( strVerbosity, c_logrMaaslin )

  # Get analysis method options
  # includes data transformations, model selection/regularization, regression models/links
  strModelSelection = tolower(strModelSelection)
  if(!strModelSelection %in% c("none","boost","penalized","forward","backward"))
  {
    logerror(paste("Received an invalid value for the selection argument, received '",strModelSelection,"'"), c_logrMaaslin)
    stop( print_help( pArgs ) )
  }
  strMethod = tolower(strMethod)
  if(!strMethod %in% c("univariate","lm","neg_binomial","quasi"))
  {
    logerror(paste("Received an invalid value for the method argument, received '",strMethod,"'"), c_logrMaaslin)
    stop( print_help( pArgs ) )
  }
  strTransform = tolower(strTransform)
  if(!strTransform %in% c("none","asinsqrt"))
  {
    logerror(paste("Received an invalid value for the transform/link argument, received '",strTransform,"'"), c_logrMaaslin)
    stop( print_help( pArgs ) )
  }

  if(!strMultTestCorrection %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY"))
  {
    logerror(paste("Received an invalid value for the multiple testing correction argument, received '",strMultTestCorrection,"'"), c_logrMaaslin)
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
  #######################################################################################################
  #  The following chages were introduced in order to be able to run under R in the following manner:
  #  source('./R/Maslin.R')
  #  main(pArgs)
  #  George Weingart   george.weingart@gmail.com   2/17/2015                    
  #######################################################################################################
  strDir = file.path( dirname( script.name )  )
  
  

  if (identical(strDir, character(0)))                   #if Running under R with no override
	{strDir <- c("./R")}

  
 
  #######################################################################################################
  #  End of changes  George Weingart   george.weingart@gmail.com   2/17/2015                    
  #######################################################################################################
 
  strSelf = basename( script.name )   
  
 
  
  
  if (identical(strSelf, character(0)))                   #if Running under R with no override - force it 
	{strSelf <- c("Maaslin.R")}
	
  for( strR in dir( strDir, pattern = "*.R$" ) )
  {
    if( strR == strSelf ) {next}
    source( file.path( strDir, strR ) )
  }

  # Get analysis modules
  afuncVariableAnalysis = funcGetAnalysisMethods(strModelSelection,strTransform,strMethod,fZeroInflated)

  # Set up parameters for variable selection
  lxParameters = list(dFreq=dSelectionFrequency, dPAlpha=dPenalizedAlpha)
  if((strMethod == "lm")||(strMethod == "univariate"))
  { lxParameters$sFamily = "gaussian"
  } else if(strMethod == "neg_binomial"){ lxParameters$sFamily = "binomial"
  } else if(strMethod == "quasi"){ lxParameters$sFamily = "poisson"}

  #Indicate start
  logdebug("Start MaAsLin", c_logrMaaslin)

  ### Output directory for the study based on the requested output file
  
if (nchar(strOutputDIR) == 0)
{print(paste("No output directory specified. Files will be logged to the current working directory"))
outputDirectory = getwd()
} else if (file.exists(strOutputDIR)){
outputDirectory = strOutputDIR
    
} else {
    dir.create(strOutputDIR)
   outputDirectory = strOutputDIR 

}

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
  ##################################################################################
  #   Modification log                                                             #
  #   Check if the input file is of type pcl - if so, transpose it                 #
  #   George Weingart  george.weingart@gmail.com  4/27/15                          #
  ##################################################################################
  
 
 
  fConvertPCLtoTSV  = FALSE
  fGenerateConfigFile  = FALSE
  OriginalstrInputTSV = strInputTSV  #  Store the original input filename
  inputFileDataSuffix  =  file_ext(strInputTSV)
  if  (inputFileDataSuffix  == "pcl")
	{
		strInputTSV = funcTransposeInputPCLtoTSV(strInputTSV, strOutputDIR)
	    fConvertPCLtoTSV  = TRUE

	}
	
 
  if  ( is.null(strInputConfig)  &&   !is.null(iLastMetadata))   #  If the User did not provide config file but provided Lastmetadata - we will build config file
	{
		if (is.null(iFirstMetadata))  # If user provided last metadata row but not first, we will try 2.....
		{
			iFirstMetadata = 2
		}
		strParmColIndeces1 = paste(toString(iFirstMetadata),"-",toString(iLastMetadata),sep="")
		funcWriteMatrixToReadConfigFile(strConfigureFileName='generated_config', strMatrixName='Metadata',  strColIndices=strParmColIndeces1 )
		strParmColIndeces2 = paste( toString(iLastMetadata + 1),'-',sep="")
		funcWriteMatrixToReadConfigFile(strConfigureFileName='generated_config', strMatrixName='Abundance',  strColIndices=strParmColIndeces2 ,fAppend=TRUE)
		strInputConfig = "generated_config"
		fGenerateConfigFile  = TRUE
	}
	
  ##################################################################################
  #   End Modification log                                                         #
  #   George Weingart  george.weingart@gmail.com  4/27/15                          #
  ##################################################################################
  
  

  
  
  inputFileData = funcReadMatrices(strInputConfig, strInputTSV, log=TRUE)
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
  if(!is.null(funcSourceScript(strInputR))){funcProcess <- get(c_strCustomProcessFunction)}

  #Clean the data and update the current data list to the cleaned data list
  funcTransformData = afuncVariableAnalysis[[c_iTransform]]
  lsQCCounts = list(aiDataCleaned = c(), aiMetadataCleaned = c())
  lsRet = list(frmeData=frmeData, aiData=aiData, aiMetadata=aiMetadata, lsQCCounts=lsQCCounts, liNaIndices=c())

  viNotTransformedDataIndices = c()
  if(!fNoQC)
  {
    c_logrMaaslin$info( "Running quality control." )
    lsRet = funcClean( frmeData=frmeData, funcDataProcess=funcProcess, aiMetadata=aiMetadata, aiData=aiData, lsQCCounts=lsData$lsQCCounts, astrNoImpute=xNoImpute, dMinSamp = dMinSamp, dMinAbd = dMinAbd, dFence=dOutlierFence, funcTransform=funcTransformData, dPOutlier=dPOutlier)

    viNotTransformedDataIndices = lsRet$viNotTransformedData

    #If using a count based model make sure all are integer (QCing can add in numeric values during interpolation for example)
    if(strMethod %in% c_vCountBasedModels)
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
  if(fPlotNA)
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
  fDoRPlot=FALSE
  #Should not occur for univariates
  if(strMethod %in% c("univariate")){ fDoRPlot=FALSE }

  #Run analysis
  alsRetBugs = funcBugs( frmeData=lsRet$frmeData, lsData=lsRet, aiMetadata=lsRet$aiMetadata, aiData=lsRet$aiData, aiNotTransformedData=viNotTransformedDataIndices, strData=strBase, dSig=dSignificanceLevel, fInvert=fInvert,
        strDirOut=outputDirectory, funcReg=afuncVariableAnalysis[[c_iSelection]], funcTransform=funcTransformData, funcUnTransform=afuncVariableAnalysis[[c_iUnTransform]], lsNonPenalizedPredictors=lsForcedParameters,
        funcAnalysis=afuncVariableAnalysis[[c_iAnalysis]], lsRandomCovariates=lsRandomCovariates, funcGetResults=afuncVariableAnalysis[[c_iResults]], fDoRPlot=fDoRPlot, fOmitLogFile=fOmitLogFile,
        fAllvAll=fAllvAll, liNaIndices=lsRet$liNaIndices, lxParameters=lxParameters, strTestingCorrection=strMultTestCorrection, 
        fIsUnivariate=afuncVariableAnalysis[[c_iIsUnivariate]], fZeroInflated=fZeroInflated )

	
	
	

		
  #Write QC files only in certain modes of verbosity
  if( c_logrMaaslin$level <= loglevels["DEBUG"] ) {
	funcWriteQCReport(strProcessFileName=file.path(strQCDir,"ProcessQC.txt"), lsQCData=alsRetBugs$lsQCCounts, liDataDim=liData, liMetadataDim=liMetaData)

    ### Write out the parameters used in the run
    unlink(file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite("Parameters used in the MaAsLin run", file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Optional input read.config file=",strInputConfig), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Optional R file=",strInputR), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("FDR threshold for pdf generation=",dSignificanceLevel), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Minimum relative abundance=",dMinAbd), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Minimum percentage of samples with measurements=",dMinSamp), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("The fence used to define outliers with a quantile based analysis. If set to 0, the Grubbs test was used=",dOutlierFence), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Ignore if the Grubbs test was not used. The significance level used as a cut-off to define outliers=",dPOutlier), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("These covariates are treated as random covariates and not fixed covariates=",strRandomCovariates), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("The type of multiple testing correction used=",strMultTestCorrection), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Zero inflated inference models were turned on=",fZeroInflated), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Feature selection step=",strModelSelection), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Statistical inference step=",strMethod), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Numeric transform used=",strTransform), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Quality control was run=",!fNoQC), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("These covariates were forced into each model=",strForcedPredictors), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("These features' data were not changed by QC processes=",strNoImpute), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Output verbosity=",strVerbosity), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Log file was generated=",!fOmitLogFile), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Data plots were inverted=",fInvert), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Ignore unless boosting was used. The threshold for the rel.inf used to select features=",dSelectionFrequency), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("All verses all inference method was used=",fAllvAll), file.path(strQCDir,"Run_Parameters.txt"))
    funcWrite(paste("Ignore unless penalized feature selection was used. Alpha to determine the type of penalty=",dPenalizedAlpha), file.path(strQCDir,"Run_Parameters.txt"))
  }
  ##################################################################################
  #   Modification log                                                             #
  #   If the input was a pcl file, notif the User we transposed it and             #
  #   move the transposed tsv file to the output directory                         #
  #   George Weingart  george.weingart@gmail.com  4/27/15                          #
  ##################################################################################
  	if (fConvertPCLtoTSV  ==  TRUE) 
		{
		    strTSVFIleNewLocation = file.path(strOutputDIR,basename(strInputTSV))  #In the output directory
			funcWrite(paste("Input pcl file: ",OriginalstrInputTSV , "transposed and converted to tsv format and stored in the output directory (with a .tsv suffix)" ),file.path(strQCDir,"Run_Parameters.txt"))
			file.rename(from=strInputTSV,to=strTSVFIleNewLocation)  # Move the new TSV file to the Output directory
		}
  if (fGenerateConfigFile  == TRUE)  #If we generated a config file - notify and move it to Output directory
				{
					strConfigFIleNewLocation = file.path(strOutputDIR,basename(strInputConfig))  #In the output directory
					file.rename(from=strInputConfig,to=strConfigFIleNewLocation)  # Move the new TSV file to the Output directory
					funcWrite(paste("Generated config file in: ",strConfigFIleNewLocation , " - You can use as a sample and create your own and pass it using the parm strInputConfig=YourConfigFile" ),file.path(strQCDir,"Run_Parameters.txt"))

				}
  
  ##################################################################################
  #   End Modification log                                                         #
  #   If the input was a pcl file, notif the User we transposed it and             #
  #   move the transposed tsv file to the output directory                         #
  #   George Weingart  george.weingart@gmail.com  4/27/15                          #
  ##################################################################################
  
  
  ### Write summary table
  # Summarize output files based on a keyword and a significance threshold
  # Look for less than or equal to the threshold (appropriate for p-value and q-value type measurements)
  # DfSummary is sorted by the q.value when it is returned
  dfSummary = funcSummarizeDirectory(astrOutputDirectory=outputDirectory,
                       strBaseName=strBase,
                       astrSummaryFileName=file.path(outputDirectory,paste(strBase,c_sSummaryFileSuffix, sep="")), 
                       astrKeyword=c_strKeywordEvaluatedForInclusion, 
                       afSignificanceLevel=dSignificanceLevel)
}


if( identical( environment( ), globalenv( ) ) &&
	!length( grep( "^source\\(", sys.calls( ) ) ) ) {
	cArgs <- parse_args( pArgs, positional_arguments = TRUE)
        # check for the correct number of positional arguments
    if(length(cArgs$args)!= 2) {
            print_help(pArgs)
            stop("Please provide an input data file ( <data.tsv> ) and an output directory ( <outputdir> ).\n\nUsage: Maaslin.R [options] <data.tsv> <outputdir>")
        }
	Maaslin( cArgs$args[1], cArgs$args[2], 
        strInputConfig=cArgs$options$strInputConfig,
        strInputR = cArgs$options$strInputR,
        dSignificanceLevel = cArgs$options$dSignificanceLevel,	
        dMinAbd = cArgs$options$dMinAbd,
        dMinSamp = cArgs$options$dMinSamp,
        dOutlierFence = cArgs$options$dOutlierFence,
        dPOutlier = cArgs$options$dPOutlier, 
        strRandomCovariates = cArgs$options$strRandomCovariates,
        strMultTestCorrection = cArgs$options$strMultTestCorrection,
        fZeroInflated = cArgs$options$fZeroInflated,
        strModelSelection = cArgs$options$strModelSelection, 
        strMethod = cArgs$options$strMethod,
        strTransform = cArgs$options$strTransform,
        fNoQC = cArgs$options$fNoQC,
        strForcedPredictors = cArgs$options$strForcedPredictors,
        strNoImpute = cArgs$options$strNoImpute,
        strVerbosity = cArgs$options$strVerbosity,
        fOmitLogFile = cArgs$options$fOmitLogFile,
        fInvert = cArgs$options$fInvert,
        dSelectionFrequency = cArgs$options$dSelectionFrequency,
        fAllvAll = cArgs$options$fAllvAll,
        fPlotNA = cArgs$options$fPlotNA,
        dPenalizedAlpha = cArgs$options$dPenalizedAlpha,
        sAlternativeLibraryLocation = cArgs$options$sAlternativeLibraryLocation,
		iLastMetadata = cArgs$options$iLastMetadata,  #Added by GW on 2015/04/28 to support Lastmetadata if pcl file was provided
		iFirstMetadata = cArgs$options$iFirstMetadata  #Added by GW on 2015/04/28 to support Lastmetadata if pcl file was provided
		) 
}
 
