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


Maaslinfn <- function(
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
strInputTSV,
strOutputDIR,
strInputConfig=lsArgs$strInputConfig,
strInputR = lsArgs$strInputR,
dSignificanceLevel = lsArgs$dSignificanceLevel,	
dMinAbd = lsArgs$dMinAbd,
dMinSamp = lsArgs$dMinSamp,
dOutlierFence = lsArgs$dOutlierFence,
dPOutlier = lsArgs$dPOutlier,
strRandomCovariates = lsArgs$strRandomCovariates,
strMultTestCorrection = lsArgs$strMultTestCorrection,
fZeroInflated = lsArgs$fZeroInflated,
strModelSelection = lsArgs$strModelSelection,
strMethod = lsArgs$strMethod,
strTransform = lsArgs$strTransform,
fNoQC = lsArgs$fNoQC,
strForcedPredictors = lsArgs$strForcedPredictors,
strNoImpute = lsArgs$strNoImpute,
strVerbosity = lsArgs$strVerbosity,
fOmitLogFile = lsArgs$fOmitLogFile,
fInvert = lsArgs$fInvert,
dSelectionFrequency = lsArgs$dSelectionFrequency,
fAllvAll = lsArgs$fAllvAll,
fPlotNA = lsArgs$fPlotNA,
dPenalizedAlpha = lsArgs$dPenalizedAlpha,
sAlternativeLibraryLocation = lsArgs$sAlternativeLibraryLocation)
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
  strDir = file.path( dirname( script.name ), "lib" )
  
  
  #######################################################################################################
  #  The following chages were introduced in order to be able to run under R in the following manner:
  #  source('./R/Maslin.R')
  #  main(pArgs)
  #  George Weingart   george.weingart@gmail.com   1/1/2015                    
  #######################################################################################################
  if (identical(strDir, character(0)))                   #if Running under R with no override
	{strDir <- c("./R/lib")}
  #######################################################################################################
  #  End of changes  George Weingart   george.weingart@gmail.com   1/1/2015                    
  #######################################################################################################
  
  # If this does not have the lib file then go for the alt lib
  if( !file.exists(strDir) )
  {
    lsPotentialListLocations = dir( path = sAlternativeLibraryLocation, pattern = "lib", recursive = TRUE, include.dirs = TRUE)
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
          strDir = file.path( sAlternativeLibraryLocation, strSearchDir )
          break
        }
      }
    }
  }

 
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
  fDoRPlot=TRUE
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
 
 