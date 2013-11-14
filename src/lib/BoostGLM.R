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
##description<< Manages the quality control of data and the performance of analysis (univariate or multivariate), regularization, and data (response) transformation.
) { return( pArgs ) }

### Load libraries quietly
suppressMessages(library( gam, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
suppressMessages(library( gbm, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
suppressMessages(library( logging, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
suppressMessages(library( outliers, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
suppressMessages(library( robustbase, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
suppressMessages(library( pscl, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))

### Get constants
#source(file.path("input","maaslin","src","Constants.R"))
#source("Constants.R")

## Get logger
c_logrMaaslin <- getLogger( "maaslin" )

funcDoGrubbs <- function(
### Use the Grubbs Test to identify outliers
iData,
### Column index in the data frame to test
frmeData,
### The data frame holding the data
dPOutlier,
### P-value threshold to indicate an outlier is significant
lsQC
### List holding the QC info of the cleaning step. Which indices are outliers is added.
){
  adData <- frmeData[,iData]

  # Original number of NA
  viNAOrig = which(is.na(adData))

  while( TRUE )
  {
    lsTest <- try( grubbs.test( adData ), silent = TRUE )
    if( ( class( lsTest ) == "try-error" ) || is.na( lsTest$p.value ) || ( lsTest$p.value > dPOutlier ) )
    {break}
    viOutliers = outlier( adData, logical = TRUE )
    adData[viOutliers] <- NA
  }

  # Record removed data
  viNAAfter = which(is.na(adData))

  # If all were set to NA then ignore the filtering
  if(length(adData)==length(viNAAfter))
  {
    viNAAfter = viNAOrig
    adData = frmeData[,iData]
    c_logrMaaslin$info( paste("Grubbs Test:: Identifed all data as outliers so was inactived for index=",iData," data=",paste(as.vector(frmeData[,iData]),collapse=","), "number zeros=", length(which(frmeData[,iData]==0)), sep = " " ))
  } else if(mean(adData, na.rm=TRUE) == 0) {
    viNAAfter = viNAOrig
    adData = frmeData[,iData]
    c_logrMaaslin$info( paste("Grubbs Test::Removed all values but 0, ignored. Index=",iData,".",sep=" " ) )
  } else {
    # Document removal
    if( sum( is.na( adData ) ) )
    {
      c_logrMaaslin$info( "Grubbs Test::Removing %d outliers from %s", sum( is.na( adData ) ), colnames(frmeData)[iData] )
			  c_logrMaaslin$info( format( rownames( frmeData )[is.na( adData )] ))
    }
  }

  return(list(data=adData,outliers=length(viNAAfter)-length(viNAOrig),indices=setdiff(viNAAfter,viNAOrig)))
}

funcDoFenceTest <- function(
### Use a threshold based on the quartiles of the data to identify outliers
iData,
### Column index in the data frame to test
frmeData,
### The data frame holding the data
dFence
### The fence outside the first and third quartiles to use as a threshold for cutt off.
### This many times the interquartile range +/- to the 3rd/1st quartiles
){
  # Establish fence
  adData <- frmeData[,iData]
  adQ <- quantile( adData, c(0.25, 0.5, 0.75), na.rm = TRUE )

  dIQR <- adQ[3] - adQ[1]
  if(!dIQR)
  {
    dIQR = sd(adData,na.rm = TRUE)
  }
  dUF <- adQ[3] + ( dFence * dIQR )
  dLF <- adQ[1] - ( dFence * dIQR )

  # Record indices of values outside of fence to remove and remove.
  aiRemove <- c()
  for( j in 1:length( adData ) )
  {
    d <- adData[j]
    if( !is.na( d ) && ( ( d < dLF ) || ( d > dUF ) ) )
    {
      aiRemove <- c(aiRemove, j)
    }
  }

  if(length(aiRemove)==length(adData))
  {
    aiRemove = c()
    c_logrMaaslin$info( "OutliersByFence:: Identified all data as outlier so was inactivated for index=", iData,"data=", paste(as.vector(frmeData[,iData]),collapse=","), "number zeros=", length(which(frmeData[,iData]==0)), sep=" " )
  } else {
    adData[aiRemove] <- NA

    # Document to screen
    if( length( aiRemove ) )
    {
      c_logrMaaslin$info( "OutliersByFence::Removing %d outliers from %s", length( aiRemove ), colnames(frmeData)[iData] )
      c_logrMaaslin$info( format( rownames( frmeData )[aiRemove] ))
    }
  }

  return(list(data=adData,outliers=length(aiRemove),indices=aiRemove))
}

funcZerosAreUneven = function(### vdRawData,### Raw data to be checked during transformationfuncTransform,### Data transform to performvsStratificationFeatures,
### Groupings to check for unevenness
dfData
### Data frame holding the features){
  # Return indicator of unevenness  fUneven = FALSE

  # Transform the data to compare  vdTransformed = funcTransform( vdRawData )

  # Go through each stratification of data  for( sStratification in vsStratificationFeatures )  {
    # Current stratification    vFactorStrats = dfData[[ sStratification ]]

    # If the metadata is not a factor then skip
    # Only binned data can be evaluated this way.
    if( !is.factor( vFactorStrats )){ next }
    
    viZerosCountsRaw = c()    for( sLevel in levels( vFactorStrats ) )    {
      vdTest = vdRawData[ which( vFactorStrats == sLevel ) ]
      viZerosCountsRaw = c( viZerosCountsRaw, length(which(vdTest == 0)))
      vdTest = vdTransformed[ which( vFactorStrats == sLevel ) ]    }
    dExpectation = 1 / length( viZerosCountsRaw )
    dMin = dExpectation / 2
    dMax = dExpectation + dMin
    viZerosCountsRaw = viZerosCountsRaw / sum( viZerosCountsRaw )
    if( ( length( which( viZerosCountsRaw <= dMin ) ) > 0 ) || ( length( which( viZerosCountsRaw >= dMax ) ) > 0 ) )
    {
      return( TRUE )
    }
  }  return( fUneven )}
funcTransformIncreasesOutliers = function(### Checks if a data transform increases outliers in a distributionvdRawData,
### Raw data to check for outlier zerosfuncTransform){  iUnOutliers = length( boxplot( vdRawData, plot = FALSE )$out )  iTransformedOutliers = length( boxplot( funcTransform( vdRawData ), plot = FALSE )$out )
  return( iUnOutliers <= iTransformedOutliers ) }

funcClean <- function(
### Properly clean / get data ready for analysis
### Includes custom analysis from the custom R script if it exists
frmeData,
### Data frame, input data to be acted on
funcDataProcess,
### Custom script that can be given to perform specialized processing before MaAsLin does.
aiMetadata,
### Indices of columns in frmeData which are metadata for analysis.
aiData,
### Indices of column in frmeData which are (abundance) data for analysis.
lsQCCounts,
### List that will hold the quality control information which is written in the output directory.
astrNoImpute = c(),
### An array of column names of frmeData not to impute.
dMinSamp,
### Minimum number of samples
dMinAbd,
# Minimum sample abundance
dFence,
### How many quartile ranges defines the fence to define outliers.
funcTransform,
### The data transformation function or a dummy function that does not affect the data
dPOutlier = 0.05
### The significance threshold for the grubbs test to identify an outlier.
){
  # Call the custom script and set current data and indicies to the processes data and indicies.
  c_logrMaaslin$debug( "Start Clean")
  if( !is.null( funcDataProcess ) )
  {
    c_logrMaaslin$debug("Additional preprocess function attempted.")

    pTmp <- funcDataProcess( frmeData=frmeData, aiMetadata=aiMetadata, aiData=aiData)
    frmeData = pTmp$frmeData
    aiMetadata = pTmp$aiMetadata
    aiData = pTmp$aiData
    lsQCCounts$lsQCCustom = pTmp$lsQCCounts
  }
  # Set data indicies after custom QC process.
  lsQCCounts$aiAfterPreprocess = aiData

  # Remove missing data, remove any sample that has less than dMinSamp * the number of data or low abundance
  aiRemove = c()
  aiRemoveLowAbundance = c()
  for( iCol in aiData )
  {
    adCol = frmeData[,iCol]
    adCol[!is.finite( adCol )] <- NA
    if( ( sum( !is.na( adCol ) ) < ( dMinSamp * length( adCol ) ) ) ||
      ( length( unique( na.omit( adCol ) ) ) < 2 ) )
    {
        aiRemove = c(aiRemove, iCol)
    }
    if( sum(adCol > dMinAbd, na.rm=TRUE ) < (dMinSamp * length( adCol)))
    {
      aiRemoveLowAbundance = c(aiRemoveLowAbundance, iCol)
    }
  }
  # Remove and document
  aiData = setdiff( aiData, aiRemove )
  aiData = setdiff( aiData, aiRemoveLowAbundance )
  lsQCCounts$iMissingData = aiRemove
  lsQCCounts$iLowAbundanceData = aiRemoveLowAbundance
  if(length(aiRemove))
  {
    c_logrMaaslin$info( "Removing the following for data lower bound.")
    c_logrMaaslin$info( format( colnames( frmeData )[aiRemove] ))
  }
  if(length(aiRemoveLowAbundance))
  {
    c_logrMaaslin$info( "Removing the following for too many low abundance bugs.")
    c_logrMaaslin$info( format( colnames( frmeData )[aiRemoveLowAbundance] ))
  }

  #Transform data
  iTransformed = 0
  viNotTransformedData = c()
  for(aiDatum in aiData)
  {
    adValues = frmeData[,aiDatum]
#    if( ! funcTransformIncreasesOutliers( adValues, funcTransform ) )
#    {
      frmeData[,aiDatum] = funcTransform( adValues )
#      iTransformed = iTransformed + 1
#    } else {
#      viNotTransformedData = c( viNotTransformedData, aiDatum )
#    }
  }
  c_logrMaaslin$info(paste("Number of features transformed = ",iTransformed))

  # Metadata: Properly factorize all logical data and integer and number data with less than iNonFactorLevelThreshold
  # Also record which are numeric metadata
  aiNumericMetadata = c()
  for( i in aiMetadata )
  {
    if( ( class( frmeData[,i] ) %in% c("integer", "numeric", "logical") ) &&
      ( length( unique( frmeData[,i] ) ) < c_iNonFactorLevelThreshold ) ) {
      c_logrMaaslin$debug(paste("Changing metadatum from numeric/integer/logical to factor",colnames(frmeData)[i],sep="="))
      frmeData[,i] = factor( frmeData[,i] )
    } 
    if( class( frmeData[,i] ) %in% c("integer","numeric") )
    {
      aiNumericMetadata = c(aiNumericMetadata,i)
    }
  }

  # Remove outliers
  # If the dFence Value is set use the method of defining the outllier as
  # dFence * the interquartile range + or - the 3rd and first quartile respectively.
  # If not the gibbs test is used.
  lsQCCounts$aiDataSumOutlierPerDatum = c()
  lsQCCounts$aiMetadataSumOutlierPerDatum = c()
  lsQCCounts$liOutliers = list()

  if( dFence > 0.0 )
  {
    # For data
    for( iData in aiData )
    {
      lOutlierInfo <- funcDoFenceTest(iData=iData,frmeData=frmeData,dFence=dFence)
      frmeData[,iData] <- lOutlierInfo[["data"]]
      lsQCCounts$aiDataSumOutlierPerDatum <- c(lsQCCounts$aiDataSumOutlierPerDatum,lOutlierInfo[["outliers"]])
      if(lOutlierInfo[["outliers"]]>0)
      {
        lsQCCounts$liOutliers[[paste(iData,sep="")]] <- lOutlierInfo[["indices"]]
      }
    }

    # Remove outlier non-factor metadata
    for( iMetadata in aiNumericMetadata )
    {
      lOutlierInfo <- funcDoFenceTest(iData=iMetadata,frmeData=frmeData,dFence=dFence)
      frmeData[,iMetadata] <- lOutlierInfo[["data"]]
      lsQCCounts$aiMetadataSumOutlierPerDatum <- c(lsQCCounts$aiMetadataSumOutlierPerDatum,lOutlierInfo[["outliers"]])
      if(lOutlierInfo[["outliers"]]>0)
      {
        lsQCCounts$liOutliers[[paste(iMetadata,sep="")]] <- lOutlierInfo[["indices"]]
      }
    }
  #Do not use the fence, use the Grubbs test
  } else if(dPOutlier!=0.0){
    # For data
    for( iData in aiData )
    {
      lOutlierInfo <- funcDoGrubbs(iData=iData,frmeData=frmeData,dPOutlier=dPOutlier)
      frmeData[,iData] <- lOutlierInfo[["data"]]
      lsQCCounts$aiDataSumOutlierPerDatum <- c(lsQCCounts$aiDataSumOutlierPerDatum,lOutlierInfo[["outliers"]])
      if(lOutlierInfo[["outliers"]]>0)
      {
        lsQCCounts$liOutliers[[paste(iData,sep="")]] <- lOutlierInfo[["indices"]]
      }
    }
    for( iMetadata in aiNumericMetadata )
    {
      lOutlierInfo <- funcDoGrubbs(iData=iMetadata,frmeData=frmeData,dPOutlier=dPOutlier)
      frmeData[,iMetadata] <- lOutlierInfo[["data"]]
      lsQCCounts$aiMetadataSumOutlierPerDatum <- c(lsQCCounts$aiMetadataSumOutlierPerDatum,lOutlierInfo[["outliers"]])
      if(lOutlierInfo[["outliers"]]>0)
      {
        lsQCCounts$liOutliers[[paste(iMetadata,sep="")]] <- lOutlierInfo[["indices"]]
      }
    }
  }

  # Metadata: Remove missing data
  # This is defined as if there is only one non-NA value or
  # if the number of NOT NA data is less than a percentage of the data defined by dMinSamp
  aiRemove = c()
  for( iCol in c(aiMetadata) )
  {
    adCol = frmeData[,iCol]
    if( ( sum( !is.na( adCol ) ) < ( dMinSamp * length( adCol ) ) ) ||
      ( length( unique( na.omit( adCol ) ) ) < 2 ) )
    {
      aiRemove = c(aiRemove, iCol)
    }
  }

  # Remove metadata
  aiMetadata = setdiff( aiMetadata, aiRemove )

  # Update the data which was removed.
  lsQCCounts$iMissingMetadata = aiRemove
  if(length(aiRemove))
  {
    c_logrMaaslin$info("Removing the following metadata for too much missing data or only one data value outside of NA.")
    c_logrMaaslin$info(format(colnames( frmeData )[aiRemove]))
  }

  # Keep track of factor levels in a list for later use
  lslsFactors <- list()
  for( iCol in c(aiMetadata) )
  {
    aCol <- frmeData[,iCol]
    if( class( aCol ) == "factor" )
    {
      lslsFactors[[length( lslsFactors ) + 1]] <- list(iCol, levels( aCol ))
    }
  }

  # Replace missing data values by the mean of the data column.
  # Remove samples that were all NA from the cleaning and so could not be imputed.
  aiRemoveData = c()
  for( iCol in aiData )
  {
    adCol <- frmeData[,iCol]
    adCol[is.infinite( adCol )] <- NA
    adCol[is.na( adCol )] <- mean( adCol[which(adCol>0)], na.rm = TRUE )
    frmeData[,iCol] <- adCol

    if(length(which(is.na(frmeData[,iCol]))) == length(frmeData[,iCol]))
    {
      print( paste("Removing data", iCol, "for being all NA after QC"))
      aiRemoveData = c(aiRemoveData,iCol)
    }
  }

  # Remove and document
  aiData = setdiff( aiData, aiRemoveData )
  lsQCCounts$iMissingData = c(lsQCCounts$iMissingData,aiRemoveData)
  if(length(aiRemoveData))
  {
    c_logrMaaslin$info( "Removing the following for having only NAs after cleaning (maybe due to only having NA after outlier testing).")
    c_logrMaaslin$info( format( colnames( frmeData )[aiRemoveData] ))
  }

  #Use na.gam.replace to manage NA metadata
  aiTmp <- setdiff( aiMetadata, which( colnames( frmeData ) %in% astrNoImpute ) )
  # Keep tack of NAs so the may not be plotted later.
  liNaIndices = list()
  lsNames = names(frmeData)
  for( i in aiTmp)
  {
    liNaIndices[[lsNames[i]]] = which(is.na(frmeData[,i]))
  }
  frmeData[,aiTmp] <- na.gam.replace( frmeData[,aiTmp] )

  #If NA is a value in factor data, set the NA as a level.
  for( lsFactor in lslsFactors )
  {
    iCol <- lsFactor[[1]]
    aCol <- frmeData[,iCol]
    if( "NA" %in% levels( aCol ) )
    {
      if(! lsNames[iCol] %in% astrNoImpute)
      {
        liNaIndices[[lsNames[iCol]]] = union(which(is.na(frmeData[,iCol])),which(frmeData[,iCol]=="NA"))
      }
      frmeData[,iCol] <- factor( aCol, levels = c(lsFactor[[2]], "NA") )
    }
  }

  # Make sure there is a minimum number of non-0 measurements
  aiRemove = c()
  for( iCol in aiData )
  {
    adCol = frmeData[,iCol]
    if(length( which(adCol!=0)) < ( dMinSamp * length( adCol ) ) )
    {
      aiRemove = c(aiRemove, iCol)
    }
  }

  # Remove and document
  aiData = setdiff( aiData, aiRemove)
  lsQCCounts$iZeroDominantData = aiRemove
  if(length(aiRemove))
  {
    c_logrMaaslin$info( "Removing the following for having not enough non-zero measurments for analysis.")
    c_logrMaaslin$info( format( colnames( frmeData )[aiRemove] ))
  }

  c_logrMaaslin$debug("End FuncClean")
  return( list(frmeData = frmeData, aiMetadata = aiMetadata, aiData = aiData, lsQCCounts = lsQCCounts, liNaIndices=liNaIndices, viNotTransformedData = viNotTransformedData) )
  ### Return list of
  ### frmeData: The Data after cleaning
  ### aiMetadata: The indices of the metadata still being used after filtering
  ### aiData: The indices of the data still being used after filtering
  ### lsQCCOunts: QC info
}

funcBugs <- function(
### Run analysis of all data features against all metadata
frmeData,
### Cleaned data including metadata, and data
lsData,
### This list is a general container for data as the analysis occurs, think about it as a cache for the analysis
aiMetadata,
### Indices of metadata used in analysis
aiData,
### Indices of response data
aiNotTransformedData,
### Indicies of the data not transformed
strData,
### Log file name
dSig,
### Significance threshold for the qvalue cut off
fInvert=FALSE,
### Invert images to have a black background
strDirOut = NA,
### Output project directory
funcReg=NULL,
### Function for regularization
funcTransform=NULL,
### Function used to transform the data
funcUnTransform=NULL,
### If a transform is used the opposite of that transfor must be used on the residuals in the partial residual plots
lsNonPenalizedPredictors=NULL,
### These predictors will not be penalized in the feature (model) selection step
funcAnalysis=NULL,
### Function to perform association analysis
lsRandomCovariates=NULL,
### List of string names of metadata which will be treated as random covariates
funcGetResults=NULL,
### Function to unpack results from analysis
fDoRPlot=TRUE,
### Plot residuals
fOmitLogFile = FALSE,
### Stops the creation of the log file
fAllvAll=FALSE,
### Flag to turn on all against all comparisons
liNaIndices = list(),
### Indicies of imputed NA data
lxParameters=list(),
### List holds parameters for different variable selection techniques
strTestingCorrection = "BH",
### Correction for multiple testing
fIsUnivariate = FALSE,
### Indicates if the function is univariate
fZeroInflated = FALSE
### Indicates to use a zero infalted model
){
  c_logrMaaslin$debug("Start funcBugs")
  # If no output directory is indicated
  # Then make it the current directory
  if( is.na( strDirOut ) || is.null( strDirOut ) )
  {
    if( !is.na( strData ) )
    {
      strDirOut <- paste( dirname( strData ), "/", sep = "" )
    } else { strDirOut = "" }
  }

  # Make th log file and output file names based on the log file name
  strLog = NA
  strBase = ""
  if(!is.na(strData))
  {
    strBaseOut <- paste( strDirOut, sub( "\\.([^.]+)$", "", basename(strData) ), sep = "/" )
    strLog <- paste( strBaseOut,c_sLogFileSuffix, ".txt", sep = "" )
  }

  # If indicated, stop the creation of the log file
  # Otherwise delete the log file if it exists and log
  if(fOmitLogFile){ strLog = NA }
  if(!is.na(strLog))
  {
    c_logrMaaslin$info( "Outputting to: %s", strLog )
    unlink( strLog )
  }
 
  # Will contain pvalues
  adP = c()
  adPAdj = c()

  # List of lists with association information
  lsSig <- list()
  # Go through each data that was not previously removed and perform inference
  for( iTaxon in aiData )
  {
    # Log to screen progress per 10 associations.
    # Can be thown off if iTaxon is missing a mod 10 value
    # So the taxons may not be logged every 10 but not a big deal
    if( !( iTaxon %% 10 ) )
    {
      c_logrMaaslin$info( "Taxon %d/%d", iTaxon, max( aiData ) )
    }

    # Call analysis method
    lsOne <- funcBugHybrid( iTaxon=iTaxon, frmeData=frmeData, lsData=lsData, aiMetadata=aiMetadata, dSig=dSig, adP=adP, lsSig=lsSig, funcTransform=funcTransform, funcUnTransform=funcUnTransform, strLog=strLog, funcReg=funcReg, lsNonPenalizedPredictors=lsNonPenalizedPredictors, funcAnalysis=funcAnalysis, lsRandomCovariates=lsRandomCovariates, funcGetResult=funcGetResults, fAllvAll=fAllvAll, fIsUnivariate=fIsUnivariate, lxParameters=lxParameters, fZeroInflated=fZeroInflated, fIsTransformed= ! iTaxon %in% aiNotTransformedData )

    # If you get a NA (happens when the lmm gets all random covariates) move on
    if( is.na( lsOne ) ){ next }

    # The updating of the following happens in the inference method call in the funcBugHybrid call
    # New pvalue array
    adP <- lsOne$adP
    # New lsSig contains data about significant feature v metadata comparisons
    lsSig <- lsOne$lsSig
    # New qc data
    lsData$lsQCCounts = lsOne$lsQCCounts
  }

  # Log the QC info
  c_logrMaaslin$debug("lsData$lsQCCounts")
  c_logrMaaslin$debug(format(lsData$lsQCCounts))

  if( is.null( adP ) ) { return( NULL ) }

  # Perform bonferonni corrections on factor data (for levels), calculate the number of tests performed, and FDR adjust for multiple hypotheses
  # Perform Bonferonni adjustment on factor data
  for( iADIndex in 1:length( adP ) )
  {
    # Only perform on factor data
    if( is.factor( lsSig[[ iADIndex ]]$metadata ) )
    {
      adPAdj = c( adPAdj, funcBonferonniCorrectFactorData( dPvalue = adP[ iADIndex ], vsFactors = lsSig[[ iADIndex ]]$metadata, fIgnoreNAs = length(liNaIndices)>0) )
    } else {
      adPAdj = c( adPAdj, adP[ iADIndex ] )
    }
  }

  iTests = funcCalculateTestCounts(iDataCount = length(aiData), asMetadata = intersect( lsData$astrMetadata, colnames( frmeData )[aiMetadata] ), asForced = lsNonPenalizedPredictors, asRandom = lsRandomCovariates, fAllvAll = fAllvAll)

  #Get indices of sorted data after the factor correction but before the multiple hypothesis corrections.
  aiSig <- sort.list( adPAdj )

  # Perform FDR BH
  adQ = p.adjust(adPAdj, method=strTestingCorrection, n=max(length(adPAdj), iTests))

  # Find all covariates that had significant associations
  astrNames <- c()
  for( i in 1:length( lsSig ) )
  {
    astrNames <- c(astrNames, lsSig[[i]]$name)
  }
  astrNames <- unique( astrNames )

  # Sets up named label return for global plotting
  lsReturnTaxa <- list()
  for( j in aiSig )
  {
    if( adQ[j] > dSig ) { next }
    strTaxon <- lsSig[[j]]$taxon
    if(strTaxon %in% names(lsReturnTaxa))
    {
      lsReturnTaxa[[strTaxon]] = min(lsReturnTaxa[[strTaxon]],adQ[j])
    } else { lsReturnTaxa[[strTaxon]] = adQ[j]}
  }

  # For each covariate with significant associations
  # Write out a file with association information
  for( strName in astrNames )
  {
    strFileTXT <- NA
    strFilePDF <- NA
    for( j in aiSig )
    {
      lsCur		<- lsSig[[j]]
      strCur		<- lsCur$name

      if( strCur != strName ) { next }

      strTaxon		<- lsCur$taxon
      adData		<- lsCur$data
      astrFactors	<- lsCur$factors
      adCur		<- lsCur$metadata
      adY <- adData

      if( is.na( strData ) ) { next }

      ## If the text file output is not written to yet
      ## make the file names, and delete any previous file output 
      if( is.na( strFileTXT ) )
      {
        strFileTXT <- sprintf( "%s-%s.txt", strBaseOut, strName )
        unlink(strFileTXT)
        funcWrite( c("Variable", "Feature", "Value", "Coefficient", "N", "N not 0", "P-value", "Q-value"), strFileTXT )
      }

      ## Write text output
      funcWrite( c(strName, strTaxon, lsCur$orig, lsCur$value, length( adData ), sum( adData > 0 ), adP[j], adQ[j]), strFileTXT )

      ## If the significance meets the threshold
      ## Write PDF file output
      if( adQ[j] > dSig ) { next }

      # Do not make residuals plots if univariate is selected
      strFilePDF = funcPDF( frmeTmp=frmeData, lsCur=lsCur, curPValue=adP[j], curQValue=adQ[j], strFilePDF=strFilePDF, strBaseOut=strBaseOut, strName=strName, funcUnTransform= funcUnTransform, fDoResidualPlot=fDoRPlot, fInvert=fInvert, liNaIndices=liNaIndices )
   }
    if( dev.cur( ) != 1 ) { dev.off( ) }
  }
  aiTmp <- aiData

  logdebug("End funcBugs", c_logMaaslin)
  return(list(lsReturnBugs=lsReturnTaxa, lsQCCounts=lsData$lsQCCounts))
  ### List of data features successfully associated without error and quality control data
}

#Lightly Tested
### Performs analysis for 1 feature
### iTaxon: integer Taxon index to be associated with data
### frmeData: Data frame The full data
### lsData: List of all associated data
### aiMetadata: Numeric vector of indices
### dSig: Numeric significance threshold for q-value cut off
### adP: List of pvalues from associations
### lsSig: List which serves as a cache of data about significant associations
### strLog: String file to log to
funcBugHybrid <- function(
### Performs analysis for 1 feature
iTaxon,
### integer Taxon index to be associated with data
frmeData,
### Data frame, the full data
lsData,
### List of all associated data
aiMetadata,
### Numeric vector of indices
dSig,
### Numeric significance threshold for q-value cut off
adP,
### List of pvalues from associations
lsSig,
### List which serves as a cache of data about significant associations
funcTransform,
### The tranform used on the data
funcUnTransform,
### The reverse transform on the data
strLog = NA,
### String, file to which to log
funcReg=NULL,
### Function to perform regularization
lsNonPenalizedPredictors=NULL,
### These predictors will not be penalized in the feature (model) selection step
funcAnalysis=NULL,
### Function to perform association analysis
lsRandomCovariates=NULL,
### List of string names of metadata which will be treated as random covariates
funcGetResult=NULL,
### Function to unpack results from analysis
fAllvAll=FALSE,
### Flag to turn on all against all comparisons
fIsUnivariate = FALSE,
### Indicates the analysis function is univariate
lxParameters=list(),
### List holds parameters for different variable selection techniques
fZeroInflated = FALSE,
### Indicates if to use a zero infalted model
fIsTransformed = TRUE
### Indicates that the bug is transformed
){
#dTime00 <- proc.time()[3]
  #Get metadata column names
  astrMetadata = intersect( lsData$astrMetadata, colnames( frmeData )[aiMetadata] )

  #Get data measurements that are not NA
  aiRows <- which( !is.na( frmeData[,iTaxon] ) )

  #Get the dataframe of non-na data measurements
  frmeTmp <- frmeData[aiRows,]

  #Set the min boosting selection frequency to a default if not given
  if( is.na( lxParameters$dFreq ) )
  {
    lxParameters$dFreq <- 0.5 / length( c(astrMetadata) )
  }

  # Get the full data for the bug feature
  adCur = frmeTmp[,iTaxon]
  lxParameters$sBugName = names(frmeTmp[iTaxon])

  # This can run multiple models so some of the results are held in lists and some are not
  llmod = list()
  liTaxon = list()
  lastrTerms = list()

  # Build formula for simple mixed effects models
  # Removes random covariates from variable selection
  astrMetadata  = setdiff(astrMetadata, lsRandomCovariates)
  strFormula <- paste( "adCur ~", paste( sprintf( "`%s`", astrMetadata ), collapse = " + " ), sep = " " )

  # Document the model
  funcWrite( c("#taxon", colnames( frmeTmp )[iTaxon]), strLog )
  funcWrite( c("#metadata", astrMetadata), strLog )
  funcWrite( c("#samples", rownames( frmeTmp )), strLog )

  #Model terms
  astrTerms <- c()

  # Attempt feature (model) selection
  if(!is.na(funcReg))
  {
    #Count model selection method attempts
    lsData$lsQCCounts$iBoosts = lsData$lsQCCounts$iBoosts + 1
    #Perform model selection
    astrTerms <- funcReg(strFormula=strFormula, frmeTmp=frmeTmp, adCur=adCur, lsParameters=lxParameters, lsForcedParameters=lsNonPenalizedPredictors, strLog=strLog)
    #If the feature selection function is set to None, set all terms of the model to all the metadata
  } else { astrTerms = astrMetadata }

  # Get look through the boosting results to get a model
  # Holds the predictors in the predictors in the model that were selected by the boosting
  if(is.null( astrTerms )){lsData$lsQCCounts$iBoostErrors = lsData$lsQCCounts$iBoostErrors + 1}

  # Get the indices that are transformed
  # Of those indices check for uneven metadata
  # Untransform any of the metadata that failed
  # Failed means true for uneven occurences of zeros
#  if( fIsTransformed )
#  {
#    vdUnevenZeroCheck = funcUnTransform( frmeData[[ iTaxon ]] )
#    if( funcZerosAreUneven( vdRawData=vdUnevenZeroCheck, funcTransform=funcTransform, vsStratificationFeatures=astrTerms, dfData=frmeData ) )
#    {
#      frmeData[[ iTaxon ]] = vdUnevenZeroCheck
#      c_logrMaaslin$debug( paste( "Taxon transformation reversed due to unevenness of zero distribution.", iTaxon ) )
#    }
#  }

  # Run association analysis if predictors exist and an analysis function is specified
  # Run analysis
  if(!is.na(funcAnalysis) )
  {
    #If there are selected and forced fixed covariates
    if( length( astrTerms ) )
    {
      #Count the association attempt
      lsData$lsQCCounts$iLms = lsData$lsQCCounts$iLms + 1

      #Make the lm formula
      #Build formula for simple mixed effects models using random covariates
      strRandomCovariatesFormula = NULL
      #Random covariates are forced
      if(length(lsRandomCovariates)>0)
      {
        #Format for lme
        #Needed for changes to not allowing random covariates through the boosting process
        strRandomCovariatesFormula <- paste( "adCur ~ ", paste( sprintf( "1|`%s`", lsRandomCovariates), collapse = " + " ))
      }

      #Set up a list of formula containing selected fixed variables changing and the forced fixed covariates constant
      vstrFormula = c()
      #Set up suppressing forced covariates in a all v all scenario only
      asSuppress = c()
      #Enable all against all comparisons
      if(fAllvAll && !fIsUnivariate)
      {
        lsVaryingCovariates = setdiff(astrTerms,lsNonPenalizedPredictors)
        lsConstantCovariates = setdiff(lsNonPenalizedPredictors,lsRandomCovariates)
        strConstantFormula = paste( sprintf( "`%s`", lsConstantCovariates ), collapse = " + " )
        asSuppress = lsConstantCovariates

        if(length(lsVaryingCovariates)==0L)
        {
          vstrFormula <- c( paste( "adCur ~ ", paste( sprintf( "`%s`", lsConstantCovariates ), collapse = " + " )) )
        } else {
          for( sVarCov in lsVaryingCovariates )
          {
            strTempFormula = paste( "adCur ~ `", sVarCov,"`",sep="")
            if(length(lsConstantCovariates)>0){ strTempFormula = paste(strTempFormula,strConstantFormula,sep=" + ") }
            vstrFormula <- c( vstrFormula, strTempFormula )
          }
        }
      } else {
        #This is either the multivariate case formula for all covariates in an lm or fixed covariates in the lmm
        vstrFormula <- c( paste( "adCur ~ ", paste( sprintf( "`%s`", astrTerms ), collapse = " + " )) )
      }

      #Run the association
      for( strAnalysisFormula in vstrFormula )
      {
        i = length(llmod)+1

        llmod[[i]] = funcAnalysis(strFormula=strAnalysisFormula, frmeTmp=frmeTmp, iTaxon=iTaxon, lsHistory=list(adP=adP, lsSig=lsSig, lsQCCounts=lsData$lsQCCounts), strRandomFormula=strRandomCovariatesFormula, fZeroInflated=fZeroInflated)

        liTaxon[[i]] = iTaxon
        lastrTerms[[i]] = funcFormulaStrToList(strAnalysisFormula)
      }
    } else {
      #If there are no selected or forced fixed covariates
      lsData$lsQCCounts$iNoTerms = lsData$lsQCCounts$iNoTerms + 1
      return(list(adP=adP, lsSig=lsSig, lsQCCounts=lsData$lsQCCounts))
    }
  }

  #Call funcBugResults and return it's return
  if(!is.na(funcGetResult))
  {
    #Format the results to a consistent expected result.
    return( funcGetResult( llmod=llmod, frmeData=frmeData, liTaxon=liTaxon, dSig=dSig, adP=adP, lsSig=lsSig, strLog=strLog, lsQCCounts=lsData$lsQCCounts, lastrCols=lastrTerms, asSuppressCovariates=asSuppress ) )
  } else {
    return(list(adP=adP, lsSig=lsSig, lsQCCounts=lsData$lsQCCounts))
  }
  ### List containing a list of pvalues, a list of significant data per association, and a list of QC data
}
