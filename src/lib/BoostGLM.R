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

#strTagX = "PC9"
#strTagY = "k__Bacteria|p__Bacteroidetes|c__Bacteroidia|o__Bacteroidales|f__Prevotellaceae|g__Prevotella"

## Get logger
c_logrMaaslin <- getLogger( "maaslin" )

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

  # Remove missing data, remove any sample has less than dMinSamp * the number of data or low abundance
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
  for(aiDatum in aiData)
  {
    frmeData[,aiDatum] = funcTransform(frmeData[,aiDatum])
  }

  # Set data indicies after custom QC process.
  lsQCCounts$aiAfterPreprocess = aiData

  # Metadata: Properly factorize all logical data and integer and number data with less than iNonFactorLevelThreshold
  for( i in aiMetadata )
  {
    if( ( class( frmeData[,i] ) %in% c("integer", "numeric", "logical") ) &&
      ( length( unique( frmeData[,i] ) ) < c_iNonFactorLevelThreshold ) ) {
      c_logrMaaslin$debug(paste("Changing metadatum from numeric/integer/logical to factor",colnames(frmeData)[i],sep="="))
      frmeData[,i] = factor( frmeData[,i] )
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

  # Remove outliers
  # If the dFence Value is set use the method of defining the outllier as
  # dFence * the interquartile range + or - the 3rd and first quartile respectively.
  # If not the gibbs test is used.
  aiSumOutlierPerDatum = c()
  if( dFence > 0.0 )
  {
    # For each sample measurements
    for( iData in aiData )
    {
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

      # Record indices of values outside of fence to remove.
      aiRemove <- c()
      for( j in 1:length( adData ) )
      {
        d <- adData[j]
        if( !is.na( d ) && ( ( d < dLF ) || ( d > dUF ) ) )
        {
          aiRemove <- c(aiRemove, j)
        }
      }

      # Document
      if( length( aiRemove ) )
      {
        c_logrMaaslin$info( "OutliersByFence::Removing %d outliers from %s", length( aiRemove ), colnames(frmeData)[iData] )
        c_logrMaaslin$info( format( rownames( frmeData )[aiRemove] ))
      }
      adData[aiRemove] <- NA
      frmeData[,iData] <- adData
      aiSumOutlierPerDatum = c(aiSumOutlierPerDatum,length(aiRemove))
    }
  #Do not use the fence, use the Grubbs test
  } else if(dPOutlier!=0.0){
    for( iData in aiData )
    {
      adData <- frmeData[,iData]
      while( TRUE )
      {
        lsTest <- try( grubbs.test( adData ), silent = TRUE )
        if( ( class( lsTest ) == "try-error" ) || is.na( lsTest$p.value ) || ( lsTest$p.value > dPOutlier ) )
        {break}
        adData[outlier( adData, logical = TRUE )] <- NA
      }

      # Document removal
      if( sum( is.na( adData ) ) )
      {
        c_logrMaaslin$info( "Grubbs Test::Removing %d outliers from %s", sum( is.na( adData ) ), colnames(frmeData)[iData] )
			  c_logrMaaslin$info( format( rownames( frmeData )[is.na( adData )] ))
      }
      aiSumOutlierPerDatum = c(aiSumOutlierPerDatum,length(is.na(adData)))
      frmeData[,iData] <- adData
    }
  }
  lsQCCounts$aiSumOutlierPerDatum = aiSumOutlierPerDatum

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
    adCol[is.na( adCol )] <- mean( adCol, na.rm = TRUE )
    frmeData[,iCol] <- adCol

    if(length(which(is.na(frmeData[,iCol]))) == length(frmeData[,iCol]))
    {
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
  return( list(frmeData = frmeData, aiMetadata = aiMetadata, aiData = aiData, lsQCCounts = lsQCCounts, liNaIndices=liNaIndices) )
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
strData,
### Log file name
dSig,
### Significance threshold for the qvalue cut off
dMinSamp,
### Minimum number of samples
fInvert=FALSE,
### Invert images to have a black background
strDirOut = NA,
### Output project directory
astrScreen = c(),
### 
funcReg=NULL,
### Function for regularization
funcUnTransform=NULL,
### If a transform is used the opporite of that transfor must be used on the residuals in the partial residual plots
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
fIsUnivariate = FALSE
### Indicates if the function is univariate
){
  c_logrMaaslin$debug("Start funcBugs")
  if( is.na( strDirOut )||is.null(strDirOut))
  {
    if(!is.na(strData))
    {
      strDirOut <- paste( dirname( strData ), "/", sep = "" )
    } else { strDirOut = "" }
  }

  strLog = NA
  strBase = ""
  if(!is.na(strData))
  {
    strBaseOut <- paste( strDirOut, sub( "\\.([^.]+)$", "", basename(strData) ), sep = "/" )
    strLog <- paste( strBaseOut,c_sLogFileSuffix, ".txt", sep = "" )
  }
  if(fOmitLogFile){ strLog = NA }
  if(!is.na(strLog))
  {
    c_logrMaaslin$info( "Outputting to: %s", strLog )
    unlink( strLog )
  }

  #Will contain pvalues
  #Will contain objects associated with significance
  adP = c()
  lsSig <- list()
  for( iTaxon in aiData )
  {
    if( !( iTaxon %% 10 ) )
    {
      c_logrMaaslin$info( "Taxon %d/%d", iTaxon, max( aiData ) )
    }
    #Call analysis method
    lsOne <- funcBugHybrid( iTaxon, frmeData, lsData, aiMetadata, dSig, dMinSamp, adP, lsSig, strLog, funcReg, lsNonPenalizedPredictors, funcAnalysis, lsRandomCovariates, funcGetResults, fAllvAll, fIsUnivariate, lxParameters )

    #TODO Check#If you get a NA (happens when the lmm gets all random covariates) move on
    if(is.na(lsOne)){next}

    #Update pvalue array
    adP <- lsOne$adP
    #lsSig contains data about significant feature v metadata comparisons
    lsSig <- lsOne$lsSig
    #Update the qc data
    lsData$lsQCCounts = lsOne$lsQCCounts
  }
  c_logrMaaslin$debug("lsData$lsQCCounts")
  c_logrMaaslin$debug(format(lsData$lsQCCounts))

  #Presort for order for FDR calculation
  if( is.null( adP ) ) { return( NULL ) }

  #Get indices of sorted data
  aiSig <- sort.list( adP )
  iTests = funcCalculateTestCounts(iDataCount = length(aiData), asMetadata = intersect( lsData$astrMetadata, colnames( frmeData )[aiMetadata] ), asForced = lsNonPenalizedPredictors, asRandom = lsRandomCovariates, fAllvAll = fAllvAll)

  #Perform FDR BH
  adQ  = p.adjust(adP,method=strTestingCorrection,n=max(length(adP),iTests))

  astrNames <- c()
  for( i in 1:length( lsSig ) )
  {
    astrNames <- c(astrNames, lsSig[[i]]$name)
  }
  astrNames <- unique( astrNames )

  # Sets up named label return for MFA
  astrRet <- c()
  for( j in aiSig )
  {
    if( adQ[j] > dSig ) { next }
    #This allows only c_iMFA count of significant data to be passed to the MFA plot
    #Or in general to be pass out of the function.
    if( length( astrRet ) >= c_iMFA ) { break }

    lsCur <- lsSig[[j]]
    astrFactors <- lsCur$factors
    strTaxon <- lsCur$taxon
	
    if( !length( astrScreen ) || sum( astrFactors %in% astrScreen ) ) {
      astrRet <- c(astrRet, strTaxon)
    }
    astrRet <- unique( astrRet )
  }
		
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
  aiReturnBugs = aiTmp[colnames( frmeData )[aiTmp] %in% astrRet]
#  return( aiTmp[colnames( frmeData )[aiTmp] %in% astrRet] )
  return(list(aiReturnBugs=aiReturnBugs, lsQCCounts=lsData$lsQCCounts))
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
dMinSamp,
### Minimum amount of samples, used in QC
adP,
### List of pvalues from associations
lsSig,
### List which serves as a cache of data about significant associations
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
lxParameters=list()
### List holds parameters for different variable selection techniques
){
#dTime00 <- proc.time()[3]
  #Get metadata column names
  astrMetadata = intersect( lsData$astrMetadata, colnames( frmeData )[aiMetadata] )
  #Get data measurements that are not NA
  aiRows <- which( !is.na( frmeData[,iTaxon] ) )

  #Get the dataframe of non-na data measurements
  frmeTmp <- frmeData[aiRows,]
  #For each metadata, check it's data and see if there are too many NA to move forward.
  astrRemove <- c()
  for( strMetadatum in astrMetadata )
  {
    aMetadatum <- frmeTmp[,strMetadatum]

    #TODO remove?
    #Find the amount of NA and if over a certain ratio, remove that metadata from analysis
    iNA = sum( is.na( aMetadatum ) ) + sum( aMetadatum == "NA", na.rm = TRUE )
	#Be kinder to genetics than other metadata
	#if( length( grep( "^((chr)|(rs)\\d+)", strMetadatum ) ) ) {
	#	dMult <- 8 }
	#else {
	dMult <- 2 # }

    if( ( iNA / length( aiRows ) ) > ( dMult * dMinSamp ) )
    {
      astrRemove <- c(astrRemove, strMetadatum)
    }
  }

  #Document removal in logging
  if(length(astrRemove))
  {
    c_logrMaaslin$debug("These metadata will be removed in func bug")
    c_logrMaaslin$debug( format(c(colnames( frmeData )[iTaxon], astrRemove) ))
  }

  #Reset metadata with removed metadata removed.
  astrMetadata <- setdiff( astrMetadata, astrRemove )
  #Set the min boosting selection frequency to a default if not given
  if( is.na( lxParameters$dFreq ) )
  {
    lxParameters$dFreq <- 0.5 / length( c(astrMetadata) )
  }

  # Get the full data for the bug feature
  adCur = frmeTmp[,iTaxon]

  # This can run multiple models so some of the results are held in lists and some are not
  llmod = list()
  liTaxon = list()
  lastrTerms = list()

  #Build formula for simple mixed effects models
  #Removes random covariates from variable selection
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

  #Get look through the boosting results to get a model
  #Holds the predictors in the predictors in the model that were selected by the boosting
  if(is.null( astrTerms )){lsData$lsQCCounts$iBoostErrors = lsData$lsQCCounts$iBoostErrors + 1}

  #Run association analysis if predictors exist and an analysis function is specified
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

        llmod[[i]] = funcAnalysis(strFormula=strAnalysisFormula, frmeTmp=frmeTmp, iTaxon=iTaxon, lsHistory=list(adP=adP, lsSig=lsSig, lsQCCounts=lsData$lsQCCounts), strRandomFormula=strRandomCovariatesFormula)

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
