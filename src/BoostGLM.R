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

funcClean <- function(
### Properly clean / get data ready for analysis
### Includes custom analysis from the custom R script if it exists
frmeData,
### Data frame, input data to be acted on
funcDataProcess,
### Custom script that can be given to perform specialized processing before MaAsLin does.
aiMetadata,
### Indices of columns in frmeData which are metadata for analysis.
aiGenetics,
### Indices of columns in frmeData which are genetics for analysis.
aiData,
### Indices of column in frmeData which are (abundance) data for analysis.
lsQCCounts,
### List that will hold the quality control information which is written in the output directory.
astrNoImpute = c(),
### An array of column names of frmeData not to impute.
dMinSamp,
### Minimum number of samples
dFence,
### How many quartile ranges defines the fence to define outliers.
funcTransform
### The data transformation function or a dummy function that does not affect the data
){
  # Call the custom script and set current data and indicies to the processes data and indicies.
  c_logrMaaslin$debug( "Start Clean")
  if( !is.null( funcDataProcess ) )
  {
    c_logrMaaslin$debug("Additional preprocess function attempted.")

    pTmp <- funcDataProcess( frmeData=frmeData, aiMetadata=aiMetadata, aiGenetics=aiGenetics, aiData=aiData, funcTransform)
    frmeData = pTmp$frmeData
    aiMetadata = pTmp$aiMetadata
    aiGenetics = pTmp$aiGenetics
    aiData = pTmp$aiData
    lsQCCounts$lsQCCustom = pTmp$lsQCCounts
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

  # Metadata and genetics: Remove missing data
  # This is defined as if there is only one non-NA value or
  # if the number of NOT NA data is less than a percentage of the data defined by dMinSamp
  aiRemove = c()
  for( iCol in c(aiMetadata, aiGenetics) )
  {
    adCol = frmeData[,iCol]
    if( ( sum( !is.na( adCol ) ) < ( dMinSamp * length( adCol ) ) ) ||
      ( length( unique( na.omit( adCol ) ) ) < 2 ) )
    {
      aiRemove = c(aiRemove, iCol)
    }
  }

  # Remove metadata and genetics data
  aiMetadata = setdiff( aiMetadata, aiRemove )
  aiGenetics = setdiff( aiGenetics, aiRemove )

  # Update the data which was removed.
  lsQCCounts$iMissingMetadata = aiRemove
  if(length(aiRemove))
  {
    c_logrMaaslin$info("Removing the following metadata/genetics for too much missing data or only one data value outside of NA.")
    c_logrMaaslin$info(format(colnames( frmeData )[aiRemove]))
  }

  # Remove genetics where the minimum number of values over one is less tan a percentage (dMinSamp) of the sample
  aiRemove = c()
  for( iCol in aiGenetics )
  {
    adCol = frmeData[,iCol]
    if( sum( adCol > 0, na.rm = TRUE ) < ( dMinSamp * length( adCol ) ) )
    {
      aiRemove = c(aiRemove, iCol)
    }
  }

  # Remove index of genetics which where filtered out
  aiGenetics <- setdiff( aiGenetics, aiRemove )
  lsQCCounts$iRemovedSNPs = aiRemove

  # Document
  if(length(aiRemove))
  {
    c_logrMaaslin$info("Removing the following genetics indicies, too sparse.")
    c_logrMaaslin$info(format(colnames( frmeData )[aiRemove]))
  }

  # Remove outliers
  # If the dFence Value is set use the method of defining the outllier as
  # dFence * the interquartile range + or - the 3rd and first quartile respectively.
  # If not the gibbs test is used.
  aiSumOutlierPerDatum = c()
  if( dFence > 0.0 )
  {
    print("OUTLIERS: Using fence")
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
        c_logrMaaslin$info( "Removing %d outliers from %s", length( aiRemove ), colnames(frmeData)[iData] )
        c_logrMaaslin$info( format( rownames( frmeData )[aiRemove] ))
      }
      adData[aiRemove] <- NA
      frmeData[,iData] <- adData
      aiSumOutlierPerDatum = c(aiSumOutlierPerDatum,length(aiRemove))
    }
  #Do not use the fence, use the Grubbs test
  } else {
    for( iData in aiData )
    {
      adData <- frmeData[,iData]

      while( TRUE )
      {
        lsTest <- try( grubbs.test( adData ), silent = TRUE )
        if( ( class( lsTest ) == "try-error" ) || is.na( lsTest$p.value ) || ( lsTest$p.value > c_dPOutlier ) )
        {break}
        adData[outlier( adData, logical = TRUE )] <- NA
      }

      # Document removal
      if( sum( is.na( adData ) ) )
      {
        c_logrMaaslin$info( "Removing %d outliers from %s", sum( is.na( adData ) ), colnames(frmeData)[iData] )
			  c_logrMaaslin$info( format( rownames( frmeData )[is.na( adData )] ))
      }
      frmeData[,iData] <- adData
    }
  }
  lsQCCounts$aiSumOutlierPerDatum = aiSumOutlierPerDatum

  # Remove missing data, remove any sample has less than dMinSamp * the number of data
  aiRemove = c()
  for( iCol in aiData )
  {
    adCol = frmeData[,iCol]
    adCol[!is.finite( adCol )] <- NA
    if( ( sum( !is.na( adCol ) ) < ( dMinSamp * length( adCol ) ) ) ||
      ( length( unique( na.omit( adCol ) ) ) < 2 ) )
    {
        aiRemove = c(aiRemove, iCol)
    }
  }

  # Remove and document
  aiData = setdiff( aiData, aiRemove )
  lsQCCounts$iMissingData = aiRemove
  if(length(aiRemove))
  {
    c_logrMaaslin$info( "Removing the following for missing data.")
    c_logrMaaslin$info( format( colnames( frmeData )[aiRemove] ))
  }

  # Keep track of factor levels in a list for later use
  lslsFactors <- list()
  for( iCol in c(aiGenetics, aiMetadata) )
  {
    aCol <- frmeData[,iCol]
    if( class( aCol ) == "factor" )
    {
      lslsFactors[[length( lslsFactors ) + 1]] <- list(iCol, levels( aCol ))
    }
  }

  # Replace missing data values by the mean of the data column.
  for( iCol in aiData )
  {
    adCol <- frmeData[,iCol]
    adCol[is.infinite( adCol )] <- NA
    adCol[is.na( adCol )] <- mean( adCol, na.rm = TRUE )
    frmeData[,iCol] <- adCol
  }

  #Use na.gam.replace to manage NA metadata and genetics
  aiTmp <- c(aiGenetics, aiMetadata)
  aiTmp <- setdiff( aiTmp, which( colnames( frmeData ) %in% astrNoImpute ) )
  frmeData[,aiTmp] <- na.gam.replace( frmeData[,aiTmp] )

  #If NA is a value in  factor data, set the NA as a level.
  for( lsFactor in lslsFactors )
  {
    iCol <- lsFactor[[1]]
    aCol <- frmeData[,iCol]
    if( "NA" %in% levels( aCol ) )
    {
      frmeData[,iCol] <- factor( aCol, levels = c(lsFactor[[2]], "NA") )
    }
  }

  c_logrMaaslin$debug("End FuncClean")
  return( list(frmeData = frmeData, aiMetadata = aiMetadata, aiGenetics = aiGenetics, aiData = aiData, lsQCCounts = lsQCCounts) )
  ### Return list of
  ### frmeData: The Data after cleaning
  ### aiMetadata: The indices of the metadata still being used after filtering
  ### aiGenetics: The indices of the genetics still being used after filtering
  ### aiData: The indices of the data still being used after filtering
  ### lsQCCOunts: QC info
}

funcBugs <- function(
### Run analysis of all data features against all metadata
frmeData,
### Cleaned data including metadata, genetics, and data
lsData,
### This list is a general container for data as the analysis occurs, think about it as a cache for the analysis
aiMetadata,
### Indices of metadata used in analysis
aiGenetics,
### Indices of genetics used in analysis
aiData,
### Indices of response data
strData,
### Log file name
dFreq,
### 
dSig,
### Signifcance threshold for the qvalue cut off
dMinSamp,
### Minimum number of samples
fInvert,
### Invert images to have a black background
strDirOut = NA,
### Output project directory
astrScreen = c(),
### 
funcReg=NULL,
### Function for regularization
funcAnalysis=NULL,
### Function to perform association analysis
funcGetResults=NULL
### Function to unpack results from analysis
){
  c_logrMaaslin$debug("Start funcBugs")
  if( is.na( strDirOut ) )
  {
    strDirOut <- paste( dirname( strData ), "/", sep = "" )
  }
  strBaseOut <- paste( strDirOut, sub( "\\.([^.]+)$", "", basename(strData) ), sep = "/" )
  strLog <- paste( strBaseOut, ".txt", sep = "" )
  c_logrMaaslin$info( "Outputting to: %s", strLog )
  unlink( strLog )

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
    lsOne <- funcBugHybrid( iTaxon, frmeData, lsData, aiMetadata, aiGenetics, dFreq, dSig, dMinSamp, adP, lsSig, strLog, funcReg, funcAnalysis, funcGetResults )

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
  adQ <- adP
  iTests <- length( intersect( lsData$astrMetadata, colnames( frmeData )[aiMetadata] ) ) * length( aiData )
  #Perform FDR BH
  for( i in 1:length( aiSig ) )
  {
    iSig <- aiSig[i]
    adQ[iSig] <- adP[iSig] * iTests / i
  }

  #
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
	
    if( length( aiGenetics ) )
    {
      if( length( intersect( astrFactors, colnames( frmeData )[aiGenetics] ) ) )
      {
      	astrRet <- c(astrRet, astrFactors)
      }
    } else if( !length( astrScreen ) || sum( astrFactors %in% astrScreen ) ) {
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
      strFilePDF = funcPDF( lsCur=lsCur, curPValue=adP[j], curQValue=adQ[j], aiGenetics=aiGenetics, strFilePDF=strFilePDF, strBaseOut=strBaseOut, strName=strName, fInvert=fInvert )
    }

    if( dev.cur( ) != 1 ) { dev.off( ) }
  }

  if( length( aiGenetics ) )
  {
    aiTmp <- aiGenetics
  } else {
    aiTmp <- aiData
  }

  logdebug("End funcBugs", c_logMaaslin)
  aiReturnBugs = aiTmp[colnames( frmeData )[aiTmp] %in% astrRet]
#  return( aiTmp[colnames( frmeData )[aiTmp] %in% astrRet] )
  return(list(aiReturnBugs=aiReturnBugs, lsQCCounts=lsData$lsQCCounts))
  ### List of data features successfully associated without error and quality control data
}

### Performs analysis for 1 feature
### iTaxon: integer Taxon index to be associated with data
### frmeData: Data frame The full data
### lsData: List of all associated data
### aiMetadata: Numeric vector of indices
### aiGenetics: Numeric vector of genetics data
### dFreq: Used to select metadata from the boosting, selected refrequency must be larger than this
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
aiGenetics,
### Numeric vector of genetics data
dFreq,
### Used to select metadata from the boosting, selected refrequency must be larger than this
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
funcAnalysis=NULL,
### Function to perform association analysis
funcGetResult=NULL
### Function to unpack results from analysis
){

#dTime00 <- proc.time()[3]
  #Get metadata and genetics column names
  astrMetadata = intersect( lsData$astrMetadata, colnames( frmeData )[aiMetadata] )
  astrGenetics = colnames( frmeData )[aiGenetics]

  #Get data measurements that are not NA
  aiRows <- which( !is.na( frmeData[,iTaxon] ) )

  #Get the dataframe of non-na data measurements
  frmeTmp <- frmeData[aiRows,]

  #For each metadata, check it's data and see if there are too many NA to move forward.
  astrRemove <- c()
  for( strMetadatum in astrMetadata )
  {
    aMetadatum <- frmeTmp[,strMetadatum]

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
  if( is.na( dFreq ) )
  {
    dFreq <- 0.5 / length( c(astrMetadata, astrGenetics) )
  }

  # Get the full data for the bug feature
  adCur = frmeTmp[,iTaxon]

  #Create a linear additive model including all metadata or genetics which were not filtered
  lmod <- NA
  strFormula <- paste( "adCur ~", paste( sprintf( "`%s`", astrMetadata ), collapse = " + " ),
  ifelse( length( astrGenetics ), "+", "" ), paste( astrGenetics, collapse = " + " ), sep = " " )

  # Document the model
  funcWrite( c("#taxon", colnames( frmeTmp )[iTaxon]), strLog )
  funcWrite( c("#metadata", astrMetadata), strLog )
  funcWrite( c("#Genetics", ifelse(length(astrGenetics),astrGenetics,"Not Genetics")), strLog )
  funcWrite( c("#samples", rownames( frmeTmp )), strLog )

  # Regularisation terms
  astrTerms <- c()

  # Attempt model selection
  lmod <- NA
  if(!is.null(funcReg))
  {
    #Count model selection method attempts
    lsData$lsQCCounts$iBoosts = lsData$lsQCCounts$iBoosts + 1
    #Perform model selection
    astrTerms <- funcReg(strFormula, frmeTmp, adCur, list(dFreq=dFreq), strLog)
    #If the boosting is set to None, set all terms of the model to all the metadata
  } else { astrTerms = astrMetadata }

  #Get look through the boosting results to get a model
  #Holds the predictors in the predictors in the model that were selected by the boosting
  if(is.na( astrTerms ))
  {lsData$lsQCCounts$iBoostErrors = lsData$lsQCCounts$iBoostErrors + 1}

  #Run association analysis if predictors exist and an analysis function is specified
  #Reset model for the glm
  lmod <- NA
  if( length( astrTerms ) && !is.null(funcAnalysis) )
  {
    #Count the association attempt
    lsData$lsQCCounts$iLms = lsData$lsQCCounts$iLms + 1
    #Make the lm formula
    strFormula <- paste( "adCur ~", paste( sprintf( "`%s`", astrTerms ), collapse = " + " ), sep = " " )
    #Run the association
    lmod <- funcAnalysis(strFormula, frmeTmp, adCur)
  } else {
    lsData$lsQCCounts$iNoTerms = lsData$lsQCCounts$iNoTerms + 1
  }

  #Call funBugResults and return it's return
  return( funcGetResult( lmod=lmod, frmeData=frmeData, iTaxon=iTaxon, dSig=dSig, adP=adP, lsSig=lsSig, strLog=strLog, lsQCCounts=lsData$lsQCCounts, astrCols=astrTerms ) )
  ### List containing a list of pvalues, a list of significant data per association, and a list of QC data
}
