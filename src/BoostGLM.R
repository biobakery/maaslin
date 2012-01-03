####################################
# Summary: Boosting + GLM
# Author: Timothy Tickle
# Start Date: 10-26-2011
####################################

#External libraries
library( gam )
library( mboost )

#Properly clean / get data ready for analysis
funcClean = function( frmeData, funcDataProcess, aiMetadata, aiGenetics, aiData, lsQCCounts, astrNoImpute = c() )
{
  print("Start Clean")

  if( !is.null( funcDataProcess ) )
  {
    print("Additional preprocess function attempted.")

    pTmp = funcDataProcess( frmeData=frmeData, aiMetadata=aiMetadata, aiGenetics=aiGenetics, aiData=aiData)
    if( class( pTmp ) == "data.frame" )
    {
      frmeData = pTmp
    } else {
      frmeData = pTmp$frmeData
      aiMetadata = pTmp$aiMetadata
      aiGenetics = pTmp$aiGenetics
      aiData = pTmp$aiData
      lsQCCounts$lsQCCustom = pTmp$lsQCCounts
    }
  }

  lsQCCounts$aiAfterPreprocess = aiData

  c_dFence = 3

  # Properly factorize all categorical features
  for( i in aiMetadata )
  {
    if( ( class( frmeData[,i] ) %in% c("integer", "numeric", "logical") ) &&
      ( length( unique( frmeData[,i] ) ) < 5 ) )
    {
      print("Changing metadatum from numeric/integer/logical to factor")
      print(colnames(frmeData)[i])
      frmeData[,i] = factor( frmeData[,i] )
    }
  }

  # Remove missing metadata
  aiRemove = c()
  for( iCol in c(aiMetadata, aiGenetics) )
  {
    adCol = frmeData[,iCol]
    if( ( sum( !is.na( adCol ) ) < ( c_dMinSamp * length( adCol ) ) ) ||
      ( length( unique( na.omit( adCol ) ) ) < 2 ) )
    {
      aiRemove = c(aiRemove, iCol)
    }
  }
  aiMetadata = setdiff( aiMetadata, aiRemove )
  aiGenetics = setdiff( aiGenetics, aiRemove )
  lsQCCounts$iMissingMetadata = aiRemove
  if(length(aiRemove))
  {
    print("Removing the following metadata/genetics indicies, too much missing data.")
    print(aiRemove)
  }

  # Remove crummy SNPs
  aiRemove = c()
  for( iCol in aiGenetics )
  {
    adCol = frmeData[,iCol]
    if( sum( adCol > 0, na.rm = TRUE ) < ( c_dMinSamp * length( adCol ) ) )
    {
      aiRemove = c(aiRemove, iCol)
    }
  }
  aiGenetics <- setdiff( aiGenetics, aiRemove )
  lsQCCounts$iRemovedSNPs = aiRemove
  if(length(aiRemove))
  {
    print("Removing the following genetics indicies, too sparse.")
    print(aiRemove)
  }

  # Remove outliers
  aiSumOutlierPerDatum = c()
  if( c_dFence )
  {
    for( iData in aiData )
    {
      adData <- frmeData[,iData]
      adQ <- quantile( adData, c(0.25, 0.5, 0.75), na.rm = TRUE )
      dIQR <- adQ[3] - adQ[1]
      if(!dIQR)
      {
        dIQR = sd(adData,na.rm = TRUE)
      }
      dUF <- adQ[3] + ( c_dFence * dIQR )
      dLF <- adQ[1] - ( c_dFence * dIQR )
      aiRemove <- c()
      for( j in 1:length( adData ) )
      {
        d <- adData[j]
        if( !is.na( d ) && ( ( d < dLF ) || ( d > dUF ) ) )
        {
          aiRemove <- c(aiRemove, j)
        }
      }
#      if(length(aiRemove))
#      {
#        print(paste("Changing the following outliers to NA for data ",colnames(frmeData)[aiData],sep=""))
#        print(aiRemove)
#      }
      adData[aiRemove] <- NA
      frmeData[,iData] <- adData
      aiSumOutlierPerDatum = c(aiSumOutlierPerDatum,length(aiRemove))
    }
  }
  lsQCCounts$aiSumOutlierPerDatum = aiSumOutlierPerDatum

  # Remove missing data
  aiRemove = c()
  for( iCol in aiData )
  {
    adCol = frmeData[,iCol]
    if( ( sum( !is.na( adCol ) ) < ( c_dMinSamp * length( adCol ) ) ) ||
      ( length( unique( na.omit( adCol ) ) ) < 2 ) )
    {
      aiRemove = c(aiRemove, iCol)
    }
  }
  aiData = setdiff( aiData, aiRemove )
  lsQCCounts$iMissingData = aiRemove
  if(length(aiRemove))
  {
    print("Removing the following data indicies, for missing data.")
    print(aiRemove)
  }

  # Keep track of factor levels
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

  for( lsFactor in lslsFactors )
  {
    iCol <- lsFactor[[1]]
    aCol <- frmeData[,iCol]
    if( "NA" %in% levels( aCol ) )
    {
      frmeData[,iCol] <- factor( aCol, levels = c(lsFactor[[2]], "NA") )
    }
  }

  print("End FuncClean")
  return( list(frmeData = frmeData, aiMetadata = aiMetadata, aiGenetics = aiGenetics, aiData = aiData, lsQCCounts = lsQCCounts) )
}

#
#iTaxon integer Taxon index being looked at
#frmeData Data frame The full data
#lsData List of all associated data
#aiMetadata Numeric vector of indices
#aiGenetics Numeric vector of genetics data
#dSig Numeric significance threshold
#adP List of pvalues
#lsSig Complex list of significant info
#strLog String file to log to
funcBugHybrid <- function( iTaxon, frmeData, lsData, aiMetadata, aiGenetics, dSig, adP, lsSig, strLog = NA )
{
  #Get metadata and genetics strings
  astrMetadata = intersect( lsData$astrMetadata, colnames( frmeData )[aiMetadata] )
  astrGenetics = colnames( frmeData )[aiGenetics]

  #Get data rows that are not NA
  aiRows <- which( !is.na( frmeData[,iTaxon] ) )

  #Get the non NA row data
  frmeTmp <- frmeData[aiRows,]

  #For each metadata, check it's data and see if there are too many NA to move forward.
  astrRemove <- c()
  for( strMetadatum in astrMetadata )
  {
    aMetadatum <- frmeTmp[,strMetadatum]

    #Find the amount of NA and if over a certain ratio, remove that metadata from analysis
    iNA = max( sum( is.na( aMetadatum ) ), sum( aMetadatum == "NA", na.rm = TRUE ) )
    if( ( iNA / length( aiRows ) ) > ( 2 * c_dMinSamp ) )
    {
      astrRemove <- c(astrRemove, strMetadatum)
    }
  }
  if(length(astrRemove))
  {
    print("These metadata will be removed in func bug")
    print( c(colnames( frmeData )[iTaxon], astrRemove) )
  }
  #Reset metadata with removed metadata removed.
  astrMetadata <- setdiff( astrMetadata, astrRemove )
	
  # For the love of anything, I simply can't get the binomial link to be stable
  # asin/sqrt is ugly, but it seems to work a lot more nicely with a plain old normal link
  #adCur <- funcTransform( frmeTmp[,iTaxon] ) #asin(sqrt()) removed
  adCur = frmeTmp[,iTaxon]

  if( !is.na( strLog ) )
  {
    funcWrite( c("#taxon", colnames( frmeData )[iTaxon]), strLog )
    funcWrite( c("#metadata", astrMetadata), strLog )
    if( length( astrGenetics ) )
    {
      funcWrite( c("#genetics", astrGenetics), strLog )
    }
    funcWrite( c("#samples", rownames( frmeTmp )), strLog )
  }

  lmod <- NA
  #Create a linear additive mode including all metadata or genetics not removed at this iteration
  strFormula <- paste( "adCur ~", paste( sprintf( "`%s`", astrMetadata ), collapse = " + " ),
    ifelse( length( astrGenetics ), "+", "" ), paste( astrGenetics, collapse = " + " ), sep = " " )

  #Count model selection
  lsData$lsQCCounts$iBoosts = lsData$lsQCCounts$iBoosts + 1

  #Boost the model for model selection
  lmod <- try( glmboost( as.formula( strFormula ), data = frmeTmp,
    control = boost_control( nu = min( 1 / length( c(astrMetadata, astrGenetics) ) ), mstop = 5000 ) ) )

  astrTerms <- c()
  if( !is.na( lmod ) && ( class( lmod ) != "try-error" ) )
  {
    lsSum <- summary( lmod )
    if( !is.na( strLog ) )
    {
      funcWrite( "#model-glmboost", strLog )
      funcWrite( lmod, strLog )
      funcWrite( "#summary-glmboost", strLog )
      funcWrite( lsSum, strLog )
    }
    #For each metadata coefficient
    #Check selprob in the model summary
    for( strMetadata in names( coefficients( lmod ) ) )
    {
      #If the selprob is less than a certain frequency, skip
      #Added#
      if( is.na(lsSum$selprob[strMetadata])) { next }
      if( lsSum$selprob[strMetadata] < c_dFreq ) { next }
      #Get the name of the metadata
      strTerm <- funcCoef2Col( strMetadata, frmeData, c(astrMetadata, astrGenetics) )

      #If you cant find the coefficient name, write
      if( is.null( strTerm ) ) { next }
      if( is.na( strTerm ) )
      {
        print( c(strMetadata, strTerm) )
        next
      }
      #Collect metadata names
      astrTerms <- c(astrTerms, strTerm)
    }
    #Collect metadata names
    astrMetadata <- astrTerms
  } else {
    lsData$lsQCCounts$iBoostErrors = lsData$lsQCCounts$iBoostErrors + 1
  }

  lmod <- NA
  if( length( astrTerms ) )
  {
    lsData$lsQCCounts$iLms = lsData$lsQCCounts$iLms + 1
    strFormula <- paste( "adCur ~", paste( sprintf( "`%s`", astrTerms ), collapse = " + " ), sep = " " )
    lmod <- try( lm( as.formula( strFormula ), data = frmeTmp ) )
  } else {
    lsData$lsQCCounts$iNoTerms = lsData$lsQCCounts$iNoTerms + 1
  }

  return( funcBugResult( lmod=lmod, frmeData=frmeData, iTaxon=iTaxon, dSig=dSig, adP=adP, lsSig=lsSig, strLog=strLog, lsQCCounts=lsData$lsQCCounts, astrCols=astrTerms ) )
}

# Run metadata
# frmeData  
# funcBug Method call for performing analysis
# lsData Lists of data used in analysis
# aiMetadata Metadata indices
# aiGenetics Genetics indicies
# aiData Data indicies (integer)
# strData 
# dSig Significane (float)
# astrScreen 
#strData Log file name. NA indicates no logging. No append to previous sessions, file is deleted if old.
funcBugs <- function( frmeData, funcBug, lsData, aiMetadata, aiGenetics, aiData, strData, dSig, fInvert, strDirOut = NA, astrScreen = c() )
{
  print("Start funcBugs")
  if( is.na( strDirOut ) ) {
	  strDirOut <- paste( dirname( strData ), "/", sep = "" ) }
  strBaseOut <- paste( strDirOut, sub( "\\.(\\S+)$", "", basename(strData) ), sep = "/" )
  strLog <- paste( strBaseOut, ".txt", sep = "" )
  unlink( strLog )

  #Will contain pvalues
  #Will contain objects associated with significance
  adP = c()
  lsSig <- list()

  for( iTaxon in aiData )
  {
    if( !( iTaxon %% 10 ) )
    {
      print( c(iTaxon, max( aiData )) )
    }
    #Call analysis method
    lsOne <- funcBug( iTaxon, frmeData, lsData, aiMetadata, aiGenetics, dSig, adP, lsSig, strLog )

    #Update pvalue array
    adP <- lsOne$adP

    #lsSig contains data about significant feature v metadata comparisons
    lsSig <- lsOne$lsSig

    #Update the qc data
    lsData$lsQCCounts = lsOne$lsQCCounts
  }
  print("lsData$lsQCCounts")
  print(lsData$lsQCCounts)
  #Presort for order for FDR calculation
  if( is.null( adP ) ) { return( NULL ) }
  #Get indices of sorted data
  aiSig <- sort.list( adP )
  adQ <- adP

#######################
# THIS IS WHAT WE JUST CHANGED
#######################
  iTests <- length( intersect( lsData$astrMetadata, colnames( frmeData )[aiMetadata] ) ) * length( aiData )
#  aiTmp <- aiMetadata[colnames( frmeData )[aiMetadata] %in% lsData$astrMetadata]
#  iTests <- 0
#  for( i in aiTmp ) {
#    if( is.factor( frmeData[,i] ) ) {
#      iAdd <- nlevels( frmeData[,i] ) - 1 }
#    else {
#      iAdd <- 1 }
#    iTests <- iTests + iAdd }
#  iTests <- iTests * length( aiData )

  #Perform FDR BH
  for( i in 1:length( aiSig ) )
  {
    iSig <- aiSig[i]
    adQ[iSig] <- adP[iSig] * iTests / i
  }
	
  astrNames <- c()
  for( i in 1:length( lsSig ) )
  {
    astrNames <- c(astrNames, lsSig[[i]]$name)
  }
  astrNames <- unique( astrNames )

  astrRet <- c()
  for( j in aiSig )
  {
    if( adQ[j] > dSig ) { next }
    if( length( astrRet ) >= c_iMFA ) { break }
    lsCur <- lsSig[[j]]
    astrFactors <- lsCur$factors
    strTaxon <- lsCur$taxon
    if( length( aiGenetics ) )
    {
      if( length( intersect( astrFactors, colnames( frmeData )[aiGenetics] ) ) )
      { astrRet <- c(astrRet, astrFactors) }
    } else if( !length( astrScreen ) || sum( astrFactors %in% astrScreen ) ){
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
      lsCur <- lsSig[[j]]
      strCur <- lsCur$name
      if( strCur != strName ) { next }
      strTaxon <- lsCur$taxon
      adData <- lsCur$data
      astrFactors <- lsCur$factors
      adCur <- lsCur$metadata
      if( is.na( strData ) ) { next }
			
      if( is.na( strFileTXT ) )
      {
        strFileTXT <- sprintf( "%s-%s.txt", strBaseOut, strName )
        unlink(strFileTXT)
        funcWrite( c("astrFactors", "Taxon", "Coefficient", "N", "N not 0", "P-value", "Q-value"), strFileTXT )
      }
      funcWrite( c(strName, strTaxon, lsCur$orig, length( adData ), sum( adData > 0 ), adP[j], adQ[j]), strFileTXT )
      if( adQ[j] > dSig ) { next }
      if( is.na( strFilePDF ) )
      {
        strFilePDF <- sprintf( "%s-%s.%s", sub( "\\.\\S+$", "", strData ), strName, "pdf" )
        pdf( strFilePDF, width = 11 )
        if( fInvert )
        {
          par( bg = "black", fg = "white", col.axis = "white", col.lab = "white", col.main = "white", col.sub = "white" )
          adColorMin <- c(1, 1, 0)
          adColorMax <- c(0, 1, 1)
          adColorMed <- c(1, 1, 1)
        } else {
          adColorMin <- c(1, 0, 0)
          adColorMax <- c(0, 1, 0)
          adColorMed <- c(0, 0, 0)
        }
      }
      strTitle <- sprintf( "%s (%.3g sd %.3g, p=%.3g, q=%.3g)", lsCur$orig, lsCur$value, lsCur$std, adP[j], adQ[j] )
      adMar <- c(5, 4, 4, 2) + 0.1
      dLine <- NA
      if( nchar( strTaxon ) > 80 )
      {
        dCEX <- 0.75
        iLen <- nchar( strTaxon )
        if( iLen > 120 )
        {
          dLine <- 2.5
          i <- round( iLen / 2 )
          strTaxon <- paste( substring( strTaxon, 0, i ), substring( strTaxon, i + 1 ), sep = "\n" )
          adMar[2] <- adMar[2] + 1
        }
      } else { dCEX = 1 }
			
      if( class( adCur ) == "factor" )
      {
        if( "NA" %in% levels( adCur ) )
        {
          afNA <- adCur == "NA"
          adData <- adData[!afNA]
          adCur <- adCur[!afNA]
          adCur <- factor( adCur, levels = setdiff( levels( adCur ), "NA" ) )
        }
        astrNames <- c()
        astrColors <- c()
        dMed <- median( adData[adCur == levels( adCur )[1]], na.rm = TRUE )
        adIQR <- quantile( adData, probs = c(0.25, 0.75), na.rm = TRUE )
        dIQR <- adIQR[2] - adIQR[1]
        if( dIQR <= 0 )
        {
          dIQR <- sd( adData, na.rm = TRUE )
        }
        dIQR <- dIQR / 2

        #Print boxplots/strip charts of raw data. Add model data to it.
        for( strLevel in levels( adCur ) )
        {
          astrNames <- c(astrNames, sprintf( "%s (%d)", strLevel, sum( adCur == strLevel, na.rm = TRUE ) ))
          astrColors <- c(astrColors, sprintf( "%sAA", funcColor( ( median( adData[adCur == strLevel], na.rm = TRUE ) - dMed ) /
            dIQR, dMax = 3, dMin = -3, adMax = adColorMin, adMin = adColorMax, adMed = adColorMed ) ))
        }
        boxplot( adData ~ adCur, notch = TRUE, names = astrNames, mar = adMar, col = astrColors,
          main = strTitle, xlab = strCur, ylab = NA, cex.lab = dCEX, outpch = 4, outcex = 0.5 )
        stripchart( adData ~ adCur, add = TRUE, col = astrColors, method = "jitter", vertical = TRUE, pch = 20 )
        title( ylab = strTaxon, cex.lab = dCEX, line = dLine )
      } else {
        fGenetics <- length( aiGenetics ) && ( class( adCur ) == "integer" ) &&
          length( intersect( astrFactors, colnames( frmeData )[aiGenetics] ) )
        if( fGenetics )
        {
          astrLabels <- c()
          for( i in 0:2 )
          {
            astrLabels <- c(astrLabels, sprintf( "%d (%d)", i, sum( adCur == i, na.rm = TRUE ) ))
          }
          adCur <- adCur + rnorm( length( adCur ), sd = 0.05 )
        }
        plot( adCur, adData, mar = adMar, main = strTitle, xlab = strCur, pch = 20,
          col = sprintf( "%s99", funcGetColor( ) ), ylab = NA, xaxt = ifelse( fGenetics, "n", "s" ) )
        if( fGenetics )
        {
          axis( 1, at = 0:2, labels = astrLabels )
        }
        title( ylab = strTaxon, cex.lab = dCEX )
        lmod <- lm( adData ~ adCur )
        dColor <- lmod$coefficients[2] * mean( adCur, na.rm = TRUE ) / mean( adData, na.rm = TRUE )
#        dColor <- lsCur$value / lsCur$std
        strColor <- sprintf( "%sDD", funcColor( dColor, adMax = adColorMin, adMin = adColorMax, adMed = adColorMed ) )
        abline( reg = lmod, col = strColor, lwd = 3 )
      }
    }
    if( dev.cur( ) != 1 ) { dev.off( ) }
  }
  if( length( aiGenetics ) )
  {
    aiTmp <- aiGenetics
  } else {
    aiTmp <- aiData
  }
  print("End funcBugs")
  aiReturnBugs = aiTmp[colnames( frmeData )[aiTmp] %in% astrRet]
#  return( aiTmp[colnames( frmeData )[aiTmp] %in% astrRet] )
  return(list(aiReturnBugs=aiReturnBugs, lsQCCounts=lsData$lsQCCounts))
}

# Writes summary information to a file
#lmod Linear model information
#frmeData
#iTaxon Integer index of taxon (column)
funcBugResult = function( lmod, frmeData, iTaxon, dSig, adP, lsSig, strLog = NA, lsQCCounts, astrCols = c() )
{
#  print("Start funcBugResult")
  #Validate parameters
  #Exclude none and errors
  if( !is.na( lmod ) && ( class( lmod ) != "try-error" ) )
  {
    #Get the column name of the iTaxon index
    strTaxon = colnames( frmeData )[iTaxon]
    #Get summary information from the linear model
    lsSum = summary( lmod )
    #Write summary information to log file
    if( !is.na( strLog ) )
    {
      funcWrite( "#model", strLog )
      funcWrite( lmod, strLog )
      funcWrite( "#summary", strLog )
      funcWrite( lsSum, strLog )
    }

    #Get the coefficients
    if( is.null( coefficients( lsSum ) ) )
    {
      adCoefs = coefficients( lmod )
    } else {
      adCoefs = coefficients( lsSum )[,1]
    }

    #Go through each coefficient
    astrRows <- names( adCoefs )
    for( iMetadata in 1:length( astrRows ) )
    {
      #Select 
      strOrig = astrRows[iMetadata]
      #Skip y interscept
      if( strOrig == "(Intercept)" ) { next }
      
      if( "mboost" %in% class( lmod ) )
      {
        if( lsSum$selprob[iMetadata] > dSig )
        {
          dP = 1e-10 * ( 1 - lsSum$selprob[iMetadata] ) / ( 1 - dSig )
        } else {
          dP = 1e-10 + ( ( dSig - lsSum$selprob[iMetadata] ) / dSig )
        }
        dStd = 0
      } else {
        dP = lsSum$coefficients[strOrig,4]
        dStd = lsSum$coefficients[strOrig,2]
      }
      if( is.na( dP ) ) { next }

      dCoef = adCoefs[iMetadata]
      if( strOrig == "aiAlleles" )
      {
        strMetadata = strOrig
        adMetadata = aiAlleles
      } else if( length( grep( ":aiAlleles", strOrig, fixed = TRUE ) ) ){
        strMetadata = "interaction"
        adMetadata = aiAlleles
      } else {
        strMetadata = funcCoef2Col( strOrig, frmeData, astrCols )
        if( is.na( strMetadata ) )
        {
          if( substring( strOrig, nchar( strOrig ) - 1 ) == "NA" ) { next }
          print( c(strOrig, strMetadata) )
        }
        if( substring( strOrig, nchar( strMetadata ) + 1 ) == "NA" ){ next }
        adMetadata <- frmeData[,strMetadata]
      }
      #Bonferonni correct the factor p-values based on the factor levels-1 comparisons
      if( class( adMetadata ) == "factor" )
      {
        dP <- dP * ( nlevels( adMetadata ) - 1 )
      }
      adP <- c(adP, dP)
      lsSig[[length( lsSig ) + 1]] <- list(
        name = strMetadata,
        orig = strOrig,
        taxon = strTaxon,
        data = frmeData[,iTaxon],
        factors = c(strMetadata),
        metadata = adMetadata,
        value = dCoef,
        std = dStd )
    }
  }
  return( list(adP=adP, lsSig=lsSig, lsQCCounts=lsQCCounts) )
}

#########
#Perform or reverse Arcsin transformation on data adX
#adX Data to be transformed
#fReverse Boolean indicating reversing (FALSE) or performing (TRUE)
#funcTransform <- function( adX, fReverse = FALSE )
#{
#  if(is.null(adX)){return( 0 )}
#  if(fReverse)
#  {
#    return(sin(adX)^2 )
#  }
#  return(asin(sqrt(adX)))
#}

#Searches for the coeffient name or, if a factor, coefficient name and factor combination
#In the dataframe column names. If found returned, if not NA is returned.
#If intercept. NULL is returned.
#StrCoef String coefficient name
#frmeData Data frame of data
#aStrCols Column names of interest (if NULL is given, all column names are inspected.
funcCoef2Col <- function( strCoef, frmeData, astrCols = c() )
{
  #If the coefficient is the intersept there is no data column to return so return null
  if( strCoef == "(Intercept)" )
  {
    return( NULL )
  }
  #Remove ` from coefficient
  strCoef <- gsub( "`", "", strCoef )

  #If the coefficient name is not in the data frame
  if( !( strCoef %in% colnames( frmeData ) ) )
  {
    fHit <- FALSE
    #If the column names are not provided, use the column names of the dataframe.
    if( is.null( astrCols ) )
    {
      astrCols <- colnames( frmeData )
    }
    #Search through the different column names (factors)
    for( strFactor in astrCols )
    {
      #Select a column, if it is not a factor or does not begin with the factor's name then skip
      adCur <- frmeData[,strFactor]
      if( ( class( adCur ) != "factor" ) ||
        ( substr( strCoef, 1, nchar( strFactor ) ) != strFactor ) )
      { next }

      #For the factors, create factor-level name combinations to read in factors
      #Then check to see the factor-level combination is the coeffient of interest
      #If it is then store that factor as the coefficient of interest
      #And break
      for( strValue in levels( adCur ) )
      {
        strCur <- paste( strFactor, strValue, sep = "" )
        if( strCur == strCoef )
        {
	  strCoef <- strFactor
          fHit <- TRUE
          break
        }
      }
      #If the factor was found, return
      if( fHit )
      {  break }
    }
  }
  #If the original coefficient or the coeficient factor combination name are in the
  #data frame, return the name. Otherwise return NA.
  return( ifelse( ( strCoef %in% colnames( frmeData ) ), strCoef, NA ) )
}
