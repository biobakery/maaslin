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


### AnalysisModules ###

notfuncGetZeroInflatedResults <- function
### Reduce the lm object return to just the data needed for further analysis
( llmod,
### The result from a linear model
frmeData,
### Data analysis is performed on
liTaxon,
### The response id
dSig,
### Significance level for q-values
adP,
### List of pvalues from all associations performed
lsSig,
### List of information from the lm containing, metadata name, metadatda name (as a factor level if existing as such), Taxon feature name, Taxon data / response, All levels, Metadata values, Current coeficient value, Standard deviation, Model coefficients
strLog,
### File to which to document logging
lsQCCounts,
### Records of counts associated with quality control
lastrCols,
### Predictors used in the association
asSuppressCovariates=c()
### Vector of covariates to suppress and not give results for
){
  ilmodIndex = 0
  for(lmod in llmod)
  {
    ilmodIndex = ilmodIndex + 1
    lmod = llmod[[ilmodIndex]]
    iTaxon = liTaxon[[ilmodIndex]]
    astrCols = lastrCols[[ilmodIndex]]

    #Exclude none and errors
    if( !is.na( lmod ) && ( class( lmod ) != "try-error" ) )
    {
      #holds the location of the pvlaues if an lm, if lmm is detected this will be changed
      iPValuePosition = 4

      #Get the column name of the iTaxon index
      #frmeTmp needs to be what?
      strTaxon = colnames( frmeData )[iTaxon]

      #Write summary information to log file
      funcWrite( "#model", strLog )
      funcWrite( lmod, strLog )

      #Get the coefficients
      #This should work for linear models
      frmeCoefs <- summary(lmod)
      frmeCoefs = frmeCoefs$coefficients$count
      funcWrite( "#Coefs", strLog )
      funcWrite( frmeCoefs, strLog )

      adCoefs = frmeCoefs[,1]
      names(adCoefs) = row.names(frmeCoefs)

      #Go through each coefficient
      astrRows <- row.names( frmeCoefs )

      for( iMetadata in 1:length( astrRows ) )
      {
        #Current coef which is being evaluated 
        strOrig = astrRows[iMetadata]

        #Skip y interscept
        if( strOrig %in% c("(Intercept)", "Intercept", "Log(theta)") ) { next }
        #Skip suppressed covariates
        if( funcCoef2Col(strOrig,frmeData) %in% asSuppressCovariates){next}

        #Extract pvalue and std in standard model
        dP = frmeCoefs[strOrig, iPValuePosition]
        if(is.nan(dP)){ next }
        dStd = frmeCoefs[strOrig,2]
        dCoef = adCoefs[iMetadata]

        #Setting adMetadata
        #Metadata values
        strMetadata = funcCoef2Col( strOrig, frmeData, astrCols )
        if( is.na( strMetadata ) )
        {
          if( substring( strOrig, nchar( strOrig ) - 1 ) == "NA" ) { next }
          c_logrMaaslin$error( "Unknown coefficient: %s", strOrig )
        }
        if( substring( strOrig, nchar( strMetadata ) + 1 ) == "NA" ) { next }
        adMetadata <- frmeData[,strMetadata]

        #Store (factor level modified) p-value
        #Store general results for each coef
        adP <- c(adP, dP)
        lsSig[[length( lsSig ) + 1]] <- list(
          #Current metadata name
          name		= strMetadata,
          #Current metadatda name (as a factor level if existing as such)
          orig		= strOrig,#
          #Taxon feature name
          taxon		= strTaxon,
          #Taxon data / response
          data		= frmeData[,iTaxon],
          #All levels
          factors	= c(strMetadata),
          #Metadata values
          metadata	= adMetadata,
          #Current coeficient value
          value		= dCoef,
          #Standard deviation
          std		= dStd,
          #Model coefficients
          allCoefs	= adCoefs)
      }
    }
  }
  return(list(adP=adP, lsSig=lsSig, lsQCCounts=lsQCCounts))
  ### List containing a list of pvalues, a list of significant data per association, and a list of QC data
}


oldfuncGetZeroInflatedResults <- function
### Reduce the lm object return to just the data needed for further analysis
( llmod,
### The result from a linear model
frmeData,
### Data analysis is performed on
liTaxon,
### The response id
dSig,
### Significance level for q-values
adP,
### List of pvalues from all associations performed
lsSig,
### List of information from the lm containing, metadata name, metadatda name (as a factor level if existing as such), Taxon feature name, Taxon data / response, All levels, Metadata values, Current coeficient value, Standard deviation, Model coefficients
strLog,
### File to which to document logging
lsQCCounts,
### Records of counts associated with quality control
lastrCols,
### Predictors used in the association
asSuppressCovariates=c()
### Vector of covariates to suppress and not give results for
){
  ilmodIndex = 0
  for(lmod in llmod)
  {
    ilmodIndex = ilmodIndex + 1
    lmod = llmod[[ilmodIndex]]
    iTaxon = liTaxon[[ilmodIndex]]
    astrCols = lastrCols[[ilmodIndex]]

    #Exclude none and errors
    if( !is.na( lmod ) && ( class( lmod ) != "try-error" ) )
    {
      #holds the location of the pvlaues if an lm, if lmm is detected this will be changed
      iPValuePosition = 4

      #Get the column name of the iTaxon index
      #frmeTmp needs to be what?
      strTaxon = colnames( frmeData )[iTaxon]

      #Write summary information to log file
      funcWrite( "#model", strLog )
      funcWrite( lmod, strLog )

      #Get the coefficients
      #This should work for linear models
      frmeCoefs <- summary(lmod)
      adCoefs = frmeCoefs[,1]
      names(adCoefs) = row.names(frmeCoefs)

      #Go through each coefficient
      astrRows <- row.names( frmeCoefs )

      for( iMetadata in 1:length( astrRows ) )
      {
        #Current coef which is being evaluated 
        strOrig = astrRows[iMetadata]
        #Skip y interscept
        if( strOrig %in% c("(Intercept)", "Intercept", "Log(theta)") ) { next }
        #Skip suppressed covariates
        if( funcCoef2Col(strOrig,frmeData) %in% asSuppressCovariates){next}

        #Extract pvalue and std in standard model
        dP = frmeCoefs[strOrig, iPValuePosition]
        dStd = frmeCoefs[strOrig,2]

        dCoef = adCoefs[iMetadata]

        #Setting adMetadata
        #Metadata values
        strMetadata = funcCoef2Col( strOrig, frmeData, astrCols )
        if( is.na( strMetadata ) )
        {
          if( substring( strOrig, nchar( strOrig ) - 1 ) == "NA" ) { next }
          c_logrMaaslin$error( "Unknown coefficient: %s", strOrig )
        }
        if( substring( strOrig, nchar( strMetadata ) + 1 ) == "NA" ) { next }
        adMetadata <- frmeData[,strMetadata]

        #Store (factor level modified) p-value
        #Store general results for each coef
        adP <- c(adP, dP)
        lsSig[[length( lsSig ) + 1]] <- list(
          #Current metadata name
          name		= strMetadata,
          #Current metadatda name (as a factor level if existing as such)
          orig		= strOrig,#
          #Taxon feature name
          taxon		= strTaxon,
          #Taxon data / response
          data		= frmeData[,iTaxon],
          #All levels
          factors	= c(strMetadata),
          #Metadata values
          metadata	= adMetadata,
          #Current coeficient value
          value		= dCoef,
          #Standard deviation
          std		= dStd,
          #Model coefficients
          allCoefs	= adCoefs)
      }
    }
  }
  return(list(adP=adP, lsSig=lsSig, lsQCCounts=lsQCCounts))
  ### List containing a list of pvalues, a list of significant data per association, and a list of QC data
}

funcSquareSin <- function(
# Transform data with square sin transformation
# Opposite of the funcArcsinSqrt
aData
### The data on which to perform the transformation
){
  return(sin(aData)^2)
  ### Transformed data
}

### Utilities ###

funcTrim=function(
### Remove whitespace at the beginning or the end of a string
tempString
### tempString String to be trimmed.
){
  return(gsub("^\\s+|\\s+$","",tempString))
}

funcRename <- function(
### Modifies labels for plotting
### If the name is not an otu collapse to the last two clades
### Otherwise use the most terminal clade
astrNames
### Names to modify for plotting
){
  astrRet <- c()
  for( strName in astrNames )
  {
    astrName <- strsplit( strName, c_cFeatureDelimRex )[[1]]
    i <- length( astrName )
    if( ( astrName[i] == c_strUnclassified ) || !is.na( as.numeric( astrName[i] ) ) )
    {
      strRet <- paste( astrName[( i - 1 ):i], collapse = c_cFeatureDelim )
    } else {
    strRet <- astrName[i]
    }
    astrRet <- c(astrRet, strRet)
  }
  return( astrRet )
  ### List of modified names
}

funcMFAValue2Col = function(
### Given a value in a column, the column name is returned.
xValue,
dfData,
aiColumnIndicesToSearch = NULL
){
  lsColumnNames = names(dfData)

  if(is.null(aiColumnIndicesToSearch))
  {
    aiColumnIndicesToSearch = c(1:dim(dfData)[2])
  }

  # Could be the column name
  if(xValue %in% lsColumnNames){return(xValue)}

  # Could be the column name and value
  iValueLength = length(xValue)
  for( iColIndex in c(1:length(lsColumnNames) ))
  {
    adCur = dfData[[lsColumnNames[iColIndex]]]
    if(is.factor(adCur))
    {
      for(strValue in levels(adCur))
      {
        strCurVersion1 <- paste( lsColumnNames[iColIndex], strValue, sep = c_sMFANameSep1 )
        strCurVersion2 <- paste( lsColumnNames[iColIndex], strValue, sep = c_sMFANameSep2 )
        if((xValue == strCurVersion1) || (xValue == strCurVersion2)){return(lsColumnNames[iColIndex])}
      }
    }
  }

  # Could be the value
  for(iColIndex in aiColumnIndicesToSearch)
  {
    if(xValue %in% dfData[[lsColumnNames[iColIndex]]]){return(lsColumnNames[iColIndex])}
  }
  return(NULL)
}

funcLMToNoNAFormula <-function(
lMod,
frmeTmp,
adCur
){
  dfCoef = coef(lMod)
  astrCoefNames = setdiff(names(dfCoef[as.vector(!is.na(dfCoef))==TRUE]),"(Intercept)")
  astrPredictors = unique(as.vector(sapply(astrCoefNames,funcCoef2Col, frmeData=frmeTmp)))
  strFormula = paste( "adCur ~", paste( sprintf( "`%s`", astrPredictors ), collapse = " + " ), sep = " " )
  return(try( lm(as.formula( strFormula ), data=frmeTmp )))
}

funcGetRandomColors=function(
#Generates a given number of random colors
tempNumberColors = 1
### Number of colors to generate
){
  adRet = c()
  return(sapply(1:tempNumberColors, function(x){
    adRGB <- ( runif( 3 ) * 0.66 ) + 0.33
    adRet <- c(adRet, rgb( adRGB[1], adRGB[2], adRGB[3] ))
  }))
}

funcFormulaListToString <- function(
# Using covariate and random covariate names, creates a lm or mixed model formula
# returns a vector of c(strLM, strMixedModel), one will be NA given the existance of random covariates.
# On error c(NA,NA) is given
astrTerms,
#Fixed covariates or all covariates if using an lm
astrRandomCovariates = NULL
#Random covariates for a mixed model
){
  strRetLMFormula = NA
  strRetMMFormula = NA

  #If no covariates return NA
  if(is.null(astrTerms)){return(c(strRetLMFormula, strRetMMFormula))}

  #Get fixed covariates
  astrFixedCovariates = setdiff(astrTerms,astrRandomCovariates)

  #If no fixed coavariates return NA
  # Can not run a model with no fixed covariate, restriction of lmm
  if(length(astrFixedCovariates)==0){return(c(strRetLMFormula, strRetMMFormula))}

  # Fixed Covariates
  strFixedCovariates = paste( sprintf( "`%s`", astrFixedCovariates ), collapse = " + " )

  #If random covariates, set up a formula for mixed models
  if(length(astrRandomCovariates)>0)
  {
    #Format for lmer
    #strRetFormula <- paste( "adCur ~ ", paste( sprintf( "(1|`%s`))", intersect(astrRandomCovariates, astrTerms)), collapse = " + " ))
    #Format for glmmpql
    strRandomCovariates = paste( sprintf( "1|`%s`", setdiff(astrRandomCovariates, astrTerms)), collapse = " + " )
    strRetMMFormula <- paste( "adCur ~ ", strFixedCovariates, " + ", strRandomCovariates, sep="")
  } else {
    #This is either the formula for all covariates in an lm or fixed covariates in the lmm
    strRetLMFormula <- paste( "adCur ~ ", strFixedCovariates, sep="")
  }
  return(c(strRetLMFormula, strRetMMFormula))
}

funcColToMFAValue = function(
### Given a column name, return the MFA values that could be associated with the name
lsColNames,
### String list of column names (as you would get from names(dataframe))
dfData
### Data frame of data the column names refer to
){
  lsMFAValues = c()

  for(sColName in lsColNames)
  {
    axCur = dfData[[sColName]]

    if(is.logical(axCur)){axCur=as.factor(axCur)}
    if(is.factor(axCur))
    {
      lsLevels = levels(axCur)
      if((length(lsLevels)==2) && (!is.na(as.numeric(lsLevels[1]))) && (!is.na(as.numeric(lsLevels[2]))))
      {
        lsMFAValues = c(lsMFAValues,paste(sColName,lsLevels[1],sep=c_sMFANameSep1),paste(sColName,lsLevels[2],sep=c_sMFANameSep1))
      }else{
        for(sLevel in levels(axCur))
        {
          lsMFAValues = c(lsMFAValues,sLevel)
        }
      }
    } else {
      lsMFAValues = c(lsMFAValues,sColName)
    }
  }
  return(setdiff(lsMFAValues,c("NA",NA)))
}

funcColors <- function(
### Generate a range of colors
dMax = 1,
### Max possible value
dMin = -1,
### Min possible value
dMed = NA,
### Central value if you don't want to be the average
adMax = c(1, 1, 0),
### Is used to generate the color for the higher values in the range, this can be changed to give different colors set to green
adMin = c(0, 0, 1),
### Is used to generate the color for the lower values in the range, this can be changed to give different colors set to red
adMed = c(0, 0, 0),
### Is used to generate the color for the central values in the range, this can be changed to give different colors set to black
iSteps = 64
### Number of intermediary colors made in the range of colors
){
  lsTmp <- funcColorHelper( dMax, dMin, dMed )
  dMax <- lsTmp$dMax
  dMin <- lsTmp$dMin
  dMed <- lsTmp$dMed
  aRet <- c ()
  for( dCur in seq( dMin, dMax, ( dMax - dMin ) / ( iSteps - 1 ) ) )
  {
    aRet <- c(aRet, funcColor( dCur, dMax, dMin, dMed, adMax, adMin, adMed ))
  }
  return( aRet )
  ### List of colors
}

