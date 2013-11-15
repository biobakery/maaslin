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
##description<< Allows one to plug in new modules to perform analysis (univariate or multivariate), regularization, and data (response) transformation.
) { return( pArgs ) }

# Libraries
suppressMessages(library( agricolae, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Needed for the pot-hoc Kruskal wallis comparisons
suppressMessages(library( penalized, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Needed for stepAIC
suppressMessages(library( MASS, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Needed for na action behavior
suppressMessages(library( gam, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Needed for boosting
suppressMessages(library( gbm, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Needed for LASSO
suppressMessages(library( glmnet, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Needed for mixed models
#suppressMessages(library( lme4, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
suppressMessages(library( nlme, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))

# Needed for zero inflated models
#suppressMessages(library( MCMCglmm, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
suppressMessages(library( pscl, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
suppressMessages(library( gamlss, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Do not use #suppressMessages(library( glmmADMB, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))

fAddBack = TRUE
dUnevenMax = .9


### Helper functions
# OK
funcMakeContrasts <- function
### Makes univariate contrasts of all predictors in the model formula with the response.
(strFormula, 
### lm style string defining reponse and predictors
strRandomFormula,
### mixed model string defining the fixed covariates
frmeTmp,
### The data frame to find predictor data in
iTaxon,
### Taxon
functionContrast,
### functionContrast The univariate test to perform
lsQCCounts
### QC info
){
  #TODO are we updating the QCCounts?
  lsSig = list()
  ### Holds all the significance results from the tests
  adP = c()
  ### Holds the p-values
  sCurDataName = names(frmeTmp)[iTaxon]
  ### The name of the taxon (data row) that is being associated (always assumed to be numeric)
  #Get test comparisons (predictor names from formula string)
  asComparisons  = unique(c(funcFormulaStrToList(strFormula),funcFormulaStrToList(strRandomFormula)))

  #Change metadata in formula to univariate comparisons
  for(sComparison in asComparisons)
  {
    # Metadata values
    vxTest = frmeTmp[[sComparison]]

    # Get the levels in the comparison
    # Can ignore the first level because it is the reference level
    asLevels = sComparison
    if(is.factor(vxTest)){asLevels = levels(vxTest)[2:length(vxTest)]}

    lasComparisonResults = functionContrast(x=sComparison, adCur=frmeTmp[[sCurDataName]], dfData=frmeTmp)
    for(asComparison in lasComparisonResults)
    {
      if( is.na( asComparison$p.value ) ) { next }
      # Get pvalue
      adP = c(adP, asComparison$p.value)
      # Get SD, if not available, give SD of covariate
      dSTD = asComparison$SD
      # TODO Is using sd on factor and binary data correct?
      if(is.na(dSTD) || is.null(dSTD)){dSTD = sd(vxTest)}

      lsSig[[length( lsSig ) + 1]] <- list(
        #Current metadata name (string column name) ok
        name = sComparison,
        #Current metadatda name (string, as a factor level if existing as such) ok
        orig = asComparison$name,
        #Taxon feature name (string) ok
        taxon = colnames( frmeTmp )[iTaxon],
        #Taxon data / response (double vector) ok
        data = frmeTmp[,iTaxon],
        #Name of column ok
        factors = sComparison,
        #Metadata values (metadata as a factor or raw numeric) ok
        metadata = vxTest,
        #Current coefficient value (named coef value with level name (from coefs) ok
        value = asComparison$coef,
        #Standard deviation (numeric) ok
        std = dSTD,
        #Model coefficients (output from coefs with intercept) ok
        allCoefs = asComparison$coef)
    }
  }
  return(list(adP=adP, lsSig=lsSig, lsQCCounts=lsQCCounts))
  ### Returns a list of p-value, standard deviation, and comparison which produced the p-value
}

#Ok
funcGetStepPredictors <- function
### Retrieve the predictors of the reduced model after stepwise model selection
(lmod,
### Linear model resulting from step-wise model selection
frmeTmp,
strLog
### File to document logging
){
  #Document
  funcWrite( "#model", strLog )
  funcWrite( lmod$fit, strLog )
#TODO  funcWrite( lmod$train.error, strLog )
#TODO  funcWrite( lmod$Terms, strLog )
  funcWrite( "#summary-gbm", strLog )
  funcWrite( summary(lmod), strLog )

  #Get Names from coefficients
  asStepCoefsFactors = coefficients(lmod)
  astrCoefNames = setdiff(names(asStepCoefsFactors[as.vector(!is.na(asStepCoefsFactors))==TRUE]),"(Intercept)")
  asStepCoefsFactors = unique(as.vector(sapply(astrCoefNames,funcCoef2Col, frmeData=frmeTmp)))

  if(length(asStepCoefsFactors)<1){ return(NA) }
  return( asStepCoefsFactors )
  ### Vector of string predictor names
}

funcGetUnivariateResults <- function
### Reduce the list of list of results to the correct format
( llmod,
### The list of univariate models
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
  adP = c()
  lsSig = list()
  for(lmod in llmod)
  {
    adP = c(adP,lmod$adP)
    lsSig = c(lsSig,lmod$lsSig)
  }
  return(list(adP=adP, lsSig=lsSig, lsQCCounts=llmod[length(llmod)]$lsQCCounts))
}

# OK
funcGetLMResults <- function
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
  for( lmod in llmod )
  {
    ilmodIndex = ilmodIndex + 1
    lmod = llmod[[ ilmodIndex ]]
    iTaxon = liTaxon[[ ilmodIndex ]]
    astrCols = lastrCols[[ ilmodIndex ]]

    #Exclude none and errors
    if( !is.na( lmod ) && ( class( lmod ) != "try-error" ) )
    {
      #holds the location of the pvlaues if an lm, if lmm is detected this will be changed
      iPValuePosition = 4

      #Get the column name of the iTaxon index
      #frmeTmp needs to be what?
      strTaxon = colnames( frmeData )[iTaxon]
      #Get summary information from the linear model
      lsSum = try( summary( lmod ) )
      #The following can actually happen when the stranger regressors return broken results
      if( class( lsSum ) == "try-error" )
      {
        next
      }

      #Write summary information to log file
      funcWrite( "#model", strLog )
      funcWrite( lmod, strLog )
      funcWrite( "#summary", strLog )
      #Unbelievably, some of the more unusual regression methods crash out when _printing_ their results 
      try( funcWrite( lsSum, strLog ) )

      #Get the coefficients
      #This should work for linear models
      frmeCoefs <- try( coefficients( lsSum ) )
      
      if( ( class(frmeCoefs ) == "try-error" ) || is.null( frmeCoefs ) )
      {
        adCoefs = try( coefficients( lmod ))
        if(class( adCoefs ) == "try-error")
        {
          adCoefs = coef(lmod)
        }
        frmeCoefs <- NA
      } else {
        if( class( frmeCoefs ) == "list" )
        {
          frmeCoefs <- frmeCoefs$count
        }
        adCoefs = frmeCoefs[,1]
      }

      #Go through each coefficient
      astrRows <- names( adCoefs )

      ##lmm
      if( is.null( astrRows ) )
      {
        astrRows = rownames( lsSum$tTable )
        frmeCoefs = lsSum$tTable
        iPValuePosition = 5
        adCoefs = frmeCoefs[,1]
      }

      for( iMetadata in 1:length( astrRows ) )
      {
        #Current coef which is being evaluated 
        strOrig = astrRows[ iMetadata ]
        #Skip y interscept
        if( strOrig %in% c("(Intercept)", "Intercept", "Log(theta)") ) { next }
        #Skip suppressed covariates
        if( funcCoef2Col( strOrig, frmeData ) %in% asSuppressCovariates){ next }

        #Extract pvalue and std in standard model
        dP = frmeCoefs[ strOrig, iPValuePosition ]
        dStd = frmeCoefs[ strOrig, 2 ]

        #Attempt to extract the pvalue and std in mixed effects summary 
        #Could not get the pvalue so skip result
        if( is.nan( dP ) || is.na( dP ) || is.null( dP ) ) { next }

        dCoef = adCoefs[ iMetadata ]

        #Setting adMetadata
        #Metadata values
        strMetadata = funcCoef2Col( strOrig, frmeData, astrCols )
        if( is.na( strMetadata ) )
        {
          if( substring( strOrig, nchar( strOrig ) - 1 ) == "NA" ) { next }
          c_logrMaaslin$error( "Unknown coefficient: %s", strOrig )
        }
        if( substring( strOrig, nchar( strMetadata ) + 1 ) == "NA" ) { next }
        adMetadata <- frmeData[, strMetadata ]

        # Store (factor level modified) p-value
        # Store general results for each coef
        adP <- c( adP, dP )

         # Append to the list of information about associations
        lsSig[[ length( lsSig ) + 1 ]] <- list(
          # Current metadata name
          name		= strMetadata,
          # Current metadatda name (as a factor level if existing as such)
          orig		= strOrig,
          # Taxon feature name
          taxon		= strTaxon,
          # Taxon data / response
          data		= frmeData[, iTaxon ],
          # All levels
          factors	= c( strMetadata ),
          # Metadata values
          metadata	= adMetadata,
          # Current coeficient value
          value		= dCoef,
          # Standard deviation
          std		= dStd,
          # Model coefficients
          allCoefs	= adCoefs )
      }
    }
  }
  return( list( adP = adP, lsSig = lsSig, lsQCCounts = lsQCCounts ) )
  ### List containing a list of pvalues, a list of significant data per association, and a list of QC data
}

funcGetZeroInflatedResults <- function
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
      if(! is.null( frmeCoefs$coefficients$count ) ) # Needed for zeroinfl
      {
        frmeCoefs = frmeCoefs$coefficients$count
      }
      adCoefs = frmeCoefs[,1]
      names(adCoefs) = row.names(frmeCoefs)

      funcWrite( "#Coefs", strLog )
      funcWrite( frmeCoefs, strLog )

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
        if(is.nan(dP)){next}
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
        print("adMetadata")
        print(adMetadata)

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


### Options for variable selection 
# OK
funcBoostModel <- function(
### Perform model selection / regularization with boosting
strFormula,
### The formula of the full model before boosting
frmeTmp,
### The data on which to perform analysis
adCur,
### The response data
lsParameters,
### User controlled parameters needed specific to boosting
lsForcedParameters = NULL,
### Force these predictors to be in the model
strLog
### File to which to document logging
){
  funcWrite( c("#Boost formula", strFormula), strLog )
  lmod = try( gbm( as.formula( strFormula ), data=frmeTmp, distribution="laplace", verbose=FALSE, n.minobsinnode=min(10, round(0.1 * nrow( frmeTmp ) ) ), n.trees=1000 ) )
#TODO#  lmod = try( gbm( as.formula( strFormula ), data=frmeTmp, distribution="gaussian", verbose=FALSE, n.minobsinnode=min(10, round(0.1 * nrow( frmeTmp ) ) ), n.trees=1000 ) )

  astrTerms <- c()
  if( !is.na( lmod ) && ( class( lmod ) != "try-error" ) )
  {
    #Get boosting summary results
    lsSum <- summary( lmod, plotit = FALSE )

    #Document
    funcWrite( "#model-gbm", strLog )
    funcWrite( lmod$fit, strLog )
    funcWrite( lmod$train.error, strLog )
    funcWrite( lmod$Terms, strLog )
    funcWrite( "#summary-gbm", strLog )
    funcWrite( lsSum, strLog )

    # Uneven metadata
    vstrUneven = c()
    # Kept metadata
    vstrKeepMetadata = c()

    #Select model predictors
    #Check the frequency of selection and skip if not selected more than set threshold dFreq
    for( strMetadata in lmod$var.names )
    {
      #Get the name of the metadata
      strTerm <- funcCoef2Col( strMetadata, frmeTmp, c(astrMetadata, astrGenetics) )

      #Add back in uneven metadata
      if(fAddBack)
      {
        ldMetadata = frmeTmp[[strMetadata]]
        if(length(which(table(ldMetadata)/length(ldMetadata)>dUnevenMax))>0)
        {
          astrTerms <- c(astrTerms, strTerm)
          vstrUneven = c(vstrUneven,strMetadata)
          next
        }
      }

      #If the selprob is less than a certain frequency, skip
      dSel <- lsSum$rel.inf[which( lsSum$var == strMetadata )] / 100
      if( is.na(dSel) || ( dSel <= lsParameters$dFreq ) ){ next }
#TODO#      if( is.na(dSel) || ( dSel < lsParameters$dFreq ) ){ next }

      #If you should ignore the metadata, continue
      if( is.null( strTerm ) ) { next }

      #If you cant find the metadata name, write
      if( is.na( strTerm ) )
      {
        c_logrMaaslin$error( "Unknown coefficient: %s", strMetadata )
        next
      }

      #Collect metadata names
      astrTerms <- c(astrTerms, strTerm)
      vstrKeepMetadata = c(vstrKeepMetadata,strTerm)
    }
  } else { astrTerms = lsForcedParameters }

#  funcBoostInfluencePlot(vdRelInf=lsSum$rel.inf, sFeature=lsParameters$sBugName, vsPredictorNames=lsSum$var, vstrKeepMetadata=vstrKeepMetadata, vstrUneven=vstrUneven)

  return(unique(c(astrTerms,lsForcedParameters)))
  ### Return a vector of predictor names to use in a reduced model
}

#Glmnet default is to standardize the variables.
#used as an example for implementation
#http://r.789695.n4.nabble.com/estimating-survival-times-with-glmnet-and-coxph-td4614225.html
funcPenalizedModel <- function(
### Perform penalized regularization for variable selection
strFormula,
### The formula of the full model before boosting
frmeTmp,
### The data on which to perform analysis
adCur,
### The response data
lsParameters,
### User controlled parameters needed specific to boosting
lsForcedParameters = NULL,
### Force these predictors to be in the model
strLog
### File to which to document logging
){
  #Convert the data frame to a model matrix
  mtrxDesign = model.matrix(as.formula(strFormula), data=frmeTmp)

  #Cross validate the lambda
  cvRet = cv.glmnet(x=mtrxDesign,y=adCur,alpha=lsParameters$dPAlpha)

  #Perform lasso
  glmnetMod = glmnet(x=mtrxDesign,y=adCur,family=lsParameters$family,alpha=lsParameters$dPAlpha,lambda=cvRet$lambda.min)

  #Get non zero beta and return column names for covariate names.
  ldBeta = glmnetMod$beta[,which.max(glmnetMod$dev.ratio)]
  ldBeta = names(ldBeta[which(abs(ldBeta)>0)])
  return(sapply(ldBeta,funcCoef2Col,frmeData=frmeTmp))
}

# OK
funcForwardModel <- function(
### Perform model selection with forward stepwise selection
strFormula,
### lm style string defining reposne and predictors 
frmeTmp,
### Data on which to perform analysis
adCur,
### Response data
lsParameters,
### User controlled parameters needed specific to boosting
lsForcedParameters = NULL,
### Force these predictors to be in the model
strLog
### File to which to document logging
){
  funcWrite( c("#Forward formula", strFormula), strLog )

  strNullFormula = "adCur ~ 1"
  if(!is.null(lsForcedParameters))
  {
    strNullFormula = paste( "adCur ~", paste( sprintf( "`%s`", lsForcedParameters ), collapse = " + " ))
  }
  lmodNull <- try( lm(as.formula( strNullFormula ), data=frmeTmp))
  lmodFull <- try( lm(as.formula( strFormula ), data=frmeTmp ))
  if(!("try-error" %in% c(class( lmodNull ),class( lmodFull ))))
  {
    lmod = stepAIC(lmodNull, scope=list(lower=lmodNull,upper=lmodFull), direction="forward", trace=0)
    return(funcGetStepPredictors(lmod, frmeTmp, strLog))
  }
  return( lsForcedParameters )
  ### Return a vector of predictor names to use in a reduced model or NA on error
}

# OK
# Select model with backwards selection
funcBackwardsModel <- function(
### Perform model selection with backwards stepwise selection
strFormula,
### lm style string defining reponse and predictors 
frmeTmp,
### Data on which to perform analysis
adCur,
### Response data
lsParameters,
### User controlled parameters needed specific to boosting
lsForcedParameters = NULL,
### Force these predictors to be in the model
strLog
### File to which to document logging
){
  funcWrite( c("#Backwards formula", strFormula), strLog )

  strNullFormula = "adCur ~ 1"
  if(!is.null(lsForcedParameters))
  {
    strNullFormula = paste( "adCur ~", paste( sprintf( "`%s`", lsForcedParameters ), collapse = " + " ))
  }

  lmodNull <- try( lm(as.formula( strNullFormula ), data=frmeTmp))
  lmodFull <- try( lm(as.formula( strFormula ), data=frmeTmp ))

  if(! class( lmodFull ) == "try-error" )
  {
    lmod = stepAIC(lmodFull, scope=list(lower=lmodNull, upper=lmodFull), direction="backward")
    return(funcGetStepPredictors(lmod, frmeTmp, strLog))
  } else { 
    return( lsForcedParameters ) }
  ### Return a vector of predictor names to use in a reduced model or NA on error
}

### Analysis methods
### Univariate options

# Sparse Dir. Model
#TODO# Implement in sfle

# Tested
# Correlation
# NOTE: Ignores the idea of random and fixed covariates
funcSpearman <- function(
### Perform multiple univariate comparisons producing spearman correlations for association
strFormula,
### lm style string defining reponse and predictors, for mixed effects models this holds the fixed variables
frmeTmp,
### Data on which to perform analysis
iTaxon,
### Index of the response data
lsQCCounts,
### List recording anything important to QC
strRandomFormula = NULL
### Has the formula for random covariates
){
  return(funcMakeContrasts(strFormula=strFormula, strRandomFormula=strRandomFormula, frmeTmp=frmeTmp, iTaxon=iTaxon,
    functionContrast=function(x,adCur,dfData)
    {
      retList = list()
      ret = cor.test(as.formula(paste("~",x,"+ adCur")), data=dfData, method="spearman", na.action=c_strNA_Action)
      #Returning rho for the coef in a named vector
      vdCoef = c()
      vdCoef[[x]]=ret$estimate
      retList[[1]]=list(p.value=ret$p.value,SD=sd(dfData[[x]]),name=x,coef=vdCoef)
      return(retList)
    }, lsQCCounts))
  ### List of contrast information, pvalue, contrast and std per univariate test
}

# Tested
# Wilcoxon (T-Test)
# NOTE: Ignores the idea of random and fixed covariates
funcWilcoxon <- function(
### Perform multiple univariate comparisons performing wilcoxon tests on discontinuous data with 2 levels
strFormula,
### lm style string defining reponse and predictors, for mixed effects models this holds the fixed variables
frmeTmp,
### Data on which to perform analysis
iTaxon,
### Index of the response data
lsQCCounts,
### List recording anything important to QC
strRandomFormula = NULL
### Has the formula for random covariates
){
  return(funcMakeContrasts(strFormula=strFormula, strRandomFormula=strRandomFormula, frmeTmp=frmeTmp, iTaxon=iTaxon,
    functionContrast=function(x,adCur,dfData)
    {
      retList = list()
      ret = wilcox.test(as.formula(paste("adCur",x,sep=" ~ ")), data=dfData, na.action=c_strNA_Action)
      #Returning NA for the coef in a named vector
      vdCoef = c()
      vdCoef[[x]]=ret$statistic
      retList[[1]]=list(p.value=ret$p.value,SD=sd(dfData[[x]]),name=x,coef=vdCoef)
      return(retList)
    }, lsQCCounts))
  ### List of contrast information, pvalue, contrast and std per univariate test
}

# Tested
# Kruskal.Wallis (Nonparameteric anova)
# NOTE: Ignores the idea of random and fixed covariates
funcKruskalWallis <- function(
### Perform multiple univariate comparisons performing Kruskal wallis rank sum tests on discontuous data with more than 2 levels
strFormula,
### lm style string defining reponse and predictors, for mixed effects models this holds the fixed variables
frmeTmp,
### Data on which to perform analysis
iTaxon,
### Index of the response data
lsQCCounts,
### List recording anything important to QC
strRandomFormula = NULL
### Has the formula for random covariates
){
  return(funcMakeContrasts(strFormula=strFormula, strRandomFormula=strRandomFormula, frmeTmp=frmeTmp, iTaxon=iTaxon,
    functionContrast=function(x,adCur,dfData)
    {
      retList = list()
      lmodKW = kruskal(adCur,dfData[[x]],group=FALSE,p.adj="holm")

      asLevels = levels(dfData[[x]])
      # The names of the generated comparisons, sometimes the control is first sometimes it is not so
      # We will just check which is in the names and use that
      asComparisons = row.names(lmodKW$comparisons)
      #Get the comparison with the control
      for(sLevel in asLevels[2:length(asLevels)])
      {
        sComparison = intersect(c(paste(asLevels[1],sLevel,sep=" - "),paste(sLevel,asLevels[1],sep=" - ")),asComparisons)
        #Returning NA for the coef in a named vector
        vdCoef = c()
        vdCoef[[paste(x,sLevel,sep="")]]=lmodKW$comparisons[sComparison,"Difference"]
#        vdCoef[[paste(x,sLevel,sep="")]]=NA
        retList[[length(retList)+1]]=list(p.value=lmodKW$comparisons[sComparison,"p.value"],SD=1.0,name=paste(x,sLevel,sep=""),coef=vdCoef)
      }
      return(retList)
    }, lsQCCounts))
  ### List of contrast information, pvalue, contrast and std per univariate test
}

# Tested
# NOTE: Ignores the idea of random and fixed covariates
funcDoUnivariate <- function(
### Perform multiple univariate comparisons producing spearman correlations for association
strFormula,
### lm style string defining reponse and predictors, for mixed effects models this holds the fixed variables
frmeTmp,
### Data on which to perform analysis
iTaxon,
### Index of the response data
lsHistory,
### List recording p-values, association information, and QC counts
strRandomFormula = NULL,
### Has the formula for random covariates
fZeroInflate = FALSE
){
  if(fZeroInflate)
  {
    throw("There are no zero-inflated univariate models to perform your analysis.")
  }

  # Get covariates
  astrCovariates = unique(c(funcFormulaStrToList(strFormula),funcFormulaStrToList(strRandomFormula)))

  # For each covariate
  for(sCovariate in astrCovariates)
  {
    ## Check to see if it is discrete
    axData = frmeTmp[[sCovariate]]
    lsRet = NA
    if(is.factor(axData) || is.logical(axData))
    {
      ## If discrete check how many levels
      lsDataLevels = levels(axData)
      ## If 2 levels do wilcoxon test
      if(length(lsDataLevels) < 3)
      {
        lsRet = funcWilcoxon(strFormula=paste("adCur",sCovariate,sep=" ~ "), frmeTmp=frmeTmp, iTaxon=iTaxon, lsQCCounts=lsHistory$lsQCCounts)
      } else {
      ## If 3 or more levels do kruskal wallis test
        lsRet = funcKruskalWallis(strFormula=paste("adCur",sCovariate,sep=" ~ "), frmeTmp=frmeTmp, iTaxon=iTaxon, lsQCCounts=lsHistory$lsQCCounts)
      }
    } else {
      ## If not discrete do spearman test (list(adP=adP, lsSig=lsSig, lsQCCounts=lsQCCounts))
      lsRet = funcSpearman(strFormula=paste("adCur",sCovariate,sep=" ~ "), frmeTmp=frmeTmp, iTaxon=iTaxon, lsQCCounts=lsHistory$lsQCCounts)
    }
    lsHistory[["adP"]] = c(lsHistory[["adP"]], lsRet[["adP"]])
    lsHistory[["lsSig"]] = c(lsHistory[["lsSig"]], lsRet[["lsSig"]])
    lsHistory[["lsQCCounts"]] = lsRet[["lsQCCounts"]]
  }
  return(lsHistory)
}

### Multivariate

# Tested
funcLM <- function(
### Perform vanilla linear regression
strFormula,
### lm style string defining reponse and predictors, for mixed effects models this holds the fixed variables
frmeTmp,
### Data on which to perform analysis
iTaxon,
### Index of the response data
lsHistory,
### List recording p-values, association information, and QC counts
strRandomFormula = NULL,
### Has the formula for random covariates
fZeroInflated = FALSE
### Turns on the zero inflated model
){
  adCur = frmeTmp[,iTaxon]
  if(fZeroInflated)
  {
    return(try(gamlss(formula=as.formula(strFormula), family=BEZI, data=frmeTmp))) # gamlss
  } else {
    if(!is.null(strRandomFormula))
    {
      return(try(glmmPQL(fixed=as.formula(strFormula), random=as.formula(strRandomFormula), family=gaussian(link="identity"), data=frmeTmp)))
      #lme4 package but does not have pvalues for the fixed variables (have to use a mcmcsamp/pvals.fnc function which are currently disabled)
      #return(try( lmer(as.formula(strFormula), data=frmeTmp, na.action=c_strNA_Action) ))
    } else {
      return(try( lm(as.formula(strFormula), data=frmeTmp, na.action=c_strNA_Action) ))
    }
  }
  ### lmod result object from lm
}

# Tested
funcBinomialMult <- function(
### Perform linear regression with negative binomial link
strFormula,
### lm style string defining reponse and predictors, for mixed effects models this holds the fixed variables
frmeTmp,
### Data on which to perform analysis
iTaxon,
### Index of the response data
lsHistory,
### List recording p-values, association information, and QC counts
strRandomFormula = NULL,
### Has the formula for random covariates
fZeroInflated = FALSE
### Turns on the zero inflated model
){
  adCur = frmeTmp[,iTaxon]
  if(fZeroInflated)
  {
    return(try(zeroinfl(as.formula(strFormula), data=frmeTmp, dist="negbin"))) # pscl
#    return(try(gamlss(as.formula(strFormula), family=ZINBI, data=frmeTmp))) # pscl
  } else {
    if(!is.null(strRandomFormula))
    {
      throw("This analysis flow is not completely developed, please choose an option other than negative bionomial with random covariates")
      #TODO need to estimate the theta
      #return(try(glmmPQL(fixed=as.formula(strFormula), random=as.formula(strRandomFormula), family=negative.binomial(theta = 2, link=log), data=frmeTmp)))
      #lme4 package but does not have pvalues for the fixed variables (have to use a mcmcsamp/pvals.fnc function which are currently disabled)
    } else {
      return(try( glm.nb(as.formula(strFormula), data=frmeTmp, na.action=c_strNA_Action )))
    }
  }
  ### lmod result object from lm
}

# Tested
funcQuasiMult <- function(
### Perform linear regression with quasi-poisson link
strFormula,
### lm style string defining reponse and predictors, for mixed effects models this holds the fixed variables
frmeTmp,
### Data on which to perform analysis
iTaxon,
### Index of the response data
lsHistory,
### List recording p-values, association information, and QC counts
strRandomFormula = NULL,
### Has the formula for random covariates
fZeroInflated = FALSE
### Turns on a zero infalted model
){
  adCur = frmeTmp[,iTaxon]
  if(fZeroInflated)
  {
#    return(try(gamlss(formula=as.formula(strFormula), family=ZIP, data=frmeTmp))) # gamlss
    return(try(zeroinfl(as.formula(strFormula), data=frmeTmp, dist="poisson"))) # pscl
  } else {
    #Check to see if | is in the model, if so use a lmm otherwise the standard glm is ok
    if(!is.null(strRandomFormula))
    {
      return(try(glmmPQL(fixed=as.formula(strFormula), random=as.formula(strRandomFormula), family= quasipoisson, data=frmeTmp)))
      #lme4 package but does not have pvalues for the fixed variables (have to use a mcmcsamp/pvals.fnc function which are currently disabled)
      #return(try ( glmer(as.formula(strFormula), data=frmeTmp, family=quasipoisson, na.action=c_strNA_Action) ))
    } else {
      return(try( glm(as.formula(strFormula), family=quasipoisson, data=frmeTmp, na.action=c_strNA_Action) ))
    }
  }
  ### lmod result object from lm
}

### Transformations
# Tested
funcArcsinSqrt <- function(
# Transform data with arcsin sqrt transformation
aData
### The data on which to perform the transformation
){
  return(asin(sqrt(aData)))
  ### Transformed data
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

# Tested
funcNoTransform <-function(
### Pass data without transform
aData
### The data on which to perform the transformation
### Only given here to preserve the pattern, not used.
){
  return(aData)
  ### Transformed data
}
