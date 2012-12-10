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
suppressMessages(library( penalized, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Need for stepAIC
suppressMessages(library( MASS, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Needed for na action behavior
suppressMessages(library( gam, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Needed for boosting
suppressMessages(library( gbm, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Needed for mixed models
suppressMessages(library( lme4, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))


### Helper functions
# OK
funcMakeContrasts <- function
### Makes univariate contrasts of all predictors in the model formula with the response.
(strFormula, 
### lm style string defining reponse and predictors 
frmeTmp,
### The data frame to find predictor data in
adCur,
### adCur Response data
iTaxon,
### Taxon
functionContrast,
### functionContrast The univariate test to perform
lsQCCounts,
### QC info
fDummy = FALSE
### Indicates if dummy variables are need is needed (tests of heterogenous data types will be encountered but the test can not handle them)
){
#  print("Start funcMakeContrasts")
#print("strFormula")
#print(strFormula)
#print("adCur")
#print(adCur)
#print("iTaxon")
#print(iTaxon)
#print("functionContrast")
#print(functionContrast)
#print("fDummy")
#print(fDummy)


  #TODO are we updating the QCCounts?
  lsSig = list()
  ### Holds all the significance results from the tests
  adP = c()
  ### Holds the p-values

  #Get test comparisons (predictor names from formula string)
  asComparisons = gsub("`","",setdiff(unlist(strsplit(unlist(strsplit(strFormula,"~"))[2]," ")),c("","+")))

  #Change metadata in formula to univariate comparisons
  for(sComparison in asComparisons)
  {
    #Removed fixed covariate formating
    lsParse = unlist(strsplit(sComparison, "[\\(\\|\\)]", perl=FALSE))
    sComparison = lsParse[length(lsParse)]
 
    viTest = as.vector(frmeTmp[[sComparison]])
    if(is.character(viTest)){viTest = as.factor(viTest)}

    #Only create dummy variables on univariate comparisons that require homogenous data.
    #I am dummy variabling it because that is equivalent to the handling in the lms
#    if(is.factor(viTest)){viTest = as.numeric(viTest)}
    if(is.factor(viTest) && fDummy)
    {
      lsLevels = levels(viTest)
      for(sLevel in lsLevels)
      {
        liLevelTest = ifelse(viTest == sLevel, 1, 0)

        sComparisonResults = functionContrast(x=liLevelTest, y=adCur)
        dPvalue = sComparisonResults$p.value
        if( is.na( dPvalue ) ) { next }

        #Bonferonni correct the factor p-values based on the factor levels-1 comparisons
        adP = c(dPvalue * ( nlevels( viTest ) - 1 ),adP)

        lsSig[[length( lsSig ) + 1]] <- list(
          #Current metadata name
          name = sComparison,
          #Current metadatda name (as a factor level if existing as such)
          orig = sLevel,
          #Taxon feature name
          taxon = colnames( frmeTmp )[iTaxon],
          #Taxon data / response
          data = frmeTmp[,iTaxon],
          #All levels
          factors = levels(liLevelTest),
          #Metadata values
          metadata = liLevelTest,
          #Current coefficient value
          value = sComparison,
          #Standard deviation
          std = sd(liLevelTest),
          #Model coefficients
          allCoefs = sComparison)
      }
    } else {
      sComparisonResults = functionContrast(x=viTest, y=adCur)
      if( is.na( sComparisonResults$p.value ) ) { next }
      adP = c(sComparisonResults$p.value,adP)

      lsSig[[length( lsSig ) + 1]] <- list(
        #Current metadata name
        name = sComparison,
        #Current metadatda name (as a factor level if existing as such)
        orig = sComparison,
        #Taxon feature name
        taxon = colnames( frmeTmp )[iTaxon],
        #Taxon data / response
        data = frmeTmp[,iTaxon],
        #All levels
        factors = levels(viTest),
        #Metadata values
        metadata = viTest,
        #Current coefficient value
        value = sComparison,
        #Standard deviation
        std = sd(viTest),
        #Model coefficients
        allCoefs = sComparison)
    }
  }
#  print("Stop funcMakeContrasts")
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
# ok
funcGetLMResults <- function
### Reduce the lm object return to just the data needed for further analysis
( lmod=lmod,
### The result from a linear model
frmeData=frmeData,
### Data analysis is perfromed on
iTaxon=iTaxon,
### The response id
dSig=dSig,
### Significance level for q-values
adP=adP,
### List of pvalues from all associations performed
lsSig=lsSig,

strLog=strLog,
### File to which to document logging
lsQCCounts=lsData$lsQCCounts,
### Records of counts associated with quality control
astrCols=astrTerms
### Predictors used in the association
){
  #TODO are we updating the QCCounts?
  #Exclude none and errors
  if( !is.na( lmod ) && ( class( lmod ) != "try-error" ) )
  {
    #Get the column name of the iTaxon index
    #frmeTmp needs to be what?
    strTaxon = colnames( frmeData )[iTaxon]
    #Get summary information from the linear model
    lsSum = try( summary( lmod ) )
    #The following can actually happen when the stranger regressors return broken results
    if( class( lsSum ) == "try-error" )
    {
      return( list(adP=adP, lsSig=lsSig, lsQCCounts=lsQCCounts) )
    }

    #Write summary information to log file
    funcWrite( "#model", strLog )
    funcWrite( lmod, strLog )
    funcWrite( "#summary", strLog )
    #Unbelievably, some of the more unusual regression methods crash out when _printing_ their results 
    try( funcWrite( lsSum, strLog ) )

    #Get the coefficients
    frmeCoefs <- try( coefficients( lsSum ) )
    if( ( class( frmeCoefs ) == "try-error" ) || is.null( frmeCoefs ) )
    {
      adCoefs = coefficients( lmod )
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
    for( iMetadata in 1:length( astrRows ) )
    {
      #Current coef which is being evaluated 
      strOrig = astrRows[iMetadata]
      #Skip y interscept
      if( strOrig %in% c("(Intercept)", "Intercept", "Log(theta)") ) { next }

      dP = frmeCoefs[strOrig,4]
      dStd = frmeCoefs[strOrig,2]

      if( is.na( dP ) ) { next }

      dCoef = adCoefs[iMetadata]

      #Setting adMetadata
      #Metadata values
      #If it is a genetics run, switch alleles listing
      #Else change coefficients to column names
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
          c_logrMaaslin$error( "Unknown coefficient: %s", strOrig )
        }
        if( substring( strOrig, nchar( strMetadata ) + 1 ) == "NA" ) { next }
        adMetadata <- frmeData[,strMetadata]
      }

      #Bonferonni correct the factor p-values based on the factor levels-1 comparisons
      if( class( adMetadata ) == "factor" )
      {
        dP <- dP * ( nlevels( adMetadata ) - 1 )
      }

      #Store (factor level modified) p-value
      #Store general results for each coef
      adP <- c(adP, dP)
      lsSig[[length( lsSig ) + 1]] <- list(
        #Current metadata name
        name		= strMetadata,
        #Current metadatda name (as a factor level if existing as such)
        orig		= strOrig,
        #Taxon feature name
        taxon		= strTaxon,
        #Taxon data / response
        data		= frmeData[,iTaxon],
        #All levels
        factors		= c(strMetadata),
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
  return(list(adP=adP, lsSig=lsSig, lsQCCounts=lsQCCounts))
  ### List containing a list of pvalues, a list of significant data per association, and a list of QC data
}

### Options for regularization 

#TODO# Add in forced into regularization
#TODO# make sure the qvuale are made of the right number, univariate case has more comparisons
# I now dummy the variable levels and have a pvalue for each. Should be ok.

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
  lmod = try( gbm( as.formula( strFormula ), data=frmeTmp, distribution="laplace", verbose=FALSE, n.minobsinnode=min(1, round(0.2 * nrow( frmeTmp ) ) ), n.trees=1000 ) )

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

    #Select model predictors
    #Check the frequency of selection and skip if not selected more than set threshold dFreq
    for( strMetadata in lmod$var.names )
    {
      #If the selprob is less than a certain frequency, skip
      dSel <- lsSum$rel.inf[which( lsSum$var == strMetadata )] / 100
      if( is.na(dSel) || ( dSel < lsParameters$dFreq ) ){ next }

      #Get the name of the metadata
      strTerm <- funcCoef2Col( strMetadata, frmeTmp, c(astrMetadata, astrGenetics) )

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
    }
  } else { astrTerms = lsForcedParameters }
  return(unique(c(astrTerms,lsForcedParameters)))
  ### Return a vector of predictor names to use in a reduced model
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

  strNULLFormula = ifelse(is.null(lsForcedParameter), "adCur ~ 1", paste( "adCur ~", paste( sprintf( "`%s`", lsForcedParameter ), collapse = " + " )))
  lmodNull <- try( lm(as.formula( strNULLFormula ), data=frmeTmp))
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

  strNULLFormula = ifelse(is.null(lsForcedParameter), "adCur ~ 1", paste( "adCur ~", paste( sprintf( "`%s`", lsForcedParameter ), collapse = " + " )))
  lmodNull <- try( lm(as.formula( strNULLFormula ), data=frmeTmp))
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

# GUnifrac
#TODO# Implemented in sfle

# OK
# Wilcoxon (T-Test)
# Does multiple lmod results integrated into maaslin?
funcWilcoxon <- function(
### Perform multiple univariate comparisons performing wilcoxon tests to measure association
strFormula,
### lm style string defining reposne and predictors 
frmeTmp,
### Data on which to perform analysis
iTaxon,
### Index of the response data
lsQCCounts
### List recording anything important to QC
){
  adCur = frmeTmp[,iTaxon]
  return(funcMakeContrasts(strFormula, frmeTmp, adCur, iTaxon=iTaxon,  functionContrast=function(x,y){wilcox.test(x=x,y=y, na.action=c_strNA_Action)}, lsQCCounts, fDummy=TRUE))
  ### List of contrast information, pvalue, contrast and std per univariate test
}

# OK
# Correlation
# Does multiple lmod results integrated into maaslin?
funcSpearman <- function(
### Perform multiple univariate comparisons producing spearman correlations for association
strFormula,
### lm style string defining reposne and predictors 
frmeTmp,
### Data on which to perform analysis
iTaxon,
### Index of the response data
lsQCCounts
### List recording anything important to QC
){
  adCur = frmeTmp[,iTaxon]
  return(funcMakeContrasts(strFormula=strFormula, frmeTmp=frmeTmp, adCur=adCur, iTaxon=iTaxon,  functionContrast=function(x,y){cor.test(x=x, y=y, method="spearman", na.action=c_strNA_Action)}, lsQCCounts, fDummy=TRUE))
  ### List of contrast information, pvalue, contrast and std per univariate test
}

### Multivariate

#TODO do I need to standardize?
funcLasso <- function(
### Perform lasso for regualrization and associations
strFormula,
### lm style string defining reposne and predictors 
frmeTmp,
### Data on which to perform analysis
iTaxon,
### Index of the response data
lsQCCounts
### List recording anything important to QC
){
  adCur = frmeTmp[,iTaxon]
  return(try(penalized(response=adCur, penalized=as.formula(strFormula), lambda1=1, data=frmeTmp, standardize=TRUE)))
  ### lmod result object from lasso lm
}

# OK
funcLM <- function(
### Perform vanilla linear regression
strFormula,
### lm style string defining reposne and predictors 
frmeTmp,
### Data on which to perform analysis
iTaxon,
### Index of the response data
lsQCCounts
### List recording anything important to QC
){
  adCur = frmeTmp[,iTaxon]
#  return(try( lmer(as.formula(strFormula), data=frmeTmp, na.action=c_strNA_Action) ))
  return(try( lm(as.formula(strFormula), data=frmeTmp, na.action=c_strNA_Action) ))
  ### lmod result object from lm
}

# Multistep maaslin in Curtis' latest code
#funcMultiStepLM <- function(strFormula, frmeTmp, adCur)
#{
#  lmod <- NA
#  if( ( sum( !adCur, na.rm = TRUE ) / sum( !is.na( adCur ) ) ) >= c_dMinSamp )
#  {
#    adCur <- round( 1e1 * adCur / min( abs( adCur[adCur != 0] ), na.rm = TRUE ) )
#    lmod <- try( zeroinfl( as.formula(strFormula), data = frmeTmp, dist = "negbin", link = "logit" ) )
#  }
#  if( is.na( lmod ) || ( class( lmod ) == "try-error" ) )
#  {
#    lmod <- try( ltsReg( as.formula(strFormula), data = frmeTmp, nsamp = "best", adjust = TRUE, alpha = 1, mcd = FALSE ) )
#  }
#  if( is.na( lmod ) || ( class( lmod ) == "try-error" ) )
#  {
#    lmod <- try( lm( as.formula(strFormula), data = frmeTmp ) )
#  }
#  return(lmod)
#}

### Link functions / Transformations

# OK
funcArcsinSqrt <- function(
# Transform data with arcsin sqrt transformation
aData
### The data on which to perform the transformation
){
  return(asin(sqrt(aData)))
  ### Transformed data
}

# OK
funcNoTransform <-function(
### Pass data without transform
aData
### The data on which to perform the transformation
### Only given here to preserve the pattern, not used.
){
  return(aData)
  ### Transformed data
}

#OK
funcBinomialMult <- function(
### Perform linear regression with binomial link
strFormula,
### lm style string defining reposne and predictors 
frmeTmp,
### Data on which to perform analysis
iTaxon,
### Index of the response data
lsQCCounts
### List recording anything important to QC
){
  adCur = frmeTmp[,iTaxon]
  return(try ( glmer(as.formula(strFormula), data=frmTmp, family=binomial(link=logit), na.action=c_strNA_Action) ))
#  return(try( glm(as.formula(strFormula), family=binomial(link=logit), data=frmeTmp, na.action=c_strNA_Action) ))
  ### lmod result object from lm
}

#OK
funcQuasiMult <- function(
### Perform linear regression with quasi-poisson link
strFormula,
### lm style string defining reposne and predictors 
frmeTmp,
### Data on which to perform analysis
iTaxon,
### Index of the response data
lsQCCounts
### List recording anything important to QC
){
  adCur = frmeTmp[,iTaxon]
  return(try ( glmer(as.formula(strFormula), data=frmTmp, family=quasipoisson, na.action=c_strNA_Action) ))
#  return(try( glm(as.formula(strFormula), family=quasipoisson, data=frmeTmp, na.action=c_strNA_Action) ))
  ### lmod result object from lm
}
