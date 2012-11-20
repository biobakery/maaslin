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


### Helper functions

funcMakeContrasts <- function
### Makes univariate contrasts of all predictors in the model formula with the response.
(strFormula, 
### lm style string defining reposne and predictors 
frmeTmp,
### The data frame to find predictor data in
adCur,
### adCur Response data
functionContrast
### functionContrast The univariate test to perfom
){
  #Get test comparisons (predictor names from formula string)
  asComparisons = gsub("`","",setdiff(unlist(strsplit(unlist(strsplit(strFormula,"~"))[2]," ")),c("","+")))
  lliTests = NULL
  lData = c()
  for(sComparison in asComparisons)
  {
    lliTests = list("1"=frmeTmp[[sComparison]])
    if(is.factor(lliTests[[1]]))
    {
      lliTests = apply(as.matrix(model.matrix(as.formula(paste("~",sComparison)),frmeTmp))[,-1],2,list)
    }

    sComparisonResults = lapply(lliTests, function(x) functionContrast(x=as.vector(unlist(x)), y=adCur))
    #Pass important data.
    #TODO finish
    lData = c(lData, list(dPValue=sComparisonResults[[1]]$p.value,dSTD=NA,sMetadata=sComparison))
  }
  return(lData)
  ### Returns a list of p-value, standard deviation, and comparison which produced the p-value
}

funcGetStepPredictors <- function
### Retrieve the predictors of the reduced model after stepwise model selection
(lmod,
### Linear model resulting from step-wise model selection
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

  # Select model predictors
  asStepCoefsFactors = lmod$anov$Step
  if(length(asStepCoefsFactors)==1){ return(NA) }
  return( setdiff(unlist(strsplit(paste(as.character(asStepCoefsFactors), collapse=" "),split=" ")),c("","+")) )
  ### Vector of string predictor names
}

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

#TODO finish

funcGetUnivariateResults <- function( 
### Reduce the lm object return to just the data needed for further analysis
mod=lmod,
### The result from a linear model
frmeData=frmeData,
### Data analysis is perfromed on
iTaxon=iTaxon,
### The response id
dSig=dSig,

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
  return(list(adP=adP, lsSig=lsSig, lsQCCounts=lsQCCounts))
  ### List containing a list of pvalues, a list of significant data per association, and a list of QC data
}


### Options for regularization 

#TODO# Add in forced into regularization

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
strLog,
### File to which to document logging
lsForcedParameters = NULL
### Force these predictors to be in the model
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
      if( is.na(dSel) || ( dSel < lsParameters$dFreq ) ) { next }
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
  } else { astrTerms = NA }
  return(astrTerms)
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
strLog,
### File to which to document logging
lsForcedParameters = NULL
### Force these predictors to be in the model
){
  astrTerms <- c()
  funcWrite( c("#Forward formula", strFormula), strLog )

  lmodNull <- try( lm(as.formula( "adCur ~ 1"), data=frmeTmp))
  lmodFull <- try( lm(as.formula( strFormula ), data=frmeTmp ))
  if(!("try-error" %in% c(class( lmodNull ),class( lmodFull ))))
  {
    lmod = stepAIC(lmodNull, scope=list(lower=lmodNull,upper=lmodFull), direction="forward", trace=0)
    return(funcGetStepPredictors(lmod,strLog))
  }
  return( NA )
  ### Return a vector of predictor names to use in a reduced model or NA on error
}

# Not done
# Select model with backwards selection
funcBackwardsModel <- function(
### Perform model selection with backwards stepwise selection
strFormula,
### lm style string defining reposne and predictors 
frmeTmp,
### Data on which to perform analysis
adCur,
### Response data
lsParameters,
### User controlled parameters needed specific to boosting
strLog,
### File to which to document logging
lsForcedParameters = NULL
### Force these predictors to be in the model
){
  astrTerms <- c()
  funcWrite( c("#Backwards formula", strFormula), strLog )

  lmodFull <- try( lm(as.formula( strFormula ), data=frmeTmp ))
  if(! class( lmodFull ) == "try-error" )
  {
    lmod = stepAIC(lmodFull, direction="backward")
    return(funcGetStepPredictors(lmod,strLog))
  } else { return( NA ) }
  return(astrTerms)
  ### Return a vector of predictor names to use in a reduced model or NA on error
}

### Analysis methods
### Univariate options

# Lefse
#TODO# Implemented in sfle

# GUnifrac
#TODO# Implemented in sfle

# Wilcoxon (T-Test)
# Does multiple lmod results integrated into maaslin?
funcWilcoxon <- function(
### Perform multiple univariate comparisons performing wilcoxon tess to measure association
strFormula,
### lm style string defining reposne and predictors 
frmeTmp,
### Data on which to perform analysis
adCur
### Response data
){
  return(funcMakeContrasts(strFormula, frmeTmp, adCur, function(x,y){wilcox.test(x=x,y=y, na.action=c_strNA_Action)}))
  ### List of contrast information, pvalue, contrast and std per univariate test
}

# Correlation
# Does multiple lmod results integrated into maaslin?
funcSpearman <- function(
### Perform multiple univariate comparisons producing spearman correlations for association
strFormula,
### lm style string defining reposne and predictors 
frmeTmp,
### Data on which to perform analysis
adCur
### Response data
){
  return(funcMakeContrasts(strFormula, frmeTmp, adCur, function(x,y){cor.test(x=x, y=y, method="spearman", na.action=c_strNA_Action)}))
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
adCur
### Response data
){
  return(try(penalized(response=adCur, penalized=as.formula(strFormula), lambda1=1, data=frmeTmp, standardize=TRUE)))
  ### lmod result object from lasso lm
}

funcLM <- function(
### Perform vanilla linear regression
strFormula,
### lm style string defining reposne and predictors 
frmeTmp,
### Data on which to perform analysis
adCur
### Response data
){
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

## Arcsin square root
#funcArcsinSqrt <- function(frmeData, aiData)
#{
#  for(aiDatum in aiData)
#  {
#    frmeData[,aiDatum] = asin(sqrt(frmeData[,aiDatum]))
#  }
#  return( frmeData=frmeData )
#}

funcArcsinSqrt <- function(
# Transform data with arcsin sqrt transformation
aData
### The data on which to perform the transformation
){
  return(asin(sqrt(aData)))
  ### Transformed data
}

funcNoTransform <-function(
### Pass data without transform
aData
### The data on which to perform the transformation
### Only given here to preserve the pattern, not used.
){
  return(aData)
  ### Transformed data
}

funcBinomialMult <- function(
### Perform linear regression with binomial link
strFormula,
### lm style string defining reposne and predictors 
frmeTmp,
### Data on which to perform analysis
adCur
### Response data
){
  return(try( glm(as.formula(strFormula), family=binomial(link=logit), data=frmeTmp, na.action=c_strNA_Action) ))
  ### lmod result object from lm
}

funcQuasiMult <- function(
### Perform linear regression with quasi-poisson link
strFormula,
### lm style string defining reposne and predictors 
frmeTmp,
### Data on which to perform analysis
adCur
### Response data
){
  return(try( glm(as.formula(strFormula), family=quasipoisson, data=frmeTmp, na.action=c_strNA_Action) ))
  ### lmod result object from lm
}
