
# Libraries
suppressMessages(library( penalized, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Need for stepAIC
suppressMessages(library( MASS, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Needed for na action behavior
suppressMessages(library( gam, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Needed for boosting
suppressMessages(library( gbm, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))


### Helper functions
funcMakeContrasts <- function(strFormula, frmeTmp, adCur, functionContrast)
{
  #Get test comparisons
  asComparisons = gsub("`","",setdiff(unlist(strsplit(strFormula," ")),c("adCur","~","+")))
  lliTests = NULL
  almod = c()
  for(sComparison in asComparisons)
  {
    lliTests = list("1"=frmeTmp[[sComparison]])
    if(is.factor(lliTests[[1]]))
    {
      lliTests = apply(as.matrix(model.matrix(as.formula(paste("~",sComparison)),frmeTmp))[,-1],2,list)
    }
    sComparisonResults = lapply(lliTests, function(x) functionContrast(x=as.vector(unlist(x)), y=adCur))
    print(sComparisonResults)
    #TODO finish
    almod = c(list("coefficients"=NA,"residuals"=NA,
    "fitted.values"=NA,"rank"=NA,"weights"=NA,
    "df.residual"=NA,"call"=NA,"terms"=NA,
    "contrasts"=NA,"xlevels"=NA,"offset"=NA,
    "y"=adCur,"x"=NA,"model"=NA,"na.action"=c_strNA_Action))
  }
  return(almod)
}

### Options for regularization 

#TODO# Add in forced into regularization

# Ok
# Gradient Boosting
funcBoostModel <- function(strFormula, frmeTmp, adCur, lsParameters, strLog)
{
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

    #For each metadata coefficient
    #Check the frequency of selection and skip if not selected more than set threshold dFreq
    for( strMetadata in lmod$var.names )
    {
      #If the selprob is less than a certain frequency, skip
      dSel <- lsSum$rel.inf[which( lsSum$var == strMetadata )] / 100
      if( is.na(dSel) || ( dSel < lsParameters$dFreq ) ) { next }
      #Get the name of the metadata
      strTerm <- funcCoef2Col( strMetadata, frmeData, c(astrMetadata, astrGenetics) )

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
}

# Select model with forward selection
funcForwardModel <- function(strFormula, frmeTmp, adCur, lsParameters, strLog)
{
  astrTerms <- c()
  funcWrite( c("#Forward formula", strFormula), strLog )
  print("#Forwards")
  lmodNull <- try( lm(as.formula( "adCur ~ 1"), data=frmeTmp))
  lmodFull <- try( lm(as.formula( strFormula ), data=frmeTmp ))
  if(!("try-error" %in% c(class( lmodNull ),class( lmodFull ))))
  {
    lmod = stepAIC(lmodNull, scope=list(lower=lmodNull,upper=lmodFull), direction="forward", trace=0)
    print("*1**")
    print(class(lmod))
    print("*2**")
    print(coef(lmod))
    print("*3**")
    print(lmod$anov)
    print("*4**")
  } else { return( NA ) }
  return(astrTerms)
}

# Select model with backwards selection
funcBackwardsModel <- function(strFormula, frmeTmp, adCur, lsParameters, strLog)
{
  astrTerms <- c()
  funcWrite( c("#Backwards formula", strFormula), strLog )
  print("#Backwards")
  lmodFull <- try( lm(as.formula( strFormula ), data=frmeTmp ))
  if(! class( lmodFull ) == "try-error" )
  {
    lmod = stepAIC(lmodFull, direction="backward")
    print("*1**")
    print(lmod)
    print("*2**")
    print(lmod$anova)
    print("*3**")
  } else { return( NA ) }
  return(astrTerms)
}

# Lasso
# Performed in analysis

### Analysis methods
### Univariate options

# Lefse
#TODO# Implemented in sfle

# GUnifrac
#TODO# Implemented in sfle

# Wilcoxon (T-Test)
# Does multiple lmod results integrated into maaslin?
funcWilcoxon <- function( strFormula, frmeTmp, adCur )
{
  return(funcMakeContrasts(strFormula, frmeTmp, adCur, function(x,y){wilcox.test(x=x,y=y, na.action=c_strNA_Action)}))
}

# Correlation
# Does multiple lmod results integrated into maaslin?
funcSpearman <- function( strFormula, frmeTmp, adCur )
{
  return(funcMakeContrasts(strFormula, frmeTmp, adCur, function(x,y){cor.test(x=x, y=y, method="spearman", na.action=c_strNA_Action)}))
}

### Multivariate

# Lasso
#TODO do I need to standardize?
funcLasso <- function(strFormula, frmeTmp, adCur)
{
  return(try(penalized(response=adCur, penalized=as.formula(strFormula), lambda1=1, data=frmeTmp, standardize=TRUE)))
}

# Standard LM
funcLM <- function(strFormula, frmeTmp, adCur)
{
  return( try( lm(as.formula(strFormula), data=frmeTmp, na.action=c_strNA_Action) ))
}

# Multistep maaslin in Curtis' latest code
funcMultiStepLM <- function(strFormula, frmeTmp, adCur)
{
  lmod <- NA
  if( ( sum( !adCur, na.rm = TRUE ) / sum( !is.na( adCur ) ) ) >= c_dMinSamp )
  {
    adCur <- round( 1e1 * adCur / min( abs( adCur[adCur != 0] ), na.rm = TRUE ) )
    lmod <- try( zeroinfl( as.formula(strFormula), data = frmeTmp, dist = "negbin", link = "logit" ) )
  }
  if( is.na( lmod ) || ( class( lmod ) == "try-error" ) )
  {
    lmod <- try( ltsReg( as.formula(strFormula), data = frmeTmp, nsamp = "best", adjust = TRUE, alpha = 1, mcd = FALSE ) )
  }
  if( is.na( lmod ) || ( class( lmod ) == "try-error" ) )
  {
    lmod <- try( lm( as.formula(strFormula), data = frmeTmp ) )
  }
  return(lmod)
}

### Link functions / Transformations

# Arcsin square root
funcArcsinSqrt <- function(frmeData, aiData)
{
  for(aiDatum in aiData)
  {
    frmeData[,aiDatum] = asin(sqrt(frmeData[,aiDatum]))
  }
  return( frmeData=frmeData )
}

# Transform data with arcsin sqrt transformation
funcArcsinSqrt <- function(aData)
{
  return(asin(sqrt(aData)))
}

# Pass data without transform
funcNoTransform <-function(aData)
{
  return(aData)
}

# Binomial
funcBinomialMult <- function(strFormula, frmeTmp, adCur)
{
  return(try( glm(as.formula(strFormula), family=binomial(link=logit), data=frmeTmp, na.action=c_strNA_Action) ))
}

# Quasi-poisson
funcQuasiMult <- function(strFormula, frmeTmp, adCur)
{
  return (try( glm(as.formula(strFormula), family=quasipoisson, data=frmeTmp, na.action=c_strNA_Action) ))
}

# To add a new method insert an entry in the switch for either the selection, transform, or method
# Insert them by using the pattern optparse_keyword_without_quotes = function_in_AnalysisModules
# Order in the return listy is curretly set and expected to be selection, transforms/links, analsis method
# none returns null
funcGetAnalysisMethods <- function(sModelSelectionKey,sTransformKey,sMethodKey)
{
  lRetMethods = list()
  #Insert selection methods here
  lRetMethods[[c_iSelection]] = switch(sModelSelectionKey,
    boost = funcBoostModel,
    forward = funcForwardModel,
    backward = funcBackwardsModel,
    none = NULL)

  #Insert transforms
  lRetMethods[[c_iTransform]] = switch(sTransformKey,
    asinsqrt = funcArcsinSqrt,
    none = funcNoTransform)

  #Insert analysis
  lRetMethods[[c_iAnalysis]] = switch(sMethodKey,
    neg_binomial = funcBinomialMult,
    quasi = funcQuasiMult,
    spearman = funcSpearman,
    wilcoxon = funcWilcoxon,
    lasso = funcLasso,
    lm = funcLM,
    none = NULL)

  return(lRetMethods)
}

