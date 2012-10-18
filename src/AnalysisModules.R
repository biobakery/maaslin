
# Libraries
#suppressMessages(library( glmnet, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Need for stepAIC
suppressMessages(library( MASS, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Needed for na action behavior
suppressMessages(library( gam, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
# Needed for boosting
suppressMessages(library( gbm, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))

### Options for regularization 

#TODO# Add in forced in Boost GLM

# Input data and formula, output model?
# Gradient Boosting
funcBoostModel <- function(strFormula, frmeTmp, adCur, lsForced=NA)
{
  return( try( gbm( as.formula( strFormula ), data=frmeTmp, distribution="laplace", verbose=FALSE, n.minobsinnode=min(1, round(0.2 * nrow( frmeTmp ) ) ), n.trees=1000 ) ))
}

# Select model with forward selection
#funcForwardModel <- function(lsFormula, frmeTmp, adCur, lsforced)
#{
#  lmod <- try( lm(as.formula( strFormula ), data=frmeTmp ))
#  if(! class( lmod ) == "try-error" )
#  { lmod = stepAIC(lmod, direction="forward") }
#  return(lmod)
#}

# Select model with backwards selection
#funcBackwardsModel <- function(lsFormula, frmeTmp, adCur, lsForced)
#{
#  lmod <- try( lm(as.formula( strFormula ), data=frmeTmp ))
#  if(! class( lmod ) == "try-error" )
#  { lmod = stepAIC(lmod, direction="backward") }
#  return(lmod)
#}

# Lasso
# Performed in analysis

### Analysis methods
### Univariate options

# Lefse
#TODO# Implemented in sfle

# GUnifrac
#TODO# Implemented in sfle

# Wilcoxon (T-Test)
#funcWilcox <- function( strFormula, frmeTmp )
#{
#  wilcox.text(x=,y=, na.action="na.omit")
#}

# Correlation
#funcSpearman <- function( strFormula, frmeTmp )
#{
#  cor.test(x=, y=, method="spearman", na.action="na.omit")
#}

### Multivariate

# Lasso
# x = observation matrix
# y = response matrix
# Alpha = 1 is lasso, = 0 is ridge
#TODO do I need to standardize?
#funcLasso <- function(strFormula, frmeTmp)
#{
#  glmnet(x=, y=, family="gaussian", alpha=1, standardize=TRUE)
#}

# Standard LM
funcLM <- function(strFormula, frmeTmp, adCur)
{
  return( try( lm(as.formula(strFormula), data=frmeTmp, na.action="na.omit") ))
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
#funcArcsinSqrt <- function(frmeData, aiData)
#{
#  for(aiDatum in aiData)
#  {
#    frmeData[,aiDatum] = asin(sqrt(frmeData[,aiDatum]))
#  }
#  return( frmeData=frmeData )
#}

# Binomial
#funcBinomialMult <- function(strFormula, frmeTmp)
#{
#  return(try( glm(as.formula(strFormula), family=binomial(link=logit), data=frmeTmp, na.action="na.omit") ))
#}

# Quasi-poisson
#funcQuasiMult <- function(strFormula, frmeTmp)
#{
#  return (try( glm(as.formula(strFormula), family=quasipoisson, data=frmeTmp, na.action="na.omit") ))
#}

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
    none = NULL)

  #Insert transforms/links
  lRetMethods[[c_iTransform]] = switch(sTransformKey,
    none = NULL)

  #Insert analysis
  lRetMethods[[c_iAnalysis]] = switch(sMethodKey,
    lm = funcMultiStepLM,
#    lm = funcLM,
    none = NULL)

  return(lRetMethods)
}

