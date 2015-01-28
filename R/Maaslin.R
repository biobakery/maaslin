#!/usr/bin/env Rscript
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
##description<< Main driver script. Should be called to perform MaAsLin Analysis.
) { return( pArgs ) }


### Load packages
vDepLibrary = c("agricolae", "gam", "gamlss", "gbm", "glmnet", "inlinedocs", "logging", "MASS", "nlme", "optparse", "outliers", "penalized", "pscl", "robustbase")
for(sDepLibrary in vDepLibrary)
{
  if(! require(sDepLibrary, character.only=TRUE) )
  {
    stop(paste("Please install the required package:",sDepLibrary,sep=" "))
  }
}

### Logging class
suppressMessages(library( logging, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
### Class for commandline argument processing
suppressMessages(library( optparse, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))


### Create command line argument parser
pArgs <- OptionParser( usage = "%prog [options] <data.tsv> <outputdir>" )

# Input files for MaAsLin
## Data configuration file
pArgs <- add_option( pArgs, c("-i", "--input_config"), type="character", action="store", dest="strInputConfig", metavar="data.read.config", help="Optional configuration file describing data input format.")
## Data manipulation/normalization file
pArgs <- add_option( pArgs, c("-I", "--input_process"), type="character", action="store", dest="strInputR", metavar="data.R", help="Optional configuration script normalizing or processing data.")

# Settings for MaAsLin
## Maximum false discovery rate
pArgs <- add_option( pArgs, c("-d", "--fdr"), type="double", action="store", dest="dSignificanceLevel", default=0.25, metavar="significance", help="The threshold to use for significance for the generated q-values (BH FDR). Anything equal to or lower than this is significant.  [Default %default]")
## Minimum feature relative abundance filtering
pArgs <- add_option( pArgs, c("-r", "--minRelativeAbundance"), type="double", action="store", dest="dMinAbd", default=0.0001, metavar="minRelativeAbundance", help="The minimum relative abundance allowed in the data. Values below this are removed and imputed as the median of the sample data.  [Default %default]")
## Minimum feature prevalence filtering
pArgs <- add_option( pArgs, c("-p", "--minPrevalence"), type="double", action="store", dest="dMinSamp", default=0.1, metavar="minPrevalence", help="The minimum percentage of samples a feature can have abundance in before being removed. Also is the minimum percentage of samples a metadata can have that are not NA before being removed.  [Default %default]")
## Fence for outlier, if not set Grubbs test is used
pArgs <- add_option( pArgs, c("-o", "--outlierFence"), type="double", action="store", dest="dOutlierFence", default=0, metavar="outlierFence", help="Outliers are defined as this number times the interquartile range added/subtracted from the 3rd/1st quartiles respectively. If set to 0 (default), outliers are defined by the Grubbs test.  [Default %default]")
## Significance for Grubbs test
pArgs <- add_option(pArgs, c("-G","--grubbsSig"), type="double", action="store", dest="dPOutlier", default=0.05, metavar="grubbsAlpha", help="This is the significance cuttoff used to indicate an outlier or not. The closer to zero, the more significant an outlier must be to be removed.  [Default %default]")
## Fixed (not random) covariates
pArgs <- add_option( pArgs, c("-R","--random"), type="character", action="store", dest="strRandomCovariates", default=NULL, metavar="fixed", help="These metadata will be treated as random covariates. Comma delimited data feature names. These features must be listed in the read.config file. Example '-R RandomMetadata1,RandomMetadata2'.  [Default %default]")
## Change the type of correction fo rmultiple corrections
pArgs <- add_option( pArgs, c("-T","--testingCorrection"), type="character", action="store", dest="strMultTestCorrection", default="BH", metavar="multipleTestingCorrection", help="This indicates which multiple hypothesis testing method will be used, available are holm, hochberg, hommel, bonferroni, BH, BY.  [Default %default]")
## Use a zero inflated model of the inference method indicate in -m
pArgs <- add_option( pArgs, c("-z","--doZeroInfated"), type="logical", action="store_true", default = FALSE, dest="fZeroInflated", metavar="fZeroInflated", help="If true, the zero inflated version of the inference model indicated in -m is used. For instance if using lm, zero-inflated regression on a gaussian distribution is used.  [Default %default].")

# Arguments used in validation of MaAsLin
## Model selection (enumerate) c("none","boost","penalized","forward","backward")
pArgs <- add_option( pArgs, c("-s", "--selection"), type="character", action="store", dest="strModelSelection", default="boost", metavar="model_selection", help="Indicates which of the variable selection techniques to use.  [Default %default]")
## Argument indicating which method should be ran (enumerate) c("univariate","lm","neg_binomial","quasi")
pArgs <- add_option( pArgs, c("-m", "--method"), type="character", action="store", dest="strMethod", default="lm", metavar="analysis_method", help="Indicates which of the statistical inference methods to run.  [Default %default]")
## Argument indicating which link function is used c("none","asinsqrt")
pArgs <- add_option( pArgs, c("-l", "--link"), type="character", action="store", dest="strTransform", default="asinsqrt", metavar="transform_method", help="Indicates which link or transformation to use with a glm, if glm is not selected this argument will be set to none.  [Default %default]")
pArgs <- add_option( pArgs, c("-Q","--NoQC"), type="logical", action="store_true", default=FALSE, dest="fNoQC", metavar="Do_Not_Run_QC", help="Indicates if the quality control will be ran on the metadata/data. Default is true.  [Default %default]")

# Arguments to suppress MaAsLin actions on certain data
## Do not perform model selection on the following data
pArgs <- add_option( pArgs, c("-F","--forced"), type="character", action="store", dest="strForcedPredictors", default=NULL, metavar="forced_predictors", help="Metadata features that will be forced into the model seperated by commas. These features must be listed in the read.config file. Example '-F Metadata2,Metadata6,Metadata10'.  [Default %default]")
## Do not impute the following
pArgs <- add_option( pArgs, c("-n","--noImpute"), type="character", action="store", dest="strNoImpute", default=NULL, metavar="no_impute", help="These data will not be imputed. Comma delimited data feature names. Example '-n Feature1,Feature4,Feature6'.  [Default %default]")

#Miscellaneouse arguments
### Argument to control logging (enumerate)
strDefaultLogging = "DEBUG"
pArgs <- add_option( pArgs, c("-v", "--verbosity"), type="character", action="store", dest="strVerbosity", default=strDefaultLogging, metavar="verbosity", help="Logging verbosity  [Default %default]")
### Run maaslin without creating a log file
pArgs <- add_option( pArgs, c("-O","--omitLogFile"), type="logical", action="store_true", default=FALSE, dest="fOmitLogFile", metavar="omitlogfile",help="Including this flag will stop the creation of the output log file.  [Default %default]")
### Argument for inverting background to black
pArgs <- add_option( pArgs, c("-t", "--invert"), type="logical", action="store_true", dest="fInvert", default=FALSE, metavar="invert", help="When given, flag indicates to invert the background of figures to black.  [Default %default]")
### Selection Frequency
pArgs <- add_option( pArgs, c("-f","--selectionFrequency"), type="double", action="store", dest="dSelectionFrequency", default=NA, metavar="selectionFrequency", help="Selection Frequency for boosting (max 1 will remove almost everything). Interpreted as requiring boosting to select metadata 100% percent of the time (or less if given a number that is less). Value should be between 1 (100%) and 0 (0%), NA (default is determined by data size).")
### All v All
pArgs <- add_option( pArgs, c("-a","--allvall"), type="logical", action="store_true", dest="fAllvAll", default=FALSE, metavar="compare_all", help="When given, the flag indicates that each fixed covariate that is not indicated as Forced is compared once at a time per data feature (bug). Made to be used with the -F option to specify one part of the model while allowing the other to cycle through a group of covariates. Does not affect Random covariates, which are always included when specified.  [Default %default]")
pArgs <- add_option( pArgs, c("-N","--PlotNA"), type="logical", action="store_true", default=FALSE, dest="fPlotNA", metavar="plotNAs",help="Plot data that was originally NA, by default they are not plotted.  [Default %default]")
### Alternative methodology settings
pArgs <- add_option( pArgs, c("-A","--pAlpha"), type="double", action="store", dest="dPenalizedAlpha", default=0.95, metavar="PenalizedAlpha",help="The alpha for penalization (1.0=L1 regularization, LASSO; 0.0=L2 regularization, ridge regression.  [Default %default]")
### Pass an alternative library dir
pArgs <- add_option( pArgs, c("-L", "--libdir"), action="store", dest="sAlternativeLibraryLocation", default=file.path( "","usr","share","biobakery" ), metavar="AlternativeLibraryDirectory", help="An alternative location to find the lib directory. This dir and children will be searched for the first maaslin/src/lib dir.")

#pArgs <- add_option( pArgs, c("-c","--MFAFeatureCount"), type="integer", action="store", dest="iMFAMaxFeatures", default=3, metavar="maxMFAFeature", help="Number of features or number of bugs to plot (default=3; 3 metadata and 3 data).")

 
#source('~/Maaslin/R/Maaslinfn.R')
source('./R/Maaslinfn.R')


if( identical( environment( ), globalenv( ) ) &&
	!length( grep( "^source\\(", sys.calls( ) ) ) ) {
	lsArgs <- parse_args( pArgs, positional_arguments = TRUE )
	Maaslinfn( lsArgs ) }
