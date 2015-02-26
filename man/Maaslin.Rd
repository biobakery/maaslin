\name{Maaslin}
\alias{Maaslin}
\title{
MaAsLin is a multivariate statistical framework that finds associations between clinical metadata and potentially high-dimensional experimental data.
}
\description{
MaAsLin performs boosted additive general linear models between one group of data (metadata/the predictors) and another group (in our case relative taxonomic abundances/the response). 
Used to discover associations between clinical metadata and microbial community relative abundance or function }
\usage{
Maaslin(strInputTSV,
	 strOutputDIR, 
	 strInputConfig = NULL, 
	 strInputR = NULL, 
	 dSignificanceLevel = 0.25,
	 dMinAbd = 0.0001,
	 dMinSamp = 0.1,
	 dOutlierFence = 0, 
	 dPOutlier = 0.05, 
	 strRandomCovariates = NULL, 
	 strMultTestCorrection = "BH",
	 fZeroInflated = FALSE, 
	 strModelSelection = "boost",
	 strMethod = "lm",
	 strTransform = "asinsqrt",
	 fNoQC = FALSE,
	 strForcedPredictors = NULL, 
	 strNoImpute = NULL,
	 strVerbosity = "DEBUG",
	 fOmitLogFile = FALSE,
	 fInvert = FALSE,
	 dSelectionFrequency = NA, 
	 fAllvAll = FALSE,
	 fPlotNA = FALSE, 
	 dPenalizedAlpha = 0.95, 
	 sAlternativeLibraryLocation = NULL)
}
\arguments{
  \item{strInputTSV}{
  The main INPUT file: The sample file is maaslin_demo2.tsv
}
  \item{strOutputDIR}{
  Output Directory   
}
  \item{strInputConfig}{
  Input Config file: The sample is located in data/maaslin_demo2.read.config
}
  \item{strInputR}{
 Optional configuration script normalizing or processing data
}
  \item{dSignificanceLevel}{
  Threshold to use for significance for the generated q-values (BH FDR). Anything equal to or lower than this is significant.
}
  \item{dMinAbd}{
  Minimum relative abundance allowed in the data. Values below this are removed and imputed as the median of the sample data.  
}
  \item{dMinSamp}{
   Minimum percentage of samples in which a feature must have the minimum relative abundance in order not to be removed. Also this is the maximum percentage of samples for which a metadata can have NAs before being removed. 
}
  \item{dOutlierFence}{
  Outliers are defined as this number times the interquartile range added/subtracted from the 3rd/1st quartiles respectively. If set to 0 (default), outliers are defined by the Grubbs test.
}
  \item{dPOutlier}{
   This is the significance cuttoff used to indicate an outlier or not. The closer to zero, the more significant an outlier must be to be removed.  
}
  \item{strRandomCovariates}{
    These metadata will be treated as random covariates. Comma delimited data feature names. These features must be listed in the read.config file. Example '-R RandomMetadata1,RandomMetadata2'.
  }
  \item{strMultTestCorrection}{
  This indicates which multiple hypothesis testing method will be used, available are holm, hochberg, hommel, bonferroni, BH, BY.}
  \item{fZeroInflated}{
   If true, the zero inflated version of the inference model indicated in -m is used. For instance if using lm, zero-inflated regression on a gaussian distribution is used.
}
  \item{strModelSelection}{
 Indicates which of the variable selection techniques to use. Default=boost
}
  \item{strMethod}{
  Indicates which of the statistical inference methods to run. Default=lm
}
  \item{strTransform}{
     Indicates which link or transformation to use with a glm, if glm is not selected this argument will be set to none.  Default=asinsqrt
}
  \item{fNoQC}{
  Indicates if the quality control will be ran on the metadata/data. Default is FALSE
}
  \item{strForcedPredictors}{
  Metadata features that will be forced into the model seperated by commas. These features must be listed in the read.config file. Example '-F Metadata2,Metadata6,Metadata10'. 
}
  \item{strNoImpute}{
   These data will not be imputed. Comma delimited data feature names. Example '-n Feature1,Feature4,Feature6'. 
}
  \item{strVerbosity}{
  Debug level
}
  \item{fOmitLogFile}{
  Including this flag will stop the creation of the output log file.  Default=FALSE
}
  \item{fInvert}{
  When given, flag indicates to invert the background of figures to black. Defaule = FALSE
}
  \item{dSelectionFrequency}{
  Selection Frequency for boosting (max 1 will remove almost everything). Interpreted as requiring boosting to select metadata 100% percent of the time (or less if given a number that is less). Value should be between 1 (100%) and 0 (0%), NA (default is determined by data size)
}
  \item{fAllvAll}{
  When given, the flag indicates that each fixed covariate that is not indicated as Forced is compared once at a time per data feature (bug). Made to be used with the -F option to specify one part of the model while allowing the other to cycle through a group of covariates. Does not affect Random covariates, which are always included when specified.
}
  \item{fPlotNA}{
   Plot data that was originally NA, by default they are not plotted.  Default=FALSE
}
  \item{dPenalizedAlpha}{
   The alpha for penalization (1.0=L1 regularization, LASSO; 0.0=L2 regularization, ridge regression.
}
  \item{sAlternativeLibraryLocation}{
    An alternative location to find the lib directory. This dir and children will be searched for the first maaslin/src/lib dir.
 }
}
\author{
Timothy Tickle<ttickle@hsph.harvard.edu>,\cr
Curtis Huttenhower <chuttenh@hsph.harvard.edu>\cr
Maintainers: Ayshwarya Subramanian<subraman@broadinstitute.org>,\cr 
Lauren McIver<lauren.j.mciver@gmail.com>,\cr
George Weingart<george.weingart@gmail.com>
}
\examples{
getwd()
list.files()
Maaslin('Maaslin/extdata/maaslin_demo2.tsv','output',strInputConfig='Maaslin/extdata/maaslin_demo2.read.config')
}

