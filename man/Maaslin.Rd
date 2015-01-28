\name{Maaslin}
\alias{Maaslin}
\title{
MaAsLin: Multivariate statistical framework main function 
}
\description{
*****************************************************************************\cr
*   Maaslin is invoked via the OS command as follows:\cr
*   "./R/Maaslin.R -i ./data/maaslin_demo2.read.config  maaslin_demo2.tsv  outputdir"\cr 
*****************************************************************************\cr
The main function manages the following:\cr
1. Optparse arguments are checked\cr
2. A logger is created if requested in the optional arguments\cr
3. The custom R script is sourced. This is the input *.R script named\cr
the same as the input *.pcl file. This script contains custom formating\cr
of data and function calls to the MFA visualization.\cr
4. Matrices are written to the project folder as they are read in seperately as metadata and data and merged together.\cr
5. Data is cleaned with custom filtering if supplied in the *.R script.\cr
6. Transformations occur if indicated by the optional arguments\cr
7. Standard quality control is performed on data\cr
8. Cleaned metadata and data are written to output project for documentation.\cr
9. A regularization method is ran (boosting by default).\cr
10. An analysis method is performed on the model (optionally boosted model).\cr
11. Data is summarized and PDFs are created for significant associations\cr
(those whose q-values {BH FDR correction} are <= the threshold given in the optional arguments.\cr
}
\usage{
Maaslin(pArgs)
}
\arguments{
  \item{pArgs}{
  
************************************************************************\cr 
* Other Parameters that can be used when invoking Maaslin \cr
* (The defaults are the following):\cr
* ./R/Maaslin.R -i ./data/maaslin_demo2.read.config  maaslin_demo2.tsv  outputdir\cr 
*************************************************************************\cr
  
  pArgs contains the parameters entered in optparse and are as follows: \cr
  
     -i DATA.READ.CONFIG \cr
                Optional configuration file describing data input format.\cr
\cr
        -I DATA.R \cr
                Optional configuration script normalizing or processing data. \cr \cr
        -d SIGNIFICANCE \cr\cr
                The threshold to use for significance for the generated q-values (BH FDR). Anything equal to or lower than this is significant.  [Default 0.25]\cr \cr
        -r MINRELATIVEABUNDANCE\cr 
                The minimum relative abundance allowed in the data. Values below this are removed and imputed as the median of the sample data.  [Default 1e-04]\cr\cr
\cr
        -p MINPREVALENCE\cr
                The minimum percentage of samples a feature can have abundance in before being removed. Also is the minimum percentage of samples a metadata can have that are not NA before being removed.  [Default 0.1]\cr\cr
\cr
        -o OUTLIERFENCE\cr
                Outliers are defined as this number times the interquartile range added/subtracted from the 3rd/1st quartiles respectively. If set to 0 (default), outliers are defined by the Grubbs test.  [Default 0]\cr\cr
\cr
        -G GRUBBSALPHA\cr
                This is the significance cuttoff used to indicate an outlier or not. The closer to zero, the more significant an outlier must be to be removed.  [Default 0.05]\cr\cr
\cr
        -R FIXED\cr
                These metadata will be treated as random covariates. Comma delimited data feature names. These features must be listed in the read.config file. Example '-R RandomMetadata1,RandomMetadata2'.  [Default NULL]\cr\cr
\cr
        -T MULTIPLETESTINGCORRECTION\cr
                This indicates which multiple hypothesis testing method will be used, available are holm, hochberg, hommel, bonferroni, BH, BY.  [Default BH]\cr\cr
\cr
        -z,\cr
                If true, the zero inflated version of the inference model indicated in -m is used. For instance if using lm, zero-inflated regression on a gaussian distribution is used.  [Default FALSE].\cr\cr
\cr
        -s MODEL_SELECTION\cr
                Indicates which of the variable selection techniques to use.  [Default boost]\cr\cr
\cr
        -m ANALYSIS_METHOD\cr
                Indicates which of the statistical inference methods to run.  [Default lm]\cr\cr
\cr
        -l TRANSFORM_METHOD\cr
                Indicates which link or transformation to use with a glm, if glm is not selected this argument will be set to none.  [Default asinsqrt]\cr\cr
\cr
        -Q, \cr
                Indicates if the quality control will be ran on the metadata/data. Default is true.  [Default FALSE]\cr\cr
\cr
        -F FORCED_PREDICTORS\cr
                Metadata features that will be forced into the model seperated by commas. These features must be listed in the read.config file. Example '-F Metadata2,Metadata6,Metadata10'.  [Default NULL]\cr\cr
\cr
        -n NO_IMPUTE\cr
                These data will not be imputed. Comma delimited data feature names. Example '-n Feature1,Feature4,Feature6'.  [Default NULL]\cr\cr
\cr
        -v VERBOSITY\cr
                Logging verbosity  [Default DEBUG]\cr\cr
\cr
        -O,  \cr
                Including this flag will stop the creation of the output log file.  [Default FALSE]\cr\cr
\cr
        -t, \cr
                When given, flag indicates to invert the background of figures to black.  [Default FALSE]\cr\cr
\cr
        -f SELECTIONFREQUENCY\cr\cr
                Selection Frequency for boosting (max 1 will remove almost everything). Interpreted as requiring boosting to select metadata 100% percent of the time (or less if given a number that is less). Value should be between 1 (100%) and 0 (0%), NA (default is determined by data size).\cr\cr
\cr
        -a, \cr\cr
                When given, the flag indicates that each fixed covariate that is not indicated as Forced is compared once at a time per data feature (bug). Made to be used with the -F option to specify one part of the model while allowing the other to cycle through a group of covariates. Does not affect Random covariates, which are always included when specified.  [Default FALSE]\cr\cr
\cr
        -N,\cr\cr
                Plot data that was originally NA, by default they are not plotted.  [Default FALSE]\cr\cr
\cr
        -A PENALIZEDALPHA\cr\cr
                The alpha for penalization (1.0=L1 regularization, LASSO; 0.0=L2 regularization, ridge regression.  [Default 0.95]\cr\cr

        -L ALTERNATIVELIBRARYDIRECTORY\cr\cr
                An alternative location to find the lib directory. This dir and children will be searched for the first maaslin/src/lib dir.\cr\cr
\cr

}
}
\author{
Timothy Tickle<ttickle@hsph.harvard.edu>, Curtis Huttenhower <chuttenh@hsph.harvard.edu>\cr
}
\examples{
setwd('~/Maaslin')
system("./R/Maaslin.R -i ./data/maaslin_demo2.read.config maaslin_demo2.tsv outputdir", intern = TRUE)
}


