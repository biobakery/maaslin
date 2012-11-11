#!/usr/bin/env Rscript
#######################################################################################
# This file is provided under the Creative Commons Attribution 3.0 license.
#
# You are free to share, copy, distribute, transmit, or adapt this work
# PROVIDED THAT you attribute the work to the authors listed below.
# For more information, please see the following web page:
# http://creativecommons.org/licenses/by/3.0/
#
# This file is a component of the MaAsLin (Multivariate Associations Using Linear Models), 
# authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Timothy Tickle, ttickle@hsph.harvard.edu).
#######################################################################################

inlinedocs <- function(
##author<< Curtis Huttenhower <chuttenh@hsph.harvard.edu> and Timothy Tickle <ttickle@hsph.harvard.edu>
##description<< Allows read config files to be created.
) { return( pArgs ) }

### Logging class
suppressMessages(library( logging, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))
### Class for commandline argument processing
suppressMessages(library( optparse, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))

### Source the IO.R for the script
source(file.path("src","lib","IO.R"))
source(file.path("src","lib","Constants.R"))

### Create command line argument parser
### The TSV (tab seperated value (column major, samples are rows) file that will be read in
### The column that is the last metadata name
### The read.config file that will be used to read in the TSV file
pArgs <- OptionParser( usage = "%prog [optional] <strOutputRC> <strMatrixName>" )
# Settings for Read config
## row indices
pArgs <- add_option( pArgs, c("-r", "--rows"), type="character", action="store", dest="strRows", default=NA, metavar="row_indices", help="Rows to read by index starting with 1.")
## column indices
pArgs <- add_option( pArgs, c("-c", "--columns"), type="character", action="store", dest="strColumns", default=NA, metavar="column_indices", help="Columns to read in by index starting with 1.")
## delimiter
pArgs <- add_option( pArgs, c("-d", "--delimiter"), type="character", action="store", dest="charDelimiter", default="\t", metavar="delimiter", help="Delimiter to read the matrix.")
## append to current file
pArgs <- add_option( pArgs, c("-a", "--append"), type="logical", action="store_true", dest="fAppend", default=FALSE, metavar="append", help="Append to existing data. Default no append.")
### Parse arguments
lsArgs <- parse_args( pArgs, positional_arguments = TRUE )

#Get positional arguments
if( !(length( lsArgs$args ) == 2) ) { stop( print_help( pArgs ) ) }

### Write to file the read config script
funcWriteMatrixToReadConfigFile(strConfigureFileName=lsArgs$args[1], strMatrixName=lsArgs$args[2], strRowIndices=lsArgs$options$strRows,
  strColIndices=lsArgs$options$strColumns,acharDelimiter=lsArgs$options$charDelimiter,fAppend=lsArgs$options$fAppend)
