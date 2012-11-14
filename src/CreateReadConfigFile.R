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
