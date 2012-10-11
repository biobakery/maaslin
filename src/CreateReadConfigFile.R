#!/usr/bin/env Rscript

inlinedocs <- function(
##author<< Timothy Tickle <ttickle@hsph.harvard.edu>
) { return( pArgs ) }

### Logging class
suppressMessages(library( logging, warn.conflicts=False, quietly=TRUE, verbose=FALSE))
### Class for commandline argument processing
suppressMessages(library( optparse, warn.conflicts=False, quietly=TRUE, verbose=FALSE))

### Source the IO.R for the script
source("IO.R")

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
## Character rows
pArgs <- add_option( pArgs, c("-s", "--stringdata"), type="character", action="store", dest="strCharData", default=NA, metavar="char_data", help="Forces these columns to be character data.")
## Factor data rows
pArgs <- add_option( pArgs, c("-f", "--factordata"), type="character", action="store", dest="strFactorData", default=NA, metavar="factor_data", help="Forces these columns to be factor data.")
## Integer data rows
pArgs <- add_option( pArgs, c("-i", "--integerdata"), type="character", action="store", dest="strIntegerData", default=NA, metavar="integer_data", help="Forces these columns to be integer data.")
## logical data rows
pArgs <- add_option( pArgs, c("-l", "--logicaldata"), type="character", action="store", dest="strLogicalData", default=NA, metavar="logical_data", help="Forces these columns to be logical data.")
## numeric data rows
pArgs <- add_option( pArgs, c("-mn", "--numericdata"), type="character", action="store", dest="strNumericData", default=NA, metavar="numeric_data", help="Forces these columns to be numeric data.")
## orderd data rows
pArgs <- add_option( pArgs, c("-o", "--orderdata"), type="character", action="store", dest="strOrderData", default=NA, metavar="order_data", help="Forces these columns to be order data.")
## delimiter
pArgs <- add_option( pArgs, c("-d", "--delimiter"), type="character", action="store", dest="charDelimiter", default="\t", metavar="delimiter", help="Delimiter to read the matrix.")
## append to current file
pArgs <- add_option( pArgs, c("-a", "--append"), type="character", action="store_true", dest="fAppend", default=FALSE, metavar="append", help="Append to existing data. Default no append.")
## append to current file
pArgs <- add_option( pArgs, c("-t", "--tsv"), type="character", action="store", dest="strTSV", default=NA, metavar="TSVfile", help="TSV file to give.")

### Parse arguments
lsArgs <- parse_args( pArgs, positional_arguments = TRUE )

#Get positional arguments
if( !length( lsArgs$args ) == 3 ) { stop( print_help( pArgs ) ) }

### Write to file the read config script
funcWriteMatrixToReadConfigFile(strConfigureFileName=lsArgs$args[1], strMatrixFile=lsArgs$options$strTSV, strMatrixName=lsArgs$args[2], strRowIndices=lsArgs$options$strRows, strColIndices=lsArgs$options$strColumns,
  strDtCharacter=lsArgs$options$strCharData, strDtFactoral=lsArgs$options$strFactorData, strDtInteger=lsArgs$options$srtIntegerData, strDtLogical=lsArgs$options$strLogicalData, strDtNumeric=lsArgs$options$strNumericData, strDtOrdered=lsArgs$options$strOrderData,
  acharDelimiter=lsArgs$options$charDelimiter, fAppend=lsArgs$options$fAppend)