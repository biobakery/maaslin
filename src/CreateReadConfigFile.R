#!/usr/bin/env Rscript

inlinedocs <- function(
##author<< Timothy Tickle <ttickle@hsph.harvard.edu>
) { return( pArgs ) }

### Logging class
library( logging )
### Class for commandline argument processing
library( optparse )

### Source the IO.R for the script
source("IO.R")

### Create command line argument parser
### The TSV (tab seperated value (column major, samples are rows) file that will be read in
### The column that is the last metadata name
### The read.config file that will be used to read in the TSV file
pArgs <- OptionParser( usage = "%prog [optional] <strOutputRC> <strInputTsv> <strMatrixName>" )

# Settings for Read config
## row indices
pArgs <- add_option( pArgs, c("-r", "--rows"), type="character", action="store", dest="strRows", default="-", metavar="row_indices", help="Rows to read by index starting with 1.")
## column indices
pArgs <- add_option( pArgs, c("-c", "--columns"), type="character", action="store", dest="strColumns", default="-", metavar="column_indices", help="Columns to read in by index starting with 1.")
## Character rows
pArgs <- add_option( pArgs, c("-s", "--stringdata"), type="character", action="store", dest="strCharData", default="", metavar="char_data", help="Forces these columns to be character data.")
## Factor data rows
pArgs <- add_option( pArgs, c("-f", "--factordata"), type="character", action="store", dest="strFactorData", default="", metavar="factor_data", help="Forces these columns to be factor data.")
## Integer data rows
pArgs <- add_option( pArgs, c("-i", "--integerdata"), type="character", action="store", dest="srtIntegerData", default="", metavar="integer_data", help="Forces these columns to be integer data.")
## logical data rows
pArgs <- add_option( pArgs, c("-l", "--logicaldata"), type="character", action="store", dest="strLogicalData", default="", metavar="logical_data", help="Forces these columns to be logical data.")
## numeric data rows
pArgs <- add_option( pArgs, c("-mn", "--numericdata"), type="character", action="store", dest="strNumericData", default="", metavar="numeric_data", help="Forces these columns to be numeric data.")
## orderd data rows
pArgs <- add_option( pArgs, c("-o", "--orderdata"), type="character", action="store", dest="strOrderData", default="", metavar="order_data", help="Forces these columns to be order data.")
## delimiter
pArgs <- add_option( pArgs, c("-d", "--delimiter"), type="character", action="store", dest="charDelimiter", default="\t", metavar="delimiter", help="Delimiter to read the matrix.")
## append to current file
pArgs <- add_option( pArgs, c("-a", "--append"), type="character", action="store_true", dest="fAppend", default=FALSE, metavar="append", help="Append to existing data. Default no append.")

### Parse arguments
lsArgs <- parse_args( pArgs, positional_arguments = TRUE )

#Get positional arguments
if( !length( lsArgs$args ) == 3 ) { stop( print_help( pArgs ) ) }

### Write to file the read config script
funcWriteMatrixToReadConfigFile = function(strConfigureFileName=lsArgs$args[1], strMatrixFile=lsArgs$args[3], strMatrixName=lsArgs$args[2], strRowIndices=lsArgs$strRows, strColIndices=lsArgs$strColumns,
  strDtCharacter=lsArgs$strCharData, strDtFactoral=lsArgs$strFactorData, strDtInteger=lsArgs$srtIntegerData, strDtLogical=lsArgs$strLogicalData, strDtNumeric=lsArgs$strNumericData, strDtOrdered=lsArgs$strOrderData,
  acharDelimiter=lsArgs$charDelimiter, fAppend=lsArgs$fAppend)