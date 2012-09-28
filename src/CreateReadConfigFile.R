#!/usr/bin/env Rscript

inlinedocs <- function(
##author<< Timothy Tickle <ttickle@hsph.harvard.edu>
) { return( pArgs ) }

### Logging class
library( logging )
### Class for commandline argument processing
library( optparse )

### Source the IO.R forthe script
source("IO.R")

### Create command line argument parser
### The TSV (tab seperated value (column major, samples are rows) file that will be read in
### The column that is the last metadata name
### The read.config file that will be used to read in the TSV file
pArgs <- OptionParser( usage = "%prog <strInputTsv> <strLastMetadata> <strOutputRC>" )


### Parse arguments
lsArgs <- parse_args( pArgs, positional_arguments = TRUE )

#Get positional arguments
if( !length( lsArgs$args ) == 3 ) { stop( print_help( pArgs ) ) }

### Write to file the read config script
funcWriteMatrixToReadConfigFile = function(strConfigureFileName, strMatrixName, strRowIndices = "-", strColIndices "-",
  strDtCharacter="", strDtFactoral="", strDtInteger="", strDtLogical="", strDtNumeric="", strDtOrdered="",
  acharDelimiter="\t", fAppend=FALSE)


### Ren needs to call this for multiple microbiome of the same file
### Based on the input settings we can get this info
### Need to expose it all

### George needs to call this for the tsv file that is created.
### 
funcWriteMatrixToReadConfigFile = function(strConfigureFileName, strMatrixName="Metadata", strColIndices "-", fAppend=FALSE)

