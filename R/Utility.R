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
##description<< Collection of minor utility scripts
) { return( pArgs ) }

#source("Constants.R")

funcBonferonniCorrectFactorData <- function
### Bonferroni correct for factor data
(dPvalue,
### P-value to correct
vsFactors,
### Factors of the data to correct
fIgnoreNAs = TRUE
){
  vsUniqueFactors = unique( vsFactors )
  if( fIgnoreNAs ){ vsUniqueFactors = setdiff( vsUniqueFactors, c("NA","na","Na","nA") ) }
  return( dPvalue * max( 1, ( length( vsUniqueFactors ) - 1 ) ) )
  ### Numeric p-value that is correct for levels (excluding NA levels)
}

funcCalculateTestCounts <- function(
### Calculates the number of tests used in inference
iDataCount,
asMetadata,
asForced,
asRandom,
fAllvAll
){
  iMetadata = length(asMetadata)
  iForced = length(setdiff(intersect( asForced, asMetadata ), asRandom))
  iRandom = length(intersect( asRandom, asMetadata ))
  if(fAllvAll)
  {
    #AllvAll flow formula
    return((iMetadata-iForced-iRandom) * iDataCount)
  }

  #Normal flow formula
  return((iMetadata-iRandom) * iDataCount)
}

funcCoef2Col <- function(
### Searches through a dataframe and looks for a column that would match the coefficient
### by the name of the column or the column name and level appended together.
strCoef,
### String coefficient name
frmeData,
### Data frame of data
astrCols = c()
### Column names of interest (if NULL is given, all column names are inspected).
){
  #If the coefficient is the intercept there is no data column to return so return null
  if( strCoef %in% c("(Intercept)", "Intercept") ) { return( NULL ) }
  #Remove ` from coefficient
  strCoef <- gsub( "`", "", strCoef )

  #If the coefficient name is not in the data frame
  if( !( strCoef %in% colnames( frmeData ) ) )
  {
    fHit <- FALSE
    #If the column names are not provided, use the column names of the dataframe.
    if( is.null( astrCols ) ){astrCols <- colnames( frmeData )}

    #Search through the different column names (factors)
    for( strFactor in astrCols )
    {
      #Select a column, if it is not a factor or does not begin with the factor's name then skip
      adCur <- frmeData[,strFactor]
      if( ( class( adCur ) != "factor" ) ||
        ( substr( strCoef, 1, nchar( strFactor ) ) != strFactor ) ) { next }

      #For the factors, create factor-level name combinations to read in factors
      #Then check to see the factor-level combination is the coefficient of interest
      #If it is then store that factor as the coefficient of interest
      #And break
      for( strValue in levels( adCur ) )
      {
        strCur <- paste( strFactor, strValue, sep = c_sFactorNameSep )
        if( strCur == strCoef )
        {
          strCoef <- strFactor
          fHit <- TRUE
          break
        }
      }

      #If the factor was found, return
      if( fHit ){break }
    }
  }

  #If the original coefficient or the coefficient factor combination name are in the
  #data frame, return the name. Otherwise return NA.
  return( ifelse( ( strCoef %in% colnames( frmeData ) ), strCoef, NA ) )
  ### Coefficient name
}

funcColorHelper <- function(
### Makes sure the max is max and the min is min, and dmed is average
dMax = 1,
### Max number
dMin = -1,
### Min number
dMed = NA
### Average value
){
  #Make sure max is max and min is min
  vSort = sort(c(dMin,dMax))
  return( list( dMin = vSort[1], dMax = vSort[2], dMed = ifelse((is.na(dMed)), (dMin+dMax)/2.0, dMed ) ))
  ### List of min, max and med numbers
}

funcColor <- function(
### Generate a color based on a number that is forced to be between a min and max range.
### The color is based on how far the number is from the center of the given range
### From red to green (high) are produced with default settings
dX,
### Number from which to generate the color
dMax = 1,
### Max possible value
dMin = -1,
### Min possible value
dMed = NA,
### Central value if you don't want to be the average
adMax = c(1, 1, 0),
### Is used to generate the color for the higher values in the range, this can be changed to give different colors set to green
adMin = c(0, 0, 1),
### Is used to generate the color for the lower values in the range, this can be changed to give different colors set to red
adMed = c(0, 0, 0)
### Is used to generate the color for the central values in the range, this can be changed to give different colors set to black
){
  lsTmp <- funcColorHelper( dMax, dMin, dMed )
  dMax <- lsTmp$dMax
  dMin <- lsTmp$dMin
  dMed <- lsTmp$dMed
  if( is.na( dX ) )
  {
    dX <- dMed
  }
  if( dX > dMax )
  {
    dX <- dMax
  } else if( dX < dMin )
  {
    dX <- dMin }
  if( dX < dMed )
  {
    d <- ( dMed - dX ) / ( dMed - dMin )
    adCur <- ( adMed * ( 1 - d ) ) + ( adMin * d )
  } else {
    d <- ( dMax - dX ) / ( dMax - dMed )
    adCur <- ( adMed * d ) + ( adMax * ( 1 - d ) )
  }
  return( rgb( adCur[1], adCur[2], adCur[3] ) )
  ### RGB object
}

funcGetColor <- function(
### Get a color based on col parameter
) {
  adCol <- col2rgb( par( "col" ) )
  return( sprintf( "#%02X%02X%02X", adCol[1], adCol[2], adCol[3] ) )
  ### Return hexadecimal color
}


funcWrite <- function(
### Write a string or a table of data
### This transposes a table before it is written
pOut,
### String or table to write
strFile
### File to which to write
){
  if(!is.na(strFile))
  {
   if( length( intersect( class( pOut ), c("character", "numeric") ) ) )
    {
      write.table( t(pOut), strFile, quote = FALSE, sep = c_cTableDelimiter, col.names = FALSE, row.names = FALSE, na = "", append = TRUE )
    } else {
      capture.output( print( pOut ), file = strFile, append = TRUE )
    }
  }
}

funcWriteTable <- function(
### Log a table to a file
frmeTable,
### Table to write
strFile,
### File to which to write
fAppend = FALSE
### Append when writing
){
  if(!is.na(strFile))
  {
    write.table( frmeTable, strFile, quote = FALSE, sep = c_cTableDelimiter, na = "", col.names = NA, append = fAppend )
  }
}

funcWriteQCReport <- function(
### Write out the quality control report
strProcessFileName,
###  File name
lsQCData,
### List of QC data generated by maaslin to be written
liDataDim,
### Dimensions of the data matrix
liMetadataDim
### Dimensions of the metadata matrix
){
  unlink(strProcessFileName)
  funcWrite( paste("Initial Metadata Matrix Size: Rows ",liMetadataDim[1],"  Columns ",liMetadataDim[2],sep=""), strProcessFileName )
  funcWrite( paste("Initial Data Matrix Size: Rows ",liDataDim[1],"  Columns ",liDataDim[2],sep=""), strProcessFileName )
  funcWrite( paste("\nInitial Data Count: ",length(lsQCData$aiDataInitial),sep=""), strProcessFileName )
  funcWrite( paste("Initial Metadata Count: ",length(lsQCData$aiMetadataInitial),sep=""), strProcessFileName )
  funcWrite( paste("Data Count after preprocess: ",length(lsQCData$aiAfterPreprocess),sep=""), strProcessFileName )
  funcWrite( paste("Removed for missing metadata: ",length(lsQCData$iMissingMetadata),sep=""), strProcessFileName )
  funcWrite( paste("Removed for missing data: ",length(lsQCData$iMissingData),sep=""), strProcessFileName )
  funcWrite( paste("Number of data with outliers: ",length(which(lsQCData$aiDataSumOutlierPerDatum>0)),sep=""), strProcessFileName )
  funcWrite( paste("Number of metadata with outliers: ",length(which(lsQCData$aiMetadataSumOutlierPerDatum>0)),sep=""), strProcessFileName )
  funcWrite( paste("Metadata count which survived clean: ",length(lsQCData$aiMetadataCleaned),sep=""), strProcessFileName )
  funcWrite( paste("Data count which survived clean: ",length(lsQCData$aiDataCleaned),sep=""), strProcessFileName )
  funcWrite( paste("\nBoostings: ",lsQCData$iBoosts,sep=""), strProcessFileName )
  funcWrite( paste("Boosting Errors: ",lsQCData$iBoostErrors,sep=""), strProcessFileName )
  funcWrite( paste("LMs with no terms suriving boosting: ",lsQCData$iNoTerms,sep=""), strProcessFileName )
  funcWrite( paste("LMs performed: ",lsQCData$iLms,sep=""), strProcessFileName )
  if(!is.null(lsQCData$lsQCCustom))
  {
    funcWrite("Custom preprocess QC data: ", strProcessFileName )
    funcWrite(lsQCData$lsQCCustom, strProcessFileName )
  } else {
    funcWrite("No custom preprocess QC data.", strProcessFileName )
  }
  funcWrite( "\n#Details###########################", strProcessFileName )
  funcWrite("\nInitial Data Count: ", strProcessFileName )
  funcWrite(lsQCData$aiDataInitial, strProcessFileName )
  funcWrite("\nInitial Metadata Count: ", strProcessFileName )
  funcWrite(lsQCData$aiMetadataInitial, strProcessFileName )
  funcWrite("\nData Count after preprocess: ", strProcessFileName )
  funcWrite(lsQCData$aiAfterPreprocess, strProcessFileName )
  funcWrite("\nRemoved for missing metadata: ", strProcessFileName )
  funcWrite(lsQCData$iMissingMetadata, strProcessFileName )
  funcWrite("\nRemoved for missing data: ", strProcessFileName )
  funcWrite(lsQCData$iMissingData, strProcessFileName )
  funcWrite("\nDetailed outlier indices: ", strProcessFileName )
  for(sFeature in names(lsQCData$liOutliers))
  {
    funcWrite(paste("Feature",sFeature,"Outlier indice(s):", paste(lsQCData$liOutliers[[sFeature]],collapse=",")), strProcessFileName )
  }
  funcWrite("\nMetadata which survived clean: ", strProcessFileName )
  funcWrite(lsQCData$aiMetadataCleaned, strProcessFileName )
  funcWrite("\nData which survived clean: ", strProcessFileName )
  funcWrite(lsQCData$aiDataCleaned, strProcessFileName )
}

funcFormulaStrToList <- function(
#Takes a lm or mixed model formula and returns a list of covariate names in the formula
strFormula
#Formula to extract covariates from
){
  #Return list
  lsRetComparisons = c()

  #If you get a null or na just return
  if(is.null(strFormula)||is.na(strFormula)){return(lsRetComparisons)}

  #Get test comparisons (predictor names from formula string)
  asComparisons = gsub("`","",setdiff(unlist(strsplit(unlist(strsplit(strFormula,"~"))[2]," ")),c("","+")))

  #Change metadata in formula to univariate comparisons
  for(sComparison in asComparisons)
  {
    #Removed random covariate formating
    lsParse = unlist(strsplit(sComparison, "[\\(\\|\\)]", perl=FALSE))
    lsRetComparisons = c(lsRetComparisons,lsParse[length(lsParse)])
  }
  return(lsRetComparisons)
}

