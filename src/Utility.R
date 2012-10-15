#source(file.path("input","maaslin","src","Constants.R"))

source("Constants.R")

### Modifies labels for plotting
### If the name is not an otu collapse to the last two clades
### Otherwise use the most terminal clade
### astrNames Names to modify for plotting
funcRename <- function( astrNames )
{
  astrRet <- c()
  for( strName in astrNames )
  {
    astrName <- strsplit( strName, c_cFeatureDelimREx )[[1]]
    i <- length( astrName )
    if( ( astrName[i] == c_strUnclassified ) || !is.na( as.numeric( astrName[i] ) ) )
    {
      strRet <- paste( astrName[( i - 1 ):i], collapse = c_cFeatureDelim )
    } else {
    strRet <- astrName[i]
    }
    astrRet <- c(astrRet, strRet)
  }
  return( astrRet )
}

#0 Tests
#Generates a given number of random colors
#tempNumberColors Number of colors to generate
funcGetRandomColors=function(tempNumberColors=1)
{
  adRet = c()
  return(sapply(1:tempNumberColors, function(x){
    adRGB <- ( runif( 3 ) * 0.66 ) + 0.33
    adRet <- c(adRet, rgb( adRGB[1], adRGB[2], adRGB[3] ))
  }))
}

### Searches through a dataframe and looks for a column that would match the coefficient
### by the name of the column or the column name and level appended together.
### StrCoef String coefficient name
### frmeData Data frame of data
### aStrCols Column names of interest (if NULL is given, all column names are inspected.
funcCoef2Col <- function( strCoef, frmeData, astrCols = c() )
{
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
}

### Makes sure the max is max and the min is min, and dmed is average
### dMax Max number
### dMin Min number
### dMed Average value
funcColorHelper <- function( dMax = 1, dMin = -1, dMed = NA )
{
  #Make sure max is max and min is min
  vSort = sort(c(dMin,dMax))

  return( list( dMin = vSort[1], dMax = vSort[2], dMed = ifelse((is.na(dMed)), (dMin+dMax)/2.0, dMed ) )) }

#0 TestCases
### Generate a color based on a number that is forced to be between a min and max range.
### The color is based on how far the number is from the center of the given range
### From red to green (high) are produced with default settings
### dX Number to generate color from
### dMax Max value the number can be
### dMin Min value the number can be
### dMed central value if you don't want to to be the average
### adMax Is used to generate the color for the higher values in the range, this can be changed to give different colors set to green
### adMin Is used to generate the color for the lower values in the range, this can be changed to give different colors set to red
### adMed Is used to generate the color for the central values in the range, this can be changed to give different colors set to black
funcColor <- function( dX, dMax = 1, dMin = -1, dMed = NA, adMax = c(1, 1, 0), adMin = c(0, 0, 1), adMed = c(0, 0, 0) )
{
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
}

funcColors <- function( dMax = 1, dMin = -1, dMed = NA, adMax = c(1, 1, 0), adMin = c(0, 0, 1), adMed = c(0, 0, 0), iSteps = 64 )
{
  lsTmp <- funcColorHelper( dMax, dMin, dMed )
  dMax <- lsTmp$dMax
  dMin <- lsTmp$dMin
  dMed <- lsTmp$dMed
  aRet <- c ()
  for( dCur in seq( dMin, dMax, ( dMax - dMin ) / ( iSteps - 1 ) ) )
  {
    aRet <- c(aRet, funcColor( dCur, dMax, dMin, dMed, adMax, adMin, adMed ))
  }
  return( aRet )
}

#0 Testcases
### Get a color based on col parameter
funcGetColor <- function( ) 
{
  adCol <- col2rgb( par( "col" ) )
  return( sprintf( "#%02X%02X%02X", adCol[1], adCol[2], adCol[3] ) )
}

#Remove whitespace at the beginning or the end of a string
#tempString String to be trimmed.
funcTrim=function(tempString)
{
  return(gsub("^\\s+|\\s+$","",tempString))
}

### Write a string or a table of data
### This transposes a table before it is written
### Always appends
### pOut String or table to write to file
### strFile String name of file
funcWrite <- function( pOut, strFile )
{
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

#Log a table to a file
#frmeTable Table to write
#strFile Log file
#fAppend Should appending occur.
funcWriteTable <- function( frmeTable, strFile, fAppend = FALSE )
{
  write.table( frmeTable, strFile, quote = FALSE, sep = c_cTableDelimiter, na = "", col.names = NA, append = fAppend )
}

#
#funcFormat <- function( pValue ) {
#	c_iMax	<- 8192
#
#	strRet <- paste( format( pValue ), collapse = " " )
#	if( nchar( strRet ) > c_iMax ) {
#	return( strRet )
#}

### Write out the quality control report
### strProcessFilename File name
### lsQCData List of QC data generated by maaslin to be written
### liDataDim Dimensions of the data matrix
### liMetadataDim Dimensions of the metadata matrix
#funcWriteQCReport <- function(strProcessFileName, lsQCData, liDataDim, liMetadataDim)
#{
#  unlink(strProcessFileName)
#  funcWrite( paste("Initial Metadata Matrix Size: Rows ",liMetadataDim[1],"  Columns ",liMetadataDim[2],sep=""), strProcessFileName )
#  funcWrite( paste("Initial Data Matrix Size: Rows ",liDataDim[1],"  Columns ",liDataDim[2],sep=""), strProcessFileName )
#  funcWrite( paste("\nInitial Data Count: ",length(lsQCData$aiDataInitial),sep=""), strProcessFileName )
#  funcWrite( paste("Initial Metadata Count: ",length(lsQCData$aiMetadataInitial),sep=""), strProcessFileName )
#  funcWrite( paste("Data Count after preprocess: ",length(lsQCData$aiAfterPreprocess),sep=""), strProcessFileName )
#  funcWrite( paste("Removed for missing metadata: ",length(lsQCData$iMissingMetadata),sep=""), strProcessFileName )
#  funcWrite( paste("Removed for missing data: ",length(lsQCData$iMissingData),sep=""), strProcessFileName )
#  funcWrite( paste("Data with outliers: ",length(lsQCData$aiSumOutlierPerDatum[lsQCData$aiSumOutlierPerDatum>0]),sep=""), strProcessFileName )
#  funcWrite( paste("Metadata count which survived clean: ",length(lsQCData$aiMetadataCleaned),sep=""), strProcessFileName )
#  funcWrite( paste("Data count which survived clean: ",length(lsQCData$aiDataCleaned),sep=""), strProcessFileName )
#  funcWrite( paste("\nBoostings: ",lsQCData$iBoosts,sep=""), strProcessFileName )
#  funcWrite( paste("Boosting Errors: ",lsQCData$iBoostErrors,sep=""), strProcessFileName )
#  funcWrite( paste("LMs with no terms suriving boosting: ",lsQCData$iNoTerms,sep=""), strProcessFileName )
#  funcWrite( paste("LMs performed: ",lsQCData$iLms,sep=""), strProcessFileName )
#  if(!is.null(lsQCData$lsQCCustom))
#  {
#    funcWrite("Custom preprocess QC data: ", strProcessFileName )
#    funcWrite(lsQCData$lsQCCustom, strProcessFileName )
#  } else {
#    funcWrite("No custom preprocess QC data.", strProcessFileName )
#  }
#  funcWrite( "\n#Details###########################", strProcessFileName )
#  funcWrite("\nInitial Data Count: ", strProcessFileName )
#  funcWrite(lsQCData$aiDataInitial, strProcessFileName )
#  funcWrite("\nInitial Metadata Count: ", strProcessFileName )
#  funcWrite(lsQCData$aiMetadataInitial, strProcessFileName )
#  funcWrite("\nData Count after preprocess: ", strProcessFileName )
#  funcWrite(lsQCData$aiAfterPreprocess, strProcessFileName )
#  funcWrite("\nRemoved for missing metadata: ", strProcessFileName )
#  funcWrite(lsQCData$iMissingMetadata, strProcessFileName )
#  funcWrite("\nRemoved for missing data: ", strProcessFileName )
#  funcWrite(lsQCData$iMissingData, strProcessFileName )
#  funcWrite("\nOutlier Count per Datum: ", strProcessFileName )
#  funcWrite(lsQCData$aiSumOutlierPerDatum, strProcessFileName )
#  funcWrite("\nMetadata which survived clean: ", strProcessFileName )
#  funcWrite(lsQCData$aiMetadataCleaned, strProcessFileName )
#  funcWrite("\nData which survived clean: ", strProcessFileName )
#  funcWrite(lsQCData$aiDataCleaned, strProcessFileName )
#}
#
#funcLogMatrixRead = function()
#{
#  #Give feed back to user
#  if(tempLog)
#  {
#    print("Generating FeedBack.")
#    #Feedback from file
#    print(paste("Successfully loaded Data frame for matrix, details are:",sep=""))
#    print(paste("Matrix Name: ", tempMatrixName, sep=""))
#    print(paste("File Name: ", tempFileName, sep=""))
#    delimiterString = tempDelimiter
#    if(delimiterString == "\t"){delimiterString = "TAB"
#    }else if(delimiterString == " "){delimiterString = "SPACE"
#    }else if(delimiterString == "\r"){delimiterString = "RETURN"
#    }else if(delimiterString == "\n"){delimiterString = "ENDLINE"}
#    print(paste("Delimiter: ",delimiterString, sep=""))
#    if(!funcIsValid(tempIdRow))
#    {
#      print(paste("No IDs for rows given.", sep=""))
#    } else {
#      print(paste("Row ids found in the col index ",paste(tempIdRow,sep="",collapse=","),".",sep=""))
#      print(paste("Row ids are ",paste(row.names(dataMatrix),sep="",collapse=","),".",sep=""))
#    }
#    if(!funcIsValid(tempIdCol))
#    {
#      print(paste("No IDs for columns given.", sep=""))
#    } else {
#      print(paste("Column ids found in the row index ",paste(tempIdCol,sep="",collapse=","),".",sep=""))
#      print(paste("Column ids are ",paste(colnames(dataMatrix),sep="",collapse=","),".",sep=""))
#    }
#    if(!funcIsValid(tempRows))
#    {
#      print(paste("No column indicies for data given so all column data kept.", sep=""))
#    } else {
#      print(paste("Data column indices are ",paste(tempColumns,sep="",collapse=","),".", sep=""))
#    }
#    if(!funcIsValid(tempColumns))
#    {
#      print(paste("No row indicies for data given so all column data kept.", sep=""))
#    } else {
#      print(paste("Data row indices are ",paste(tempRows,sep="",collapse=","),".", sep=""))
#    }
#    if(funcIsValid(tempDtCharacter))
#    {
#      print(paste("Character data are found at column indices ",paste(tempDtCharacter,sep="",collapse=","),".", sep=""))
#    }
#    if(funcIsValid(tempDtFactor))
#    {
#      print(paste("Factor data are found at column indices ",paste(tempDtFactor,sep="",collapse=","),".", sep=""))
#    }
#    if(funcIsValid(tempDtInteger))
#    {
#      print(paste("Integer data are found at column indices ",paste(tempDtInteger,sep="",collapse=","),".", sep=""))
#    }
#    if(funcIsValid(tempDtLogical))
#    {
#      print(paste("Logical data are found at column indices ",paste(tempDtLogical,sep="",collapse=","),".", sep=""))
#    }
#    if(funcIsValid(tempDtNumeric))
#    {
#      print(paste("Numeric data are found at column indices ",paste(tempDtNumeric,sep="",collapse=","),".", sep=""))
#    }
#    if(funcIsValid(tempDtOrderedFactor))
#    {
#      print(paste("Ordered Factor data are found at column indices ",paste(tempDtOrderedFactor,sep="",collapse=","),".", sep=""))
#    }
#    if((!funcIsValid(tempDtOrderedFactor))&&(!funcIsValid(tempDtNumeric))&&(!funcIsValid(tempDtLogical))&&(!funcIsValid(tempDtInteger))&&(!funcIsValid(tempDtFactor))&&(!funcIsValid(tempDtCharacter)))
#    {
#      print(paste("No data type information was given, all data was assumed to be numeric.", sep=""))
#    }
#    #Feedback from object
#    print(paste("The shape of the date is ",paste(dim(dataMatrix),collapse=","),".", sep=""))
#    print(paste("The column names of the matrix are ",paste(colnames(dataMatrix),collapse=","),".", sep=""))
#    print(paste("The row names of the matrix are ",paste(row.names(dataMatrix),collapse=","),".", sep=""))
#    if(modeError == TRUE){ print("Please note errors occured on converting data modes of some rows. Please check output, unsuccessful conversion leaves data as default (character mode).")}
#  }
#}