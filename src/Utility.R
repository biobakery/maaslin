#Modifies taxa names for plotting
#@params tempTaxaNames A list of string taxa names
funcRenameTaxa = function(tempTaxaNames)
{
  returnList=c()
  newName=""
  for( name in tempTaxaNames )
  {
    modifiedName=strsplit(name, "\\|")[[1]]
    i=length(modifiedName)
    if((modifiedName[i] == "unclassified") || !is.na(as.numeric(modifiedName[i])))
    {
      newName=paste(modifiedName[(i-1):i], collapse = "|")
    }else{
      newName=modifiedName[i]
    }
    returnList=c(returnList, newName)
  }
  return(returnList)
}

#Generates a given number of random colors
#@params tempNumberColors Number of colors to generate
funcGetRandomColors=function(tempNumberColors=1)
{
  adRet <- c()
  for(i in 1:tempNumberColors)
  {
    adRGB <- ( runif( 3 ) * 0.66 ) + 0.33
    adRet <- c(adRet, rgb( adRGB[1], adRGB[2], adRGB[3] ))
  }
  return(adRet)
}

funcColorHelper <- function( dMax = 1, dMin = -1, dMed = NA )
{
  if( dMax < dMin )
  {
    d <- dMax
    dMax <- dMin
    dMin <- d
  }
  if( is.na( dMed ) )
  {
    dMed <- ( dMin + dMax ) / 2
  }
  return( list( dMin = dMin, dMax = dMax, dMed = dMed ) ) }
	
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

#Write a string or a table of data
#pOut String or table to write to file
#strFile String name of file
funcWrite <- function( pOut, strFile )
{
#  write( pOut, strFile, ncolumns = length( astrOut ), append = TRUE ) }
  if( length( intersect( class( pOut ), c("character", "numeric") ) ) )
  {
    write.table( t(pOut), strFile, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE, na = "", append = TRUE )
  } else {
    capture.output( print( pOut ), file = strFile, append = TRUE )
  }
}

#Log a table to a file
#frmeTable Table to write
#strFile Log file
#fAppend Should appending occur.
funcWriteTable <- function( frmeTable, strFile, fAppend = FALSE )
{
  write.table( frmeTable, strFile, quote = FALSE, sep = "\t", na = "", col.names = NA, append = fAppend )
}
