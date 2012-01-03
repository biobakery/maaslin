####################################
# Summary: Multifactor Analysis
# Author: Timothy Tickle
# Start Date: 11-21-2011
####################################

#External libraries
library( FactoMineR )

#@params 
funcMFA = function( frmeData, aiMetadata, aiBugs, aiGenes = c() )
{
  #Check / Select numeric data
  #Eventually holds numeric rows to be extracted from the data frame
  aiRows = c()
  #For each row in the data frame
  for( iRow in 1:nrow( frmeData ) )
  {

    #Indicates if the row passed and is good
    rowOK = TRUE
    #Look through all the data rows of each metadata, bug, or gene of interest
    for( aiCols in list(aiMetadata, aiBugs, aiGenes) )
    {

      #Attempt to convert to numeric
      aData = as.numeric( frmeData[iRow, aiCols] )

      #Check to see if is not NA
      naSumOne = sum( !is.na( aData ) )
      naSumTwo = sum( aData != "NA" )
      iNNA = 0
      if(is.na(naSumOne))
      {
        if(is.na(naSumTwo))
        {
          #Both NA
          iNNA = 0
        }else{
          #naSumOne NA
          iNNA = naSumTwo
        }
      }else{
        if(is.na(naSumTwo))
        {
          #naSumTwo NA
          iNNA = naSumOne
        }else{
          #Both not NA
          iNNA = max( naSumOne, naSumTwo )
        }
      }

      #If there is data and the amount of NA is under a certain threshold
      if( length( aiCols ) && ( iNNA < ( c_dMinSamp * length( aiCols ) ) ) )
      {
        rowOK = FALSE
        break
      }
    }
    #if the row is ok, then add it to the rows that will be used
    if(rowOK)
    {
      aiRows = c(aiRows, iRow)
    }
  }
  #Reduce the data frame to rows
  frmeMFA = frmeData[aiRows,]

  #Check / Select factor data
  aiMDN <- c()
  aiMDC <- c()
  for( iCol in c(aiMetadata, aiBugs, aiGenes) )
  {
    adCol <- frmeMFA[,iCol]
    if( class( adCol ) == "factor" )
    {
      aCol <- frmeMFA[,iCol]
      aCol[aCol == "NA"] <- NA
      frmeMFA[,iCol] <- factor( aCol, levels = setdiff( levels( aCol ), "NA" ) )
      aiMDN <- c(aiMDN, iCol)
      next
    }
    if( iCol %in% aiGenes ) { next }
    if( iCol %in% aiBugs )
    {
      next
#      adCol <- funcTransform( adCol )
    } else {
      aiMDC <- c(aiMDC, iCol)
    }
    dSD <- sd( adCol, na.rm = TRUE )
    if( is.na( dSD ) )
    {
      dSD=1
    }
    frmeMFA[,iCol] <- ( adCol - mean( adCol, na.rm = TRUE ) ) / dSD
  }

  for( iCol in c(aiMDC, aiBugs, aiGenes) )
  {
    adCol <- frmeMFA[,iCol]
    if( class( adCol ) != "factor" )
    {
      frmeMFA[,iCol][is.na( adCol )] <- 0 
    }
  }
  if( length( aiGenes ) )
  {
    aiCols <- c(aiMDN, aiMDC, aiGenes, aiBugs)
    aiLengths <- c(length( aiMDN ), length( aiMDC ), length( aiGenes ), length( aiBugs ))
    astrTypes <- c("n", "c", ifelse( class( frmeMFA[,aiGenes[1]] ) == "factor", "n", "s" ), "c")
    astrNames <- c("metadata_nom", "metadata_cont", "genetics", "taxa")
  } else {
    aiCols <- c(aiMDN, aiMDC, aiBugs)
    aiLengths <- c(length( aiMDN ), length( aiMDC ), length( aiBugs ))
    astrTypes <- c("n", "c", "c")
    astrNames <- c("metadata_nom", "metadata_cont", "taxa")
  }
  #MFA requires
  lsRet <- try( MFA( frmeMFA[,aiCols], group = aiLengths, type = astrTypes, name.group = astrNames, graph = FALSE ) )
  return( lsRet )
}

#Plot Multiple Factoral Analysis
#@params lsMFA MFA output object
#@params tempSaveFileName File to save as a pdf (.pdf extension will be added
#@params tempWidth Width of pdf graphics region in inches
#@params tempHeight Heigth of pdf graphics region in inches
#@parmas tempPCH Numeric pch symbol
#@params tempBackgroundColor Pdf image background color
#@params tempForeground Pdf image foreground color
#@params tempAxesColor Pdf image axes color
#@params tempLabelColor Pdf image label color
#@params tempTitleColor Pdf image title color
#@params tempSubtitleColor Pdf image subtitle color
funcPlotMFA <- function(lsMFA, tempSaveFileName="MFA", tempWidth=(tempHeight*1.5), tempHeight=6, tempPCH=20, tempBackgroundColor="white", tempForeground="black", tempAxesColor="black", tempLabelColor="black", tempTitleColor="black", tempSubtitleColor="black")
{
  #Set pdf settings
  pdf(paste(tempSaveFileName,".pdf",sep=""), width=tempWidth, height=tempHeight)
  par(bg=tempBackgroundColor, fg=tempForeground, col.axis=tempAxesColor, col.lab=tempLabelColor, col.main=tempTitleColor, col.sub=tempSubtitleColor)

  #Get MFA pca data and set dimensions
  lsPCA = lsMFA$global.pca
  pointColors = funcGetRandomColors( length( lsPCA$ind$coord ) )
  dX1 = max( lsPCA$ind$coord[,1] ) / max( lsPCA$var$coord[,1] )
  dX2 = min( lsPCA$ind$coord[,1] ) / min( lsPCA$var$coord[,1] )
  dY1 = max( lsPCA$ind$coord[,2] ) / max( lsPCA$var$coord[,2] )
  dY2 = min( lsPCA$ind$coord[,2] ) / min( lsPCA$var$coord[,2] )
  dScale = 0.95 * min( abs( c(dX1, dX2, dY1, dY2) ) )

  #Set X and Y coordinates and plot points
  tempXTitle = sprintf( "Dimension 1 (%.2f%%)", lsPCA$eig$`percentage of variance`[1] )
  tempYTitle = sprintf( "Dimension 2 (%.2f%%)", lsPCA$eig$`percentage of variance`[2] )
  plot( lsPCA$ind$coord, pch=tempPCH, col=pointColors, xlab=tempXTitle, ylab=tempYTitle )

  #Plot metadata
  metadata = rownames(lsPCA$var$coord) %in% lsMFA$summary.quali$modalite
  plot( lsPCA$ind$coord, pch=tempPCH, col = pointColors, xlab = tempXTitle, ylab = tempYTitle )
  text( lsPCA$var$coord[metadata,] * dScale, labels = funcRenameTaxa( rownames( lsPCA$var$coord )[metadata] ), font = 2 )

  #Plot factors
  plot( lsPCA$ind$coord, pch=tempPCH, col = pointColors, xlab = tempXTitle, ylab = tempYTitle )
  text( lsPCA$var$coord[!metadata,] * dScale, labels = funcRenameTaxa( rownames( lsPCA$var$coord )[!metadata] ), cex = 0.75, font = 3 )

  plot( lsPCA$ind$coord, pch=tempPCH, col = pointColors, xlab = tempXTitle, ylab = tempYTitle )
  text( lsPCA$var$coord * dScale, labels = funcRenameTaxa( rownames( lsPCA$var$coord ) ), cex = ifelse(metadata, 1, 0.75 ), font = ifelse(metadata, 2, 3 ) )
  dev.off()
}

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
