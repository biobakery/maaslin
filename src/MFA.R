####################################
# Summary: Multifactor Analysis
# Author: Timothy Tickle
# Start Date: 11-21-2011 current
####################################

#Import library
suppressMessages(library( FactoMineR, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))

### 
### frmeData: Data
### aiMetadata: Indicies or metadata to use
### aiBugs: Bugs to 
funcMFA = function( frmeData, aiMetadata, aiBugs, aiGenes = c() )
{
  #Update aiBugs with factor and gene data that is being plotted
  #If a custom plotting function is given, then use, otherwise give defaults
  lsFeatures = c()
  if(exists("funcPlotFeatures",mode="function")){lsFeatures = funcPlotFeatures()}
  lsasMetadata = ifelse(exists("funcPlotMetadata",mode="function"), funcPlotMetadata(), c(list(asNames=c(),asLabels=c())))

  #Get indicies of data to plot and add to those given as parameters
  liFeatures = which(colnames(frmeData) %in% lsFeatures)
  liMetadata = which(colnames(frmeData) %in% lsasMetadata$asNames)
  aiMetadata = union(aiMetadata,liMetadata)
  aiBugs = union(aiBugs,liFeatures)

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

  #Needing to set if the data is numeric or categorical for the MFA
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

  aiCols <- c()
  aiLengths <- c()
  astrTypes <- c()
  astrNames <- c()
  for( lsCur in list(
                  list(aiCols = aiMDN, strType = "n", strName = "metadata_nom"),
                  list(aiCols = aiMDC, strType = "c", strName = "metadata_cont"),
                  list(aiCols = aiGenetics, strType = "n", strName = "genetics"),
                  list(aiCols = aiBugs, strType = "c", strName = "taxa")) )
  {
    aiCur <- lsCur$aiCols
    if( !length( aiCur ) )
    {
      next
    }
    aiCols <- c(aiCols, aiCur)
    aiLengths <- c(aiLengths, length( aiCur ))
    astrTypes <- c(astrTypes, lsCur$strType)
    astrNames <- c(astrNames, lsCur$strName)
  }

  #MFA requires
  lsRet <- try( MFA( frmeMFA[,aiCols], group = aiLengths, type = astrTypes, name.group = astrNames, graph = FALSE ) )
  return( lsRet )
}

### Plot Multiple Factoral Analysis
### lsMFA: MFA output object
### tempSaveFileName: File to save as a pdf (.pdf extension will be added
### tempWidth: Width of pdf graphics region in inches
### tempHeight: Heigth of pdf graphics region in inches
### tempPCH: Numeric pch symbol
### tempBackgroundColor: Pdf image background color
### tempForeground: Pdf image foreground color
### tempAxesColor: Pdf image axes color
### tempLabelColor: Pdf image label color
### tempTitleColor: Pdf image title color
### tempSubtitleColor: Pdf image subtitle color
funcPlotMFA <- function(lsMFA, fInvert = FALSE, tempSaveFileName="MFA", funcPlotColors = NULL, funcPlotPoints = NULL, funcPlotLegend = NULL, tempWidth=(tempHeight*1.5), tempHeight=6, tempPCH=20)
{
  #Set pdf settings
  pdf(paste(tempSaveFileName,".pdf",sep=""), width = c_dHeight * 1.5, height = c_dHeight, useDingbats=FALSE )
  if( fInvert )
  {
    par( bg = "black", fg = "white", col.axis = "white", col.lab = "white", col.main = "white", col.sub = "white" )
  }
  #Get MFA pca data and set dimensions
  lsPCA = lsMFA$global.pca

  #Plot colors
  astrCols = ifelse(exists("funcPlotColors",mode="function"),funcPlotColors( frmeData ),funcGetRandomColors( length( lsPCA$ind$coord ) ))

  #Plot points
  aiPoints = ifelse(exists("funcPlotPoints",mode="function"),funcPlotPoints( frmeData ),16)

  dX1 = max( lsPCA$ind$coord[,1] ) / max( lsPCA$var$coord[,1] )
  dX2 = min( lsPCA$ind$coord[,1] ) / min( lsPCA$var$coord[,1] )
  dY1 = max( lsPCA$ind$coord[,2] ) / max( lsPCA$var$coord[,2] )
  dY2 = min( lsPCA$ind$coord[,2] ) / min( lsPCA$var$coord[,2] )

  #Scale the metadate labels so they are viewable
  dScaleFactor = ifelse(exists("funcGetMetadataScale",mode="function"),funcMetadataScale(),1)

  #Scale feature and metadata text locations so they are on the page
  dScale = dScaleFactor * min( abs( c(dX1, dX2, dY1, dY2) ) )
  dBugScaleFactor = ifelse(exists("funcGetFeatureScale",mode="function"),funcGetFeatureScale(),1)
  dBugScale = dBugScaleFactor * min( abs( c(dX1, dX2, dY1, dY2) ) )

  #Set X and Y coordinates and plot points
  strX <- sprintf( "Dimension 1 (%.2f%%)", lsPCA$eig$`percentage of variance`[1] )
  strY <- sprintf( "Dimension 2 (%.2f%%)", lsPCA$eig$`percentage of variance`[2] )

  #Get metadata names
  lsasMetadata = ifelse(exists("funcPlotMetadata",mode="function"),funcPlotMetadata(),
    list(asNames=lsMFA$summary.quali$modalite, asLabels=lsMFA$summary.quali$modalite))
  afMetadata <- rownames(lsPCA$var$coord) %in% lsasMetadata$asNames
  afLabels = lsasMetadata$asLabels

  #Get feature names
  lsFeaturesToPlot = ifelse(exists("funcPlotFeatures",mode="function"),funcPlotFeatures(),lsMFA$summary.quanti$variable)
  afFeatures = rownames(lsPCA$var$coord) %in% lsFeaturesToPlot

  #Order the metadata
  asOrderedPlotLabels = c()
  for(sMetadata in rownames(lsPCA$var$coord)[afMetadata])
  {
    asOrderedPlotLabels=c(asOrderedPlotLabels,afLabels[which(lsasMetadata$asNames==sMetadata)])
  }
  if(length(asOrderedPlotLabels) != length(afMetadata[afMetadata==TRUE]))
  {
    asOrderedPlotLabels = afMetadata
  }

  #Do plots
  #Plot just points 
  funcPlotMFAPage(coordinatesPlot=lsPCA$ind$coord, coordinatesText=NA, strX=strX, strY=strY, aiPoints=aiPoints, astrCols=astrCols)

  if( sum( afMetadata ) )
  {
    funcPlotMFAPage(coordinatesPlot=lsPCA$ind$coord, coordinatesText=lsPCA$var$coord, strX=strX, strY=strY, aiPoints=aiPoints,
      astrCols=astrCols, afMetadata=afMetadata, dScale=dScale, lsMetadataLabels=asOrderedPlotLabels, sLegendLoc="topright")
  }

  #Plot features
  if( sum( afFeatures ) )
  {
    funcPlotMFAPage(coordinatesPlot=lsPCA$ind$coord, coordinatesText=lsPCA$var$coord, strX=strX, strY=strY, aiPoints=aiPoints,
      astrCols=astrCols, afFeatures=afFeatures, dBugScale=dBugScale, lsFeatureLabels=funcRename( rownames( lsPCA$var$coord )[afFeatures] ),
      sLegendLoc="topright")

    #Plot metadata and features
    if( sum( afMetadata ) )
    {
      funcPlotMFAPage(coordinatesPlot=lsPCA$ind$coord, coordinatesText=lsPCA$var$coord, strX=strX, strY=strY, aiPoints=aiPoints,
        astrCols=astrCols, afMetadata=afMetadata, dScale=dScale, lsMetadataLabels=asOrderedPlotLabels, afFeatures=afFeatures,
        dBugScale=dBugScale, lsFeatureLabels=funcRename( rownames( lsPCA$var$coord )[afFeatures] ), sLegendLoc="topright")
    }
  }

  dev.off( )
}

funcPlotMFAPage <- function(coordinatesPlot, coordinatesText, strX, strY, aiPoints, astrCols, afMetadata=NA, dScale=NA, lsMetadataLabels=NA, afFeatures=NA, dBugScale=NA, lsFeatureLabels=NA, sLegendLoc="topright")
{
  plot( coordinatesPlot, pch = aiPoints, col = astrCols, xlab = strX, ylab = strY, cex.axis=1.5, cex.lab=1.5, cex=1.5)
  if( !is.na(afMetadata) && !is.na(dScale) && !is.na(lsMetadataLabels) )
  {
    text( coordinatesText[afMetadata,] * dScale, labels=lsMetadataLabels, cex=1.8, font = 2 )
  }
  if( !is.na( afFeatures ) && !is.na(dBugScale) && !is.na(lsFeatureLabels) )
  {
    text( coordinatesText[afFeatures,] * dBugScale, labels=lsFeatureLabels, cex=1.6, font = 3 )
  }
  if(exists("funcPlotLegend",mode="function"))
  {
    funcPlotLegend(sLegendLoc, NULL)
  }
}