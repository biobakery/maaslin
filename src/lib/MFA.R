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
##description<< Performs MultiFactorial Analysis and plots ordination.
) { return( pArgs ) }

#Import library
suppressMessages(library( FactoMineR, warn.conflicts=FALSE, quietly=TRUE, verbose=FALSE))

# Get the most influencial variables or individuals in the MFA output
# A value for iCount less than 1 give you all variables or individuals
funcGetHighestContribution = function(MFAOut,dfOriginal,fReturnValues=FALSE,iCount=-1,iOnlyInDimension=-1,lsNames=NULL)
{
  # Get and combine the first two components of the contributions unless a specific component is indicated
  adfFeatures = NULL
  if(iOnlyInDimension<1)
  {
    adfFeatures = rbind(MFAOut$global.pca$var$contrib[,1], MFAOut$global.pca$var$contrib[,2])
  } else {
    adfFeatures = MFAOut$global.pca$var$contrib[,iOnlyInDimension]
  }

  #Currently the data frame has two row entries for each column
  #I only want the greatest instance of contribution
  #Update the first row to the max value of the two rows in each column
  if(!is.null(dim(adfFeatures)))
  {
    for(i in 1:dim(adfFeatures)[2])
    {
      adfFeatures[1,i] = max(adfFeatures[,i])
    }
    #Get sorted names
    lsTargetGroupNames = names(adfFeatures[1,order(adfFeatures[1,],decreasing=TRUE)])
  } else {
    #Get sorted names
    lsTargetGroupNames = names(adfFeatures[order(adfFeatures,decreasing=TRUE)])
  }

  #Filter to features or variables
  #Change the values to column names because values are being used in mfa but names are given as covariates or column names
  if(is.null(lsNames)){lsNames=names(dfOriginal)}

  lsTargetGroupColumnNames = c()
  lsReducedValuesList = c()
  # Go through each of the values that contribute
  for( iTargetNameIndex in 1:length(lsTargetGroupNames))
  {
    #Change the value to a data column
    sCurrentColumnName = funcMFAValue2Col(lsTargetGroupNames[iTargetNameIndex],dfOriginal)

    #If the column name is in the covariates to keep then either keep the value or all values associated with the covariate.
    if(sCurrentColumnName %in% lsNames)
    {
      lsReducedValuesList = unique(c(lsReducedValuesList, as.character(dfOriginal[[sCurrentColumnName]])))
      lsTargetGroupColumnNames = c(lsTargetGroupColumnNames, sCurrentColumnName)
      if(length(lsTargetGroupColumnNames) == iCount){break}
    }
  }

  #This process can get values that were not influential (or in the components of interest, in that case remove those values)
  if(fReturnValues)
  {
    return(unique(lsReducedValuesList))
  } else {
    return(unique(lsTargetGroupColumnNames))
  }
}

#  if(is.null(iCount) || (iCount< 1) || (iCount>length(lsTargetGroupNames))){iCount = length(lsTargetGroupNames)}
#  lsTargetGroupNames = lsTargetGroupNames[1:iCount]
#  return(lsTargetGroupNames[!is.na(lsTargetGroupNames)])
#}

# Returns the first of a given data mode (discontinous or continuous)
funcGetFirstOfDataMode = function(lsSortedNames,dfData,fContinuous=TRUE)
{
  sMarkerMetadata = NULL

  for(i in 1:length(lsSortedNames))
  {
    # Get the column name from the mfa value which can be the column name or the value or a combination of both
    strColumnName = funcMFAValue2Col(lsSortedNames[i], dfData)

    sCurData = dfData[[strColumnName]]

    if(fContinuous)
    {
      if(is.numeric(sCurData) || is.integer(sCurData))
      {
        sMarkerMetadata = strColumnName
        break
      }
    } else {
      if(is.factor(sCurData) || is.ordered(sCurData) || is.logical(sCurData) || is.character(sCurData))
      {
        sMarkerMetadata = strColumnName
        break
      }
    }
  }
  # Return NA if nothing is found
  return(sMarkerMetadata)
}

# Get the most influencial variables or individuals in the MFA output
funcGetMarkerValues = function(MFAOut,dfOriginal,lsMetadata=NULL,strMFAColorCovariate=NULL,strMFAShapeCovariate=NULL,fPlotNA=FALSE)
{
  # Default marker values
  viMarkers = as.integer(c(18:15,0:14,"#","@"))

  # Returned values
  liMarkerShapes = NULL
  lsMarkerColors = NULL

  # Features to base color and shape on
  sFeatureContinuous = strMFAColorCovariate
  sFeature = strMFAShapeCovariate

  # Legend information
  lLegendInfo = list(sLocation="topright")

  if(is.null(sFeatureContinuous) || is.null(sFeature))
  {
    # Get all contributors in the first and second dim of the contributions in sorted order of contribution (as dataframe column names)
    lsSortedMetadataDim1 = funcGetHighestContribution(MFAOut,dfOriginal=dfOriginal,iOnlyInDimension=1,lsNames=lsMetadata)
    lsSortedMetadataDim2 = funcGetHighestContribution(MFAOut,dfOriginal=dfOriginal,iOnlyInDimension=2,lsNames=lsMetadata)

    # Get the column name from the mfa value which can be the column name or the value or a combination of both
    strColumnName = funcMFAValue2Col(lsSortedMetadataDim1[1], dfOriginal)

    # See if the first variable is continous or factor
    adMostInfluential = dfOriginal[[strColumnName]]

    if(is.logical(adMostInfluential) || is.factor(adMostInfluential) || is.character(adMostInfluential))
    {
      if(is.null(sFeature))
      {
        sFeature = lsSortedMetadataDim1[1]
        lLegendInfo[["sShapeMetadata"]] = sFeature
      }
    } else {
      if(is.null(sFeatureContinuous))
      {
        sFeatureContinuous = lsSortedMetadataDim1[1]
        lLegendInfo[["sColorMetatadata"]]  = sFeatureContinuous
      }
    }
  }

  if(is.null(sFeature))
  {
    sFeature = funcGetFirstOfDataMode(lsSortedNames=lsSortedMetadataDim2,dfData=dfOriginal,fContinuous=FALSE)
    if(is.null(sFeature)){liMarkerShapes = rep(20,dim(dfOriginal)[2])}
    lLegendInfo[["sShapeMetadata"]] = sFeature
  }
  if(is.null(sFeatureContinuous))
  {
    sFeatureContinuous = funcGetFirstOfDataMode(lsSortedNames =lsSortedMetadataDim2,dfData=dfOriginal,fContinuous=TRUE)
    if(is.null(sFeatureContinuous)){lsMarkerColors = funcGetRandomColors(dim(dfOriginal)[2])}
    lLegendInfo[["sColorMetatadata"]] = sFeatureContinuous
  }

  # Get colors
  if(is.null(lsMarkerColors))
  {
    xData = dfOriginal[[sFeatureContinuous]]

    # Get colors based on the scaled data
    ldUniqueData = unique(xData)
    ldSortedData = ldUniqueData[order(ldUniqueData)]
    lsColors = funcColors(adMax=c(1,.5,.25),adMin=c(0,1,1),iSteps=length(ldSortedData))

    # Update legend information
    lLegendInfo[["lsText"]] = c(lLegendInfo$sText,paste(sFeatureContinuous,"(High)"),paste(sFeatureContinuous,"(Low)"))
    lLegendInfo[["lsColor"]] = c(lLegendInfo$lsColor, lsColors[length(lsColors)], lsColors[1])
    lLegendInfo[["lsMarker"]] = c(lLegendInfo$lsMarker,20,20)

    #Initialize to first color
    lsMarkerColors = rep(lsColors[1],length(xData))
    
    for(i in 1:length(lsColors))
    {
      lsMarkerColors[which(xData==ldSortedData[i])] = lsColors[i]
    }
  }

  # Get shapes
  if(is.null(liMarkerShapes))
  {
    # Get data levels
    xData = dfOriginal[[sFeature]]
    sLevels = levels(xData)
    asLevels = as.character(sLevels)

    # If there are more levels than given in the default then return all points,
    # They will have to use the custom script
    if(length(sLevels) > length(viMarkers)){return(list(markers=rep(20,dim(dfOriginal)[2]), levels=NA, pairs=NA))}

    for(iMarkerIndex in 1:length(asLevels))
    {
      if(!fPlotNA && (is.na(asLevels[iMarkerIndex]) || (tolower(as.character(asLevels[iMarkerIndex]))=="na"))){ next }
      # Update legend information
      lLegendInfo[["lsText"]] = c(lLegendInfo$lsText,paste(asLevels[iMarkerIndex]," (", sFeature,")",sep=""))
      lLegendInfo[["lsColor"]] = c(lLegendInfo$lsColor, "#333333")
      lLegendInfo[["lsMarker"]] = c(lLegendInfo$lsMarker,viMarkers[iMarkerIndex])
    } 

    # Return markers and levels
    liMarkerShapes = rep(20, length(xData))
    for(i in 1:length(sLevels))
    {
      liMarkerShapes[xData==sLevels[i]]=viMarkers[i]
    }
  }

  return(list(shapes=liMarkerShapes, colors=lsMarkerColors, lLegendInfo=lLegendInfo))
}

funcMFA = function(
### Performs multi-factorial analysis
frmeData,
### Data frame of data
dMinSamp,
## Minimum samples
aiMetadata,
### Indicies of metadata to use
aiBugs,
### Indicies of features
aiGenes = c()
### Indicies of genetics
){
  #Update aiBugs with factor and gene data that is being plotted
  #If a custom plotting function is given, then use, otherwise give defaults
  lsFeatures = c()
  if(exists("funcPlotFeatures",mode="function")){lsFeatures = funcPlotFeatures()}
  asMetadata = NA
  if(exists("funcPlotMetadata",mode="function"))
  {
    asMetadata  = funcPlotMetadata()
  } else {
    asMetadata  = c()
  }

  #Get indicies of data to plot and add to those given as parameters
  liFeatures = which(colnames(frmeData) %in% lsFeatures)
  liMetadata = which(colnames(frmeData) %in% asMetadata)
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
      if( length( aiCols ) && ( iNNA < ( dMinSamp * length( aiCols ) ) ) )
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
  ### Multifactorial analysis results
}

#Todo Not done
funcPlotMFA <- function(
### Plot Multiple Factoral Analysis
lsMFA,
### MFA output object
frmeData,
lsMetadata,
lsFeatures,
### Analysis data
iMaxFeatures = 10,
### Make number of bugs and max number of metadata seperately shown
strMFAColorCovariate=NULL,
strMFAShapeCovariate=NULL,
dMFAMetadataScale=NULL,
dMFADataScale=NULL,
lsPlotFeatures=NULL,
fPlotNA=FALSE,
fInvert = FALSE,
### Invert figure to make background black
tempSaveFileName="MFA",
### File to save as a pdf (.pdf extension will be added
funcPlotColors = NULL,
### Function to control plotting colors, read from custom *.R script input file
funcPlotPoints = NULL,
### Function to control plotting point shapes, read from custom *.R script input file
funcPlotLegend = NULL,
### Function to control plotting legend, read from custom *.R script input file
tempWidth=(tempHeight*1.5),
### Width of plot
tempHeight=6,
### Height of plot
tempPCH=20
### Point size
){
  #Set pdf settings
  pdf(paste(tempSaveFileName,".pdf",sep=""), width = c_dHeight * 1.5, height = c_dHeight, useDingbats=FALSE )
  if( fInvert )
  {
    par( bg = "black", fg = "white", col.axis = "white", col.lab = "white", col.main = "white", col.sub = "white" )
  }

  #Get MFA pca data and set dimensions
  lsPCA = lsMFA$global.pca

  # Get markers shapes for plotting (find the highest contributing dim one discontinuous data)
  # Also get colors for plotting (find the highest contributing dim one continuous data)
  llMarkerInfo = funcGetMarkerValues(MFAOut=lsMFA,dfOriginal=frmeData,lsMetadata=lsMetadata,strMFAColorCovariate=strMFAColorCovariate, strMFAShapeCovariate=strMFAShapeCovariate,fPlotNA=fPlotNA)

  #Use default derived colors unless a function has been specified to make the colors. This is over-ridden by setting the covariate to base color on through commandline.
  astrCols = llMarkerInfo$colors
  if(exists("funcPlotColors",mode="function") && is.null(strMFAColorCovariate)){astrCols = funcPlotColors( frmeData )}

  #Use default derived points unless a function has been specified to make the points. This is over-ridden by setting the covariate to base points on through commandline.
  aiPoints =  llMarkerInfo$shapes
  if(exists("funcPlotPoints",mode="function") && is.null(strMFAShapeCovariate)){aiPoints = funcPlotPoints( frmeData )}

  #Get metadata names (in MFA format, could be column names, values or names_values)
  lsMetadataToPlot = funcColToMFAValue(funcGetHighestContribution(MFAOut=lsMFA,dfOriginal=frmeData,iCount=iMaxFeatures,lsNames=lsMetadata),frmeData)
  if(exists("funcPlotMetadata",mode="function")){lsMetadataToPlot = funcPlotMetadata()}
  if(!is.null(lsPlotFeatures))
  {
    lsMetadataToPlot = funcColToMFAValue(intersect(lsPlotFeatures,lsMetadata),frmeData)
  }
  afMetadata <- rownames(lsPCA$var$coord) %in% lsMetadataToPlot

  #Get feature names (in MFA format, could be column names, values or names_values)
  lsFeaturesToPlot = funcColToMFAValue(funcGetHighestContribution(MFAOut=lsMFA,dfOriginal=frmeData,iCount=iMaxFeatures,lsNames=lsFeatures),frmeData)
  if(exists("funcPlotFeatures",mode="function")){lsFeaturesToPlot = funcPlotFeatures()}
  if(!is.null(lsPlotFeatures))
  {
    lsFeaturesToPlot = funcColToMFAValue(intersect(lsPlotFeatures, lsFeatures),frmeData)
  }
  afFeatures = rownames(lsPCA$var$coord) %in% lsFeaturesToPlot

  dScale = dMFAMetadataScale
  if(exists("funcGetMetadataScale",mode="function") && is.null(dScale)){dScale = funcGetMetadataScale()}
  tryCatch(if(is.null(dScale)){dScale = funcGetScale(lsMFA, lsMetadataToPlot)})
  if(is.null(dScale)){dScale = c_dDefaultScale }

  dBugScale = dMFADataScale
  if(exists("funcGetFeatureScale",mode="function") && is.null(dBugScale)){dBugScale = funcGetFeatureScale()}
  if(length(lsFeaturesToPlot) && is.null(dBugScale)){dBugScale = funcGetScale(lsMFA, lsFeaturesToPlot)}

  #Set X and Y coordinates and plot points
  strX <- sprintf( "Dimension 1 (%.2f%%)", lsPCA$eig$`percentage of variance`[1] )
  strY <- sprintf( "Dimension 2 (%.2f%%)", lsPCA$eig$`percentage of variance`[2] )

  #Get coordinates
  lsPointCoordinates = lsPCA$ind$coord

  #Reduce coordinates, markers and color lists to no NA unless specified to plot
  if(!fPlotNA)
  {
    aiNotNA = intersect(which(!is.na(frmeData[[llMarkerInfo$lLegendInfo$sShapeMetadata]])),which(!tolower(frmeData[[llMarkerInfo$lLegendInfo$sShapeMetadata]])=="na"))
    lsPointCoordinates = lsPointCoordinates[aiNotNA,]
    aiPoints = aiPoints[aiNotNA]
    astrCols = astrCols[aiNotNA]
  }

  #Do plots
  #Plot just points 
  funcPlotMFAPage(coordinatesPlot=lsPointCoordinates, coordinatesText=NA, strX=strX, strY=strY, aiPoints=aiPoints, astrCols=astrCols)

  if( sum( afMetadata ) )
  {
    funcPlotMFAPage(coordinatesPlot=lsPointCoordinates, coordinatesText=lsPCA$var$coord, strX=strX, strY=strY, aiPoints=aiPoints,
      astrCols=astrCols, afMetadata=afMetadata, dScale=dScale, lsMetadataLabels=rownames( lsPCA$var$coord )[afMetadata], lLegendLoc=llMarkerInfo$lLegendInfo)
  }

  #Plot features
  if( sum( afFeatures ) )
  {
    funcPlotMFAPage(coordinatesPlot=lsPointCoordinates, coordinatesText=lsPCA$var$coord, strX=strX, strY=strY, aiPoints=aiPoints,
      astrCols=astrCols, afFeatures=afFeatures, dBugScale=dBugScale, lsFeatureLabels=funcRename( rownames( lsPCA$var$coord )[afFeatures] ), lLegendLoc=llMarkerInfo$lLegendInfo)

    #Plot metadata and features
    if( sum( afMetadata ) )
    {
      funcPlotMFAPage(coordinatesPlot=lsPointCoordinates, coordinatesText=lsPCA$var$coord, strX=strX, strY=strY, aiPoints=aiPoints,
        astrCols=astrCols, afMetadata=afMetadata, dScale=dScale, lsMetadataLabels=rownames( lsPCA$var$coord )[afMetadata], afFeatures=afFeatures,
        dBugScale=dBugScale, lsFeatureLabels=funcRename( rownames( lsPCA$var$coord )[afFeatures] ), lLegendLoc=llMarkerInfo$lLegendInfo)
    }
  }
  dev.off( )
}

funcPlotMFAPage <- function(
### Plot a MFA ordination
coordinatesPlot,
### The coordinates of the points from the MFA
coordinatesText,
### The coordinates of the labels from the MFA
strX,
### X axis label
strY,
### Y axis label
aiPoints,
### Shapes to plot points with
astrCols,
### Colors to plot points with
afMetadata=NULL,
### The metadata to plot (indices)
dScale=NULL,
### Scale factor for metadata ordination
lsMetadataLabels=NULL,
### Labels for the metadata to be plotted
afFeatures=NULL,
### Indicies of features to plot
dBugScale=NULL,
### Scale factor for the data feature ordination
lsFeatureLabels=NULL,
### Names of features to plot
lLegendLoc=NULL
### Location to place the legend
){
  plot( coordinatesPlot, pch = aiPoints, col = astrCols, xlab = strX, ylab = strY, cex.axis=1.5, cex.lab=1.5, cex=1.5)
  if( !is.null(afMetadata) && !is.null(dScale) && !is.null(lsMetadataLabels) )
  {
    vCoordinates = coordinatesText[afMetadata,c(1,2)]
    if((length(vCoordinates)==2) && is.null(rownames(vCoordinates))){vCoordinates = data.frame(Dim.1=vCoordinates[1],Dim.2=vCoordinates[2])}
    text( vCoordinates * dScale, labels=lsMetadataLabels, cex=1.8, font = 2 )
  }
  if( !is.null( afFeatures ) && !is.null(dBugScale) && !is.null(lsFeatureLabels) )
  {
    vCoordinates = coordinatesText[afFeatures,c(1,2)]
    if((length(vCoordinates)==2) && is.null(rownames(vCoordinates))){vCoordinates = data.frame(Dim.1=vCoordinates[1],Dim.2=vCoordinates[2])}
    text( vCoordinates * dBugScale, labels=lsFeatureLabels, cex=1.6, font = 3 )
  }
  if(exists("funcPlotLegend",mode="function"))
  {
    funcPlotLegend(sLegendLoc, NULL)
  } else if(!is.null(lLegendLoc)){
    legend(lLegendLoc$sLocation, legend=lLegendLoc$lsText, col=lLegendLoc$lsColor, pch=lLegendLoc$lsMarker, cex=1.5, pt.cex=1.5)
  }
}

funcGetScale <- function(
lsMFA,
lsToPlot
){
  #Scale the metadate labels so they are viewable
  # Max coordinates for samples (dim1,dim2)
  sdDim1Max = max(lsMFA$global.pca$ind$coord[,1])
  sdDim1Min = min(lsMFA$global.pca$ind$coord[,1])
  sdDim2Max = max(lsMFA$global.pca$ind$coord[,2])
  sdDim2Min = min(lsMFA$global.pca$ind$coord[,2])

  # Max x and y coordinates for the select metadata
  # Guard against division by zero
  dMetadataMaxDim1 = max(lsMFA$global.pca$var$coord[lsToPlot,1])
  if(dMetadataMaxDim1==0){dMetadataMaxDim1=0.00000001}
  dMetadataMinDim1 = min(lsMFA$global.pca$var$coord[lsToPlot,1])
  if(dMetadataMinDim1==0){dMetadataMinDim1=0.00000001}
  dMetadataMaxDim2 = max(lsMFA$global.pca$var$coord[lsToPlot,2])
  if(dMetadataMaxDim2==0){dMetadataMaxDim2=0.00000001}
  dMetadataMinDim2 = min(lsMFA$global.pca$var$coord[lsToPlot,2])
  if(dMetadataMinDim2==0){dMetadataMinDim2=0.00000001}

  # Get the scale for metadata
  # See if there are entries less both positive and negative in an axis,
  # You have to use either the pos or neg valuses to determine the scaling depending on
  # The maginitude of their difference from the lowest or highest axis value
  # This is because the MFA plots are not square.
  ldScales = c(sdDim1Max/dMetadataMaxDim1,sdDim2Max/dMetadataMaxDim2,sdDim1Min/dMetadataMinDim1,sdDim2Min/dMetadataMinDim2)
  return(min(abs(ldScales)))
}