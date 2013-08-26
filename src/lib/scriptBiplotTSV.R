#!/usr/bin/env Rscript

library(vegan)
library(optparse)

funcGetCentroidForMetadatum <- function(
### Given a binary metadatum, calculate the centroid of the samples associated with the metadata value of 1
# 1. Get all samples that have the metadata value of 1
# 2. Get the x and y coordinates of the selected samples
# 3. Get the median value for the x and ys
# 4. Return those coordinates as the centroid's X and Y value
vfMetadata,
### Logical or integer (0,1) vector, TRUE or 1 values indicate correspoinding samples in the
### mSamplePoints which will be used to define the centroid
mSamplePoints
### Coordinates (columns;n=2) of samples (rows) corresponding to the vfMetadata
){
  # Check the lengths which should be equal
  if(length(vfMetadata)!=nrow(mSamplePoints))
  {
    print(paste("funcGetCentroidForMetadata::Error: Should have received metadata and samples of the same length, received metadata length ",length(vfMetadata)," and sample ",nrow(mSamplePoints)," length.",sep=""))
    return( FALSE )
  }

  # Get all the samples that have the metadata value of 1
  viMetadataSamples = which(as.integer(vfMetadata)==1)

  # Get the x and y coordinates for the selected samples
  mSelectedPoints = mSamplePoints[viMetadataSamples,]

  # Get the median value for the x and the ys
  if(!is.null(nrow(mSelectedPoints)))
  {
    return( list(x=median(mSelectedPoints[,1],na.rm = TRUE),y=median(mSelectedPoints[,2],na.rm = TRUE)) )
  } else {
    return( list(x=mSelectedPoints[1],y=mSelectedPoints[2]) )
  }
}

funcGetMaximumForMetadatum <- function(
### Given a continuous metadata
### 1. Use the x and ys from mSamplePoints for coordinates and the metadata value as a height (z)
### 2. Use lowess to smooth the landscape
### 3. Take the maximum of the landscape
### 4. Return the coordiantes for the maximum as the centroid
vdMetadata,
### Continuous (numeric or integer) metadata
mSamplePoints
### Coordinates (columns;n=2) of samples (rows) corresponding to the vfMetadata
){
  # Work with data frame
  if(class(mSamplePoints)=="matrix")
  {
    mSamplePoints = data.frame(mSamplePoints)
  }
  # Check the lengths of the dataframes and the metadata
  if(length(vdMetadata)!=nrow(mSamplePoints))
  {
    print(paste("funcGetMaximumForMetadatum::Error: Should have received metadata and samples of the same length, received metadata length ",length(vdMetadata)," and sample ",nrow(mSamplePoints)," length.",sep=""))
    return( FALSE )
  }

  # Add the metadata value to the points
  mSamplePoints[3] = vdMetadata
  names(mSamplePoints) = c("x","y","z") 

  # Create lowess to smooth the surface
  # And calculate the fitted heights
  # x = sample coordinate 1
  # y = sample coordinate 2
  # z = metadata value
  loessSamples = loess(z~x*y, data=mSamplePoints, degree = 1, normalize = FALSE, na.action=na.omit)

  # Naively get the max
  vdCoordinates = loessSamples$x[which(loessSamples$y==max(loessSamples$y)),]
  return(list(lsmod = loessSamples, x=vdCoordinates[1],y=vdCoordinates[2]))
}

funcMakeShapes <- function(
### Takes care of defining shapes for the plot
dfInput,
### Data frame of metadata measurements
sShapeBy,
### The metadata to shape by
sShapes,
### List of custom metadata (per level if factor).
### Should correspond to the number of levels in shapeBy; the format is level:shape,level:shape for example HighLuminosity:14,LowLuminosity:2,HighPH:10,LowPH:18 
cDefaultShape
### Shape to default to if custom shapes are not used
){
  lShapes = list()
  vsShapeValues = c()
  vsShapeShapes = c()
  vsShapes = c()
  sMetadataId = sShapeBy

  # Set default shape, color, and color ranges 
  if(!is.null(cDefaultShape))
  {
    # Default shape should be an int for the int pch options
    if(!is.na(as.integer(cDefaultShape)))
    {
      cDefaultShape = as.integer(cDefaultShape)
    }
  } else {
    cDefaultShape = 16
  }

  # Make shapes
  vsShapes = rep(cDefaultShape,nrow(dfInput))

  if(!is.null(sMetadataId))
  {
    if(is.null(sShapes))
    {
      vsShapeValues = unique(dfInput[[sMetadataId]])
      vsShapeShapes = 1:length(vsShapeValues)
    } else {
      # Put the markers in the order of the values)
      vsShapeBy = unlist(strsplit(sShapes,","))
      for(sShapeBy in vsShapeBy)
      {
        vsShapeByPieces = unlist(strsplit(sShapeBy,":"))
        lShapes[vsShapeByPieces[1]] = as.integer(vsShapeByPieces[2])
      }
      vsShapeValues = names(lShapes)
   }

    # Shapes in the correct order
    if(!is.null(sShapes))
    {
      vsShapeShapes = unlist(lapply(vsShapeValues,function(x) lShapes[[x]]))
    }
    vsShapeValues = paste(vsShapeValues)

    # Make the list of shapes
    for(iShape in 1:length(vsShapeValues))
    {
      vsShapes[which(paste(dfInput[[sMetadataId]])==vsShapeValues[iShape])]=vsShapeShapes[iShape]
    }

    # If they are all numeric characters, make numeric
    viIntNas = which(is.na(as.integer(vsShapes)))
    viNas = which(is.na(vsShapes))
    if(length(setdiff(viIntNas,viNas))==0)
    {
      vsShapes = as.integer(vsShapes)
    } else {
      print("funcMakeShapes::Error: Please supply numbers 1-25 for shape in the -y,--shapeBy option")
      vsShapeValues = c()
      vsShapeShapes = c()
    }
  }
  return(list(PlotShapes=vsShapes,Values=vsShapeValues,Shapes=vsShapeShapes,ID=sMetadataId,DefaultShape=cDefaultShape))
}

### Global defaults
c_sDefaultColorBy = NULL
c_sDefaultColorRange = "orange,cyan"
c_sDefaultTextColor = "black"
c_sDefaultArrowColor = "cyan"
c_sDefaultArrowTextColor = "Blue"
c_sDefaultNAColor = "grey"
c_sDefaultShapeBy = NULL
c_sDefaultShapes = NULL
c_sDefaultMarker = "16"
c_sDefaultRotateByMetadata = NULL
c_sDefaultResizeArrow = 1
c_sDefaultTitle = "Custom Biplot of Bugs and Samples - Metadata Plotted with Centroids"
c_sDefaultOutputFile = NULL

### Create command line argument parser
pArgs <- OptionParser( usage = "%prog last_metadata input.tsv" )

# Selecting features to plot
pArgs <- add_option( pArgs, c("-b", "--bugs"), type="character", action="store", default=NULL, dest="sBugs", metavar="BugsToPlot", help="Comma delimited list of data to plot as text. Bug|1,Bug|2")
pArgs <- add_option( pArgs, c("-m", "--metadata"), type="character", action="store", default=NULL, dest="sMetadata", metavar="MetadataToPlot", help="Comma delimited list of metadata to plot as arrows. metadata1,metadata2,metadata3")

# Colors
pArgs <- add_option( pArgs, c("-c", "--colorBy"), type="character", action="store", default=c_sDefaultColorBy, dest="sColorBy", metavar="MetadataToColorBy", help="The id of the metadatum to use to make the marker colors. Expected to be a continuous metadata.")
pArgs <- add_option( pArgs, c("-r", "--colorRange"), type="character", action="store", default=c_sDefaultColorRange, dest="sColorRange", metavar="ColorRange", help=paste("Colors used to color the samples; a gradient will be formed between the color.Default=", c_sDefaultColorRange))
pArgs <- add_option( pArgs, c("-t", "--textColor"), type="character", action="store", default=c_sDefaultTextColor, dest="sTextColor", metavar="TextColor", help=paste("The color bug features will be plotted with as text. Default =", c_sDefaultTextColor))
pArgs <- add_option( pArgs, c("-a", "--arrowColor"), type="character", action="store", default=c_sDefaultArrowColor, dest="sArrowColor", metavar="ArrowColor", help=paste("The color metadata features will be plotted with as an arrow and text. Default", c_sDefaultArrowColor))
pArgs <- add_option( pArgs, c("-w", "--arrowTextColor"), type="character", action="store", default=c_sDefaultArrowTextColor, dest="sArrowTextColor", metavar="ArrowTextColor", help=paste("The color for the metadata text ploted by the head of the metadata arrow. Default", c_sDefaultArrowTextColor))
pArgs <- add_option(pArgs, c("-n","--plotNAColor"), type="character", action="store", default=c_sDefaultNAColor, dest="sPlotNAColor", metavar="PlotNAColor", help=paste("Plot NA values as this color. Example -n", c_sDefaultNAColor))

# Shapes
pArgs <- add_option( pArgs, c("-y", "--shapeby"), type="character", action="store", default=c_sDefaultShapeBy, dest="sShapeBy", metavar="MetadataToShapeBy", help="The metadata to use to make marker shapes. Expected to be a discrete metadatum. An example would be -y Environment")
pArgs <- add_option( pArgs, c("-s", "--shapes"), type="character", action="store", default=c_sDefaultShapes, dest="sShapes", metavar="ShapesForPlotting", help="This is to be used to specify the shapes to use for plotting. Can use numbers recognized by R as shapes (see pch). Should correspond to the number of levels in shapeBy; the format is level:shape,level:shape for example HighLuminosity:14,LowLuminosity:2,HighPH:10,LowPH:18 . Need to specify -y/--shapeBy for this option to work.")
pArgs <- add_option( pArgs, c("-d", "--defaultMarker"), type="character", action="store", default=c_sDefaultMarker, dest="sDefaultMarker", metavar="DefaultColorMarker", help="Default shape for markers which are not otherwise indicated in --shapes, can be used for unspecified values or NA. Must not be a shape in --shapes.")

# Plot manipulations
pArgs <- add_option( pArgs, c("-e","--rotateByMetadata"), type="character", action="store", default=c_sDefaultRotateByMetadata, dest="sRotateByMetadata", metavar="RotateByMetadata", help="Rotate the ordination by a metadata. Give both the metadata and value to weight it by. The larger the weight, the more the ordination is influenced by the metadata. If the metadata is continuous, use the metadata id; if the metadata is discrete, the ordination will be by one of the levels so use the metadata ID and level seperated by a '_'. Discrete example -e Environment_HighLumninosity,100 ; Continuous example -e Environment,100 .")
pArgs <- add_option( pArgs, c("-z","--resizeArrow"), type="numeric", action="store", default=c_sDefaultResizeArrow, dest="dResizeArrow", metavar="ArrowScaleFactor", help="A constant to multiple the length of the arrow to expand or shorten all arrows together. This will not change the angle of the arrow nor the relative length of arrows to each other.")

# Misc
pArgs <- add_option( pArgs, c("-i", "--title"), type="character", action="store", default=c_sDefaultTitle, dest="sTitle", metavar="Title", help="This is the title text to add to the plot.")
pArgs <- add_option( pArgs, c("-o", "--outputfile"), type="character", action="store", default=c_sDefaultOutputFile, dest="sOutputFileName", metavar="OutputFile", help="This is the name for the output pdf file. If an output file is not given, an output file name is made based on the input file name.")

funcDoBiplot <- function(
### Perform biplot. Samples are markers, bugs are text, and metadata are text with arrows. Markers and bugs are dtermined usiing NMDS and Bray-Curtis dissimilarity. Metadata are placed on the ordination in one of two ways: 1. Factor data - for each level take the ordination points for the samples that have that level and plot the metadata text at the average orindation point. 2. For continuous data - make a landscape (x and y form ordination of the points) and z (height) as the metadata value. Use a lowess line to get the fitted values for z and take the max of the landscape. Plot the metadata text at that smoothed max.
sBugs,
### Comma delimited list of data to plot as text. Bug|1,Bug|2
sMetadata,
### Comma delimited list of metadata to plot as arrows. metadata1,metadata2,metadata3.
sColorBy = c_sDefaultColorBy,
### The id of the metadatum to use to make the marker colors. Expected to be a continuous metadata.
sColorRange = c_sDefaultColorRange,
### Colors used to color the samples; a gradient will be formed between the color. Example orange,cyan
sTextColor = c_sDefaultTextColor,
### The color bug features will be plotted with as text. Example black
sArrowColor = c_sDefaultArrowColor,
### The color metadata features will be plotted with as an arrow and text. Example cyan
sArrowTextColor = c_sDefaultArrowTextColor,
### The color for the metadata text ploted by the head of the metadat arrow. Example Blue
sPlotNAColor = c_sDefaultNAColor,
### Plot NA values as this color. Example grey
sShapeBy = c_sDefaultShapeBy,
### The metadata to use to make marker shapes. Expected to be a discrete metadatum.
sShapes = c_sDefaultShapes,
### This is to be used to specify the shapes to use for plotting. Can use numbers recognized by R as shapes (see pch). Should correspond to the number of levels in shapeBy; the format is level:shape,level:shape for example HighLuminosity:14,LowLuminosity:2,HighPH:10,LowPH:18 .  Works with sShapesBy.
sDefaultMarker = c_sDefaultMarker,
### The default marker shape to use if shapes are not otherwise indicated.
sRotateByMetadata = c_sDefaultRotateByMetadata,
### Metadata and value to rotate by. example Environment_HighLumninosity,100
dResizeArrow = c_sDefaultResizeArrow,
### Scale factor to resize tthe metadata arrows
sTitle = c_sDefaultTitle,
### The title for the figure.
sInputFileName,
### File to input (tsv file: tab seperated, row = sample file)
sLastMetadata,
### Last metadata that seperates data and metadata
sOutputFileName = c_sDefaultOutputFile
### The file name to save the figure.
){
  # Define the colors
  vsColorRange = c("blue","orange")
  if(!is.null(sColorRange))
  {
    vsColorRange = unlist(strsplit(sColorRange,","))
  }
  cDefaultColor = "grey"
  if(!is.null(vsColorRange) && length(vsColorRange)>0)
  {
    cDefaultColor = vsColorRange[1]
  }

  # List of bugs to plot
  # If there is a list it needs to be more than one.
  vsBugsToPlot = c()
  if(!is.null(sBugs))
  {
    vsBugsToPlot = unlist(strsplit(sBugs,","))
  }

  # Metadata to plot
  vsMetadata = c()
  if(!is.null(sMetadata))
  {
    vsMetadata = unlist(strsplit(sMetadata,","))
  }

  ### Load table
  if(class(sInputFileName)=="character")
  {
    dfInput = read.table(sInputFileName, sep = "\t", header=TRUE)
    names(dfInput) = unlist(lapply(names(dfInput),function(x) gsub(".","|",x,fixed=TRUE)))
    row.names(dfInput) = dfInput[,1]
    dfInput = dfInput[-1]
  } else {dfInput = sInputFileName}

  iLastMetadata = which(names(dfInput)==sLastMetadata)
  viMetadata = 1:iLastMetadata
  viData = (iLastMetadata+1):ncol(dfInput)

  ### Dummy the metadata if discontinuous
  ### Leave the continous metadata alone but include
  listMetadata = list()
  vsRowNames = c()
  viContinuousMetadata = c()
  for(i in viMetadata)
  {
    vCurMetadata = unlist(dfInput[i])
    if(is.numeric(vCurMetadata)||is.integer(vCurMetadata))
    {
      vCurMetadata[which(is.na(vCurMetadata))] = mean(vCurMetadata,na.rm=TRUE)
      listMetadata[[length(listMetadata)+1]] = vCurMetadata
      vsRowNames = c(vsRowNames,names(dfInput)[i])
      viContinuousMetadata = c(viContinuousMetadata,length(listMetadata))
    } else {
      vCurMetadata = as.factor(vCurMetadata)
      vsLevels = levels(vCurMetadata)
      for(sLevel in vsLevels)
      { 
        vNewMetadata = rep(0,length(vCurMetadata))
        vNewMetadata[which(vCurMetadata == sLevel)] = 1
        listMetadata[[length(listMetadata)+1]] = vNewMetadata
        vsRowNames = c(vsRowNames,paste(names(dfInput)[i],sLevel,sep="_"))
      }
    }
  }

  # Convert to data frame
  dfDummyMetadata = as.data.frame(sapply(listMetadata,rbind))
  names(dfDummyMetadata) = vsRowNames
  iNumberMetadata = ncol(dfDummyMetadata)

  # Data to use in ordination in NMDS
  dfData = dfInput[viData]

  # If rotating the ordination by a metadata
  # 1. Add in the metadata as a bug
  # 2. Multiply the bug by the weight
  # 3. Push this through the NMDS
  if(!is.null(sRotateByMetadata))
  {
    vsRotateMetadata = unlist(strsplit(sRotateByMetadata,","))
    sMetadata = vsRotateMetadata[1]
    dWeight = as.numeric(vsRotateMetadata[2])
    sOrdinationMetadata = dfDummyMetadata[sMetadata]*dWeight
    dfData[sMetadata] = sOrdinationMetadata
  }

  # Run NMDS on bug data (Default B-C)
  # Will have species and points because working off of raw data
  mNMDSData = metaMDS(dfData,k=2)

  ## Make shapes
  # Defines thes shapes and the metadata they are based on
  # Metadata to use as shapes
  lShapeInfo = funcMakeShapes(dfInput=dfInput, sShapeBy=sShapeBy, sShapes=sShapes, cDefaultShape=sDefaultMarker)

  sMetadataShape = lShapeInfo[["ID"]]
  vsShapeValues = lShapeInfo[["Values"]]
  vsShapeShapes = lShapeInfo[["Shapes"]]
  vsShapes = lShapeInfo[["PlotShapes"]]
  cDefaultShape = lShapeInfo[["DefaultShape"]]

  # Colors
  vsColors = rep(cDefaultColor,nrow(dfInput))
  vsColorValues = c()
  vsColorRBG = c()
  if(!is.null(sColorBy))
  {
    vsColorValues = paste(sort(unique(unlist(dfInput[[sColorBy]])),na.last=TRUE))
    iLengthColorValues = length(vsColorValues)

    vsColorRBG = lapply(1:iLengthColorValues/iLengthColorValues,colorRamp(vsColorRange))
    vsColorRBG = unlist(lapply(vsColorRBG, function(x) rgb(x[1]/255,x[2]/255,x[3]/255)))

    for(iColor in 1:length(vsColorRBG))
    {
      vsColors[which(paste(dfInput[[sColorBy]])==vsColorValues[iColor])]=vsColorRBG[iColor]
    }

    #If NAs are seperately given color, then color here
    if(!is.null(sPlotNAColor))
    {
      vsColors[which(is.na(dfInput[[sColorBy]]))] = sPlotNAColor
      vsColorRBG[which(vsColorValues=="NA")] = sPlotNAColor
    }
  }

  # Reduce the bugs down to the ones in the list to be plotted
  viBugsToPlot = which(row.names(mNMDSData$species) %in% vsBugsToPlot)
  viMetadataDummy = which(names(dfDummyMetadata) %in% vsMetadata)

  # Build the matrix of metadata coordinates
  mMetadataCoordinates = matrix(rep(NA, iNumberMetadata*2),nrow=iNumberMetadata)
  for( i in 1:iNumberMetadata )
  {
    lxReturn = NA
    if( i %in% viContinuousMetadata )
    {
      lxReturn = funcGetMaximumForMetadatum(dfDummyMetadata[[i]],mNMDSData$points)
    } else {
      lxReturn = funcGetCentroidForMetadatum(dfDummyMetadata[[i]],mNMDSData$points)
    }
    mMetadataCoordinates[i,] = c(lxReturn$x,lxReturn$y)
  }
  row.names(mMetadataCoordinates) = vsRowNames

  # Plot the biplot with the centroid constructed metadata coordinates
  if(length(viMetadataDummy)==0)
  {
    viMetadataDummy = 1:nrow(mMetadataCoordinates)
  }

  # Plot samples
  # Make output name
  if(is.null(sOutputFileName))
  {
    viPeriods = which(sInputFileName==".")
    if(length(viPeriods)>0)
    {
      sOutputFileName = paste(OutputFileName[1:viPeriods[length(viPeriods)]],"pdf",sep=".")
    } else {
      sOutputFileName = paste(sInputFileName,"pdf",sep=".")
    }
  }

  pdf(sOutputFileName, useDingbats=FALSE)
  plot(mNMDSData$points, xlab=paste("NMDS1","Stress=",mNMDSData$stress), ylab="NMDS2", pch=vsShapes, col=vsColors)
  title(sTitle,line=3)
  # Plot Bugs
  mPlotBugs = mNMDSData$species[viBugsToPlot,]
  if(length(viBugsToPlot)==1)
  {
    text(x=mPlotBugs[1],y=mPlotBugs[2],labels=row.names(mNMDSData$species)[viBugsToPlot],col=sTextColor)
  } else if(length(viBugsToPlot)>1){
    text(x=mPlotBugs[,1],y=mPlotBugs[,2],labels=row.names(mNMDSData$species)[viBugsToPlot],col=sTextColor)
  }

  # Add alternative axes
  axis(3, col=sArrowColor)
  axis(4, col=sArrowColor)
  box(col = "black")

  # Plot Metadata
  if(length(viMetadataDummy)>0)
  {
    for(i in viMetadataDummy)
    {
      curCoordinates = mMetadataCoordinates[i,]
      curCoordinates = curCoordinates * dResizeArrow
      # Plot Arrow
      arrows(0,0, curCoordinates[1] * 0.8, curCoordinates[2] * 0.8, col=sArrowColor, length=0.1 )
    }
    # Plot text
    if(length(viMetadataDummy)==1)
    {
      text(x=mMetadataCoordinates[viMetadataDummy,][1]*dResizeArrow*0.8, y=mMetadataCoordinates[viMetadataDummy,][2]*dResizeArrow*0.8, labels=row.names(mMetadataCoordinates)[viMetadataDummy],col=sArrowTextColor)
    } else {
      text(x=mMetadataCoordinates[viMetadataDummy,1]*dResizeArrow*0.8, y=mMetadataCoordinates[viMetadataDummy,2]*dResizeArrow*0.8, labels=row.names(mMetadataCoordinates)[viMetadataDummy],col=sArrowTextColor)
    }
  }

  sLegendText = c(paste(vsColorValues,sColorBy,sep="_"),paste(vsShapeValues,sMetadataShape,sep="_"))
  sLegendShapes = c(rep(cDefaultShape,length(vsColorValues)),vsShapeShapes)
  sLegendColors = c(vsColorRBG,rep(cDefaultColor,length(vsShapeValues)))
  if(length(sLegendText)>0)
  {
    legend("topright",legend=sLegendText,pch=sLegendShapes,col=sLegendColors)
  }

  # Original biplot call if you want to check the custom ploting of the script
  # There will be one difference where the biplot call scales an axis, this one does not. In relation to the axes, the points, text and arrows should still match.
  # Axes to the top and right are for the arrow, otherse are for markers and bug names.
  #biplot(mNMDSData$points,mMetadataCoordinates[viMetadataDummy,],xlabs=vsShapes,xlab=paste("MDS1","Stress=",mNMDSData$stress),main="Biplot function Bugs and Sampes - Metadata Plotted with Centroids")
  dev.off()
}

# This is the equivalent of __name__ == "__main__" in Python.
# That is, if it's true we're being called as a command line script;
# if it's false, we're being sourced or otherwise included, such as for
# library or inlinedocs.
if( identical( environment( ), globalenv( ) ) &&
	!length( grep( "^source\\(", sys.calls( ) ) ) )
{
  lsArgs <- parse_args( pArgs, positional_arguments=TRUE )

  funcDoBiplot(
    sBugs = lsArgs$options$sBugs,
    sMetadata = lsArgs$options$sMetadata,
    sColorBy = lsArgs$options$sColorBy,
    sColorRange = lsArgs$options$sColorRange,
    sTextColor = lsArgs$options$sTextColor,
    sArrowColor = lsArgs$options$sArrowColor,
    sArrowTextColor = lsArgs$options$sArrowTextColor,
    sPlotNAColor = lsArgs$options$sPlotNAColor,
    sShapeBy = lsArgs$options$sShapeBy,
    sShapes = lsArgs$options$sShapes,
    sDefaultMarker = lsArgs$options$sDefaultMarker,
    sRotateByMetadata = lsArgs$options$sRotateByMetadata,
    dResizeArrow = lsArgs$options$dResizeArrow,
    sTitle = lsArgs$options$sTitle,
    sInputFileName = lsArgs$args[2],
    sLastMetadata = lsArgs$args[1],
    sOutputFileName = lsArgs$options$sOutputFileName)
}
