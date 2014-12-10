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
##description<< Collection of functions centered on custom reading of data and some IO services.
) { return( pArgs ) }

#Project Constants

c_astrNA <- c(""," ","  ","NA","na")

#Do not report warnings
options(warn=-1)

funcWriteMatrixToReadConfigFile = function(
### Writes a read config file. Will write over a file by default
strConfigureFileName,
### Matrix that will be read
strMatrixName,
### Name of matrix that will be read
strRowIndices=NA,
### Rows which will be Read (TSV) by default all will be read
strColIndices=NA,
### Cols which will be Read (TSV) by default all will be read
acharDelimiter=c_strDefaultMatrixDelimiter,
### Delimiter for the matrix that will be read in\
fAppend=FALSE
### Append to a current read config file
){
  #If no append delete previous file
  if(!fAppend){unlink(strConfigureFileName)}

  #Make delimiter readable
  switch(acharDelimiter,
    "\t" = {acharDelimiter = "TAB"},
    " " = {acharDelimiter = "SPACE"},
    "\r" = {acharDelimiter = "RETURN"},
    "\n" = {acharDelimiter = "ENDLINE"})
    
  #Manage NAs
  if(is.na(strRowIndices)){strRowIndices="-"}
  if(is.na(strColIndices)){strColIndices="-"}

  #Required output
  lsDataLines = c(paste(c_MATRIX_NAME,strMatrixName,sep=" "),
    paste(c_DELIMITER,acharDelimiter,sep=" "),
    paste(c_ID_ROW,"1",sep=" "),
    paste(c_ID_COLUMN,"1",sep=" "),
    paste(c_TSVROWS,strRowIndices,sep=" "),
    paste(c_TSVCOLUMNS,strColIndices,sep=" "))

  lsDataLines = c(lsDataLines,"\n")

  #Output to file
  lapply(lsDataLines, cat, file=strConfigureFileName, sep="\n", append=TRUE)
}

funcWriteMatrices = function(
### Write data frame data files with config files
dataFrameList,
### A named list of data frames (what you get directly from the read function)
saveFileList,
### File names to save the data matrices in (one name per data frame)
configureFileName,
### Name of the configure file to be written which will direct the reading of these data
acharDelimiter=c_strDefaultMatrixDelimiter,
### Matrix delimiter
log = FALSE
### Indicates if logging should occur
){
  #Get names
  dataFrameNames = names(dataFrameList)

  #Get length of dataFrameList
  dataFrameListLength = length(dataFrameList)

  #Get length of save file list
  saveFileListLength = length(saveFileList)

  #If the save file list length and data frame list length are not equal, abort
  if(!saveFileListLength == dataFrameListLength)
  {stop(paste("Received a length of save files (",saveFileListLength,") that are different from the count of data frames (",dataFrameListLength,"). Stopped and returned false."),sep="")}

  #Delete the old config file
  unlink(configureFileName)

  #For each data save
  for (dataIndex in c(1:dataFrameListLength))
  {
    #Current data frame
    data = dataFrameList[[dataIndex]]

    #Get column count
    columnCount = ncol(data)

    #Get row and column names
    rowNames = row.names(data)
    rowNamesString = paste(rowNames,sep="",collapse=",")
    if(length(rowNamesString)==0){rowNamesString = NA}

    columnNamesString = paste(colnames(data),sep="",collapse=",")
    if(length(columnNamesString)==0){columnNamesString = NA}

    #Get row indices
    rowStart = 1
    if(!is.na(rowNamesString)){rowStart = 2}
    rowEnd = nrow(data)+rowStart - 1
    rowIndices = paste(c(rowStart:rowEnd),sep="",collapse=",")

    #Get col indices
    colStart = 1
    if(!is.na(columnNamesString)){  colStart = 2}
    colEnd = columnCount+colStart - 1
    colIndices = paste(c(colStart:colEnd),sep="",collapse=",")

    #Write Data to file
    write.table(data, saveFileList[dataIndex], quote = FALSE, sep = acharDelimiter, col.names = NA, row.names = rowNames, na = "NA", append = FALSE)

    #Write the read config file
    funcWriteMatrixToReadConfigFile(strConfigureFileName=configureFileName, strMatrixName=dataFrameNames[dataIndex],
      strRowIndices=rowIndices, strColIndices=colIndices, acharDelimiter=acharDelimiter, fAppend=TRUE)
  }
  return(TRUE)
}

funcReadMatrices = function(
### Dynamically Read a Matrix/Matrices from a configure file
configureFile,
### Read config file to guide reading in data
defaultFile = NA,
### Default data file to read
log = FALSE
){
  #Named vector to return data frames read
  returnFrames = list()
  #Holds the names of the frames as they are being added
  returnFrameNames = c()
  returnFramesIndex = 1

  #Read in config file info
  #Read each data block extracted from the config file
  lsDataBlocks <- funcReadConfigFile(configureFile, defaultFile)
  if(!length(lsDataBlocks)) {
	  astrMetadata <- NULL
	  astrMetadata[2] <- defaultFile
	  astrMetadata[5] <- "2"
	  astrData <- NULL
	  astrData[2] <- defaultFile
	  astrData[5] <- "3-"
	  lsDataBlocks <- list(astrMetadata, astrData)
  }
  for(dataBlock in lsDataBlocks)
  {
    #Read in matrix
    returnFrames[[returnFramesIndex]] = funcReadMatrix(tempMatrixName=dataBlock[1], tempFileName=dataBlock[2], tempDelimiter=dataBlock[3], tempColumns=dataBlock[5], tempRows=dataBlock[4], tempLog=log)
    returnFrameNames = c(returnFrameNames,dataBlock[1])
    returnFramesIndex = returnFramesIndex + 1
  }
  names(returnFrames) = returnFrameNames
  return(returnFrames)
}

funcReadMatrix = function(
### Read one matrix
### The name to give the block of data read in from file
tempMatrixName,
### ID rows and columns are assumed to be 1
tempFileName=NA,
### Data file to read
tempDelimiter=NA,
### Data matrix delimiter
tempColumns=NA,
### Data columns to read
tempRows=NA,
### Data rows to read
tempLog=FALSE
### Indicator to log
){
  if(is.na(tempDelimiter)){tempDelimiter <- c_strDefaultMatrixDelimiter}
  if(is.na(tempColumns)){tempColumns <- c_strDefaultReadCols}
  if(is.na(tempRows)){tempRows <- c_strDefaultReadRows}
  
  #Check parameter and make sure not NA
  if(is.na(tempMatrixName)){tempMatrixName <- ""}
  if(!funcIsValid(tempMatrixName)){stop(paste("Did not receive a valid matrix name, received ",tempMatrixName,"."))}

  #Check to make sure there is a file name for the matrix
  if(! funcIsValidFileName(tempFileName))
  {stop(paste("No valid file name is given for the matrix ",tempMatrixName," from file: ",tempFileName,". Please add a valid file name to read the matrix from.", sep=""))}

  #Read in superset matrix and give names if indicated 
  #Read in matrix
  dataMatrix = read.table(tempFileName, sep = tempDelimiter, as.is = TRUE, na.strings=c_astrNA, quote = "", comment.char = "")
  dataFrameDimension = dim(dataMatrix)

  #Remove special characters from column names
  tempColumns=gsub("[^a-zA-Z0-9_|,-]|@|\\?|\\]|\\[|\\^","_",tempColumns)
  for(i in seq_along(dataMatrix[1,]))
  {
  	#First remove all special characters, then remove ",-" as these are
        #part of the read.config organization
	dataMatrix[1,i]=gsub("[^a-zA-Z0-9_|,-]|@|\\?|\\]|\\[|\\^","_",dataMatrix[1,i])
	newColumnName=gsub("-|,","_",dataMatrix[1,i])
	if(newColumnName!=dataMatrix[1,i])
	{
    	tempColumns=gsub(dataMatrix[1,i],newColumnName,tempColumns)
    	dataMatrix[1,i]=newColumnName
    }
  }
  
  #Get column names
  columnNameList = as.matrix(dataMatrix[1,])
  rowNameList = dataMatrix[1][[1]]
  
  #Convert characters to vectors of indices
  tempColumns = funcParseIndexSlices(ifelse(is.na(tempColumns),"-",tempColumns), columnNameList)
  tempRows = funcParseIndexSlices(ifelse(is.na(tempRows),"-", tempRows), rowNameList)

  #Check indices
  #Check to make sure valid id col/rows and data col/rows
  if((!funcIsValid(tempColumns)) || (!funcIsValid(tempRows)))
  {stop(paste("Received invalid row or col. Rows=",tempRows," Cols=", tempColumns))}

  #Check to make sure only 1 row id is given and it is not repeated in the data rows
  if(length(intersect(1,tempColumns)) == 1)
  {stop(paste("Index indicated as an id row but was found in the data row indices, can not be both. Index=1 Data indices=",tempColumns,sep=""))}

  #Check to make sure only one col id is given and it is not repeated in the data columns
  #Id row/col should not be in data row/col
  if(length(intersect(1, tempRows)) == 1)
  {stop(paste("Index indicated as an id column but was found in the data column indices, can not be both. ID Index=1 Data Indices=", tempRows,".",sep=""))}

  #If the row names have the same length as the column count and has column names 
  #it is assumed that the tempIdCol index item is associated with the column names.
  #Visa versa for rows, either way it is removed
  #Remove ids from name vector
  rowNameList = rowNameList[(-1)]
  #Remove ids from data
  dataMatrix = dataMatrix[(-1)]
  #Adjust row ids given the removal of the id row
  tempColumns=(tempColumns-1)

  ## Remove id rows/columns and set row/col names
  #Remove ids from vector
  columnNameList = columnNameList[(-1)]
  #Remove ids from data
  dataMatrix = dataMatrix[(-1),]
  #Adjust column ids given the removal of the id column
  tempRows =(tempRows-1)
  #Add row and column names
  row.names(dataMatrix) = as.character(rowNameList)
  colnames(dataMatrix) = as.character(columnNameList)

  #Reduce matrix
  #Account for when both column ranges and row ranges are given or just a column or a row range is given
  dataMatrix = dataMatrix[tempRows, tempColumns, drop=FALSE]

  #Set all columns data types to R guessed default
  for(i in 1:ncol(dataMatrix)){
    dataMatrix[,i] <- type.convert(dataMatrix[,i], na.strings = c_astrNA)}

  #Return matrix
  return(dataMatrix)
}

funcReadConfigFile = function(
### Reads in configure file and extracts the pieces needed for reading in a matrix
configureFile,
### Configure file = string path to configure file
defaultFile = NA
### Used to set a default data file
){
  #Read configure file
  fileDataList <- list()
  if(!is.null( configureFile ) ) {
    fileDataList <- scan( file = configureFile, what = character(), sep="\n", quiet=TRUE) }
  newList = list()
  for(sLine in fileDataList)
  {
    sLine = gsub("\\s","",sLine)
    vUnits = unlist(strsplit(sLine,":"))
    if(length(vUnits)>1)
    {
      vUnits[1] = paste(vUnits[1],":",sep="")
      newList[[length(newList)+1]] = vUnits
    }
  }
  fileDataList = unlist(newList)

  matrixName <- NA
  fileName <- defaultFile

  #Hold information on matrices to be read
  matrixInformationList = list()
  matrixInformationListCount = 1

  for(textIndex in c(1:length(fileDataList)))
  {
    if(textIndex > length(fileDataList)) {break}
    #Start at the Matrix name
    #Keep this if statement first so that you scan through until you find a matrix block
    if(fileDataList[textIndex] == c_MATRIX_NAME)
    {
      #If the file name is not NA then that is sufficient for a matrix, store
      #Either way reset
      if(funcIsValid(fileName)&&funcIsValid(matrixName))
      {
        matrixInformationList[[matrixInformationListCount]] = c(matrixName,fileName,delimiter,rows,columns)
        matrixInformationListCount = matrixInformationListCount + 1
      }

      #Get the matrix name and store
      matrixName = fileDataList[textIndex + 1]

      fileName = defaultFile
      delimiter = "\t"
      rows = NA
      columns = NA
      #If is not matrix name and no matrix name is known skip until you find the matrix name
      #If matrix name is known, continue to collect information about that matrix
    } else if(is.na(matrixName)){next}

    #Parse different keywords
    strParseKey = fileDataList[textIndex]
    if(strParseKey == c_FILE_NAME){fileName=fileDataList[textIndex+1]}
    else if(strParseKey==c_FILE_NAME){fileName=fileDataList[textIndex+1]}
    else if(strParseKey %in% c(c_TSVROWS,c_PCLCOLUMNS,c_ROWS)){rows=fileDataList[textIndex+1]}
    else if(strParseKey %in% c(c_TSVCOLUMNS,c_PCLROWS,c_COLUMNS)){columns=fileDataList[textIndex+1]}
    else if(strParseKey==c_DELIMITER)
    {
        switch(fileDataList[textIndex + 1],
        "TAB" = {delimiter = "\t"},
        "SPACE" = {delimiter = " "},
        "RETURN" = {delimiter = "\r"},
        "ENDLINE" = {delimiter = "\n"})
    }
  }
  #If there is matrix information left
  if((!is.na(matrixName)) && (!is.na(fileName)))
  {
    matrixInformationList[[matrixInformationListCount]] = c(matrixName,fileName,delimiter,rows,columns)
    matrixInformationListCount = matrixInformationListCount + 1
  }

  return(matrixInformationList)
}

funcParseIndexSlices = function(
### Take a string of comma or dash seperated integer strings and convert into a vector
### of integers to use in index slicing
strIndexString,
### String to be parsed into indicies vector
cstrNames
### Column names of the data so names can be resolved to indicies
){
  #If the slices are NA then return
  if(is.na(strIndexString)){return(strIndexString)}

  #List of indices to return
  viRetIndicies = c()

  #Split on commas
  lIndexString = sapply(strsplit(strIndexString, c_COMMA),function(x) return(x))
  for(strIndexItem in lIndexString)
  {
    #Handle the - case
    if(strIndexItem=="-"){strIndexItem = paste("2-",length(cstrNames),sep="")}

    #Split on dash and make sure it makes sense
    lItemElement = strsplit(strIndexItem, c_DASH)[[1]]
    if(length(lItemElement)>2){stop("Error in index, too many dashes, only one is allowed. Index = ",strIndexItem,sep="")}

    #Switch names to numbers
    aiIndices = which(is.na(as.numeric(lItemElement)))
    for( iIndex in aiIndices )
    {
      lItemElement[iIndex] = which(cstrNames==lItemElement[iIndex])[1]
    }

    #Make numeric
    liItemElement = unlist(lapply(lItemElement, as.numeric))

    #If dash is at the end or the beginning add on the correct number
    if(substr(strIndexItem,1,1)==c_DASH){liItemElement[1]=2}
    if(substr(strIndexItem,nchar(strIndexItem),nchar(strIndexItem))==c_DASH){liItemElement[2]=length(cstrNames)}

    #If multiple numbers turn to a slice
    if(length(liItemElement)==2){liItemElement = c(liItemElement[1]:liItemElement[2])}

    #Update indices
    viRetIndicies = c(viRetIndicies, liItemElement)
  }
  if(length(viRetIndicies)==0){return(NA)}
  return(sort(unique(viRetIndicies)))
  ### Sorted indicies vector
}
