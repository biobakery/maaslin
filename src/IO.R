####################################
# Summary: IO
# Author: Timothy Tickle
# Start Date: 11-01-2011
####################################

#Project Constants

c_astrNA <- c(""," ","  ","NA","na")

#Do not report warnings
options(warn=-1)

#4 Test cases
### Writes a read config file. Will write over a file by default
### strConfigureFileName Read Config file to write to 
### strMatrixFile File that will be read
### strMatrixName Name of matrix that will be read
### strRowIndices Rows which will be Read (TSV) by default all will be read
### strColIndices Cols which will be Read (TSV) by default all will be read
### strDtCharacter Data columns which will be forced to character data
### strDtFactoral Data columns which will be forced to factor data
### strDtInteger Data columns which will be forced to integer data
### strDtLogical Data columns which will be forced to logical data
### strDtNumeric Data columns which will be forced to numeric data
### strDtOrdered Data columns which will be forced to ordered data
### acharDelimiter Delimiter for the matrix that will be read in\
### fAppend Append to a current read config file
funcWriteMatrixToReadConfigFile = function(strConfigureFileName, strMatrixFile, strMatrixName, strRowIndices=NA, strColIndices=NA,
  strDtCharacter=NA, strDtFactoral=NA, strDtInteger=NA, strDtLogical=NA, strDtNumeric=NA, strDtOrdered=NA, acharDelimiter=c_strDefaultMatrixDelimiter, fAppend=FALSE)
{
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

  #Optional output
  if(!is.na(strMatrixFile) && !is.null(strMatrixFile)){lsDataLines=c(lsDataLines,paste(c_FILE_NAME,strMatrixFile,sep=" "))}
  if(!is.na(strDtCharacter) && !is.null(strDtCharacter)){lsDataLines=c(lsDataLines,paste(c_CHARACTER_DATA_TYPE,strDtCharacter,sep=" "))}
  if(!is.na(strDtFactoral) && !is.null(strDtFactoral)){lsDataLines=c(lsDataLines,paste(c_FACTOR_DATA_TYPE,strDtFactoral,sep=" "))}
  if(!is.na(strDtInteger) && !is.null(strDtInteger)){lsDataLines=c(lsDataLines,paste(c_INTEGER_DATA_TYPE,strDtInteger,sep=" "))}
  if(!is.na(strDtLogical) && !is.null(strDtLogical)){lsDataLines=c(lsDataLines,paste(c_LOGICAL_DATA_TYPE,strDtLogical,sep=" "))}
  if(!is.na(strDtNumeric) && !is.null(strDtNumeric)){lsDataLines=c(lsDataLines,paste(c_NUMERIC_DATA_TYPE,strDtNumeric,sep=" "))}
  if(!is.na(strDtOrdered) && !is.null(strDtOrdered)){lsDataLines=c(lsDataLines,paste(c_ORDEREDFACTOR_DATA_TYPE,strDtOrdered,sep=" "))}
  lsDataLines = c(lsDataLines,"\n")

  #Output to file
  lapply(lsDataLines, cat, file=strConfigureFileName, sep="\n", append=TRUE)
}

#Write data frame data files with config files
#dataFrameList is a named list of data frames (what you get directly from the read function)
#saveFileList File names to save the data matrices in (one name per data frame)
#configureFileName Name of the configure file to be written which will direct the reading of these data
funcWriteMatrices = function(dataFrameList, saveFileList, configureFileName, acharDelimiter=c_strDefaultMatrixDelimiter, log = FALSE)
{
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
    colEnd = columnCount+ colStart - 1
    colIndices = paste(c(colStart:colEnd),sep="",collapse=",")

    #Data modes
    dtCharacter = c()
    dtFactoral = c()
    dtInteger = c()
    dtLogical = c()
    dtNumeric = c()
    dtOrdered = c()

    #Move through data frame and collect data
    for(colIndex in c(1:columnCount))
    {
      currentCol = data[[colIndex]]

      if(is.character(currentCol))
      {
        dtCharacter = c(dtCharacter,colIndex)
        next
      }
      if(is.integer(currentCol))
      {
        dtInteger = c(dtInteger,colIndex)
        next
      }
      if(is.logical(currentCol))
      {
        dtLogical = c(dtLogical,colIndex)
        next
      }
      if(is.numeric(currentCol))
      {
        dtNumeric = c(dtNumeric,colIndex)
        next
      }
      if(is.ordered(currentCol))
      {
        dtOrdered = c(dtOrdered,colIndex)
        next
      }
      if(is.factor(currentCol))
      {
        dtFactoral = c(dtFactoral,colIndex)
        next
      }
    }

    #Adjust for column of row names and turn into comma delimited list
    if(length(dtCharacter)>0)
    {dtCharacter=paste(dtCharacter+1,sep="",collapse=",")
    } else {dtCharacter=""}
    
    if(length(dtFactoral)>0)
    {dtFactoral=paste(dtFactoral+1,sep="",collapse=",")
    } else {dtFactoral=""}
    
    if(length(dtInteger)>0)
    {dtInteger=paste(dtInteger+1,sep="",collapse=",")
    } else {dtInteger=""}
    
    if(length(dtLogical)>0)
    {dtLogical=paste(dtLogical+1,sep="",collapse=",")
    } else {dtLogical=""}
    
    if(length(dtNumeric)>0)
    {dtNumeric=paste(dtNumeric+1,sep="",collapse=",")
    } else {dtNumeric=""}
    
    if(length(dtOrdered)>0)
    {dtOrdered=paste(dtOrdered+1,sep="",collapse=",")
    } else {dtOrdered=""}

    #Write Data to file
    write.table(data, saveFileList[dataIndex], quote = FALSE, sep = acharDelimiter, col.names = NA, row.names = rowNames, na = "NA", append = FALSE)

    #Write the read config file
    funcWriteMatrixToReadConfigFile(strConfigureFileName=configureFileName, strMatrixFile=saveFileList[dataIndex], strMatrixName=dataFrameNames[dataIndex],
      strRowIndices=rowIndices, strColIndices=colIndices, strDtCharacter=dtCharacter, strDtFactoral=dtFactoral, strDtInteger=dtInteger,
      strDtLogical=dtLogical, strDtNumeric=dtNumeric, strDtOrdered=dtOrdered, acharDelimiter=acharDelimiter, fAppend=TRUE)
  }
  return(TRUE)
}

#Dynamically Read a Matrix/Matrices from a configure file
#TODO If there is no row or column values then read the full matrix
funcReadMatrices = function( configureFile , defaultFile = NA, log = FALSE)
{
  #Named vector to return data frames read
  returnFrames = list()
  #Holds the names of the frames as they are being added
  returnFrameNames = c()
  returnFramesIndex = 1

  #Read in config file info
  #Read each data block extracted from the config file
  for(dataBlock in funcReadConfigFile(configureFile, defaultFile))
  {
    #Read in matrix
    returnFrames[[returnFramesIndex]] = funcReadMatrix(tempMatrixName=dataBlock[1], tempFileName=dataBlock[2], tempDelimiter=dataBlock[3], tempColumns=dataBlock[5], tempRows=dataBlock[4], tempDtCharacter=dataBlock[6], tempDtFactor=dataBlock[7], tempDtInteger=dataBlock[8], tempDtLogical=dataBlock[9], tempDtNumeric=dataBlock[10], tempDtOrderedFactor=dataBlock[11], tempLog=log)
    returnFrameNames = c(returnFrameNames,dataBlock[1])
    returnFramesIndex = returnFramesIndex + 1
  }
  names(returnFrames) = returnFrameNames
  return(returnFrames)
}

#Read one matrix
#ID rows and columns are assumed to be 1
funcReadMatrix = function(tempMatrixName, tempFileName, tempDelimiter=c_strDefaultMatrixDelimiter, tempColumns=c_strDefaultReadCols, tempRows=c_strDefaultReadRows, tempDtCharacter=NA, tempDtFactor=NA, tempDtInteger=NA, tempDtLogical=NA, tempDtNumeric=NA, tempDtOrderedFactor=NA, tempLog=FALSE)
{
  #Check parameter and make sure not NA
  if(!funcIsValid(tempMatrixName)){stop(paste("Did not receive a valid matrix name, received ",tempMatrixName,"."))}

  #Check to make sure there is a file name for the matrix
  if(! funcIsValidFileName(tempFileName))
  {stop(paste("No valid file name is given for the matrix ",tempMatrixName,". Please add a valid file name to read the matrix from.", sep=""))}

  #Read in superset matrix and give names if indicated 
  #Read in matrix
  dataMatrix = read.table(tempFileName, sep = tempDelimiter, as.is = TRUE, na.strings=c_astrNA, quote = "", comment.char = "")
  dataFrameDimension = dim(dataMatrix)

  #Get column names
  columnNameList = columnNameList = as.matrix(dataMatrix[1,])
  rowNameList = dataMatrix[1][[1]]

  #Convert characters to vectors of indices
  tempColumns = funcParseIndexSlices(ifelse(is.na(tempColumns),"-",tempColumns), columnNameList)
  tempRows = funcParseIndexSlices(ifelse(is.na(tempRows),"-", tempRows), rowNameList)
  
  tempDtCharacter = funcParseIndexSlices(tempDtCharacter, columnNameList)
  tempDtFactor = funcParseIndexSlices(tempDtFactor, columnNameList)
  tempDtOrderedFactor = funcParseIndexSlices(tempDtOrderedFactor, columnNameList)
  tempDtInteger = funcParseIndexSlices(tempDtInteger, columnNameList)
  tempDtLogical = funcParseIndexSlices(tempDtLogical, columnNameList)
  tempDtNumeric = funcParseIndexSlices(tempDtNumeric, columnNameList)

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
  if(funcIsValid(tempDtCharacter)){ tempDtCharacter=(tempDtCharacter-1)}
  if(funcIsValid(tempDtFactor)){ tempDtFactor=(tempDtFactor-1)}
  if(funcIsValid(tempDtInteger)){ tempDtInteger=(tempDtInteger-1)}
  if(funcIsValid(tempDtLogical)){ tempDtLogical=(tempDtLogical-1)}
  if(funcIsValid(tempDtNumeric)){ tempDtNumeric=(tempDtNumeric-1)}
  if(funcIsValid(tempDtOrderedFactor)){ tempDtOrderedFactor=(tempDtOrderedFactor-1)}

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

  #Set all columns data types to R guessed default
  for(i in 1:ncol(dataMatrix)){
    dataMatrix[,i] <- type.convert(dataMatrix[,i], na.strings = c_astrNA)}

  #Check to make sure there are no duplicates in data mode indices (remove NAs before checking uniqueness, we done care about duplicate Nas)
  if((funcIsValid(tempDtNumeric)) || (funcIsValid(tempDtFactor)) || (funcIsValid(tempDtCharacter)) || (funcIsValid(tempDtInteger)) || (funcIsValid(tempDtLogical)) || (funcIsValid(tempDtOrderedFactor)))
  {
    allModes = c(tempDtNumeric,tempDtFactor,tempDtCharacter,tempDtInteger,tempDtLogical,tempDtOrderedFactor)
    allModes = allModes[!is.na(allModes)]

    if(anyDuplicated(allModes))
    {
      stop("Data types (modes) indices appeared in more than one category. Make sure the indices in all combined data types (mode) categories are unique.")
    }
    
    #Set modes of data
    if (funcIsValid(tempDtNumeric))
    {
      #Set numerics
      for(numericIndex in tempDtNumeric)
      {
        tryCatch({
         dataMatrix[numericIndex] = as.numeric(as.character(dataMatrix[[numericIndex]]))
        }, error = function(ex) {
          stop(paste("IO::Index ",numericIndex," could not be converted to a numeric vector.",dataMatrix[numericIndex], sep=""))}, finally={})
      }
    }
    if (funcIsValid(tempDtFactor))
    {
      #Set factors
      for(factorIndex in tempDtFactor)
      {
        tryCatch({
        dataMatrix[factorIndex] = as.factor(as.character(dataMatrix[[factorIndex]]))
        }, error = function(ex) {
          stop(paste("IO::Index ",factorIndex," could not be converted to a factor vector.",dataMatrix[factorIndex], sep=""))}, finally={})
      }
    }
    if (funcIsValid(tempDtInteger))
    {
      #Set integers
      for(integerIndex in tempDtInteger)
      {
        tryCatch({
        dataMatrix[integerIndex] = as.integer(as.character(dataMatrix[[integerIndex]]))
        }, error = function(ex) {
          stop(paste("IO::Index ",integerIndex," could not be converted to an integer vector.",dataMatrix[integerIndex], sep=""))}, finally={})
      }
    }
    if (funcIsValid(tempDtLogical))
    {
      #Set logical
      for(logicalIndex in tempDtLogical)
      {
        tryCatch({
        dataMatrix[logicalIndex] = sapply(dataMatrix[[logicalIndex]],
          function(x){ifelse(is.logical(x), as.logical(as.character(x)), as.logical(as.numeric(as.character(x))))})
        }, error = function(ex) {
          stop(paste("IO::Index ",logicalIndex," could not be converted to a logical vector.",dataMatrix[logicalIndex], sep=""))}, finally={})
      }
    }
    if (funcIsValid(tempDtOrderedFactor))
    {
      #Set orderedFactor
      for(orderedIndex in tempDtOrderedFactor)
      {
        tryCatch({
        dataMatrix[orderedIndex] = as.ordered(as.character(dataMatrix[[orderedIndex]]))
        }, error = function(ex) {
          stop(paste("IO::Index ",orderedIndex," could not be converted to an ordered factor vector.",dataMatrix[orderedIndex], sep=""))}, finally={})
      }
    }
    if (funcIsValid(tempDtCharacter))
    {
      #Set Character
      for(characterIndex in tempDtCharacter)
      {
        tryCatch({
        dataMatrix[characterIndex] = as.character(dataMatrix[[characterIndex]])
        }, warning={}, error = function(ex) {
          stop(paste("IO::Index ",characterIndex," could not be converted to a character vector.",dataMatrix[characterIndex], sep=""))}, finally={})
      }
    }
  }

  #Reduce matrix
  #Account for when both column ranges and row ranges are given or just a column or a row range is given
  dataMatrix = dataMatrix[tempRows, tempColumns, drop=FALSE]

#  if(tempLog){funcLogMatrixRead()}
  
  #Return matrix
  return(dataMatrix)
}

#2 Testcases
#Reads in configure file and extracts the pieces needed for reading in a matrix
#Configure file = string path to configure file
#TODO Make sure that commented is commented out
funcReadConfigFile = function(configureFile, defaultFile = NA)
{
  #Read configure file
  fileDataList <- scan( file = configureFile, what = character(), quiet=TRUE)
  matrixName <- NA
  fileName <- defaultFile

  #Hold information on matrices to be read
  matrixInformationList = list()
  matrixInformationListCount = 1

  for(textIndex in c(1:length(fileDataList)))
  {
    #Start at the Matrix name
    #Keep this if statement first so that you scan through until you find a matrix block
    if(fileDataList[textIndex] == c_MATRIX_NAME)
    {
      #If the file name is not NA then that is sufficient for a matrix, store
      #Either way reset
      if(funcIsValid(fileName)&&funcIsValid(matrixName))
      {
        matrixInformationList[[matrixInformationListCount]] = c(matrixName,fileName,delimiter,rows,columns,dtCharacter,dtFactor,dtInteger,dtLogical,dtNumeric,dtOrderedFactor)
        matrixInformationListCount = matrixInformationListCount + 1
      }

      #Get the matrix name and store
      matrixName = fileDataList[textIndex + 1]

      fileName = defaultFile
      delimiter = "\t"
      rows = NA
      columns = NA
      dtCharacter = NA
      dtFactor = NA
      dtInteger = NA
      dtLogical = NA
      dtNumeric = NA
      dtOrderedFactor = NA 
      #If is not matrix name and no matrix name is known skip until you find the matrix name
      #If matrix name is known, continue to collect information about that matrix
    } else if(is.na(matrixName)){next}

    #Parse different keywords
    strParseKey = fileDataList[textIndex]
    if(strParseKey == c_FILE_NAME){fileName=fileDataList[textIndex+1]}
    else if(strParseKey==c_FILE_NAME){fileName=fileDataList[textIndex+1]}
    else if(strParseKey %in% c(c_TSVROWS,c_PCLCOLUMNS,c_ROWS)){rows=fileDataList[textIndex+1]}
    else if(strParseKey %in% c(c_TSVCOLUMNS,c_PCLROWS,c_COLUMNS)){columns=fileDataList[textIndex+1]}
    else if(strParseKey==c_CHARACTER_DATA_TYPE){dtCharacter=fileDataList[textIndex+1]}
    else if(strParseKey==c_FACTOR_DATA_TYPE){dtFactor=fileDataList[textIndex+1]}
    else if(strParseKey==c_INTEGER_DATA_TYPE){dtInteger=fileDataList[textIndex+1]}
    else if(strParseKey==c_LOGICAL_DATA_TYPE){dtLogical=fileDataList[textIndex+1]}
    else if(strParseKey==c_NUMERIC_DATA_TYPE){dtNumeric=fileDataList[textIndex+1]}
    else if(strParseKey==c_ORDEREDFACTOR_DATA_TYPE){dtOrderedFactor=fileDataList[textIndex+1]}
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
    matrixInformationList[[matrixInformationListCount]] = c(matrixName,fileName,delimiter,rows,columns,dtCharacter,dtFactor,dtInteger,dtLogical,dtNumeric,dtOrderedFactor)
    matrixInformationListCount = matrixInformationListCount + 1
  }
  return(matrixInformationList)
}

#16 Test cases
#Take a string of comma or dash seperated integer strings and convert into a vector
#of integers to use in index slicing
#TODO Should writing the read config be in words
funcParseIndexSlices = function(strIndexString,cstrNames)
{
  #If the slices are NA then return
  if(is.na(strIndexString)){return(strIndexString)}

  #List of indices to return
  viRetIndicies = c()

  #Split on commas
  lIndexString = strsplit(strIndexString, c_COMMA)
  for(strIndexItem in lIndexString[[1]])
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
}
