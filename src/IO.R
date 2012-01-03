####################################
# Summary: IO
# Author: Timothy Tickle
# Start Date: 11-01-2011
####################################

#Project Constants
#source("./src/Constants.R")

c_astrNA <- c(""," ","  ","NA","na")

#Do not report warnings
options(warn=-1)

#Write data frame data files with config files
#dataFrameList is a named list of data frames (what you get directly from the read function)
#saveFileList File names to save the data matrices in (one name per data frame)
#configureFileName Name of the configure file to be written which will direct the reading of these data
funcWriteMatrices = function(dataFrameList = NA, saveFileList = NA, configureFileName = NA, acharDelimiter=",", log = FALSE)
{
  #Validate input
  if(is.na(dataFrameList))
  {
    print("Received NA for the data frame list, returned false and produced no output.")
    return(FALSE)
  }
  if(is.na(saveFileList))
  {
    print("Received NA for the save file list, returned false and produced no output.")
    return(FALSE)
  }

  #Get names
  dataFrameNames = names(dataFrameList)

  #Get length of dataFrameList
  dataFrameListLength = length(dataFrameList)

  #Get length of save file list
  saveFileListLength = length(saveFileList)

  #If the save file list length and data frame list length are not equal, abort
  if(!saveFileListLength == dataFrameListLength)
  {
    print(paste("Received a length of save files (",saveFileListLength,") that are different from the count of data frames (",dataFrameListLength,"). Stopped and returned false."),sep="")
    return(FALSE)
  }

  #For each data save
  for (dataIndex in c(1:dataFrameListLength))
  {
    #Current data frame
    data = dataFrameList[[dataIndex]]

    #Get row and column count
    dataDimension = dim(data)
    rowCount = dataDimension[1]
    columnCount = dataDimension[2]

    #Get row and column names
    rowNames = row.names(data)
    rowNamesString = paste(rowNames,sep="",collapse=",")
    if(length(rowNamesString)==0)
    {
      rowNamesString = NA
    }
    columnNames = colnames(data)
    columnNamesString = paste(columnNames,sep="",collapse=",")
    if(length(columnNamesString)==0)
    {
      columnNamesString = NA
    }

    #Get row indices
    rowStart = 1
    if(!is.na(rowNamesString)){  rowStart = 2}
    rowEnd = rowCount+rowStart - 1
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

    #Adjust for column of row names
    if(length(dtCharacter)>0)
    {
      dtCharacter = dtCharacter+1
    }
    if(length(dtFactoral)>0)
    {
      dtFactoral = dtFactoral+1
    }
    if(length(dtInteger)>0)
    {
      dtInteger = dtInteger+1
    }
    if(length(dtLogical)>0)
    {
      dtLogical = dtLogical+1
    }
    if(length(dtNumeric)>0)
    {
      dtNumeric = dtNumeric+1
    }
    if(length(dtOrdered)>0)
    {
      dtOrdered = dtOrdered+1
    }

    #Convert mode data to string
    if(length(dtCharacter)==0)
    {
      dtCharacter = ""
    } else {
      dtCharacter = paste(dtCharacter,sep="",collapse=",")
    }
    if(length(dtFactoral)==0)
    {
      dtFactoral = ""
    } else {
      dtFactoral = paste(dtFactoral,sep="",collapse=",")
    }
    if(length(dtInteger)==0)
    {
      dtInteger = ""
    } else {
      dtInteger = paste(dtInteger,sep="",collapse=",")
    }
    if(length(dtLogical)==0)
    {
      dtLogical = ""
    } else {
      dtLogical = paste(dtLogical,sep="",collapse=",")
    }
    if(length(dtNumeric)==0)
    {
      dtNumeric = ""
    } else {
      dtNumeric = paste(dtNumeric,sep="",collapse=",")
    }
    if(length(dtOrdered)==0)
    {
      dtOrdered = ""
    } else {
      dtOrdered = paste(dtOrdered,sep="",collapse=",")
    }

    #Write Data to file
    write.table(data, saveFileList[dataIndex], quote = FALSE, sep = acharDelimiter, col.names = NA, row.names = rowNames, na = "NA", append = FALSE)

    #Write read config file
    #Make delimiter readable
    if(acharDelimiter == "\t"){
      acharDelimiter = "TAB"
    }else if(acharDelimiter == " "){
      acharDelimiter = "SPACE"
    }else if(acharDelimiter == "\r"){
      acharDelimiter = "RETURN"
    }else if(acharDelimiter == "\n"){
      acharDelimiter = "ENDLINE"
    }

    cat(paste(c_MATRIX_NAME," ",dataFrameNames[dataIndex],sep = ""),"\n",file=configureFileName, append=TRUE)
    cat(paste(c_FILE_NAME," ",saveFileList[dataIndex],sep = ""),"\n",file=configureFileName, append=TRUE)
    cat(paste(c_DELIMITER," ",acharDelimiter,sep = ""),"\n",file=configureFileName, append=TRUE)
    cat(paste(c_ID_ROW," ","1",sep = ""),"\n",file=configureFileName, append=TRUE)
    cat(paste(c_ID_COLUMN," ","1",sep = ""),"\n",file=configureFileName, append=TRUE)
    cat(paste(c_ROWS," ",rowIndices,sep = ""),"\n",file=configureFileName, append=TRUE)
    cat(paste(c_COLUMNS," ",colIndices,sep = ""),"\n",file=configureFileName, append=TRUE)
    cat(paste(c_CHARACTER_DATA_TYPE," ",dtCharacter,sep = ""),"\n",file=configureFileName, append=TRUE)
    cat(paste(c_FACTOR_DATA_TYPE," ",dtFactoral,sep = ""),"\n",file=configureFileName, append=TRUE)
    cat(paste(c_INTEGER_DATA_TYPE," ",dtInteger,sep = ""),"\n",file=configureFileName, append=TRUE)
    cat(paste(c_LOGICAL_DATA_TYPE," ",dtLogical,sep = ""),"\n",file=configureFileName, append=TRUE)
    cat(paste(c_NUMERIC_DATA_TYPE," ",dtNumeric,sep = ""),"\n",file=configureFileName, append=TRUE)
    cat(paste(c_ORDEREDFACTOR_DATA_TYPE," ",dtOrdered,sep = ""),"\n",file=configureFileName, append=TRUE)
    cat("\n",file=configureFileName, append=TRUE)
  }
  return(TRUE)
}

#Dynamically Read a Matrix/Matrices from a configure file
funcReadMatrices = function( configureFile , defaultFile = NA, log = FALSE)
{
  #Named vector to return data frames read
  returnFrames = list()
  #Holds the names of the frames as they are being added
  returnFrameNames = c()
  #Current Index adding to next
  returnFramesIndex = 1

  #Read in config file info
  configData = readConfigFile(configureFile, defaultFile)
  
  #Read each data block extracted from the config file
  for(dataBlock in configData)
  {
    #Read in matrix
	tempFileName = dataBlock[2]
    dataMatrix = funcReadMatrix(tempMatrixName=dataBlock[1], tempFileName=tempFileName, tempDelimiter=dataBlock[3], tempIdRow=dataBlock[4], tempIdCol=dataBlock[5], tempRows=dataBlock[6], tempColumns=dataBlock[7], tempDtCharacter=dataBlock[8], tempDtFactor=dataBlock[9], tempDtInteger=dataBlock[10], tempDtLogical=dataBlock[11], tempDtNumeric=dataBlock[12], tempDtOrderedFactor=dataBlock[13], tempLog=log)
    if(class(dataMatrix)=="data.frame")
    {
      #Add data frame
      returnFrames[[returnFramesIndex]] = dataMatrix
      returnFramesIndex = returnFramesIndex + 1
      returnFrameNames = c(returnFrameNames,dataBlock[1])
    } else {
      print(paste("Reading from file named ",tempFileName," did not produce a valid data frame.",sep=""))
      print("Received")
      print(dataMatrix)
    }
  }
  if((length(returnFrames)>0)&&(length(returnFrameNames)>0))
  {
    names(returnFrames) = returnFrameNames
  }
#  print(returnFrameNames)
  return(returnFrames)
}

#Read one matrix
funcReadMatrix = function(tempMatrixName = NA, tempFileName = NA, tempDelimiter = NA, tempIdRow = NA, tempIdCol = NA, tempRows = NA, tempColumns = NA, tempDtCharacter = NA, tempDtFactor = NA, tempDtInteger = NA, tempDtLogical = NA, tempDtNumeric = NA, tempDtOrderedFactor = NA, tempLog = FALSE)
{
  print("Start funcReadMatrix.")
  #Check parameter and make sure not NA
  if(!funcIsValid(tempMatrixName)){return(FALSE)}

  #Check to make sure there is a file name for the matrix
  if(!funcIsValid(tempFileName))
  {
    print(paste("No file name is given for the matrix ",tempMatrixName,". Please add a file name to read the matrix from.", sep=""))
    return(FALSE)
  }

  #Convert characters to vectors of indices
  if(funcIsValid(tempIdRow))
  {
    tempIdRow = funcParseIndexSlices(tempIdRow)
    if((tempIdRow == FALSE)||(length(tempIdRow) == 0)){  tempIdRow = NA}
  }
  if(funcIsValid(tempIdCol))
  {
    tempIdCol = funcParseIndexSlices(tempIdCol)
    if((tempIdCol == FALSE)||(length(tempIdCol) == 0)){  tempIdCol = NA}
  }
  if(funcIsValid(tempRows))
  {
    tempRows = funcParseIndexSlices(tempRows)
    if((tempRows == FALSE)||(length(tempRows) == 0)){  tempRows = NA}
  }
  if(funcIsValid(tempColumns))
  {
    tempColumns = funcParseIndexSlices(tempColumns)
    if((tempColumns == FALSE)||(length(tempColumns) == 0)){  tempColumns = NA}
  }
  if(funcIsValid(tempDtCharacter))
  {
    tempDtCharacter = funcParseIndexSlices(tempDtCharacter)
    if((tempDtCharacter == FALSE)||(length(tempDtCharacter) == 0)){  tempDtCharacter = NA}
  }
  if(funcIsValid(tempDtFactor))
  {
    tempDtFactor = funcParseIndexSlices(tempDtFactor)
    if((tempDtFactor == FALSE)||(length(tempDtFactor) == 0)){  tempDtFactor = NA}
  }
  if(funcIsValid(tempDtOrderedFactor))
  {
    tempDtOrderedFactor = funcParseIndexSlices(tempDtOrderedFactor)
    if((tempDtOrderedFactor == FALSE)||(length(tempDtOrderedFactor) == 0)){  tempDtOrderedFactor = NA}
  }
  if(funcIsValid(tempDtInteger))
  {
    tempDtInteger = funcParseIndexSlices(tempDtInteger)
    if((tempDtInteger == FALSE)||(length(tempDtInteger) == 0)){  tempDtInteger = NA}
  }
  if(funcIsValid(tempDtLogical))
  {
    tempDtLogical = funcParseIndexSlices(tempDtLogical)
    if((tempDtLogical == FALSE)||(length(tempDtLogical) == 0)){  tempDtLogical = NA}
  }
  if(funcIsValid(tempDtNumeric))
  {
    tempDtNumeric = funcParseIndexSlices(tempDtNumeric)
    if((tempDtNumeric == FALSE)||(length(tempDtNumeric) == 0)){  tempDtNumeric = NA}
  }
  print("Initial indices check.")
  #Check indices
  errorFound = FALSE
  #Id row/col if given should be 1 index
  if(funcIsValid(tempIdRow))
  {
    if(!length(tempIdRow) == 1)
    {
      print(paste("The indices for the row containing the ids for the column should be of length 1, received=",tempIdRow,sep=""))
      errorFound = TRUE
    }
    if(length(tempRows)>1)
    {
      if(length(intersect(tempIdRow,tempRows)) == 1)
      {
        print(paste("Index indicated as an id row but was found in the data row indices, can not be both. Index=",tempIdRow," Data indices=",tempRows,sep=""))
        errorFound = TRUE
      }
    }
  }
  if(funcIsValid(tempIdCol))
  {
    if(!length(tempIdCol) == 1)
    {
      print(paste("The indices for the col containing the ids for the rows should be of length 1, received=",tempIdCol,sep=""))
      errorFound = TRUE
    }
    if(funcIsValid(tempColumns))
    {
      #Id row/col should not be in data row/col
      if(length(intersect(tempIdCol,tempColumns)) == 1)
      {
        print(paste("Index indicated as an id column but was found in the data column indices, can not be both. ID Index=",tempIdCol," Data Indices=",tempColumns,".",sep=""))
        errorFound = TRUE
      }
    }
  }
  if(errorFound){ return(FALSE)}

  #Manage elements not entered but associated with initial read
  #Variable to indicate a initial skipping of rows if not all the initial rows 
  #are needed (where to start reading). 0 indicates no skipping
  #Make sure it is non-negative and less than the idRow and row data indices
  startSkip = 0
  if(funcIsValid(tempRows))
  {
    startSkip = (tempRows[1] -1)
  }
  #If there is a id row before where data will start reading, read before this row
  #Additionally read in rows will be removed later in the process
  if(funcIsValid(tempIdRow)){if(tempIdRow<=startSkip){startSkip = (tempIdRow-1)}}
  if(startSkip < 0){startSkip == 0}
  if(!funcIsValid(tempDelimiter)){ tempDelimiter = "\t" }

  #Read in superset matrix and give names if indicated
  print("Read data.")
  if(!file.exists(tempFileName))
  {
    print(paste("The file given does not exist, could not read data. File:",tempFileName,sep=""))
    return(-1)
  }

  dataMatrix = read.table(tempFileName, header = FALSE, sep = tempDelimiter, as.is = TRUE, skip = startSkip, na.strings=c_astrNA)
  
  #Get column names
  print("Reading ID column/row.")
  columnNameList = list()
  rowNameList = list()
  #Read first row and column for row and column names
  if(funcIsValid(tempIdCol))
  {
    columnNameList = as.matrix(dataMatrix[tempIdRow,])
  }
  if(funcIsValid(tempIdRow))
  {
    rowNameList = dataMatrix[tempIdCol][[1]]
  }

  dataFrameDimension = dim(dataMatrix)
  #If the row names have the same length as the column count and has column names 
  #it is assumed that the tempIdCol index item is associated with the column names.
  #Visa versa for rows, either way it is removed
  print("Removing Row ID Column.")
  if(length(rowNameList)>0)
  {
    if(funcIsValid(tempIdCol)&&(length(rowNameList)==dataFrameDimension[1]))
    {
      #Remove ids from name vector
      rowNameList = rowNameList[(-1*tempIdCol)]
      #Remove ids from data
      dataMatrix = dataMatrix[(-1*tempIdRow)]
      #Adjust row ids given the removal of the id row
      if(funcIsValid(tempRows)){ tempRows[tempRows>tempIdRow]=(tempRows-1)}
      if(funcIsValid(tempDtCharacter)){ tempDtCharacter[tempDtCharacter>tempIdRow]=(tempDtCharacter-1)}
      if(funcIsValid(tempDtFactor)){ tempDtFactor[tempDtFactor>tempIdRow]=(tempDtFactor-1)}
      if(funcIsValid(tempDtInteger)){ tempDtInteger[tempDtInteger>tempIdRow]=(tempDtInteger-1)}
      if(funcIsValid(tempDtLogical)){ tempDtLogical[tempDtLogical>tempIdRow]=(tempDtLogical-1)}
      if(funcIsValid(tempDtNumeric)){ tempDtNumeric[tempDtNumeric>tempIdRow]=(tempDtNumeric-1)}
      if(funcIsValid(tempDtOrderedFactor)){ tempDtOrderedFactor[tempDtOrderedFactor>tempIdRow]=(tempDtOrderedFactor-1)}
    }
  }

  print("Removing Column ID Row.")
  if(length(columnNameList)>0)
  {
    if(funcIsValid(tempIdRow)&&(length(columnNameList)==dataFrameDimension[2]))
    {
      #Remove ids from vector
      columnNameList = columnNameList[(-1*tempIdRow)]
      #Remove ids from data
      dataMatrix = dataMatrix[(-1*tempIdCol),]
      #Adjust column ids given the removal of the id column
      if(funcIsValid(tempColumns)){ tempColumns[tempColumns>tempIdCol]=(tempColumns-1)}
    }
  }

  print("Adding name adornments to data frame.")
  #Add row and column names
  if(funcIsValid(tempIdRow)){  row.names(dataMatrix) = as.character(rowNameList)}
  if(funcIsValid(tempIdCol)){  colnames(dataMatrix) = as.character(columnNameList)}

  for(i in 1:ncol(dataMatrix)){
    dataMatrix[,i] <- type.convert(dataMatrix[,i], na.strings = c_astrNA)}

  modeError = FALSE
  #Check to make sure there are no duplicates in data mode indices (remove NAs before checking uniqueness, we done care about duplicate Nas)
  if((funcIsValid(tempDtNumeric)) || (funcIsValid(tempDtFactor)) || (funcIsValid(tempDtCharacter)) || (funcIsValid(tempDtInteger)) || (funcIsValid(tempDtLogical)) || (funcIsValid(tempDtOrderedFactor)))
  {
    allModes = c(tempDtNumeric,tempDtFactor,tempDtCharacter,tempDtInteger,tempDtLogical,tempDtOrderedFactor)
    allModes = allModes[!is.na(allModes)]

    if(anyDuplicated(allModes))
    {
      print("Data types (modes) indices appeared in more than one category. Make sure the indices in each data type (mode) category are unique to it's group. The current data frame is default mode (character data only).")
      print("Not unique:")
      print(allModes[duplicated(allModes)])
      print("Data mode numeric:")
      print(tempDtNumeric)
      print("Data mode factor:")
      print(tempDtFactor)
      print("Data mode character:")
      print(tempDtCharacter)
      print("Data mode integer:")
      print(tempDtInteger)
      print("Data mode logical:")
      print(tempDtLogical)
      print("Data mode ordered factor:")
      print(tempDtOrderedFactor)
    
      return(-1)
    } else {
      #Set modes of data
      if (funcIsValid(tempDtNumeric))
      {
        #Set numerics
        for(numericIndex in tempDtNumeric)
        {
          tryCatch({
           data=dataMatrix[[numericIndex]]
           dataMatrix[numericIndex] = as.numeric(as.character(dataMatrix[[numericIndex]]))
          }, error = function(ex) { modeError = TRUE
            print(paste("IO::Index ",numericIndex," could not be converted to a numeric vector.", sep=""))
            print(dataMatrix)
            return(-1)}, finally={})
        }
      }
      if (funcIsValid(tempDtFactor))
      {
        #Set factors
        for(factorIndex in tempDtFactor)
        {
          tryCatch({
#          print(dataMatrix[factorIndex])
#          print(as.factor(as.character(dataMatrix[[factorIndex]])))
          dataMatrix[factorIndex] = as.factor(as.character(dataMatrix[[factorIndex]]))
#          print(dataMatrix[factorIndex])
          }, error = function(ex) { modeError = TRUE
            print(paste("IO::Index ",factorIndex," could not be converted to a factor vector.", sep=""))}, finally={})
        }
      }
      if (funcIsValid(tempDtInteger))
      {
        #Set integers
        for(integerIndex in tempDtInteger)
        {
          tryCatch({
          dataMatrix[integerIndex] = as.integer(as.character(dataMatrix[[integerIndex]]))
          }, error = function(ex) { modeError = TRUE
            print(paste("IO::Index ",integerIndex," could not be converted to an integer vector.", sep=""))}, finally={})
        }
      }
      if (funcIsValid(tempDtLogical))
      {
        #Set logical
        for(logicalIndex in tempDtLogical)
        {
          tryCatch({
          toBeLogicalElements = dataMatrix[[logicalIndex]]
          convertedElements = c()
          for(convertingElement in toBeLogicalElements)
          {
            characterAttempt = as.logical(as.character(convertingElement))
            if(is.na(characterAttempt))
            {
              characterAttempt = as.logical(as.numeric(as.character(convertingElement)))
            }
            convertedElements = c(convertedElements,characterAttempt)
          }
          dataMatrix[logicalIndex] = convertedElements
          }, error = function(ex) { modeError = TRUE
            print(paste("IO::Index ",logicalIndex," could not be converted to a logical vector.", sep=""))}, finally={})
        }
      }
      if (funcIsValid(tempDtOrderedFactor))
      {
        #Set orderedFactor
        for(orderedIndex in tempDtOrderedFactor)
        {
          tryCatch({
#          print(dataMatrix[orderedIndex])
#          print(as.ordered(as.character(dataMatrix[[orderedIndex]])))
          dataMatrix[orderedIndex] = as.ordered(as.character(dataMatrix[[orderedIndex]]))
#          print(dataMatrix[orderedIndex])
          }, error = function(ex) { modeError = TRUE
            print(paste("IO::Index ",orderedIndex," could not be converted to an ordered factor vector.", sep=""))}, finally={})
        }
      }
      if (funcIsValid(tempDtCharacter))
      {
        #Set Character
        for(characterIndex in tempDtCharacter)
        {
          tryCatch({
#          print(dataMatrix[characterIndex])
#          print(as.character(dataMatrix[[characterIndex]]))
          dataMatrix[characterIndex] = as.character(dataMatrix[[characterIndex]])
#          print(dataMatrix[characterIndex])
          }, warning={}, error = function(ex) { modeError = TRUE
            print(paste("IO::Index ",characterIndex," could not be converted to a character vector.", sep=""))}, finally={})
        }
      }
    }
  }

  print("Reduce Matrix.")
  #Reduce matrix
  #Account for when both column ranges and row ranges are given or just a column or a row range is given
  if(funcIsValid(tempColumns))
  {
    if(funcIsValid(tempRows))
    {
      dataMatrix = dataMatrix[tempRows,tempColumns, drop=FALSE]
    } else {
      dataMatrix = dataMatrix[,tempColumns, drop=FALSE]
    }
  } else {
    if(funcIsValid(tempRows))
    {
      dataMatrix = dataMatrix[tempRows,,drop=FALSE]
    }
  }

  #Give feed back to user
  if(tempLog)
  {
    print("Generating FeedBack.")
    #Feedback from file
    print(paste("Successfully loaded Data frame for matrix, details are:",sep=""))
    print(paste("Matrix Name: ", tempMatrixName, sep=""))
    print(paste("File Name: ", tempFileName, sep=""))
    delimiterString = tempDelimiter
    if(delimiterString == "\t"){delimiterString = "TAB"
    }else if(delimiterString == " "){delimiterString = "SPACE"
    }else if(delimiterString == "\r"){delimiterString = "RETURN"
    }else if(delimiterString == "\n"){delimiterString = "ENDLINE"}
    print(paste("Delimiter: ",delimiterString, sep=""))
    if(!funcIsValid(tempIdRow))
    {
      print(paste("No IDs for rows given.", sep=""))
    } else {
      print(paste("Row ids found in the col index ",paste(tempIdRow,sep="",collapse=","),".",sep=""))
      print(paste("Row ids are ",paste(row.names(dataMatrix),sep="",collapse=","),".",sep=""))
    }
    if(!funcIsValid(tempIdCol))
    {
      print(paste("No IDs for columns given.", sep=""))
    } else {
      print(paste("Column ids found in the row index ",paste(tempIdCol,sep="",collapse=","),".",sep=""))
      print(paste("Column ids are ",paste(colnames(dataMatrix),sep="",collapse=","),".",sep=""))
    }
    if(!funcIsValid(tempRows))
    {
      print(paste("No column indicies for data given so all column data kept.", sep=""))
    } else {
      print(paste("Data column indices are ",paste(tempColumns,sep="",collapse=","),".", sep=""))
    }
    if(!funcIsValid(tempColumns))
    {
      print(paste("No row indicies for data given so all column data kept.", sep=""))
    } else {
      print(paste("Data row indices are ",paste(tempRows,sep="",collapse=","),".", sep=""))
    }
    if(funcIsValid(tempDtCharacter))
    {
      print(paste("Character data are found at column indices ",paste(tempDtCharacter,sep="",collapse=","),".", sep=""))
    }
    if(funcIsValid(tempDtFactor))
    {
      print(paste("Factor data are found at column indices ",paste(tempDtFactor,sep="",collapse=","),".", sep=""))
    }
    if(funcIsValid(tempDtInteger))
    {
      print(paste("Integer data are found at column indices ",paste(tempDtInteger,sep="",collapse=","),".", sep=""))
    }
    if(funcIsValid(tempDtLogical))
    {
      print(paste("Logical data are found at column indices ",paste(tempDtLogical,sep="",collapse=","),".", sep=""))
    }
    if(funcIsValid(tempDtNumeric))
    {
      print(paste("Numeric data are found at column indices ",paste(tempDtNumeric,sep="",collapse=","),".", sep=""))
    }
    if(funcIsValid(tempDtOrderedFactor))
    {
      print(paste("Ordered Factor data are found at column indices ",paste(tempDtOrderedFactor,sep="",collapse=","),".", sep=""))
    }
    if((!funcIsValid(tempDtOrderedFactor))&&(!funcIsValid(tempDtNumeric))&&(!funcIsValid(tempDtLogical))&&(!funcIsValid(tempDtInteger))&&(!funcIsValid(tempDtFactor))&&(!funcIsValid(tempDtCharacter)))
    {
      print(paste("No data type information was given, all data was assumed to be numeric.", sep=""))
    }
    #Feedback from object
    print(paste("The shape of the date is ",paste(dim(dataMatrix),collapse=","),".", sep=""))
    print(paste("The column names of the matrix are ",paste(colnames(dataMatrix),collapse=","),".", sep=""))
    print(paste("The row names of the matrix are ",paste(row.names(dataMatrix),collapse=","),".", sep=""))
    if(modeError == TRUE){ print("Please note errors occured on converting data modes of some rows. Please check output, unsuccessful conversion leaves data as default (character mode).")}
  }
  print("Completed reading matrix.")
  #Return matrix
  return(dataMatrix)
}

#Reads in configure file and extracts the pieces needed for reading in a matrix
#Configure file = string path to configure file
readConfigFile = function(configureFile, defaultFile = NA)
{
  #Read configure file
  fileDataList <- scan( file = configureFile, what = character())
  matrixName <- NA
  fileName <- defaultFile

  #Hold information on matrices to be read
  matrixInformationList = list()
  matrixInformationListCount = 1

  for(textIndex in c(1:length(fileDataList)))
  {
    #Start at the Matrix name
    #Keep this if statement first so that you scan through until you find a matrix block
#    print(fileDataList[textIndex])
    if(fileDataList[textIndex] == c_MATRIX_NAME)
    {
      #If the file name is not NA then that is sufficient for a matrix, store
      #Either way reset
      if(funcIsValid(fileName)&&funcIsValid(matrixName))
      {
        matrixInformationList[[matrixInformationListCount]] = c(matrixName,fileName,delimiter,idRow,idCol,rows,columns,dtCharacter,dtFactor,dtInteger,dtLogical,dtNumeric,dtOrderedFactor)
        matrixInformationListCount = matrixInformationListCount + 1
      }

      #Get the matrix name and store
      textIndex = textIndex + 1
      matrixName = fileDataList[textIndex]

      fileName = defaultFile
      delimiter = NA
      idRow = 1
      idCol = 1
      rows = NA
      columns = NA
      dtCharacter = NA
      dtFactor = NA
      dtInteger = NA
      dtLogical = NA
      dtNumeric = NA
      dtOrderedFactor = NA 
      
#      print(paste(c("matrixName = ",matrixName)))
    } else if(is.na(matrixName)){
      #If is not matrix name and no matrix name is known skip until you find the matrix name
      #IF matrix name is known, continue to collect information about that matrix
#      print("next")
      next }

    #Find file name
    if(fileDataList[textIndex] == c_FILE_NAME)
    {
      textIndex = textIndex + 1
      fileName = fileDataList[textIndex]
#      print(paste(c("fileName = ",fileName)))
      next
    }

    #Find the row that has the ids for the columns
    if(fileDataList[textIndex] == c_DELIMITER)
    {
      textIndex = textIndex + 1
      delimiter = fileDataList[textIndex]
      if(delimiter == "TAB"){delimiter = "\t"
      }else if(delimiter == "SPACE"){delimiter = " "
      }else if(delimiter == "RETURN"){delimiter = "\r"
      }else if(delimiter == "ENDLINE"){delimiter = "\n"}
#      print(paste(c("Delimiter = ",delimiter)))
      next
    }

    #Find the row that has the ids for the columns
    if(fileDataList[textIndex] == c_ID_ROW)
    {
      textIndex = textIndex + 1
      idRow = fileDataList[textIndex]
#      print(paste(c("ID row = ",idRow)))
      next
    }

    #Find the column that has the ids for the rows
    if(fileDataList[textIndex] == c_ID_COLUMN)
    {
      textIndex = textIndex + 1
      idCol = fileDataList[textIndex]
#      print(paste(c("ID col = ",idCol)))
      next
    }

    #Find the rows for reading
    if(fileDataList[textIndex] == c_ROWS)
    {
      textIndex = textIndex + 1
      rows = fileDataList[textIndex]
#      print(paste(c("Rows = ",rows)))
      next
    }

    #Find the columns for reading
    if(fileDataList[textIndex] == c_COLUMNS)
    {
      textIndex = textIndex + 1
      columns = fileDataList[textIndex]
#      print(paste(c("Columns = ",columns)))
      next
    }

    #Find rows for character data types
    if(fileDataList[textIndex] == c_CHARACTER_DATA_TYPE)
    {
      textIndex = textIndex + 1
      dtCharacter = fileDataList[textIndex]
#      print(paste(c("Character Data Type = ",dtCharacter)))
      next
    }

    #Find rows for factor data types
    if(fileDataList[textIndex] == c_FACTOR_DATA_TYPE)
    {
      textIndex = textIndex + 1
      dtFactor = fileDataList[textIndex]
#      print(paste(c("Factor Data Type = ",dtFactor)))
      next
    }

    #Find rows for integer data types
    if(fileDataList[textIndex] == c_INTEGER_DATA_TYPE)
    {
      textIndex = textIndex + 1
      dtInteger = fileDataList[textIndex]
#      print(paste(c("Integer Data Type = ",dtInteger)))
      next
    }

    #Find rows for logical data types
    if(fileDataList[textIndex] == c_LOGICAL_DATA_TYPE)
    {
      textIndex = textIndex + 1
      dtLogical = fileDataList[textIndex]
#      print(paste(c("Logical Data Type = ",dtLogical)))
      next
    }

    #Find rows for Numeric data types
    if(fileDataList[textIndex] == c_NUMERIC_DATA_TYPE)
    {
      textIndex = textIndex + 1
      dtNumeric = fileDataList[textIndex]
#      print(paste(c("Numeric Data Type = ",dtNumeric)))
      next
    }

    #Find rows for ordered data types
    if(fileDataList[textIndex] == c_ORDEREDFACTOR_DATA_TYPE)
    {
      textIndex = textIndex + 1
      dtOrderedFactor = fileDataList[textIndex]
#      print(paste(c("Ordered FACTOR Data Type = ",dtOrderedFactor)))
      next
    }
  }
  #If there is matrix information left
  if((!is.na(matrixName)) && (!is.na(fileName)))
  {
    matrixInformationList[[matrixInformationListCount]] = c(matrixName,fileName,delimiter,idRow,idCol,rows,columns,dtCharacter,dtFactor,dtInteger,dtLogical,dtNumeric,dtOrderedFactor)
    matrixInformationListCount = matrixInformationListCount + 1
  }
  return(matrixInformationList)
}

#Take a string of comma or dash seperated integer strings and convert into a vector
#of integers to use in index slicing
funcParseIndexSlices = function(indexString)
{
  #Validate string parameter
  if(!funcIsValidString(indexString))
  { 
    print(paste("Received invalid string to parse as an index, recieved = ",indexString,".",sep=""))
    return(FALSE)
  }

  #List of indices to return
  indices = c()

  #Split on commas
  indexString = strsplit(indexString, c_COMMA)
  for(indexItem in indexString[[1]])
  {
    indexElement = strsplit(indexItem, c_DASH)[[1]]
    indexElementLength = length(indexElement)
    #If only one index is recieved
    if(indexElementLength == 1)
    {
      #Change index to numeric and add to index list
      #OR skip if it is not numeric
      index = as.numeric(indexElement)
      if(is.na(index))
      {
        print(paste("Index can not be changed to a number, skipped. Index =",index,".",sep = ""))
      } else {
        if(index > 0)
        {
          indices = c(indices,index)
        }
      }
    #If 2 indices are recieved with a dash
    #Try to change both to numeric
    #If they cant be change, skip
    } else if(indexElementLength == 2) {
        numericOne = as.numeric(indexElement[1])
        numericTwo = as.numeric(indexElement[2])
        if(is.na(numericOne) || is.na(numericTwo))
        {
          print(paste("One of the two indices can not be changed to a number, skipped. Before Indices = ",indexElementLength[1]," and ",indexElementLength[2],". After Indices =",numericOne," and ",numericTwo,".",sep = ""))
        } else {
          indices = c(indices,c(numericOne:numericTwo))
        }
    #If more than 2 indices are found then this is nonsensical and is ignored.
    } else if(indexElementLength > 2) {
      print(paste("Error in index, too many dashes, only one is allowed. Index = ",indexItem,sep=""))
    }
  }

  #Make indices unique
  indices = unique(indices)
  #Sort Indices
  indices = sort(indices)
  #Return sorted, unique indices
  return(indices)
}

#Test
#readConfigureFile = "./Test/TestReadMatrices.read.config"
#writeConfigureFile = "./Test/TestWriteMatrices.read.config"
#writeDataFiles = c("./Test/TestWriteMatrices.pcl")
#data = funcReadMatrices(configureFile = readConfigureFile, log = FALSE)
#funcWriteMatrices(dataFrameList = data, saveFileList = writeDataFiles, configureFileName = writeConfigureFile, log = FALSE)
