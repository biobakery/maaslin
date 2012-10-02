library(testthat)
source("Constants.R")
source("IO.R")
source("ValidateData.R")

context("Test funcParseIndexSlices")

cNames = c("One","Two","Three","Four","Five","Six","Seven","Eight","Nine","Ten","Eleven",
  "Twelve","Thirteen","Fourteen","Fifteen")

test_that("Just Numerics are parsed",{
  expect_equal(funcParseIndexSlices("1",cNames),c(1))
  expect_equal(funcParseIndexSlices("8,10",cNames),c(8,10))
  expect_equal(funcParseIndexSlices("2-6",cNames), c(2,3,4,5,6))
  expect_equal(funcParseIndexSlices("3,7,10-12",cNames), c(3,7,10,11,12))
})

test_that("Missing numbers are parsed",{
  expect_equal(funcParseIndexSlices("-",cNames), c(1:15))
  expect_equal(funcParseIndexSlices("-4",cNames), c(1,2,3,4))
  expect_equal(funcParseIndexSlices("3-",cNames), c(3:15))
})

test_that("Words are parsed correctly",{
  expect_equal(funcParseIndexSlices("One",cNames), c(1))
  expect_equal(funcParseIndexSlices("Eight,Ten",cNames), c(8,10))
  expect_equal(funcParseIndexSlices("Two-Six",cNames), c(2,3,4,5,6))
  expect_equal(funcParseIndexSlices("Three,Seven,Ten-Twelve",cNames), c(3,7,10,11,12)) 
})

test_that("Missing words are parsed",{
  expect_equal(funcParseIndexSlices("-Four",cNames), c(1:4))
  expect_equal(funcParseIndexSlices("Three-",cNames), c(3:15))
})

test_that("Words and numbers are parsed correctly",{
  expect_equal(funcParseIndexSlices("Eight,10",cNames), c(8,10))
  expect_equal(funcParseIndexSlices("2-Six",cNames), c(2,3,4,5,6))
  expect_equal(funcParseIndexSlices("Three,7,10-Twelve",cNames), c(3,7,10,11,12))
})

context("Test funcWriteMatrixToReadConfigFile")
# File to temporarily write to
strWriteMatrixRCTestFile = file.path(c_strTestingDirectory,c_strTemporaryFiles,"FuncWriteMatrixToReadConfigFileTemp.read.config")
# Files that hold answers
strFileSimpleRCFileAnswer = file.path(c_strTestingDirectory,c_strCorrectAnswers,"FuncWriteMatrixToReadConfigFile_SimpleAnswer.read.config")
strFileUseAllRCFileAnswer = file.path(c_strTestingDirectory,c_strCorrectAnswers,"FuncWriteMatrixToReadConfigFile_AllAnswer.read.config")
strFileAppendTrueRCFileAnswer = file.path(c_strTestingDirectory,c_strCorrectAnswers,"FuncWriteMatrixToReadConfigFile_AppendAnswer.read.config")
#Input matrix file
strFileMatrix = file.path(c_strTestingDirectory,c_strTestingInput,"TestMatrix.tsv")

#Get read config files in different scenarios
funcWriteMatrixToReadConfigFile(strWriteMatrixRCTestFile,strFileMatrix,"SimpleMatrix")
strSimpleInterface = readLines(strWriteMatrixRCTestFile)

funcWriteMatrixToReadConfigFile(strWriteMatrixRCTestFile,strFileMatrix,"AllMatrix",strRowIndices="1,2,3,4,5",
  strColIndices="10,11,12",strDtCharacter="1", strDtFactoral="2", strDtInteger="3", strDtLogical="4", strDtNumeric="5", strDtOrdered="6",
  acharDelimiter=" ")
strUseAllParametersInterface = readLines(strWriteMatrixRCTestFile)

funcWriteMatrixToReadConfigFile(strWriteMatrixRCTestFile,strFileMatrix,"SimpleMatrix")
funcWriteMatrixToReadConfigFile(strWriteMatrixRCTestFile,strFileMatrix,"SimpleMatrix")
strAppendFalseInterface = readLines(strWriteMatrixRCTestFile)

funcWriteMatrixToReadConfigFile(strWriteMatrixRCTestFile,strFileMatrix,"SimpleMatrix")
funcWriteMatrixToReadConfigFile(strWriteMatrixRCTestFile,strFileMatrix,"SimpleMatrix",fAppend=TRUE)
strAppendTrueInterface = readLines(strWriteMatrixRCTestFile)

test_that("Correct config file is written",{
  expect_equal(strSimpleInterface,readLines(strFileSimpleRCFileAnswer))
  expect_equal(strUseAllParametersInterface,readLines(strFileUseAllRCFileAnswer))
  expect_equal(strAppendFalseInterface,readLines(strFileSimpleRCFileAnswer))
  expect_equal(strAppendTrueInterface,readLines(strFileAppendTrueRCFileAnswer))
})

context("Test readConfigFile")
lsSimpleRC = funcReadConfigFile(strFileSimpleRCFileAnswer)
lsAllRC = funcReadConfigFile(strFileUseAllRCFileAnswer)

lsSimpleListAnswer = list()
lsSimpleListAnswer[[1]]=c("SimpleMatrix","testing/input/TestMatrix.tsv","\t","1","1","-","-",NA,NA,NA,NA,NA,NA)
lsAllListAnswer = list()
lsAllListAnswer[[1]]=c("AllMatrix","testing/input/TestMatrix.tsv"," ","1","1","1,2,3,4,5","10,11,12","1","2","3","4","5","6")

test_that("Test readConfigFile reads in files correctly.",{
  expect_equal(lsSimpleRC,lsSimpleListAnswer)
  expect_equal(lsAllRC,lsAllListAnswer)
})

#contest("Test funcReadMatrix")
#tempMatrixName = NA, tempFileName = NA, tempDelimiter = NA, tempIdRow = NA, tempIdCol = NA, tempRows = NA, tempColumns = NA, tempDtCharacter = NA, tempDtFactor = NA, tempDtInteger = NA, tempDtLogical = NA, tempDtNumeric = NA, tempDtOrderedFactor = NA, tempLog = FALSE

#Read in config files
#dfSimpleRead = funcReadMatrix("SimpleMatrix",strFileMatrix,1,1,"7,3,5","2,4,6","7","3","5")
#dfUseAllParametersRead = funcReadMatrix("AllMatrix",strFileMatrix,1,1,"7,3,5","2,4,6",NA,NA,NA,"7","3","5")

#dfSimpleReadCorrect = as.data.frame()
#dfUseAllReadCorrect = as.data.frame()

#test_that("Matrix file is read correctly.",{
#  expect_equal(dfSimpleRead,dfSimpleReadCorrect)
#  expect_equal(dfUseAllParametersRead,dfUseAllREadCorrect)
#})

#context("Test funcWriteMatrices")
#context("Test funcReadMatrices")
