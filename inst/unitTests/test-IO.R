c_strDir <- file.path(getwd( ),"..")

source(file.path(c_strDir,"lib","Constants.R"))
source(file.path(c_strDir,"lib","ValidateData.R"))
strTestingDirectory = file.path(c_strDir,c_strTestingDirectory)

expect_equal(funcParseIndexSlices("1",cNames),c(1))

cNames = c("One","Two","Three","Four","Five","Six","Seven","Eight","Nine","Ten","Eleven",
  "Twelve","Thirteen","Fourteen","Fifteen")

test_that("Just Numerics are parsed",{
  expect_equal(funcParseIndexSlices("1",cNames),c(1))
  expect_equal(funcParseIndexSlices("8,10",cNames),c(8,10))
  expect_equal(funcParseIndexSlices("2-6",cNames), c(2,3,4,5,6))
  expect_equal(funcParseIndexSlices("3,7,10-12",cNames), c(3,7,10,11,12))
})

test_that("Missing numbers are parsed",{
  expect_equal(funcParseIndexSlices("-",cNames), c(2:15))
  expect_equal(funcParseIndexSlices("-4",cNames), c(2,3,4))
  expect_equal(funcParseIndexSlices("3-",cNames), c(3:15))
})

test_that("Words are parsed correctly",{
  expect_equal(funcParseIndexSlices("One",cNames), c(1))
  expect_equal(funcParseIndexSlices("Eight,Ten",cNames), c(8,10))
  expect_equal(funcParseIndexSlices("Two-Six",cNames), c(2,3,4,5,6))
  expect_equal(funcParseIndexSlices("Three,Seven,Ten-Twelve",cNames), c(3,7,10,11,12)) 
})

test_that("Missing words are parsed",{
  expect_equal(funcParseIndexSlices("-Four",cNames), c(2:4))
  expect_equal(funcParseIndexSlices("Three-",cNames), c(3:15))
})

test_that("Words and numbers are parsed correctly",{
  expect_equal(funcParseIndexSlices("Eight,10",cNames), c(8,10))
  expect_equal(funcParseIndexSlices("2-Six",cNames), c(2,3,4,5,6))
  expect_equal(funcParseIndexSlices("Three,7,10-Twelve",cNames), c(3,7,10,11,12))
})


context("Test funcWriteMatrixToReadConfigFile")
# File to temporarily write to
strWriteMatrixRCTestFile = file.path(strTestingDirectory,c_strTemporaryFiles,"FuncWriteMatrixToReadConfigFileTemp.read.config")
# Files that hold answers
strFileSimpleRCFileAnswer = file.path(strTestingDirectory,c_strCorrectAnswers,"FuncWriteMatrixToReadConfigFile_SimpleAnswer.read.config")
strFileUseAllRCFileAnswer = file.path(strTestingDirectory,c_strCorrectAnswers,"FuncWriteMatrixToReadConfigFile_AllAnswer.read.config")
strFileAppendTrueRCFileAnswer = file.path(strTestingDirectory,c_strCorrectAnswers,"FuncWriteMatrixToReadConfigFile_AppendAnswer.read.config")
#Input matrix file
strFileMatrix = file.path(strTestingDirectory,c_strTestingInput,"TestMatrix.tsv")

#Get read config files in different scenarios
funcWriteMatrixToReadConfigFile(strWriteMatrixRCTestFile,"SimpleMatrix")
strSimpleInterface = readLines(strWriteMatrixRCTestFile)
funcWriteMatrixToReadConfigFile(strWriteMatrixRCTestFile,"AllMatrix",strRowIndices="1,2,3,4,5", strColIndices="10,11,12",acharDelimiter=" ")
strUseAllParametersInterface = readLines(strWriteMatrixRCTestFile)
funcWriteMatrixToReadConfigFile(strWriteMatrixRCTestFile,"SimpleMatrix")
funcWriteMatrixToReadConfigFile(strWriteMatrixRCTestFile,"SimpleMatrix")
strAppendFalseInterface = readLines(strWriteMatrixRCTestFile)
funcWriteMatrixToReadConfigFile(strWriteMatrixRCTestFile,"SimpleMatrix")
funcWriteMatrixToReadConfigFile(strWriteMatrixRCTestFile,"SimpleMatrix",fAppend=TRUE)
strAppendTrueInterface = readLines(strWriteMatrixRCTestFile)

test_that("Correct config file is written",{
  expect_equal(strSimpleInterface,readLines(strFileSimpleRCFileAnswer))
  expect_equal(strUseAllParametersInterface,readLines(strFileUseAllRCFileAnswer))
  expect_equal(strAppendFalseInterface,readLines(strFileSimpleRCFileAnswer))
  expect_equal(strAppendTrueInterface,readLines(strFileAppendTrueRCFileAnswer))
})

context("Test funcReadConfigFile")
lsSimpleRC = funcReadConfigFile(strFileSimpleRCFileAnswer,strFileMatrix)
lsAllRC = funcReadConfigFile(strFileUseAllRCFileAnswer,strFileMatrix)

lsSimpleListAnswer = list()
lsSimpleListAnswer[[1]]=c("SimpleMatrix",strFileMatrix,"\t","-","-")
lsAllListAnswer = list()
lsAllListAnswer[[1]]=c("AllMatrix",strFileMatrix," ","1,2,3,4,5","10,11,12")

test_that("Test readConfigFile reads in files correctly.",{
  expect_equal(lsSimpleRC,lsSimpleListAnswer)
  expect_equal(lsAllRC,lsAllListAnswer)
})


context("Test funcReadMatrix")

#Read in config files
dfSimpleRead = funcReadMatrix("SimpleMatrix",strFileMatrix,"\t","2,4,5","7,3,5")
dfUseAllParametersRead = funcReadMatrix("AllMatrix",strFileMatrix,"\t","2,3,4","6,2,4")

dfSimpleReadCorrect = as.data.frame(as.matrix(rbind(c(21,23,24),c(41,43,44),c(61,63,64))))
rownames(dfSimpleReadCorrect) = c("Feature2", "Feature4", "Feature6")
colnames(dfSimpleReadCorrect) = c("Sample1", "Sample3", "Sample4")

dfUseAllReadCorrect = as.data.frame(as.matrix(rbind(c(11,12,13),c(31,32,33),c(51,52,53))))
rownames(dfUseAllReadCorrect) = c("Feature1", "Feature3", "Feature5")
colnames(dfUseAllReadCorrect) = c("Sample1", "Sample2", "Sample3")

test_that("Matrix file is read correctly.",{
  expect_equal(dfSimpleRead,dfSimpleReadCorrect)
  expect_equal(dfUseAllParametersRead,dfUseAllReadCorrect)
})

context("Test funcReadMatrices")

sConfigureFile1Matrix = file.path(strTestingDirectory,c_strTestingInput,"1Matrix.read.config")
mtxOne = as.data.frame(as.matrix(rbind(c(11,12,13,14,15),c(21,22,23,24,25),c(31,32,33,34,35),c(41,42,43,44,45),
                                                        c(51,52,53,54,55),c(61,62,63,64,65),c(71,72,73,74,75),c(81,82,83,84,85),
                                                        c(91,92,93,94,95),c(101,102,103,104,105),c(111,112,113,114,115),c(121,122,123,124,125),
                                                        c(131,132,133,134,135),c(141,142,143,144,145),c(151,152,153,154,155))))
rownames(mtxOne) = c("Feature1","Feature2","Feature3","Feature4","Feature5","Feature6","Feature7","Feature8","Feature9","Feature10",
                     "Feature11","Feature12","Feature13","Feature14","Feature15")
colnames(mtxOne) = c("Sample1","Sample2","Sample3","Sample4","Sample5")
sConfigureFile2Matrix = file.path(strTestingDirectory,c_strTestingInput,"2Matrix.read.config")
mtxTwo = as.data.frame(as.matrix(rbind(c(11,12,13),c(21,22,23),c(31,32,33))))
rownames(mtxTwo) = c("Feature1","Feature2","Feature3")
colnames(mtxTwo) = c("Sample1","Sample2","Sample3")

sConfigureFile3Matrix = file.path(strTestingDirectory,c_strTestingInput,"3Matrix.read.config")
mtxThree = as.data.frame(as.matrix(rbind(c(11,12,14),c(21,22,24),c(31,32,34),c(41,42,44),
                                         c(51,52,54),c(61,62,64),c(71,72,74),c(81,82,84),c(91,92,94))))
rownames(mtxThree) = c("Feature1","Feature2","Feature3","Feature4","Feature5","Feature6","Feature7","Feature8","Feature9")
colnames(mtxThree) = c("Sample1","Sample2","Sample4")

#Read one matrix
ldfRet1 = funcReadMatrices(configureFile=sConfigureFile1Matrix,strFileMatrix)
ldfRet1Answer = list( "Matrix1" = mtxOne)

#Read two matrices
ldfRet2 = funcReadMatrices(configureFile=sConfigureFile2Matrix,strFileMatrix)
ldfRet2Answer = list( "Matrix1" = mtxOne,
                      "Matrix2" = mtxTwo)

#Read three matrices from two different files
ldfRet3 = funcReadMatrices(configureFile=sConfigureFile3Matrix,strFileMatrix)
ldfRet3Answer = list( "Matrix1" = mtxOne,
                      "Matrix2" = mtxTwo,
                      "Matrix3" = mtxThree)

test_that("Test funcReadMatrices read in the correct matrices not matter the number or source",{
  expect_equal(ldfRet1,ldfRet1Answer)
  expect_equal(ldfRet2,ldfRet2Answer)
  expect_equal(ldfRet3,ldfRet3Answer)
})

context("Test funcWriteMatrices")
strFuncWriteMatricesMatrix1 = file.path(strTestingDirectory,c_strTemporaryFiles,"FuncWriteMatrices1.tsv")
strFuncWriteMatricesMatrix2 = file.path(strTestingDirectory,c_strTemporaryFiles,"FuncWriteMatrices2.tsv")
strFuncWriteMatricesMatrix1Answer = file.path(strTestingDirectory, c_strCorrectAnswers,"FuncWriteMatrices1.tsv")
strFuncWriteMatricesMatrix2Answer = file.path(strTestingDirectory, c_strCorrectAnswers,"FuncWriteMatrices2.tsv")
strFuncWriteMatricesRCFile = file.path(strTestingDirectory,c_strTemporaryFiles,"FuncWriteMatrices.read.config")
strFuncWriteMatricesRCFileAnswer = file.path(strTestingDirectory, c_strCorrectAnswers,"FuncWriteMatrices.read.config")
funcWriteMatrices(list("1"=mtxOne, "2"=mtxThree),c(strFuncWriteMatricesMatrix1, strFuncWriteMatricesMatrix2), strFuncWriteMatricesRCFile)

test_that("Test that writing to a file occurs correctly, for both matrix and configure file.",{
  expect_equal(readLines(strFuncWriteMatricesMatrix1Answer),readLines(strFuncWriteMatricesMatrix1))
  expect_equal(readLines(strFuncWriteMatricesMatrix2Answer),readLines(strFuncWriteMatricesMatrix2))
  expect_equal(readLines(strFuncWriteMatricesRCFileAnswer),readLines(strFuncWriteMatricesRCFile))
})
