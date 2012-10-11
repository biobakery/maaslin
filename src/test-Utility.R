library(testthat)
source("Constants.R")
source("Utility.R")

context("Test funcRename")
test_that("Test that unclassified and none otus are represented as 2 terminal clades and others are 1",{
  expect_equal(funcRename(paste("A","B","C","D",c_strUnclassified, sep=c_cFeatureDelim)),paste("D",c_strUnclassified, sep=c_cFeatureDelim))
  expect_equal(funcRename(paste("A","B","C","D","101", sep=c_cFeatureDelim)),paste("D","101", sep=c_cFeatureDelim))
  expect_equal(funcRename(paste("A","B","C","D", sep=c_cFeatureDelim)),paste("D", sep=c_cFeatureDelim))
  expect_equal(funcRename(paste("A", sep=c_cFeatureDelim)),paste("A", sep=c_cFeatureDelim))
  expect_equal(funcRename(paste(c_strUnclassified, sep=c_cFeatureDelim)),paste(c_strUnclassified, sep=c_cFeatureDelim))
  expect_equal(funcRename(paste("101", sep=c_cFeatureDelim)),paste("101", sep=c_cFeatureDelim))
})

context("Test funcColorHelper")
test_that("Test that min is min and max is max and average is average even if given as NA",{
  expect_equal(funcColorHelper( dMax = 1, dMin = 1, dMed = NA ), list( dMin = 1, dMax = 1, dMed = 1))
  expect_equal(funcColorHelper( dMax = -3, dMin = 10, dMed = NA ), list( dMin = -3, dMax = 10, dMed = 3.5))
  expect_equal(funcColorHelper( dMax = 1, dMin = 11, dMed = NA ), list( dMin = 1, dMax = 11, dMed = 6))
  expect_equal(funcColorHelper( dMax = 4, dMin = 10, dMed = 5 ), list( dMin = 4, dMax = 10, dMed = 5))
  expect_equal(funcColorHelper( dMax = 10, dMin = 4, dMed = 5 ), list( dMin = 4, dMax = 10, dMed = 5))
})

context("Test funcTrim")
test_that("Test that white spaces at the beginning and end of s string are removed",{
  expect_equal(funcTrim("TRIM"),"TRIM")
  expect_equal(funcTrim(" TRIM"),"TRIM")
  expect_equal(funcTrim("  TRIM"),"TRIM")
  expect_equal(funcTrim(" TRIM "),"TRIM")
  expect_equal(funcTrim("TRIM "),"TRIM")
  expect_equal(funcTrim("      TRIM          "),"TRIM")
  expect_equal(funcTrim("TR IM"),"TR IM")
  expect_equal(funcTrim(" TR IM"),"TR IM")
  expect_equal(funcTrim("  TR I M"),"TR I M")
  expect_equal(funcTrim(" TR IM "),"TR IM")
  expect_equal(funcTrim("T R IM "),"T R IM")
  expect_equal(funcTrim("      T RIM          "),"T RIM")
})

#TODO currently the capture versio of this does not produce a tabbed table (or default table delim) which is not consistent with the rest of the code base.
context("Test funcWrite")
#Answer files
c_sAnswerWriteFile1 = file.path(c_strTestingDirectory,c_strCorrectAnswers,"FuncWriteTemp1.txt")
c_sAnswerWriteFile2 = file.path(c_strTestingDirectory,c_strCorrectAnswers,"FuncWriteTemp2.txt")
c_sAnswerWriteDFFile1 = file.path(c_strTestingDirectory,c_strCorrectAnswers,"FuncWriteTempDF1.txt")
c_sAnswerWriteDFFile2 = file.path(c_strTestingDirectory,c_strCorrectAnswers,"FuncWriteTempDF2.txt")

#Working files
c_sTempWriteFile1 = file.path(c_strTestingDirectory,c_strTemporaryFiles,"FuncWriteTemp1.txt")
c_sTempWriteFile2 = file.path(c_strTestingDirectory,c_strTemporaryFiles,"FuncWriteTemp2.txt")
c_sTempWriteDFFile1 = file.path(c_strTestingDirectory,c_strTemporaryFiles,"FuncWriteTempDF1.txt")
c_sTempWriteDFFile2 = file.path(c_strTestingDirectory,c_strTemporaryFiles,"FuncWriteTempDF2.txt")
dfTest = as.data.frame(as.matrix(cbind(c(1,11,111),c(2,22,222),c(3,33,333))))
sWriteString = "Testing, 1,2,3 anything but that."
unlink(c_sTempWriteFile1)
unlink(c_sTempWriteFile2)
unlink(c_sTempWriteDFFile1)
unlink(c_sTempWriteDFFile2)
funcWrite(sWriteString,c_sTempWriteFile1)
funcWrite(sWriteString,c_sTempWriteFile2)
funcWrite(sWriteString,c_sTempWriteFile2)
funcWrite(dfTest,c_sTempWriteDFFile1)
funcWrite(dfTest,c_sTempWriteDFFile2)
funcWrite(dfTest,c_sTempWriteDFFile2)

test_that("Test that a test file is written and appended to for strings and dataframes.",{
  expect_equal(readLines(c_sTempWriteFile1),readLines(c_sAnswerWriteFile1))
  expect_equal(readLines(c_sTempWriteFile2),readLines(c_sAnswerWriteFile2))
  expect_equal(readLines(c_sTempWriteDFFile1),readLines(c_sAnswerWriteDFFile1))
  expect_equal(readLines(c_sTempWriteDFFile2),readLines(c_sAnswerWriteDFFile2))
})

context("Test funcWriteTable")
#Answer files
c_sAnswerWriteDFFile1 = file.path(c_strTestingDirectory,c_strCorrectAnswers,"FuncWriteTableTempDF1.txt")
c_sAnswerWriteDFFile2 = file.path(c_strTestingDirectory,c_strCorrectAnswers,"FuncWriteTableTempDF2.txt")

#Working files
c_sTempWriteDFFile1 = file.path(c_strTestingDirectory,c_strTemporaryFiles,"FuncWriteTableTempDF1.txt")
c_sTempWriteDFFile2 = file.path(c_strTestingDirectory,c_strTemporaryFiles,"FuncWriteTableTempDF2.txt")
unlink(c_sTempWriteDFFile1)
unlink(c_sTempWriteDFFile2)
funcWriteTable(dfTest,c_sTempWriteDFFile1)
funcWriteTable(dfTest,c_sTempWriteDFFile2, fAppend=TRUE)
funcWriteTable(dfTest,c_sTempWriteDFFile2, fAppend=TRUE)

test_that("Test that a test file is written and appended to for dataframes.",{
  expect_equal(readLines(c_sTempWriteDFFile1),readLines(c_sAnswerWriteDFFile1))
  expect_equal(readLines(c_sTempWriteDFFile2),readLines(c_sAnswerWriteDFFile2))
})

context("Test funcCoef2Col")
dfTestWithFactors = as.data.frame(as.matrix(cbind(c(1,11,111),c(2,22,222),c(3,33,333))))
colnames(dfTestWithFactors)=c("A","B","C")
dfTestWithFactors["B"]=as.factor(as.character(dfTestWithFactors[["B"]]))
test_that("Test that a coefficients are found or not given if they exist",{
  expect_equal(funcCoef2Col(strCoef="C",frmeData=dfTestWithFactors,astrCols=c()),"C")
  expect_equal(funcCoef2Col(strCoef="A",frmeData=dfTestWithFactors,astrCols=c()),"A")
  expect_equal(funcCoef2Col(strCoef=paste("B","2",sep=c_sFactorNameSep),frmeData=dfTestWithFactors,astrCols=c()),"B")
  expect_equal(funcCoef2Col(strCoef=paste("B","22",sep=c_sFactorNameSep),frmeData=dfTestWithFactors,astrCols=c()),"B")
  expect_equal(funcCoef2Col(strCoef=paste("B","222",sep=c_sFactorNameSep),frmeData=dfTestWithFactors,astrCols=c()),"B")
  expect_equal(funcCoef2Col(strCoef="C",frmeData=dfTestWithFactors,astrCols=c("A","B","C")),"C")
  expect_equal(funcCoef2Col(strCoef=paste("B","2",sep=c_sFactorNameSep),frmeData=dfTestWithFactors,astrCols=c("A","B")),"B")
  expect_equal(funcCoef2Col(strCoef=paste("B","22",sep=c_sFactorNameSep),frmeData=dfTestWithFactors,astrCols=c("B","C")),"B")
  expect_equal(funcCoef2Col(strCoef=paste("B","222",sep=c_sFactorNameSep),frmeData=dfTestWithFactors,astrCols=c("B")),"B")
})
