c_strDir <- file.path(getwd( ),"..")

source(file.path(c_strDir,"lib","Constants.R"))
source(file.path(c_strDir,"lib","SummarizeMaaslin.R"))
source(file.path(c_strDir,"lib","Utility.R"))

context("Test funcSummarizeDirectory")
strDirectoryNone = file.path(c_strDir,c_strTestingDirectory,c_strTestingInput,"funcSummarizeDirectory","None")
strDirectory1 = file.path(c_strDir,c_strTestingDirectory,c_strTestingInput,"funcSummarizeDirectory","1")
strDirectory3 = file.path(c_strDir,c_strTestingDirectory,c_strTestingInput,"funcSummarizeDirectory","3")
strFileBase1 = "FileBase1.txt"
strFileBase2 = "FileBase2.txt"

sKeyword = "Q.value"
sAltKeyword = "P.value"
sAltSignificance = "0.35"

sBaseName = "FuncSummarizeDirectory"

#Output and answer files
sNoFileResult = file.path(c_strDir,c_strTestingDirectory,c_strTemporaryFiles,"FuncSummarizeDirectory-NoFileResult.txt")
sNoFileResultAltKeyword = file.path(c_strDir,c_strTestingDirectory,c_strTemporaryFiles,"FuncSummarizeDirectory-NoFileAltKeyResult.txt")
sNoFileResultAltSig = file.path(c_strDir,c_strTestingDirectory,c_strTemporaryFiles,"FuncSummarizeDirectory-NoFileAltSigResult.txt")
sNoFileResultAnswer = file.path(c_strDir,c_strTestingDirectory,c_strCorrectAnswers,"FuncSummarizeDirectory-NoFileAnswer.txt")
sNoFileResultAnswerAltKeyword = file.path(c_strDir,c_strTestingDirectory,c_strCorrectAnswers,"FuncSummarizeDirectory-NoFileAltKeyAnswer.txt")
sNoFileResultAnswerAltSig = file.path(c_strDir,c_strTestingDirectory,c_strCorrectAnswers,"FuncSummarizeDirectory-NoFileAltSigAnswer.txt")
unlink(sNoFileResult)
sCorrectResults1File = file.path(c_strDir,c_strTestingDirectory,c_strTemporaryFiles,"FuncSummarizeDirectory-1FileResult.txt")
sCorrectResults1FileAltKeyword = file.path(c_strDir,c_strTestingDirectory,c_strTemporaryFiles,"FuncSummarizeDirectory-1FileAltKeyResult.txt")
sCorrectResults1FileAltSig = file.path(c_strDir,c_strTestingDirectory,c_strTemporaryFiles,"FuncSummarizeDirectory-1FileAltSigResult.txt")
sCorrectResults1FileAnswer = file.path(c_strDir,c_strTestingDirectory,c_strCorrectAnswers,"FuncSummarizeDirectory-1FileResult.txt")
sCorrectResults1FileAnswerAltKeyword = file.path(c_strDir,c_strTestingDirectory,c_strCorrectAnswers,"FuncSummarizeDirectory-1FileAltKeyResult.txt")
sCorrectResults1FileAnswerAltSig = file.path(c_strDir,c_strTestingDirectory,c_strCorrectAnswers,"FuncSummarizeDirectory-1FileAltSigResult.txt")
unlink(sCorrectResults1File)
sCorrectResults3Files = file.path(c_strDir,c_strTestingDirectory,c_strTemporaryFiles,"FuncSummarizeDirectory-3FileResult.txt")
sCorrectResults3FilesAltKeyword = file.path(c_strDir,c_strTestingDirectory,c_strTemporaryFiles,"FuncSummarizeDirectory-3FileAltKeyResult.txt")
sCorrectResults3FilesAltSig = file.path(c_strDir,c_strTestingDirectory,c_strTemporaryFiles,"FuncSummarizeDirectory-3FileAltSigResult.txt")
sCorrectResults3FilesAnswer = file.path(c_strDir,c_strTestingDirectory,c_strCorrectAnswers,"FuncSummarizeDirectory-3FileResult.txt")
sCorrectResults3FilesAnswerAltKeyword = file.path(c_strDir,c_strTestingDirectory,c_strCorrectAnswers,"FuncSummarizeDirectory-3FileAltKeyResult.txt")
sCorrectResults3FilesAnswerAltSig = file.path(c_strDir,c_strTestingDirectory,c_strCorrectAnswers,"FuncSummarizeDirectory-3FileAltSigResult.txt")
unlink(sCorrectResults3Files)

#Run tests
funcSummarizeDirectory(astrOutputDirectory=strDirectoryNone, strBaseName=sBaseName, astrSummaryFileName=sNoFileResult, astrKeyword=sKeyword, afSignificanceLevel="0.25")
funcSummarizeDirectory(astrOutputDirectory=strDirectory1, strBaseName=sBaseName, astrSummaryFileName=sCorrectResults1File, astrKeyword=sKeyword, afSignificanceLevel="0.25")
funcSummarizeDirectory(astrOutputDirectory=strDirectory3, strBaseName=sBaseName, astrSummaryFileName=sCorrectResults3Files, astrKeyword=sKeyword, afSignificanceLevel="0.25")

funcSummarizeDirectory(astrOutputDirectory=strDirectoryNone, strBaseName=sBaseName, astrSummaryFileName=sNoFileResultAltKeyword, astrKeyword=sAltKeyword, afSignificanceLevel="0.25")
funcSummarizeDirectory(astrOutputDirectory=strDirectory1, strBaseName=sBaseName, astrSummaryFileName=sCorrectResults1FileAltKeyword, astrKeyword=sAltKeyword, afSignificanceLevel="0.25")
funcSummarizeDirectory(astrOutputDirectory=strDirectory3, strBaseName=sBaseName, astrSummaryFileName=sCorrectResults3FilesAltKeyword, astrKeyword=sAltKeyword, afSignificanceLevel="0.25")

funcSummarizeDirectory(astrOutputDirectory=strDirectoryNone, strBaseName=sBaseName, astrSummaryFileName=sNoFileResultAltSig, astrKeyword= sKeyword, afSignificanceLevel=sAltSignificance)
funcSummarizeDirectory(astrOutputDirectory=strDirectory1, strBaseName=sBaseName, astrSummaryFileName=sCorrectResults1FileAltSig, astrKeyword= sKeyword, afSignificanceLevel=sAltSignificance)
funcSummarizeDirectory(astrOutputDirectory=strDirectory3, strBaseName=sBaseName, astrSummaryFileName=sCorrectResults3FilesAltSig, astrKeyword= sKeyword, afSignificanceLevel=sAltSignificance)

test_that("Check the cases where no, and real summary files exist.",{
  expect_equal(readLines(sNoFileResult),readLines(sNoFileResultAnswer))
  expect_equal(readLines(sCorrectResults1File),readLines(sCorrectResults1FileAnswer))
  expect_equal(readLines(sCorrectResults3Files),readLines(sCorrectResults3FilesAnswer))
})

test_that("Check changing the keyword.",{
  expect_equal(readLines(sNoFileResultAltKeyword),readLines(sNoFileResultAnswerAltKeyword))
  expect_equal(readLines(sCorrectResults1FileAltKeyword),readLines(sCorrectResults1FileAnswerAltKeyword))
  expect_equal(readLines(sCorrectResults3FilesAltKeyword),readLines(sCorrectResults3FilesAnswerAltKeyword))
})

test_that("Check that changing the significance threshold effects inclusion.",{
  expect_equal(readLines(sNoFileResultAltSig),readLines(sNoFileResultAnswerAltSig))
  expect_equal(readLines(sCorrectResults1FileAltSig),readLines(sCorrectResults1FileAnswerAltSig))
  expect_equal(readLines(sCorrectResults3FilesAltSig),readLines(sCorrectResults3FilesAnswerAltSig))
})
