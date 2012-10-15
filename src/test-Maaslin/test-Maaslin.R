library(testthat)
source(file.path("input","maaslin","src","Constants.R"))

#source("Constants.R")

context("Test Run From Commandline")

#Input files and directories
#Scripts
#sScriptBoostGLM = file.path("input","maaslin","src","BoostGLM.R")
#sScriptConstants = file.path("input","maaslin","src","Constants.R")
#sScriptIO = file.path("input","maaslin","src","IO.R")
#sScriptMaaslin = file.path("input","maaslin","src","Maaslin.R")
#sScriptMFA = file.path("input","maaslin","src","MFA.R")
#sScriptSummarizeMaaslin = file.path("input","maaslin","src","SummarizeMaaslin.R")
#sScriptUtility = file.path("input","maaslin","src","Utility.R")
#sScriptValidateData = file.path("input","maaslin","src","ValidateData.R")

sScriptBoostGLM = "BoostGLM.R"
sScriptConstants = "Constants.R"
sScriptIO = "IO.R"
sScriptMaaslin = "./Maaslin.R"
sScriptMFA = "MFA.R"
sScriptSummarizeMaaslin = "SummarizeMaaslin.R"
sScriptUtility = "Utility.R"
sScriptValidateData = "ValidateData.R"
#Input Files
sTestReadConfig = file.path(c_strTestingDirectory, c_strTestingInput, "TestMaaslin.read.config")
sTestCustomR = file.path(c_strTestingDirectory, c_strTestingInput, "TestMaaslin.R")
sTestMaaslinDirectory = file.path(c_strTestingDirectory, c_strTemporaryFiles, "testMaaslin")
sTestOutput = file.path(sTestMaaslinDirectory,"output.txt")
sTestTSV = file.path(c_strTestingDirectory, c_strTestingInput, "TestMaaslin.tsv")
#Test file answers
sTestOutputAnswer = file.path(c_strTestingDirectory, c_strCorrectAnswers, "TestMaaslin.tsv")

#Delete Test MaAsLin output
unlink(sTestMaaslinDirectory, recursive=TRUE)
#Make neccessary directories
dir.create(sTestMaaslinDirectory)
dir.create(file.path(sTestMaaslinDirectory,"QC"))

sCommand = paste(sScriptMaaslin, "-v", "ERROR", "-d", "0.25", "-r", "0.0001", "-p", "0.1", sTestOutput, sTestTSV, sTestReadConfig, sTestCustomR, sScriptBoostGLM, sScriptConstants, sScriptIO, sScriptMFA, sScriptSummarizeMaaslin, sScriptUtility, sScriptValidateData, sep=" ")
print(sCommand)
system(sCommand)

sExpectedTitle = "\tVariable\tFeature\tValue\tCoefficient\tN\tN.not.0\tP.value\tQ.value"
iExpectedNumberOfLines = 3
lsOutputSummaryFile = readLines(sTestOutput)

test_that("Make sure that the summary output file is what is expected (generally).",{
  expect_equal(lsOutputSummaryFile[1], sExpectedTitle)
  expect_equal(length(lsOutputSummaryFile),iExpectedNumberOfLines)
})

lsDirectoryStructure = list.files(sTestMaaslinDirectory)
lsDirectoryStructureAnswer = c("output.txt","QC","TestMaaslin-age.pdf","TestMaaslin-age.txt","TestMaaslin-dx.txt","TestMaaslin.pdf","TestMaaslin.txt")
test_that("Make sure the expected directory structure is created.",{
  expect_equal(lsDirectoryStructure, lsDirectoryStructureAnswer)
})
