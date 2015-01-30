c_strCWD <- getwd()
c_strDir <- file.path(c_strCWD,"..","..")

source(file.path(c_strDir,"R","lib","Constants.R"))
strTestingDirectory = c_strCWD
sScriptMaaslin = file.path( c_strDir, "R","Maaslin.R" )

context("Test Run From Commandline")

#Input Files
sTestReadConfig = file.path(strTestingDirectory, c_strTestingInput, "TestMaaslin.read.config")
sTestCustomR = file.path(strTestingDirectory, c_strTestingInput, "TestMaaslin.R")
sTestMaaslinDirectory = file.path(strTestingDirectory, c_strTemporaryFiles, "testMaaslin")
sTestOutput = file.path(sTestMaaslinDirectory,"TestMaaslin_Summary.txt")
sTestTSV = file.path(strTestingDirectory, c_strTestingInput, "TestMaaslin.tsv")
#Test file answers
sTestOutputAnswer = file.path(strTestingDirectory, c_strCorrectAnswers, "TestMaaslin.tsv")

#Delete Test MaAsLin output
unlink(sTestMaaslinDirectory, recursive=TRUE)
#Make neccessary directories
dir.create(sTestMaaslinDirectory)
dir.create(file.path(sTestMaaslinDirectory,"QC"))

sCommand = paste(sScriptMaaslin, "-v", "ERROR", "-d", "0.25", "-r", "0.0001", "-p", "0.1", sTestOutput, sTestTSV, sTestReadConfig, sTestCustomR, sep=" ")
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
lsDirectoryStructureAnswer = c(basename(sTestOutput),"QC","TestMaaslin-age.pdf","TestMaaslin-age.txt","TestMaaslin-dx.txt","TestMaaslin.pdf","TestMaaslin.txt")
test_that("Make sure the expected directory structure is created.",{
  expect_equal(sort(lsDirectoryStructure), sort(lsDirectoryStructureAnswer))
})
