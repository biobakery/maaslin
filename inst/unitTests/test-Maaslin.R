c_strCWD <- getwd()
c_strDir <- file.path(c_strCWD,"..","..")

source(file.path(c_strDir,"R","lib","Constants.R"))
strTestingDirectory = c_strCWD
sScriptMaaslin = file.path( c_strDir, "R","Maaslin.R" )

context("Test Run From Commandline")

#Input Files
sTestReadConfig = file.path(strTestingDirectory, c_strTestingInput, "TestMaaslin.read.config")
sTestTmpDirectory = file.path(strTestingDirectory, c_strTemporaryFiles)
sTestMaaslinDirectory = file.path(sTestTmpDirectory, "testMaaslin")
sTestOutput = file.path(sTestMaaslinDirectory,"TestMaaslin.txt")
sTestOutputDir = file.path(sTestMaaslinDirectory)
sTestTSV = file.path(strTestingDirectory, c_strTestingInput, "TestMaaslin.tsv")
#Test file answers
sTestOutputAnswer = file.path(strTestingDirectory, c_strCorrectAnswers, "TestMaaslin.tsv")

#Delete Test MaAsLin output
unlink(sTestMaaslinDirectory, recursive=TRUE)
#Make neccessary directories
if (!file.exists(sTestTmpDirectory)) {
    dir.create(sTestTmpDirectory)
}
dir.create(sTestMaaslinDirectory)
dir.create(file.path(sTestMaaslinDirectory,"QC"))

setwd(c_strDir)
sCommand = paste(sScriptMaaslin, "-v", "ERROR", "-d", "0.25", "-r", "0.0001", "-p", "0.1", sTestTSV,"-i",sTestReadConfig, sTestOutputDir, sep=" ")
print(sCommand)
system(sCommand)
setwd(c_strCWD)

sExpectedTitle = "\tVariable\tFeature\tValue\tCoefficient\tN\tN.not.0\tP.value\tQ.value"
iExpectedNumberOfLines = 3
lsOutputSummaryFile = readLines(sTestOutput)

test_that("Make sure that the summary output file is what is expected (generally).",{
  expect_equal(lsOutputSummaryFile[1], sExpectedTitle)
  expect_equal(length(lsOutputSummaryFile),iExpectedNumberOfLines)
})

lsDirectoryStructure = list.files(sTestMaaslinDirectory)
lsDirectoryStructureAnswer = c("QC","TestMaaslin-age.pdf","TestMaaslin-age.txt","TestMaaslin-dx.txt","TestMaaslin_log.txt","TestMaaslin.txt")
test_that("Make sure the expected directory structure is created.",{
  expect_equal(sort(lsDirectoryStructure), sort(lsDirectoryStructureAnswer))
})
