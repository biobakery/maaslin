import sfle
import csv

Import( "*" )
pE = DefaultEnvironment( )

# Extensions
sGraphlanAnnotationFileExtension = "-ann.txt"
sGraphlanCoreAnnotFileExtension = "-ann-core.txt"
sGraphlanCoreGenesFileExtension = "-core.txt"
sGraphlanFigureExtension = "-graphlan.pdf"
sMaaslinDataFileExtension = ".txt"
sMaaslinReadConfigFileExtension = ".read.config"
sMaaslinSummaryFileExtension = ".txt"

sCustomRScriptExtension = ".R"
sPCLExtension = ".pcl"
sTransposeExtension = ".tsv"

# Files
strMaaslinGraphlanSettings = "Graphlan_settings.txt"

# Script
sScriptGraphlan = File(os.path.join("..","graphlan","graphlan.py"))
sScriptGraphlanAnnotate = File(os.path.join("..","graphlan","graphlan_annotate.py"))
sScriptMaaslinSummaryToGraphlanAnnotation = File(sfle.d(fileDirSrc,"MaaslinToGraphlanAnnotation.py"))
sScriptPCLToCoreGene = File(sfle.d(fileDirSrc,"PCLToGraphlanCoreGene.py"))

sProgMaaslin = sfle.d(fileDirSrc,"Maaslin.R")

# Settings
iGraphlanDPI = 150
iGraphlanFigureSize = 4
iGraphlanPad = 0.2
strGraphlanDirectory = "graphlan"

c_fileDirLib = sfle.d( fileDirSrc, "lib" )
c_fileInputMaaslinR = sfle.d( pE, fileDirSrc, "Maaslin.R" )
c_afileTestsR = [sfle.d( pE, c_fileDirLib, s ) for s in
						("AnalysisModules.R", "IO.R", "SummarizeMaaslin.R", "Utility.R", "ValidateData.R")]

c_afileDocsR = c_afileTestsR + [sfle.d( pE, c_fileDirLib, s ) for s in
						("BoostGLM.R", "Constants.R", "MaaslinPlots.R", "MFA.R")]

##Test scripts
for fileInputR in c_afileTestsR:
  strBase = sfle.rebase( fileInputR, True )
  #Testing summary file
  fileTestingSummary = sfle.d( pE, fileDirOutput, strBase +"-TestReport.txt" )
  dirTestingR = Dir( sfle.d( fileDirSrc, "test-" + strBase ) )
  Default( sfle.testthat( pE, fileInputR, dirTestingR, fileTestingSummary ) )

##Inline doc
for fileProg in c_afileDocsR:
  filePDF = sfle.d( pE, fileDirOutput, sfle.rebase( fileProg, sfle.c_strSufR, sfle.c_strSufPDF ) )
  Default( sfle.inlinedocs( pE, fileProg, filePDF, fileDirTmp ) )

##Start regression suite
execfile( "SConscript_maaslin.py" )

##Input pcl files
lsMaaslinInputFiles = Glob( sfle.d( fileDirInput, "*" + sfle.c_strSufPCL ) )

## Run MaAsLin and generate output
for strPCLFile in lsMaaslinInputFiles:
  Default( MaAsLin( strPCLFile ))
