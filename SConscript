import sfle

Import( "*" )

sRExtention = ".R"
lsTestedScriptNames = ["SummarizeMaaslin","ValidateData","Utility","IO","Maaslin"]
lsInlineDocs = [sfle.d(fileDirSrc,sFile) for sFile in ["Maaslin.R","AnalysisModules.R","BoostGLM.R","IO.R","MaaslinPlots.R","MFA.R","SummarizeMaaslin.R","Utility.R","ValidateData.R"]]

#Test scripts
for sScriptName in lsTestedScriptNames:
  #Testing summary file
  sTestingSummary = sfle.d(fileDirOutput, sScriptName +"-TestReport.txt")
  sfle.testthat(DefaultEnvironment(), sfle.d(fileDirSrc,sScriptName+sRExtention), sfle.d(fileDirSrc,"test-"+sScriptName), sTestingSummary)
  Default( sTestingSummary )

#Inline doc
for fileProg in lsInlineDocs:
	filePDF = sfle.d( pE, fileDirOutput, sfle.rebase( fileProg, sRExtention, sfle.c_strSufPDF ) )
	sfle.inlinedocs( pE, fileProg, filePDF, fileDirTmp )
	Default( filePDF )

#Start regression suite
execfile( "SConscript_maaslin.py" )

for filePCL in Glob( sfle.d( fileDirInput, "*" + sfle.c_strSufPCL ) ):
	Default( MaAsLin( filePCL ) )
