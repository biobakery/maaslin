import sfle

Import( "*" )

sRExtention = ".R"
lsTestedScriptNames = []#"SummarizeMaaslin"]#,"ValidateData","Utility","IO","Maaslin"]

for sScriptName in lsTestedScriptNames:
  #Testing summary file
  sTestingSummary = sfle.d(fileDirOutput, sScriptName +"-TestReport.txt")
  sfle.testthat(DefaultEnvironment(), sfle.d(fileDirSrc,sScriptName+sRExtention), sfle.d(fileDirSrc,"test-"+sScriptName), sTestingSummary)
  Default( sTestingSummary ) 

#Start regression suite
execfile( "SConscript_maaslin.py" )

for filePCL in Glob( sfle.d( fileDirInput, "*" + sfle.c_strSufPCL ) ):
	Default( MaAsLin( filePCL ) )
