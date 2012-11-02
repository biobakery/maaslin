import sfle

Import( "*" )
pE = DefaultEnvironment( )

c_fileDirLib		= sfle.d( fileDirSrc, "lib" )
c_fileInputMaaslinR	= sfle.d( pE, fileDirSrc, "Maaslin.R" )
c_afileTestsR		= [sfle.d( pE, c_fileDirLib, s ) for s in
						("IO.R", "SummarizeMaaslin.R", "Utility.R", "ValidateData.R")]

c_afileDocsR		= c_afileTestsR + [sfle.d( pE, c_fileDirLib, s ) for s in
						("AnalysisModules.R", "BoostGLM.R", "MaaslinPlots.R", "MFA.R")]

#Test scripts
for fileInputR in c_afileTestsR:
	strBase = sfle.rebase( fileInputR, True )
	#Testing summary file
	fileTestingSummary = sfle.d( pE, fileDirOutput, strBase +"-TestReport.txt" )
	dirTestingR = Dir( sfle.d( fileDirSrc, "test-" + strBase ) )
	Default( sfle.testthat( pE, fileInputR, dirTestingR, fileTestingSummary ) )

#Inline doc
for fileProg in c_afileDocsR:
	filePDF = sfle.d( pE, fileDirOutput, sfle.rebase( fileProg, sfle.c_strSufR, sfle.c_strSufPDF ) )
	Default( sfle.inlinedocs( pE, fileProg, filePDF, fileDirTmp ) )

#Start regression suite
execfile( "SConscript_maaslin.py" )

for filePCL in Glob( sfle.d( fileDirInput, "*" + sfle.c_strSufPCL ) ):
	Default( MaAsLin( filePCL ) )
