import os
import sfle
import sys

c_strSufRC		= ".read.config"

# Ugly!  Makes src relative to current file, not SCons/SConscript
c_fileDirSrc		= Dir( sfle.d( os.path.dirname( sfle.current_file( ) ), sfle.c_strDirSrc ) )
c_fileProgMaaslin	= File( sfle.d( c_fileDirSrc, "Maaslin.R" ) )
c_afileProgRs		= [File( sfle.d( c_fileDirSrc, s ) ) for s in
	("AnalysisModules.R","BoostGLM.R", "Constants.R", "IO.R", "MaaslinPlots.R", "MFA.R", "SummarizeMaaslin.R", "Utility.R", "ValidateData.R")]

def funcMaaslin( target, source, env ):
	strT, astrSs = sfle.ts( target, source )
	strProg, astrArgs = astrSs[0], astrSs[1:]
	return sfle.ex( [strProg, strT] + astrArgs )

def MaAsLin( filePCL, fileDirOut = fileDirOutput, fileDirInt = fileDirTmp ):
	strBase = str(filePCL).replace( sfle.c_strSufPCL, "" )
	fileR, fileRC = (File( strBase + s ) for s in (sfle.c_strSufR, c_strSufRC))
	strBase = os.path.basename( strBase )
	strDir = sfle.d( fileDirOut, strBase )
	fileTXT = File( sfle.d( strDir, strBase + sfle.c_strSufTXT ) )

	fileDirInt = fileDirInt or strDir
	pProc = sfle.CProcessor( None,
		sfle.CTarget( sfle.c_strSufPCL ),
		sfle.CCommand( c_fileProgTranspose ),
		sfle.CTarget( sfle.c_strSufTSV, fileDirInt ) )
	fileTSV = sfle.CProcessor.pipeline( DefaultEnvironment( ), pProc, filePCL )[0]
	return ( Command( fileTXT, [c_fileProgMaaslin, fileTSV, fileRC, fileR] + c_afileProgRs, funcMaaslin ) +
		[fileTSV] )
