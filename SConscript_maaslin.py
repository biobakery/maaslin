import os
import sfle
import sys

c_strSufRC		= ".read.config"

# Ugly!  Makes src relative to current file, not SCons/SConscript
c_fileDirSrc		= Dir( sfle.d( os.path.dirname( sfle.current_file( ) ), sfle.c_strDirSrc ) )
c_fileProgMaaslin	= File( sfle.d( c_fileDirSrc, "Maaslin.R" ) )

def MaAsLin( filePCL, strModel=None, strAnalysis=None, strLink=None, strForcedCovariates=None, strRandomCovariates=None, fileDirOut = fileDirOutput, fileDirInt = fileDirTmp, fOmitLog=False ):
	strBase = str(filePCL).replace( sfle.c_strSufPCL, "" )
	strR, strRC = (( strBase + s ) for s in (sfle.c_strSufR, c_strSufRC))
	fileR, fileRC = (( File( s ) if os.path.isfile( s ) else "" ) for s in (strR, strRC))

	# Ugly - make this use sfle.op or similar instead
	def funcMaaslin( target, source, env, fileR = fileR, fileRC = fileRC, strRandomCovariates=strRandomCovariates ):
		strT, astrSs = sfle.ts( target, source )
		strProg, strTSV = astrSs[:2]
		return sfle.ex( [strProg] + 
                                ( ["-i", fileRC] if fileRC else [] ) +
			        ( ["-I", fileR] if fileR else [] ) + 
                                ( ["-s", strModel] if strModel else []) + 
                                ( ["-m", strAnalysis] if strAnalysis else []) + 
                                ( ["-l", strLink] if strLink else []) + 
                                ( ["-F", strForcedCovariates] if strForcedCovariates else []) + 
                                ( ["-R", strRandomCovariates] if strRandomCovariates else []) + 
                                ( ["-O"] if fOmitLog else []) + 
                                  [strT, strTSV])
	
	strBase = os.path.basename( strBase )
	strDir = sfle.d( fileDirOut, strBase )
	fileTXT = File( sfle.d( strDir, strBase + sfle.c_strSufTXT ) )

	fileDirInt = fileDirInt or strDir
	pProc = sfle.CProcessor( None,
		sfle.CTarget( sfle.c_strSufPCL ),
		sfle.CCommand( c_fileProgTranspose ),
		sfle.CTarget( sfle.c_strSufTSV, fileDirInt ) )
	fileTSV = sfle.CProcessor.pipeline( DefaultEnvironment( ), pProc, filePCL )[0]
	return ( Command( fileTXT, [c_fileProgMaaslin, fileTSV] +
		( [fileRC] if fileRC else [] ) + ( [fileR] if fileR else [] ), funcMaaslin ) + [fileTSV] )
