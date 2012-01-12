import math
import os
import sfle
import sys

Import( "*" )

c_strSufPCL		= ".pcl"
c_strSufR		= ".R"
c_strSufRC		= ".read.config"
c_strSufTSV		= ".tsv"
c_strSufTXT		= ".txt"

c_fileProgMaaslin	= File( sfle.d( fileDirSrc, "Maaslin.R" ) )
c_afileProgRs		= [File( sfle.d( fileDirSrc, s ) ) for s in
	("BoostGLM.R", "Constants.R", "IO.R", "MFA.R", "SummarizeMaaslin.R", "Utility.R", "ValidateData.R")]

pE = DefaultEnvironment( )

pProc = sfle.CProcessor( "trn",
	sfle.CTarget( c_strSufPCL ),
	sfle.CCommand( c_fileProgTranspose ),
	sfle.CTarget( c_strSufTSV, fileDirTmp ) )
afileTSVs = sfle.CProcessor.pipeline( pE, pProc, Glob( sfle.d( fileDirInput, "*" ) ) )

def funcMaaslin( target, source, env ):
	strT, astrSs = sfle.ts( target, source )
	strProg, astrArgs = astrSs[0], astrSs[1:]
	return sfle.ex( [strProg, strT] + astrArgs )

for fileTSV in afileTSVs:
	strBase = sfle.rebase( fileTSV, c_strSufTSV )
	fileR, fileRC = (File( sfle.d( fileDirInput, strBase + s ) ) for s in (c_strSufR, c_strSufRC))
	fileTXT = File( sfle.d( fileDirOutput, strBase, strBase + c_strSufTXT ) )
	Default( Command( fileTXT, [c_fileProgMaaslin, fileTSV, fileRC, fileR] + c_afileProgRs, funcMaaslin ) )
