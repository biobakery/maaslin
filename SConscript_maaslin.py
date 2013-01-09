#!/usr/bin/env python
"""
Authors: Timothy Tickle and Curtis Huttenhower
Description: Find associations in two matrices of data.
"""

__author__ = "Timothy Tickle and Curtis Huttenhower"
__copyright__ = "Copyright 2012"
__credits__ = ["Timothy Tickle","Curtis Huttenhower"]
__maintainer__ = "Timothy Tickle"
__email__ = "ttickle@hsph.harvard.edu"

import argparse
import os
import sfle
import sys

c_strSufRC		= ".read.config"

c_fileDirSrc		= Dir( sfle.d( os.path.dirname( sfle.current_file( ) ), sfle.c_strDirSrc ) )
c_fileProgMaaslin	= File( sfle.d( c_fileDirSrc, "Maaslin.R" ) )

def MaAsLin( filePCL, dSignificanceThreshold=None,
                      dMinRelativeAbundance=None,
                      dMinPrevalence=None,
                      dOutlierFence=None,
                      strRandomCovariates=None,
                      strModel=None,
                      strAnalysis=None,
                      strLink=None,
                      strForcedCovariates=None,
                      strNoImpute=None,
                      strVerbosity=None,
                      fOmitLogFile=False,
                      fInvert=False,
                      dSelectionFrequency=None,
                      iMFAFeatureCount=None,
                      dMFAMetadataScale=None,
                      dMFADataScale=None,
                      strMFAColorBy=None,
                      strMFAShapeBy=None,
                      lsMFAPlotFeatures=None,
                      fileDirOut=fileDirOutput,
                      fileDirIn=fileDirTmp):
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
			        ( ["-d", dSignificanceThreshold] if dSignificanceThreshold else [] ) + 
			        ( ["-r", dMinRelativeAbundance] if dMinRelativeAbundance else [] ) + 
			        ( ["-p", dMinPrevalence] if dMinPrevalence else [] ) + 
			        ( ["-o", dOutlierFence] if dOutlierFence else [] ) + 
                                ( ["-R", strRandomCovariates] if strRandomCovariates else []) + 
                                ( ["-s", strModel] if strModel else []) + 
                                ( ["-m", strAnalysis] if strAnalysis else []) + 
                                ( ["-l", strLink] if strLink else []) +
                                ( ["-F", strForcedCovariates] if strForcedCovariates else []) + 
                                ( ["-n", strNoImpute] if strAnalysis else []) + 
                                ( ["-v", strVerbosity] if strVerbosity else []) + 
                                ( ["-O", fOmitLogFile] if fOmitLogFile else []) + 
                                ( ["-t", fInvert] if fInvert else []) + 
                                ( ["-f", dSelectionFrequency] if dSelectionFrequency else []) + 
                                ( ["-c", iMFAFeatureCount] if iMFAFeatureCount else []) + 
                                ( ["-M", dMFAMetadataScale] if dMFAMetadataScale else []) + 
                                ( ["-D", dMFADataScale] if dMFADataScale else []) + 
                                ( ["-C", strMFAColorBy] if strMFAColorBy else []) + 
                                ( ["-S", strMFAShapeBy] if strMFAShapeBy else []) + 
                                ( ["-P", lsMFAPlotFeatures] if lsMFAPlotFeatures else []) +  
                                  [strT, strTSV])
	
	strBase = os.path.basename( strBase )
	strDir = sfle.d( fileDirOut, strBase )
	fileTXT = File( sfle.d( strDir, strBase + sfle.c_strSufTXT ) )

	fileDirIn = fileDirIn or strDir
	pProc = sfle.CProcessor( None,
		sfle.CTarget( sfle.c_strSufPCL ),
		sfle.CCommand( c_fileProgTranspose ),
		sfle.CTarget( sfle.c_strSufTSV, fileDirIn ) )
	fileTSV = sfle.CProcessor.pipeline( DefaultEnvironment( ), pProc, filePCL )[0]
	return ( Command( fileTXT, [c_fileProgMaaslin, fileTSV] +
		( [fileRC] if fileRC else [] ) + ( [fileR] if fileR else [] ), funcMaaslin ) + [fileTSV] )