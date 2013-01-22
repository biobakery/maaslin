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
sArgsExt = ".args"
#Commandline to ignore
lsIgnore = ["-i","-I","--input_config","--input_process"]

def MaAsLin( filePCL ):
	#Build input file name if they exist or give ""
	strBase = filePCL.get_abspath().replace( sfle.c_strSufPCL, "" )
	strR, strRC, strArgs = (( strBase + s ) for s in (sfle.c_strSufR, c_strSufRC, sArgsExt))
	fileR, fileRC, fileArgs = (( File( s ) if os.path.exists( s ) else "" ) for s in (strR, strRC, strArgs))

	## Read in an args file if it exists
	lsArgs = []
	if fileArgs:
		fReader = csv.reader(open(fileArgs.get_abspath(),'r'), delimiter = " ")
		lsArgsTmp = []
		[lsArgsTmp.extend(lsLine) for lsLine in fReader]
		fSkip = False
		for s in lsArgsTmp:
			if s in lsIgnore:
				fSkip=True
				continue
			if fSkip:
				fSkip = not fSkip
				continue
			lsArgs.append(s)

	lsInputArgs = ["-I",[fileR]] if fileR else []
	lsInputArgs.extend(["-i",[fileRC]] if fileRC else [])
	lsArgs.extend(lsInputArgs)

	strBase = os.path.basename( strBase )
	fileTSVFile = File(sfle.d(fileDirTmp,sfle.rebase(filePCL,sfle.c_strSufPCL,sfle.c_strSufTSV)))
	strT = File( sfle.d( os.path.join(fileDirOutput.get_abspath(), strBase, strBase + sfle.c_strSufTXT) ) )

	#Transpose PCL
	sfle.spipe(pE, filePCL, c_fileProgTranspose, fileTSVFile)
	#Run MaAsLin
	sfle.op(pE, c_fileProgMaaslin, lsArgs+[[True,strT],[False, fileTSVFile]])
	if fileArgs: Depends(c_fileProgMaaslin, fileArgs)
	Default(strT)