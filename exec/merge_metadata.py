#!/usr/bin/env python
#####################################################################################
#Copyright (C) <2012>
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in the
#Software without restriction, including without limitation the rights to use, copy,
#modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
#and to permit persons to whom the Software is furnished to do so, subject to
#the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies
#or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# This file is a component of the MaAsLin (Multivariate Associations Using Linear Models), 
# authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Timothy Tickle, ttickle@hsph.harvard.edu).
#####################################################################################
"""
Examples
~~~~~~~~

``metadata.txt``::

	-	Y	Z
	a	1	x
	b	0	y
	c		z

``data.pcl``::

	-	a	b	c
	A|B	1	2	3
	A|C	4	5	6
	D|E	7	8	9

``Examples``::

	$ merge_metadata.py metadata.txt < data.pcl
	sample	a	b	c
	Y	1	0
	Z	x	y	z
	A	0.416667	0.466667	0.5
	A|B	0.0833333	0.133333	0.166667
	A|C	0.333333	0.333333	0.333333
	D|E	0.583333	0.533333	0.5

	$ merge_metadata.py metadata.txt -t 0 < data.pcl
	sample	a	b	c
	Y	1	0
	Z	x	y	z
	A|B	0.0833333	0.133333	0.166667
	A|C	0.333333	0.333333	0.333333
	D|E	0.583333	0.533333	0.5

	$ merge_metadata.py metadata.txt -t 1 < data.pcl
	sample	a	b	c
	Y	1	0
	Z	x	y	z
	A	0.416667	0.466667	0.5
	D	0.583333	0.533333	0.5

	$ merge_metadata.py metadata.txt -t 0 -n < data.pcl
	sample	a	b	c
	Y	1	0
	Z	x	y	z
	A|B	1	2	3
	A|C	4	5	6
	D|E	7	8	9

	$ merge_metadata.py metadata.txt -t 0 -m 0.8 -s "-" < data.pcl
	sample	b	c
	Y	0	-
	Z	y	z
	A|B	0.133333	0.166667
	A|C	0.333333	0.333333
	D|E	0.533333	0.5

	$ merge_metadata.py -t 0 < data.pcl
	sample	a	b	c
	A|B	1	2	3
	A|C	4	5	6
	D|E	7	8	9

.. testsetup::

	from merge_metadata import *
"""

import argparse
import blist
import csv
import re
import sys

c_dTarget	= 1.0
c_fRound	= False

class CClade:
	
	def __init__( self ):
		
		self.m_hashChildren = {}
		self.m_adValues = None
	
	def get( self, astrClade ):
		
		return self.m_hashChildren.setdefault(
			astrClade[0], CClade( ) ).get( astrClade[1:] ) if astrClade else self
			
	def set( self, adValues ):
		
		self.m_adValues = blist.blist( [0] ) * len( adValues )
		for i, d in enumerate( adValues ):
			if d:
				self.m_adValues[i] = d
		
	def impute( self ):
		
		if not self.m_adValues:
			for pChild in self.m_hashChildren.values( ):
				adChild = pChild.impute( )
				if self.m_adValues:
					for i in range( len( adChild or [] ) ):
						if adChild[i]:
							self.m_adValues[i] += adChild[i]
				elif adChild:
					self.m_adValues = adChild[:] 
					
		return self.m_adValues
	
	def _freeze( self, hashValues, iTarget, astrClade, iDepth, fLeaves ):
		
		fHit = ( not iTarget ) or ( ( fLeaves and ( iDepth == iTarget ) ) or ( ( not fLeaves ) and ( iDepth <= iTarget ) ) )
		iDepth += 1
		setiRet = set()
		if self.m_hashChildren:
			for strChild, pChild in self.m_hashChildren.items( ):
				setiRet |= pChild._freeze( hashValues, iTarget, astrClade + [strChild], iDepth, fLeaves )
			setiRet = set( ( i + 1 ) for i in setiRet )
		else:
			setiRet.add( 0 )
		if iTarget < 0:
			if fLeaves:
				fHit = -( iTarget + 1 ) in setiRet
			else:
				fHit = -( iTarget + 1 ) <= max( setiRet )
		if astrClade and self.m_adValues and fHit:
			hashValues["|".join( astrClade )] = self.m_adValues
		return setiRet
	
	def freeze( self, hashValues, iTarget, fLeaves ):
		
		self._freeze( hashValues, iTarget, [], 0, fLeaves )
	
	def _repr( self, strClade ):

		strRet = "<"
		if strClade:
			strRet += "%s %s" % (strClade, self.m_adValues)
			if self.m_hashChildren:
				strRet += " "
		if self.m_hashChildren:
			strRet += " ".join( p._repr( s ) for (s, p) in self.m_hashChildren.items( ) )
		
		return ( strRet + ">" )
		
	def __repr__( self ):
		
		return self._repr( "" )

"""
pTree = CClade( )
pTree.get( ("A", "B") ).set( [1, 2, 3] )
pTree.get( ("A", "C") ).set( [4, 5, 6] )
pTree.get( ("D", "E") ).set( [7, 8, 9] )
iTaxa = 0
if iTaxa:
	pTree.impute( )
hashFeatures = {}
pTree.freeze( hashFeatures, iTaxa )
print( pTree )
print( hashFeatures )
sys.exit( 0 )
#"""

def merge_metadata( aastrMetadata, aastrData, ostm, fNormalize, strMissing, astrExclude, dMin, iTaxa, fLeaves ):
	"""
	Joins and outputs a data matrix with a metadata matrix, optionally normalizing and filtering it.
	A pipe-delimited taxonomy hierarchy can also be dynamically added or removed.
	
	:param	aastrMetadata:	Split lines from which metadata are read.
	:type	aastrMetadata:	collection of string collections
	:param	aastrData:		Split lines from which data are read.
	:type	aastrData:		collection of string collections
	:param	ostm:			Output stream to which joined rows are written.
	:type	ostm:			output stream
	:param	fNormalize:		If true, divide data values by column sums.
	:type	fNormalize:		bool
	:param	strMissing:		Representation for missing metadata values.
	:type	strMissing:		str
	:param	astrExclude:	Lines from which excluded IDs are read.
	:type	astrExclude:	collection of strings
	:param	dMin:			Minimum fraction of maximum value for per-column quality control.
	:type	dMin:			bool
	:param	iTaxa:			Depth of taxonomy to be computed, -1 = leaves only, 0 = no change
	:type	iTaxa:			int
	:param	fLeaves:		Output only leaves, not complete taxonomy; ignored if taxa = 0
	:type	fLeaves:		bool

	Metadata are optional; if not provided, data will be optionally normalized or its taxonomy
	modified as requested.  Metadata are provided one row per sample, data one column per
	sample, both files tab-delimited text with one header row and one header column.
	
	Metadata IDs that do not match data IDs are discarded, and data IDs without corresponding
	metadata IDs are given missing values.  Missing data values are always treated (and output)
	as zero.
	
	Per-column quality control is performed if the requested minimum fraction is greater than
	zero.  Specifically, for each column i, the row j containing the maximum value d is
	identified.  If d is less than the minimum fraction of row j's maximum value over all columns,
	the entire column i is removed.
	
	A taxonomy hierarchy will be calculated by default if row IDs are pipe-delimited, i.e. of
	the form A|B|C.  All parent clades are computed by default, e.g. A|B and A, save when
	they would be identical to a more specific child clade.  Negative values are counted from the
	bottom (right) of the hierarchy rather than the top.  The special value of 0 deactivates
	hierarchy calculation.
	
	>>> aastrMetadata = [[t.strip( ) for t in s] for s in ("-YZ", "a1x", "b0y", "c z")]
	>>> aastrData = [s.split( ) for s in ( \
		"-	a	b	c",		\
		"A|B	1	2	3",	\
		"A|C	4	5	6",	\
		"D|E	7	8	9")]
	>>> merge_metadata( aastrMetadata, aastrData, sys.stdout, True, "", [], 0.01, -1, False ) #doctest: +NORMALIZE_WHITESPACE
	sample	a	b	c
	Y	1	0
	Z	x	y	z
	A	0.416667	0.466667	0.5
	A|B	0.0833333	0.133333	0.166667
	A|C	0.333333	0.333333	0.333333
	D|E	0.583333	0.533333	0.5

	>>> merge_metadata( aastrMetadata, aastrData, sys.stdout, True, "", [], 0.01, -1, True ) #doctest: +NORMALIZE_WHITESPACE
	sample	a	b	c
	Y	1	0
	Z	x	y	z
	A|B	0.0833333	0.133333	0.166667
	A|C	0.333333	0.333333	0.333333
	D|E	0.583333	0.533333	0.5

	>>> merge_metadata( aastrMetadata, aastrData, sys.stdout, True, "", [], 0, 0, True ) #doctest: +NORMALIZE_WHITESPACE
	sample	a	b	c
	Y	1	0
	Z	x	y	z
	A|B	0.0833333	0.133333	0.166667
	A|C	0.333333	0.333333	0.333333
	D|E	0.583333	0.533333	0.5

	>>> merge_metadata( aastrMetadata, aastrData, sys.stdout, True, "", [], 0, 1, False ) #doctest: +NORMALIZE_WHITESPACE
	sample	a	b	c
	Y	1	0
	Z	x	y	z
	A	0.416667	0.466667	0.5
	D	0.583333	0.533333	0.5

	>>> merge_metadata( aastrMetadata, aastrData, sys.stdout, True, "", [], 0, -1, True ) #doctest: +NORMALIZE_WHITESPACE
	sample	a	b	c
	Y	1	0
	Z	x	y	z
	A|B	0.0833333	0.133333	0.166667
	A|C	0.333333	0.333333	0.333333
	D|E	0.583333	0.533333	0.5

	>>> merge_metadata( aastrMetadata, aastrData, sys.stdout, False, "", [], 0, 0, True ) #doctest: +NORMALIZE_WHITESPACE
	sample	a	b	c
	Y	1	0
	Z	x	y	z
	A|B	1	2	3
	A|C	4	5	6
	D|E	7	8	9

	>>> merge_metadata( aastrMetadata, aastrData, sys.stdout, True, "-", [], 0.8, 0, True ) #doctest: +NORMALIZE_WHITESPACE
	sample	b	c
	Y	0	-
	Z	y	z
	A|B	0.133333	0.166667
	A|C	0.333333	0.333333
	D|E	0.533333	0.5

	>>> merge_metadata( None, aastrData, sys.stdout, False, "", [], 0, 0, True ) #doctest: +NORMALIZE_WHITESPACE
	sample	a	b	c
	A|B	1	2	3
	A|C	4	5	6
	D|E	7	8	9

	>>> merge_metadata( aastrMetadata, aastrData, sys.stdout, True, "", ["b"], 0.01, -1, False ) #doctest: +NORMALIZE_WHITESPACE
	sample	a	c
	Y	1
	Z	x	z
	A	0.416667	0.5
	A|B	0.0833333	0.166667
	A|C	0.333333	0.333333
	D|E	0.583333	0.5
	"""

        #Put metadata in a dictionary
        #{"First line element",["line element 2","line element 3","line element 4"]}
        #If there is no metadata then 
	astrMetadata = None
	hashMetadata = {}
	for astrLine in ( aastrMetadata or [] ):
		if astrMetadata:
			hashMetadata[astrLine[0]] = astrLine[1:]
		else:
			astrMetadata = astrLine[1:]
	
	astrHeaders = adSeqs = iCol = None
	pTree = CClade( )
	aastrRaw = []
	for astrLine in aastrData:
		if astrHeaders:
			if ( astrLine[0] == "EWEIGHT" ) or ( astrLine[0] == "total" ) or \
				( len( astrLine ) < 2 ):
				continue
			try:
				adCounts = [( float(strCur) if len( strCur.strip( ) ) else 0 ) for
					strCur in astrLine[iCol:]]
			except ValueError:
				aastrRaw.append( astrLine )
				continue
			for i in range( len( adCounts ) ):
				adSeqs[i] += adCounts[i]
			if ( iCol > 1 ) and ( astrLine[0] != astrLine[1] ):
				if astrLine[1].find( astrLine[0] ) >= 0:
					astrLine[0] = astrLine[1]
				else:
					astrLine[0] += " " + astrLine[1]
			pTree.get( astrLine[0].split( "|" ) ).set( adCounts )
		else:
			iCol = 2 if ( astrLine[1].upper( ) == "NAME" ) else 1
			astrHeaders = [strCur.replace( " ", "_" ) for strCur in astrLine[iCol:]]
			adSeqs = [0] * len( astrHeaders )
			
	if iTaxa:
		pTree.impute( )
	hashFeatures = {}
	pTree.freeze( hashFeatures, iTaxa, fLeaves )
	setstrFeatures = hashFeatures.keys( )
	
	afOmit = [False] * len( astrHeaders )
	if dMin > 0:
		aadData = list(hashFeatures.values( ))
		for i in range( len( astrHeaders ) ):
			iMax = max( range( len( aadData ) ), key = lambda j: aadData[j][i] )
			dMaxUs = aadData[iMax][i]
			dMaxThem = max( aadData[iMax][j] for j in ( range( i ) + range( i + 1, len( astrHeaders ) ) ) )
			if dMaxUs < ( dMin * dMaxThem ):
				sys.stderr.write( "Omitting: %s\n" % astrHeaders[i] )
				afOmit[i] = True
	
	if astrExclude:
		setstrExclude = set(s.strip( ) for s in astrExclude)
		for i in range( len( astrHeaders ) ):
			if ( not afOmit[i] ) and ( astrHeaders[i] in setstrExclude ):
				afOmit[i] = True
	
	adMult = [( ( c_dTarget / d ) if ( fNormalize and ( d > 0 ) ) else 1 ) for d in adSeqs]
	for strFeature, adCounts in hashFeatures.items( ):
		for i in range( len( adCounts ) ):
			if adCounts[i]:
				adCounts[i] *= adMult[i]
				if c_fRound:
					adCounts[i] = round( adCounts[i] )
	for strFeature, adCounts in hashFeatures.items( ):
		astrFeature = strFeature.strip( ).split( "|" )
		while len( astrFeature ) > 1:
			astrFeature = astrFeature[:-1]
			strParent = "|".join( astrFeature )
			adParent = hashFeatures.get( strParent )
			if adParent == adCounts:
				del hashFeatures[strParent]
				setstrFeatures.remove( strParent )
	
	if astrMetadata:
		for i in range( len( astrMetadata ) ):
			hashFeatures[astrMetadata[i]] = astrCur = []
			for strSubject in astrHeaders:
				astrSubject = hashMetadata.get( strSubject )
				if not astrSubject:
					strSubject = re.sub( '_.*$', "", strSubject )
					astrSubject = hashMetadata.get( strSubject, [] )
				astrCur.append( astrSubject[i] if ( i < len( astrSubject ) ) else "" )
	
	astrFeatures = sorted( astrMetadata or [] ) + sorted( setstrFeatures )
	aiHeaders = filter( lambda i: not afOmit[i], range( len( astrHeaders ) ) )
	csvw = csv.writer( sys.stdout, csv.excel_tab )
	csvw.writerow( ["sample"] + [astrHeaders[i] for i in aiHeaders] )
	for iFeature in range( len( astrFeatures ) ):
		strFeature = astrFeatures[iFeature]
		adFeature = hashFeatures[strFeature]
		astrValues = [adFeature[i] for i in aiHeaders]
		for i in range( len( astrValues ) ):
			strValue = astrValues[i]
			if type( strValue ) in (int, float):
				astrValues[i] = "%g" % astrValues[i]
			elif ( not strValue ) or ( ( type( strValue ) == str ) and
				( len( strValue ) == 0 ) ):
				astrValues[i] = strMissing
		csvw.writerow( [strFeature] + astrValues )

	for astrRaw in aastrRaw:
		csvw.writerow( [astrRaw[i] for i in aiHeaders] )

argp = argparse.ArgumentParser( prog = "merge_metadata.py",
	description = "Join a data matrix with a metadata matrix, optionally normalizing and filtering it.\n\n" +
	"A pipe-delimited taxonomy hierarchy can also be dynamically added or removed." )
argp.add_argument( "-n",		dest = "fNormalize",	action = "store_false",
	help = "Don't normalize data values by column sums" )
argp.add_argument( "-s",		dest = "strMissing",	metavar = "missing",
	type = str,		default = " ",
	help = "String representing missing metadata values" )
argp.add_argument( "-m",		dest = "dMin",			metavar = "min",
	type = float,	default = 0.01,
	help = "Per-column quality control, minimum fraction of maximum value" )
argp.add_argument( "-t",		dest = "iTaxa",			metavar = "taxa",
	type = int,		default = -1,
	help = "Depth of taxonomy to be computed, negative = from right, 0 = no change" )
argp.add_argument( "-l",		dest = "fLeaves",		action = "store_true",
	help = "Output only leaves, not complete taxonomy" )
argp.add_argument( "-x",		dest = "istmExclude",	metavar = "exclude.txt",
	type = file,
	help = "File from which sample IDs to exclude are read" )
argp.add_argument( "istmMetadata",	metavar = "metadata.txt",
	type = file,	nargs = "?",
	help = "File from which metadata is read" )
__doc__ = "::\n\n\t" + argp.format_help( ).replace( "\n", "\n\t" ) + __doc__

def _main( ):
	args = argp.parse_args( )
	merge_metadata( args.istmMetadata and csv.reader( args.istmMetadata, csv.excel_tab ),
		csv.reader( sys.stdin, csv.excel_tab ), sys.stdout, args.fNormalize, args.strMissing,
			args.istmExclude, args.dMin, args.iTaxa, args.fLeaves )
	
if __name__ == "__main__":
	_main( )
