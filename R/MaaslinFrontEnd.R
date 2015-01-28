MaaslinFrontEnd <-
function(
#*************************************************************************************
#*  	MaaslinFrontEnd                                                              *
#*************************************************************************************
	help = FALSE,							 
	strInputConfig = 'data/maaslin_demo2.read.config',
	dSignificanceLevel =  0.25,				 
	dMinAbd = 0.0001,
	dMinSamp = 0.1,
 	dOutlierFence = 0, 
	dPOutlier = 0.05,		 
	strMultTestCorrection = "BH", 
	fZeroInflated = FALSE, 
	strModelSelection = "boost", 
	strMethod = "lm",
	strTransform = "asinsqrt",
	fNoQC = FALSE,
	strVerbosity = "DEBUG",
	fOmitLogFile = FALSE,
	fInvert = FALSE,
	dSelectionFrequency = NA,
	fAllvAll = FALSE,
	fPlotNA = FALSE,
    dPenalizedAlpha = 0.95,
	sAlternativeLibraryLocation = "/usr/share/biobakery",
	strOutputTXT = "maaslin_demo2.tsv",
	strInputTSV = "outputdir"
	)					

{
	#**********************************************************************
	# Invoke Maaslin                                                      *
	#**********************************************************************
	lsArgs = list()
	lsArgs$opt$help = help
	lsArgs$options$strInputConfig   =  strInputConfig 
 	lsArgs$options$dSignificanceLevel   =  dSignificanceLevel 
	lsArgs$options$dMinAbd   =  dMinAbd 
	lsArgs$options$dMinSamp   =  dMinSamp 
	lsArgs$options$dOutlierFence   =  dOutlierFence 
	lsArgs$options$dPOutlier   =  dPOutlier 
	lsArgs$options$strMultTestCorrection   =  strMultTestCorrection 
	lsArgs$options$fZeroInflated   =  fZeroInflated 
	lsArgs$options$strModelSelection   =  strModelSelection 
	lsArgs$options$strMethod   =  strMethod 
	lsArgs$options$strTransform   =  strTransform 
	lsArgs$options$fNoQC   =  fNoQC 
	lsArgs$options$strVerbosity   =  strVerbosity 
	lsArgs$options$fOmitLogFile   =  fOmitLogFile 
	lsArgs$options$fInvert   =  fInvert 
	lsArgs$options$dSelectionFrequency   =  dSelectionFrequency 
	lsArgs$options$fAllvAll   =  fAllvAll 
	lsArgs$options$fPlotNA   =  fPlotNA 
    lsArgs$options$dPenalizedAlpha   =  dPenalizedAlpha 
	lsArgs$options$sAlternativeLibraryLocation = sAlternativeLibraryLocation
	lsArgs$args[1] =  strOutputTXT 
	lsArgs$args[2] =  strInputTSV 
	return(lsArgs)					

}
