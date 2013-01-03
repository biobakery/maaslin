import sfle

Import( "*" )
pE = DefaultEnvironment( )

# Extensions
sGraphlanAnnotationFileExtension = "-ann.txt"
sGraphlanCoreAnnotFileExtension = "-ann-core.txt"
sGraphlanCoreGenesFileExtension = "-core.txt"
sGraphlanFigureExtension = "-graphlan.pdf"
sMaaslinDataFileExtension = ".txt"
sMaaslinReadConfigFileExtension = ".read.config"
sMaaslinSummaryFileExtension = ".txt"

sCustomRScriptExtension = ".R"
sPCLExtension = ".pcl"
sTransposeExtension = ".tsv"

# Files
strMaaslinGraphlanSettings = "Graphlan_settings.txt"

# Script
sScriptGraphlan = File(os.path.join("..","graphlan","graphlan.py"))
sScriptGraphlanAnnotate = File(os.path.join("..","graphlan","graphlan_annotate.py"))
sScriptMaaslinSummaryToGraphlanAnnotation = File(sfle.d(fileDirSrc,"MaaslinToGraphlanAnnotation.py"))
sScriptPCLToCoreGene = File(sfle.d(fileDirSrc,"PCLToGraphlanCoreGene.py"))

sProgMaaslin = sfle.d(fileDirSrc,"Maaslin.R")

# Settings
iGraphlanDPI = 150
iGraphlanFigureSize = 4
iGraphlanPad = 0.2
strGraphlanDirectory = "graphlan"

c_fileDirLib = sfle.d( fileDirSrc, "lib" )
c_fileInputMaaslinR = sfle.d( pE, fileDirSrc, "Maaslin.R" )
c_afileTestsR = [sfle.d( pE, c_fileDirLib, s ) for s in
						("AnalysisModules.R", "BoostGLM.R", "IO.R", "SummarizeMaaslin.R", "Utility.R", "ValidateData.R")]

c_afileDocsR = c_afileTestsR + [sfle.d( pE, c_fileDirLib, s ) for s in
						("Constants.R", "MaaslinPlots.R", "MFA.R")]

##Test scripts
for fileInputR in c_afileTestsR:
  strBase = sfle.rebase( fileInputR, True )
  #Testing summary file
  fileTestingSummary = sfle.d( pE, fileDirOutput, strBase +"-TestReport.txt" )
  dirTestingR = Dir( sfle.d( fileDirSrc, "test-" + strBase ) )
  Default( sfle.testthat( pE, fileInputR, dirTestingR, fileTestingSummary ) )

##Inline doc
for fileProg in c_afileDocsR:
  filePDF = sfle.d( pE, fileDirOutput, sfle.rebase( fileProg, sfle.c_strSufR, sfle.c_strSufPDF ) )
  Default( sfle.inlinedocs( pE, fileProg, filePDF, fileDirTmp ) )

##Start regression suite
execfile( "SConscript_maaslin.py" )

##Input pcl files
lsMaaslinInputFiles = Glob( sfle.d( fileDirInput, "*" + sfle.c_strSufPCL ) )

## Run MaAsLin and generate output
for strPCLFile in lsMaaslinInputFiles:

  ## Run MaAsLin
  Default( MaAsLin( strPCLFile))

  ## Run Graphlan on all output projects
  strProjectName = os.path.splitext(os.path.split(strPCLFile.get_abspath())[1])[0]
  strMaaslinOutputDir = sfle.d(fileDirOutput,strProjectName)

  ##Get maaslin data files
  strMaaslinSummaryFile = sfle.d(os.path.join(strMaaslinOutputDir, strProjectName + sMaaslinSummaryFileExtension))

  # Make annotation file
  sAnnotationFile = File(sfle.d(strMaaslinOutputDir, os.path.join(strGraphlanDirectory,sfle.rebase(strMaaslinSummaryFile, sMaaslinSummaryFileExtension,sGraphlanAnnotationFileExtension))))
  sfle.op(pE, sScriptMaaslinSummaryToGraphlanAnnotation, [[False, strMaaslinSummaryFile],[False,File(sfle.d(fileDirSrc,strMaaslinGraphlanSettings))],[True,sAnnotationFile]])

  # Make core gene file
  sCoreGeneFile = File(sfle.d(strMaaslinOutputDir,  os.path.join(strGraphlanDirectory,sfle.rebase(strMaaslinSummaryFile, sMaaslinSummaryFileExtension,sGraphlanCoreGenesFileExtension))))
  sReadConfigFile = File(sfle.d(fileDirInput,sfle.rebase(strMaaslinSummaryFile, sMaaslinSummaryFileExtension,sMaaslinReadConfigFileExtension)))

  sfle.op(pE, sScriptPCLToCoreGene, [[False, strPCLFile],[False, sReadConfigFile],[True, sCoreGeneFile]])

  # Generate core gene annotation file names
  sCoreGeneAnnotationFile = File(sfle.d(strMaaslinOutputDir,  os.path.join(strGraphlanDirectory,sfle.rebase(strMaaslinSummaryFile, sMaaslinSummaryFileExtension,sGraphlanCoreAnnotFileExtension))))

  sfle.op(pE, sScriptGraphlanAnnotate, ["--annot",[sAnnotationFile],[False, sCoreGeneFile],[True, sCoreGeneAnnotationFile]])

  # Call graphlan
  # graphlan.py --dpi 150 --size 4 --pad 0.2 core_genes.annot.xml core_genes.png
  sGraphlanFigure = File(sfle.d(strMaaslinOutputDir,  os.path.join(strGraphlanDirectory, sfle.rebase(strMaaslinSummaryFile, sMaaslinSummaryFileExtension,sGraphlanFigureExtension))))

  sfle.op(pE, sScriptGraphlan, [[False, sCoreGeneAnnotationFile],[True, sGraphlanFigure],"--dpi",iGraphlanDPI,"--size",iGraphlanFigureSize,"--pad",iGraphlanPad])

  Default(sCoreGeneFile,sAnnotationFile,sCoreGeneAnnotationFile,sGraphlanFigure)

if False:
  # Permuting the following:
  #1# Model Selection 4 (Boost, forwards select, backwards selection, none)
  #2# Analysis 5 (LM, Neg-binomial, quasi-poisson, spearman, wilcoxon) Add Lasso
  #3# Link Functions 2 (arcsin-sqrt, none)
  #4# Forced 3 (no covariate, One covariate, Two covariates)
  #5# Random 3 (no random, One random, two random covariate)
 # Also testing:
  ## Vanilla MaAsLin run R files (with R file, without R file)
  ## Vanilla MaAsLin run read config (with read.config, without read.config)

  #Permutations
  # 1 Model: Boost, 2 Analysis: LM, 3 Link: Arcsin-sqrt, 4 Forced: None, 5 Random: None
  sPCLFile = File(sfle.d(fileDirInput,"check", "Check-Boost-LM-AS-None-None.pcl")).get_abspath()
  sTSVFile = sfle.d(fileDirTmp,sfle.rebase(sPCLFile, sPCLExtension, sTransposeExtension ))
  sRCFile = File(sfle.d(fileDirInput,os.path.join("check",sfle.rebase(sPCLFile, sPCLExtension, sMaaslinReadConfigFileExtension )))).get_abspath()
  sOutputFile = sfle.d(os.path.join(fileDirOutput.get_abspath(),"check", os.path.splitext(os.path.basename(sPCLFile))[0], sfle.rebase(sPCLFile, sPCLExtension, sMaaslinSummaryFileExtension)))
  sfle.spipe(pE, sPCLFile, c_fileProgTranspose, sTSVFile)
  Default(sfle.op( pE, sProgMaaslin, [[True, sOutputFile], [False, sTSVFile], "-i", sRCFile]))
  
  # 2 Model: Boost, 2 Analysis: LM, 3 Link: Arcsin-sqrt, 4 Forced: None, 5 Random: One
  sPCLFile = File(sfle.d(fileDirInput,"check", "Check-Boost-LM-AS-None-One.pcl")).get_abspath()
  sTSVFile = sfle.d(fileDirTmp,sfle.rebase(sPCLFile, sPCLExtension, sTransposeExtension ))
  sRCFile = File(sfle.d(fileDirInput,os.path.join("check",sfle.rebase(sPCLFile, sPCLExtension, sMaaslinReadConfigFileExtension )))).get_abspath()
  sOutputFile = sfle.d(os.path.join(fileDirOutput.get_abspath(),"check", os.path.splitext(os.path.basename(sPCLFile))[0], sfle.rebase(sPCLFile, sPCLExtension, sMaaslinSummaryFileExtension)))
  sfle.spipe(pE, sPCLFile, c_fileProgTranspose, sTSVFile)
  Default(sfle.op( pE, sProgMaaslin, [[True, sOutputFile], [False, sTSVFile], "-i", sRCFile, "-R", "Cohort"]))

  # 3 Model: Boost, 2 Analysis: LM, 3 Link: Arcsin-sqrt, 4 Forced: None, 5 Random: Two
  sPCLFile = File(sfle.d(fileDirInput,"check", "Check-Boost-LM-AS-None-Two.pcl")).get_abspath()
  sTSVFile = sfle.d(fileDirTmp,sfle.rebase(sPCLFile, sPCLExtension, sTransposeExtension ))
  sRCFile = File(sfle.d(fileDirInput,os.path.join("check",sfle.rebase(sPCLFile, sPCLExtension, sMaaslinReadConfigFileExtension )))).get_abspath()
  sOutputFile = sfle.d(os.path.join(fileDirOutput.get_abspath(),"check", os.path.splitext(os.path.basename(sPCLFile))[0], sfle.rebase(sPCLFile, sPCLExtension, sMaaslinSummaryFileExtension)))
  sfle.spipe(pE, sPCLFile, c_fileProgTranspose, sTSVFile)
  Default(sfle.op( pE, sProgMaaslin, [[True, sOutputFile], [False, sTSVFile], "-i", sRCFile, "-R", "Cohort,Sex"]))

  # 4 Model: Boost, 2 Analysis: LM, 3 Link: Arcsin-sqrt, 4 Forced: One, 5 Random: None
  sPCLFile = File(sfle.d(fileDirInput,"check", "Check-Boost-LM-AS-One-None.pcl")).get_abspath()
  sTSVFile = sfle.d(fileDirTmp,sfle.rebase(sPCLFile, sPCLExtension, sTransposeExtension ))
  sRCFile = File(sfle.d(fileDirInput,os.path.join("check",sfle.rebase(sPCLFile, sPCLExtension, sMaaslinReadConfigFileExtension )))).get_abspath()
  sOutputFile = sfle.d(os.path.join(fileDirOutput.get_abspath(),"check", os.path.splitext(os.path.basename(sPCLFile))[0], sfle.rebase(sPCLFile, sPCLExtension, sMaaslinSummaryFileExtension)))
  sfle.spipe(pE, sPCLFile, c_fileProgTranspose, sTSVFile)
  Default(sfle.op( pE, sProgMaaslin, [[True, sOutputFile], [False, sTSVFile], "-i", sRCFile, "-F", "Star_Trek_Fan"]))

  # 5 Model: Boost, 2 Analysis: LM, 3 Link: Arcsin-sqrt, 4 Forced: Two, 5 Random: None
  sPCLFile = File(sfle.d(fileDirInput,"check", "Check-Boost-LM-AS-Two-None.pcl")).get_abspath()
  sTSVFile = sfle.d(fileDirTmp,sfle.rebase(sPCLFile, sPCLExtension, sTransposeExtension ))
  sRCFile = File(sfle.d(fileDirInput,os.path.join("check",sfle.rebase(sPCLFile, sPCLExtension, sMaaslinReadConfigFileExtension )))).get_abspath()
  sOutputFile = sfle.d(os.path.join(fileDirOutput.get_abspath(),"check", os.path.splitext(os.path.basename(sPCLFile))[0], sfle.rebase(sPCLFile, sPCLExtension, sMaaslinSummaryFileExtension)))
  sfle.spipe(pE, sPCLFile, c_fileProgTranspose, sTSVFile)
  Default(sfle.op( pE, sProgMaaslin, [[True, sOutputFile], [False, sTSVFile], "-i", sRCFile, "-F", "Star_Trek_Fan,Smoking"]))

  # 6 Model: Boost, 2 Analysis: Spearman, 3 Link: Arcsin-sqrt, 4 Forced: None, 5 Random: None
  sPCLFile = File(sfle.d(fileDirInput,"check", "Check-Boost-Spearman-AS-None-None.pcl")).get_abspath()
  sTSVFile = sfle.d(fileDirTmp,sfle.rebase(sPCLFile, sPCLExtension, sTransposeExtension ))
  sRCFile = File(sfle.d(fileDirInput,os.path.join("check",sfle.rebase(sPCLFile, sPCLExtension, sMaaslinReadConfigFileExtension )))).get_abspath()
  sOutputFile = sfle.d(os.path.join(fileDirOutput.get_abspath(),"check", os.path.splitext(os.path.basename(sPCLFile))[0], sfle.rebase(sPCLFile, sPCLExtension, sMaaslinSummaryFileExtension)))
  sfle.spipe(pE, sPCLFile, c_fileProgTranspose, sTSVFile)
  Default(sfle.op( pE, sProgMaaslin, [[True, sOutputFile], [False, sTSVFile], "-i", sRCFile, "-m","spearman"]))

  # 7 Model: Boost, 2 Analysis: Wilcoxon, 3 Link: Arcsin-sqrt, 4 Forced: None, 5 Random: None
  sPCLFile = File(sfle.d(fileDirInput,"check", "Check-Boost-Wilcoxon-AS-None-None.pcl")).get_abspath()
  sTSVFile = sfle.d(fileDirTmp,sfle.rebase(sPCLFile, sPCLExtension, sTransposeExtension ))
  sRCFile = File(sfle.d(fileDirInput,os.path.join("check",sfle.rebase(sPCLFile, sPCLExtension, sMaaslinReadConfigFileExtension )))).get_abspath()
  sOutputFile = sfle.d(os.path.join(fileDirOutput.get_abspath(),"check", os.path.splitext(os.path.basename(sPCLFile))[0], sfle.rebase(sPCLFile, sPCLExtension, sMaaslinSummaryFileExtension)))
  sfle.spipe(pE, sPCLFile, c_fileProgTranspose, sTSVFile)
  Default(sfle.op( pE, sProgMaaslin, [[True, sOutputFile], [False, sTSVFile], "-i", sRCFile, "-m","wilcoxon"]))

  # 8 Model: Boost, 2 Analysis: quasi, 3 Link: Arcsin-sqrt, 4 Forced: None, 5 Random: None
  sPCLFile = File(sfle.d(fileDirInput,"check", "Check-Boost-Quasi-AS-None-None.pcl")).get_abspath()
  sTSVFile = sfle.d(fileDirTmp,sfle.rebase(sPCLFile, sPCLExtension, sTransposeExtension ))
  sRCFile = File(sfle.d(fileDirInput,os.path.join("check",sfle.rebase(sPCLFile, sPCLExtension, sMaaslinReadConfigFileExtension )))).get_abspath()
  sOutputFile = sfle.d(os.path.join(fileDirOutput.get_abspath(),"check", os.path.splitext(os.path.basename(sPCLFile))[0], sfle.rebase(sPCLFile, sPCLExtension, sMaaslinSummaryFileExtension)))
  sfle.spipe(pE, sPCLFile, c_fileProgTranspose, sTSVFile)
  Default(sfle.op( pE, sProgMaaslin, [[True, sOutputFile], [False, sTSVFile], "-i", sRCFile, "-m", "quasi"]))

  # 9 Model: Boost, 2 Analysis: neg-bin, 3 Link: Arcsin-sqrt, 4 Forced: None, 5 Random: None
  sPCLFile = File(sfle.d(fileDirInput,"check", "Check-Boost-NB-AS-None-None.pcl")).get_abspath()
  sTSVFile = sfle.d(fileDirTmp,sfle.rebase(sPCLFile, sPCLExtension, sTransposeExtension ))
  sRCFile = File(sfle.d(fileDirInput,os.path.join("check",sfle.rebase(sPCLFile, sPCLExtension, sMaaslinReadConfigFileExtension )))).get_abspath()
  sOutputFile = sfle.d(os.path.join(fileDirOutput.get_abspath(),"check", os.path.splitext(os.path.basename(sPCLFile))[0], sfle.rebase(sPCLFile, sPCLExtension, sMaaslinSummaryFileExtension)))
  sfle.spipe(pE, sPCLFile, c_fileProgTranspose, sTSVFile)
  Default(sfle.op( pE, sProgMaaslin, [[True, sOutputFile], [False, sTSVFile], "-i", sRCFile, "-m", "neg_binomial"]))

  # 10 Model: Boost, 2 Analysis: LM, 3 Link: None, 4 Forced: None, 5 Random: None
  sPCLFile = File(sfle.d(fileDirInput,"check", "Check-Boost-LM-None-None-None.pcl")).get_abspath()
  sTSVFile = sfle.d(fileDirTmp,sfle.rebase(sPCLFile, sPCLExtension, sTransposeExtension ))
  sRCFile = File(sfle.d(fileDirInput,os.path.join("check",sfle.rebase(sPCLFile, sPCLExtension, sMaaslinReadConfigFileExtension )))).get_abspath()
  sOutputFile = sfle.d(os.path.join(fileDirOutput.get_abspath(),"check", os.path.splitext(os.path.basename(sPCLFile))[0], sfle.rebase(sPCLFile, sPCLExtension, sMaaslinSummaryFileExtension)))
  sfle.spipe(pE, sPCLFile, c_fileProgTranspose, sTSVFile)
  Default(sfle.op( pE, sProgMaaslin, [[True, sOutputFile], [False, sTSVFile], "-i", sRCFile, "-l", "none"]))

  # 11 Model: None, 2 Analysis: LM, 3 Link: Arcsin-sqrt, 4 Forced: None, 5 Random: None
  sPCLFile = File(sfle.d(fileDirInput,"check", "Check-None-LM-AS-None-None.pcl")).get_abspath()
  sTSVFile = sfle.d(fileDirTmp,sfle.rebase(sPCLFile, sPCLExtension, sTransposeExtension ))
  sRCFile = File(sfle.d(fileDirInput,os.path.join("check",sfle.rebase(sPCLFile, sPCLExtension, sMaaslinReadConfigFileExtension )))).get_abspath()
  sOutputFile = sfle.d(os.path.join(fileDirOutput.get_abspath(),"check", os.path.splitext(os.path.basename(sPCLFile))[0], sfle.rebase(sPCLFile, sPCLExtension, sMaaslinSummaryFileExtension)))
  sfle.spipe(pE, sPCLFile, c_fileProgTranspose, sTSVFile)
  Default(sfle.op( pE, sProgMaaslin, [[True, sOutputFile], [False, sTSVFile], "-i", sRCFile, "-s", "none"]))

  # 12 Model: Forwards, 2 Analysis: LM, 3 Link: Arcsin-sqrt, 4 Forced: None, 5 Random: None
  sPCLFile = File(sfle.d(fileDirInput,"check", "Check-Forwards-LM-AS-None-None.pcl")).get_abspath()
  sTSVFile = sfle.d(fileDirTmp,sfle.rebase(sPCLFile, sPCLExtension, sTransposeExtension ))
  sRCFile = File(sfle.d(fileDirInput,os.path.join("check",sfle.rebase(sPCLFile, sPCLExtension, sMaaslinReadConfigFileExtension )))).get_abspath()
  sOutputFile = sfle.d(os.path.join(fileDirOutput.get_abspath(),"check", os.path.splitext(os.path.basename(sPCLFile))[0], sfle.rebase(sPCLFile, sPCLExtension, sMaaslinSummaryFileExtension)))
  sfle.spipe(pE, sPCLFile, c_fileProgTranspose, sTSVFile)
  Default(sfle.op( pE, sProgMaaslin, [[True, sOutputFile], [False, sTSVFile], "-i", sRCFile, "-s", "forward"]))

  # 13 Model: Backwards, 2 Analysis: LM, 3 Link: Arcsin-sqrt, 4 Forced: None, 5 Random: None
  sPCLFile = File(sfle.d(fileDirInput,"check", "Check-Backwards-LM-AS-None-None.pcl")).get_abspath()
  sTSVFile = sfle.d(fileDirTmp,sfle.rebase(sPCLFile, sPCLExtension, sTransposeExtension ))
  sRCFile = File(sfle.d(fileDirInput,os.path.join("check",sfle.rebase(sPCLFile, sPCLExtension, sMaaslinReadConfigFileExtension )))).get_abspath()
  sOutputFile = sfle.d(os.path.join(fileDirOutput.get_abspath(),"check", os.path.splitext(os.path.basename(sPCLFile))[0], sfle.rebase(sPCLFile, sPCLExtension, sMaaslinSummaryFileExtension)))
  sfle.spipe(pE, sPCLFile, c_fileProgTranspose, sTSVFile)
  #Default(sfle.op( pE, sProgMaaslin, [[True, sOutputFile], [False, sTSVFile], "-i", sRCFile, "-s", "backward"]))