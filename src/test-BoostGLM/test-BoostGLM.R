c_strDir <- file.path(getwd( ),"..")

source(file.path(c_strDir,"lib","Constants.R"))
source(file.path(c_strDir,"lib","Utility.R"))
source(file.path(c_strDir,"lib","AnalysisModules.R"))

# General setup
covX1 = c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
covX2 = c(144.4, 245.9, 141.9, 253.3, 144.7, 244.1, 150.7, 245.2, 160.1)
covX3 = as.factor(c(1,2,3,1,2,3,1,2,3))
covX4 = as.factor(c(1,1,1,1,2,2,2,2,2))
covX5 = as.factor(c(1,2,1,2,1,2,1,2,1))
covY = c(.26,  .31,  .25,  .50,  .36,  .40,  .52,  .28,  .38)
frmeTmp = data.frame(Covariate1=covX1, Covariate2=covX2, Covariate3=covX3, Covariate4=covX4, Covariate5=covX5, adCur= covY)
iTaxon = 6

lsCov1 = list()
lsCov1$name = "Covariate1"
lsCov1$orig = "Covariate1"
lsCov1$taxon = "adCur"
lsCov1$data = covY
lsCov1$factors = "Covariate1"
lsCov1$metadata = covX1
vdCoef = c()
vdCoef["(Intercept)"]=round(0.0345077486,5)
vdCoef["Covariate1"]= round(0.0052097355,5)
vdCoef["Covariate2"]= round(0.0005806568,5)
vdCoef["Covariate32"]=round(-0.1333421874,5)
vdCoef["Covariate33"]=round(-0.1072006419,5)
vdCoef["Covariate42"]=round(0.0849198280,5)
lsCov1$value = c(Covariate1=round(0.005209736,5))
lsCov1$std = round(0.0063781728,5)
lsCov1$allCoefs = vdCoef
lsCov2 = list()
lsCov2$name = "Covariate2"
lsCov2$orig = "Covariate2"
lsCov2$taxon = "adCur"
lsCov2$data = covY
lsCov2$factors = "Covariate2"
lsCov2$metadata = covX2
lsCov2$value = c(Covariate2=round(0.0005806568,5))
lsCov2$std = round(0.0006598436,5)
lsCov2$allCoefs = vdCoef
lsCov3 = list()
lsCov3$name = "Covariate3"
lsCov3$orig = "Covariate32"
lsCov3$taxon = "adCur"
lsCov3$data = covY
lsCov3$factors = "Covariate3"
lsCov3$metadata = covX3
lsCov3$value = c(Covariate32=round(-0.1333422,5))
lsCov3$std = round(0.0895657826,5)
lsCov3$allCoefs = vdCoef
lsCov4 = list()
lsCov4$name = "Covariate3"
lsCov4$orig = "Covariate33"
lsCov4$taxon = "adCur"
lsCov4$data = covY
lsCov4$factors = "Covariate3"
lsCov4$metadata = covX3
lsCov4$value = c(Covariate33=round(-0.1072006,5))
lsCov4$std = round(0.0792209541,5)
lsCov4$allCoefs = vdCoef
lsCov5 = list()
lsCov5$name = "Covariate4"
lsCov5$orig = "Covariate42"
lsCov5$taxon = "adCur"
lsCov5$data = covY
lsCov5$factors = "Covariate4"
lsCov5$metadata = covX4
lsCov5$value = c(Covariate42=round(0.08491983,5))
lsCov5$std = round(0.0701018621,5)
lsCov5$allCoefs = vdCoef

context("Test funcClean")

context("Test funcBugHybrid")
# multiple covariates, one call lm
aiMetadata = c(1:5)
aiData = c(iTaxon)
dFreq = 0.5 / length( aiMetadata )
dSig = 0.25
dMinSamp = 0.1
adP = c()
lsSig = list()
funcReg = NA
funcAnalysis = funcLM
funcGetResult = funcGetLMResults
lsData = list(frmeData=frmeTmp, aiMetadata=aiMetadata, aiData=aiData, lsQCCounts=list())
lsData$astrMetadata = names(frmeTmp)[aiMetadata]

adPExpected = round(c(0.4738687,0.4436566,0.4665972,0.5378693,0.3124672),5)
QCExpected = list(iLms=numeric(0))
lsSigExpected = list()
lsSigExpected[[1]] = lsCov1
lsSigExpected[[2]] = lsCov2 
lsSigExpected[[3]] = lsCov3
lsSigExpected[[4]] = lsCov4
lsSigExpected[[5]] = lsCov5
expectedReturn = list(adP=adPExpected,lsSig=lsSigExpected,lsQCCounts=QCExpected)
receivedReturn = funcBugHybrid(iTaxon=iTaxon,frmeData=frmeTmp,lsData=lsData,aiMetadata=aiMetadata,dFreq=dFreq,dSig=dSig,dMinSamp=dMinSamp,adP=adP,lsSig=lsSig, strLog=NA,funcReg=funcReg,lsNonPenalizedPredictors=NULL,funcAnalysis=funcAnalysis,lsRandomCovariates=NULL,funcGetResult=funcGetResult)
receivedReturn$adP = round(receivedReturn$adP,5)

vCoefs=receivedReturn$lsSig[[1]]$allCoefs
vCoefs[1]=round(vCoefs[1],5)
vCoefs[2]=round(vCoefs[2],5)
vCoefs[3]=round(vCoefs[3],5)
vCoefs[4]=round(vCoefs[4],5)
vCoefs[5]=round(vCoefs[5],5)
vCoefs[6]=round(vCoefs[6],5)
receivedReturn$lsSig[[1]]$allCoefs=vCoefs
receivedReturn$lsSig[[2]]$allCoefs=vCoefs
receivedReturn$lsSig[[3]]$allCoefs=vCoefs
receivedReturn$lsSig[[4]]$allCoefs=vCoefs
receivedReturn$lsSig[[5]]$allCoefs=vCoefs
vValue=c()
vValue[receivedReturn$lsSig[[1]]$orig]=round(receivedReturn$lsSig[[1]]$value[[1]],5)
receivedReturn$lsSig[[1]]$value=vValue
vValue=c()
vValue[receivedReturn$lsSig[[2]]$orig]=round(receivedReturn$lsSig[[2]]$value[[1]],5)
receivedReturn$lsSig[[2]]$value=vValue
vValue=c()
vValue[receivedReturn$lsSig[[3]]$orig]=round(receivedReturn$lsSig[[3]]$value[[1]],5)
receivedReturn$lsSig[[3]]$value=vValue
vValue=c()
vValue[receivedReturn$lsSig[[4]]$orig]=round(receivedReturn$lsSig[[4]]$value[[1]],5)
receivedReturn$lsSig[[4]]$value=vValue
vValue=c()
vValue[receivedReturn$lsSig[[5]]$orig]=round(receivedReturn$lsSig[[5]]$value[[1]],5)
receivedReturn$lsSig[[5]]$value=vValue
receivedReturn$lsSig[[1]]$std=round(receivedReturn$lsSig[[1]]$std,5)
receivedReturn$lsSig[[2]]$std=round(receivedReturn$lsSig[[2]]$std,5)
receivedReturn$lsSig[[3]]$std=round(receivedReturn$lsSig[[3]]$std,5)
receivedReturn$lsSig[[4]]$std=round(receivedReturn$lsSig[[4]]$std,5)
receivedReturn$lsSig[[5]]$std=round(receivedReturn$lsSig[[5]]$std,5)
test_that("funcBugHybrid works with the lm option with multiple covariates.",{expect_equal(receivedReturn,expectedReturn)})


# single covariate, single call lm
aiMetadata = c(1)
dFreq = 0.5 / length( aiMetadata )
lsData$astrMetadata = names(frmeTmp)[aiMetadata]
adPExpected = round(c(0.1081731),5)
QCExpected = list(iLms=numeric(0))
lsSigExpected = list()
lsSigExpected[[1]] = lsCov1
lsSigExpected[[1]]$std=round(0.005278468,5)
vdCoef = c()
vdCoef["(Intercept)"]=round(-0.102410716,5)
vdCoef["Covariate1"]= round(0.009718095,5)
lsSigExpected[[1]]$allCoefs= vdCoef
lsSigExpected[[1]]$value = c(Covariate1=round(0.009718095,5))

expectedReturn = list(adP=adPExpected,lsSig=lsSigExpected,lsQCCounts=QCExpected)
receivedReturn = funcBugHybrid(iTaxon=iTaxon,frmeData=frmeTmp,lsData=lsData,aiMetadata=aiMetadata,dFreq=dFreq,dSig=dSig,dMinSamp=dMinSamp,adP=adP,lsSig=lsSig, strLog=NA,funcReg=funcReg,lsNonPenalizedPredictors=NULL,funcAnalysis=funcAnalysis,lsRandomCovariates=NULL,funcGetResult=funcGetResult)
receivedReturn$adP = round(receivedReturn$adP,5)

vCoefs=receivedReturn$lsSig[[1]]$allCoefs
vCoefs[1]=round(vCoefs[1],5)
vCoefs[2]=round(vCoefs[2],5)
receivedReturn$lsSig[[1]]$allCoefs=vCoefs
vValue=c()
vValue[receivedReturn$lsSig[[1]]$orig]=round(receivedReturn$lsSig[[1]]$value[[1]],5)
receivedReturn$lsSig[[1]]$value=vValue
receivedReturn$lsSig[[1]]$std=round(0.005278468,5)
test_that("funcBugHybrid works with the lm option with 1 covariates.",{expect_equal(receivedReturn,expectedReturn)})


# multiple covariate, single call univariate
funcReg = NA
funcAnalysis = funcDoUnivariate
funcGetResult = NA
aiMetadata = c(3,1,2)
dFreq = 0.5 / length( aiMetadata )
lsData$astrMetadata = names(frmeTmp)[aiMetadata]
adPExpected = round(c(1.0,1.0,0.09679784,0.21252205),5)
QCExpected = list(iLms=numeric(0))
lsSigExpected = list()
lsCov1 = list()
lsCov1$name = "Covariate3"
lsCov1$orig = "Covariate32"
lsCov1$taxon = "adCur"
lsCov1$data = covY
lsCov1$factors = "Covariate3"
lsCov1$metadata = frmeTmp[["Covariate3"]]
vdCoef = c(Covariate32=NA)
lsCov1$value = vdCoef
lsCov1$std = sd(frmeTmp[["Covariate3"]])
lsCov1$allCoefs = vdCoef
lsCov2 = list()
lsCov2$name = "Covariate3"
lsCov2$orig = "Covariate33"
lsCov2$taxon = "adCur"
lsCov2$data = covY
lsCov2$factors = "Covariate3"
lsCov2$metadata = frmeTmp[["Covariate3"]]
vdCoef = c(Covariate33=NA)
lsCov2$value = vdCoef
lsCov2$std = sd(frmeTmp[["Covariate3"]])
lsCov2$allCoefs = vdCoef
lsCov3 = list()
lsCov3$name = "Covariate1"
lsCov3$orig = "Covariate1"
lsCov3$taxon = "adCur"
lsCov3$data = covY
lsCov3$factors = "Covariate1"
lsCov3$metadata = frmeTmp[["Covariate1"]]
vdCoef = c(Covariate1=0.6)
lsCov3$value = vdCoef
lsCov3$std = sd(frmeTmp[["Covariate1"]])
lsCov3$allCoefs = vdCoef
lsCov4 = list()
lsCov4$name = "Covariate2"
lsCov4$orig = "Covariate2"
lsCov4$taxon = "adCur"
lsCov4$data = covY
lsCov4$factors = "Covariate2"
lsCov4$metadata = frmeTmp[["Covariate2"]]
vdCoef = c(Covariate2=0.46666667)
lsCov4$value = vdCoef
lsCov4$std = sd(frmeTmp[["Covariate2"]])
lsCov4$allCoefs = vdCoef

lsSigExpected = list()
lsSigExpected[[1]] = lsCov1
lsSigExpected[[2]] = lsCov2
lsSigExpected[[3]] = lsCov3
lsSigExpected[[4]] = lsCov4

expectedReturn = list(adP=adPExpected,lsSig=lsSigExpected,lsQCCounts=QCExpected)
receivedReturn = funcBugHybrid(iTaxon=iTaxon,frmeData=frmeTmp,lsData=lsData,aiMetadata=aiMetadata,dFreq=dFreq,dSig=dSig,dMinSamp=dMinSamp,adP=adP,lsSig=lsSig, strLog=NA,funcReg=funcReg,lsNonPenalizedPredictors=NULL,funcAnalysis=funcAnalysis,lsRandomCovariates=NULL,funcGetResult=funcGetResult)
receivedReturn$adP = round(receivedReturn$adP,5)
test_that("funcBugHybrid works with the univariate option with 3 covariates.",{expect_equal(receivedReturn,expectedReturn)})


# single covariate, single call univariate
funcReg = NA
funcAnalysis = funcDoUnivariate
funcGetResult = NA
aiMetadata = c(1)
dFreq = 0.5 / length( aiMetadata )
lsData$astrMetadata = names(frmeTmp)[aiMetadata]
adPExpected = round(c(0.09679784),5)
QCExpected = list(iLms=numeric(0))
lsSigExpected = list()
lsSigExpected[[1]] = lsCov3

expectedReturn = list(adP=adPExpected,lsSig=lsSigExpected,lsQCCounts=QCExpected)
receivedReturn = funcBugHybrid(iTaxon=iTaxon,frmeData=frmeTmp,lsData=lsData,aiMetadata=aiMetadata,dFreq=dFreq,dSig=dSig,dMinSamp=dMinSamp,adP=adP,lsSig=lsSig, strLog=NA,funcReg=funcReg,lsNonPenalizedPredictors=NULL,funcAnalysis=funcAnalysis,lsRandomCovariates=NULL,funcGetResult=funcGetResult)
receivedReturn$adP = round(receivedReturn$adP,5)
test_that("funcBugHybrid works with the univariate option with 1 covariates.",{expect_equal(receivedReturn,expectedReturn)})


context("Test funcBugs")
#One LM run
frmeData=frmeTmp
aiMetadata=c(1)
aiData=c(iTaxon)
strData=NA
dFreq= 0.5 / length( aiMetadata )
dSig=0.25
dMinSamp=0.1
strDirOut=NA
funcReg=NA
lsNonPenalizedPredictors=NULL
lsRandomCovariates=NULL
funcAnalysis=funcLM
funcGetResults=funcGetLMResults
fDoRPlot=FALSE
lsData = list(frmeData=frmeData, aiMetadata=aiMetadata, aiData=aiData, lsQCCounts=list())
lsData$astrMetadata = names(frmeTmp)[aiMetadata]
QCExpected = list(iLms=numeric(0))

expectedReturn = list(aiReturnBugs=aiData,lsQCCounts=QCExpected)
receivedReturn = funcBugs(frmeData=frmeData, lsData=lsData, aiMetadata=aiMetadata, aiData=aiData, strData=strData, dFreq=dFreq, dSig=dSig, dMinSamp=dMinSamp,strDirOut=strDirOut, funcReg=funcReg,lsNonPenalizedPredictors=lsNonPenalizedPredictors,funcAnalysis=funcAnalysis,lsRandomCovariates=lsRandomCovariates,funcGetResults=funcGetResults,fDoRPlot=fDoRPlot)

test_that("funcBugs works with the lm option with 1 covariate.",{expect_equal(receivedReturn,expectedReturn)})

#multiple LM run
frmeData=frmeTmp
aiMetadata=c(1:5)
aiData=c(iTaxon)
strData=NA
dFreq= 0.5 / length( aiMetadata )
dSig=0.25
dMinSamp=0.1
strDirOut=NA
funcReg=NA
lsNonPenalizedPredictors=NULL
lsRandomCovariates=NULL
funcAnalysis=funcLM
funcGetResults=funcGetLMResults
fDoRPlot=FALSE
lsData = list(frmeData=frmeData, aiMetadata=aiMetadata, aiData=aiData, lsQCCounts=list())
lsData$astrMetadata = names(frmeTmp)[aiMetadata]
QCExpected = list(iLms=numeric(0))

expectedReturn = list(aiReturnBugs=aiData,lsQCCounts=QCExpected)
receivedReturn = funcBugs(frmeData=frmeData, lsData=lsData, aiMetadata=aiMetadata, aiData=aiData, strData=strData, dFreq=dFreq, dSig=dSig, dMinSamp=dMinSamp,strDirOut=strDirOut, funcReg=funcReg,lsNonPenalizedPredictors=lsNonPenalizedPredictors,funcAnalysis=funcAnalysis,lsRandomCovariates=lsRandomCovariates,funcGetResults=funcGetResults,fDoRPlot=fDoRPlot)

print("START START")
print(expectedReturn)
print("RECEIVED")
print(receivedReturn)
print("STOP STOP")

test_that("funcBugs works with the lm option with multiple covariates.",{expect_equal(receivedReturn,expectedReturn)})