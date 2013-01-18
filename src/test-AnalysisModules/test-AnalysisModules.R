c_strDir <- file.path(getwd( ),"..")

source(file.path(c_strDir,"lib","Constants.R"))
source(file.path(c_strDir,"lib","Utility.R"))

#Test Utilities
context("Test funcGetLMResults")
context("Test funcGetStepPredictors")

context("Test funcMakeContrasts")
covX1 = c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
covX2 = c(144.4, 245.9, 141.9, 253.3, 144.7, 244.1, 150.7, 245.2, 160.1)
covX3 = as.factor(c(1,2,3,1,2,3,1,2,3))
covX4 = as.factor(c(1,1,1,1,2,2,2,2,2))
covX5 = as.factor(c(1,2,1,2,1,2,1,2,1))
covY = c(.26,  .31,  .25,  .50,  .36,  .40,  .52,  .28,  .38)
frmeTmp = data.frame(Covariate1=covX1, Covariate2=covX2, Covariate3=covX3, Covariate4=covX4, Covariate5=covX5, adCur= covY)
iTaxon = 6
#Add in updating QC errors
#Add in random covariates
strFormula = "adCur ~ Covariate1"
strRandomFormula = NULL
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate1"
lsSig[[1]]$orig = "Covariate1"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = covY
lsSig[[1]]$factors = "Covariate1"
lsSig[[1]]$metadata = covX1
vdCoef = c(Covariate1=0.6)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(covX1)
lsSig[[1]]$allCoefs = vdCoef
ret1 = funcMakeContrasts(strFormula=strFormula, strRandomFormula=strRandomFormula, frmeTmp=frmeTmp, iTaxon=iTaxon,
    functionContrast=function(x,adCur,dfData)
    {
      retList = list()
      ret = cor.test(as.formula(paste("~",x,"+ adCur")), data=dfData, method="spearman", na.action=c_strNA_Action)
      #Returning rho for the coef in a named vector
      vdCoef = c()
      vdCoef[[x]]=ret$estimate
      retList[[1]]=list(p.value=ret$p.value,SD=sd(dfData[[x]]),name=x,coef=vdCoef)
      return(retList)
    }, lsQCCounts=list())
ret1$adP = round(ret1$adP,5)
test_that("1. Test that the funcMakeContrasts works on a continuous variable.",{
  expect_equal(ret1,list(adP=round(c(0.09679784),5),lsSig=lsSig,lsQCCounts=list()))})

strFormula = "adCur ~ Covariate1 + Covariate2"
strRandomFormula = NULL
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate1"
lsSig[[1]]$orig = "Covariate1"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = covY
lsSig[[1]]$factors = "Covariate1"
lsSig[[1]]$metadata = covX1
vdCoef = c(Covariate1=0.6)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(covX1)
lsSig[[1]]$allCoefs = vdCoef
lsSig[[2]] = list()
lsSig[[2]]$name = "Covariate2"
lsSig[[2]]$orig = "Covariate2"
lsSig[[2]]$taxon = "adCur"
lsSig[[2]]$data = covY
lsSig[[2]]$factors = "Covariate2"
lsSig[[2]]$metadata = covX2
vdCoef = c(Covariate2=0.46666667)
lsSig[[2]]$value = vdCoef
lsSig[[2]]$std = sd(covX2)
lsSig[[2]]$allCoefs = vdCoef
ret1 = funcMakeContrasts(strFormula=strFormula, strRandomFormula=strRandomFormula, frmeTmp=frmeTmp, iTaxon=iTaxon,
    functionContrast=function(x,adCur,dfData)
    {
      retList = list()
      ret = cor.test(as.formula(paste("~",x,"+ adCur")), data=dfData, method="spearman", na.action=c_strNA_Action)
      #Returning rho for the coef in a named vector
      vdCoef = c()
      vdCoef[[x]]=ret$estimate
      retList[[1]]=list(p.value=ret$p.value,SD=sd(dfData[[x]]),name=x,coef=vdCoef)
      return(retList)
    }, lsQCCounts=list())
ret1$adP = round(ret1$adP,5)
test_that("Test that the funcMakeContrasts works on 2 continuous variables.",{
  expect_equal(ret1,list(adP=round(c(0.09679784,0.21252205),5),lsSig=lsSig,lsQCCounts=list()))})

strFormula = "adCur ~ Covariate4"
strRandomFormula = NULL
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate4"
lsSig[[1]]$orig = "Covariate42"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = covY
lsSig[[1]]$factors = "Covariate4"
lsSig[[1]]$metadata = covX4 #update
vdCoef = c(Covariate42=NA)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(covX4) #update
lsSig[[1]]$allCoefs = vdCoef
# Get return
rets = funcMakeContrasts(strFormula=strFormula, strRandomFormula=strRandomFormula, frmeTmp=frmeTmp, iTaxon=iTaxon,
    functionContrast=function(x,adCur,dfData)
    {
      retList = list()
      lmodKW = kruskal(adCur,dfData[[x]],group=FALSE,p.adj="holm")
      asLevels = levels(dfData[[x]])
      # The names of the generated comparisons, sometimes the control is first sometimes it is not so
      # We will just check which is in the names and use that
      asComparisons = row.names(lmodKW$comparisons)
      #Get the comparison with the control
      for(sLevel in asLevels[2:length(asLevels)])
      {
        sComparison = intersect(c(paste(asLevels[1],sLevel,sep=" - "),paste(sLevel,asLevels[1],sep=" - ")),asComparisons)
        #Returning NA for the coef in a named vector
        vdCoef = c()
        vdCoef[[paste(x,sLevel,sep="")]]=NA
        retList[[length(retList)+1]]=list(p.value=lmodKW$comparisons[sComparison,"p.value"],SD=NA,name=paste(x,sLevel,sep=""),coef=vdCoef)
      }
      return(retList)
    }, lsQCCounts=list())
rets$adP=round(rets$adP,digits=5)
test_that("Test that the funcMakeContrasts works on 1 factor covariate with 2 levels.",{
  expect_equal(rets,list(adP=round(c(0.24434),5),lsSig=lsSig,lsQCCounts=list()))})

strFormula = "adCur ~ Covariate3"
strRandomFormula = NULL
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate3"
lsSig[[1]]$orig = "Covariate32"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = covY
lsSig[[1]]$factors = "Covariate3"
lsSig[[1]]$metadata = covX3 #update
vdCoef = c(Covariate32=NA)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(covX3) #update
lsSig[[1]]$allCoefs = vdCoef
lsSig[[2]] = list()
lsSig[[2]]$name = "Covariate3"
lsSig[[2]]$orig = "Covariate33"
lsSig[[2]]$taxon = "adCur"
lsSig[[2]]$data = covY
lsSig[[2]]$factors = "Covariate3"
lsSig[[2]]$metadata = covX3 #update
vdCoef = c(Covariate33=NA)
lsSig[[2]]$value = vdCoef
lsSig[[2]]$std = sd(covX3) #update
lsSig[[2]]$allCoefs = vdCoef
ret1 = funcMakeContrasts(strFormula=strFormula, strRandomFormula=strRandomFormula, frmeTmp=frmeTmp, iTaxon=iTaxon,
    functionContrast=function(x,adCur,dfData)
    {
      retList = list()
      lmodKW = kruskal(adCur,dfData[[x]],group=FALSE,p.adj="holm")
      asLevels = levels(dfData[[x]])
      # The names of the generated comparisons, sometimes the control is first sometimes it is not so
      # We will just check which is in the names and use that
      asComparisons = row.names(lmodKW$comparisons)
      #Get the comparison with the control
      for(sLevel in asLevels[2:length(asLevels)])
      {
        sComparison = intersect(c(paste(asLevels[1],sLevel,sep=" - "),paste(sLevel,asLevels[1],sep=" - ")),asComparisons)
        #Returning NA for the coef in a named vector
        vdCoef = c()
        vdCoef[[paste(x,sLevel,sep="")]]=NA
        retList[[length(retList)+1]]=list(p.value=lmodKW$comparisons[sComparison,"p.value"],SD=NA,name=paste(x,sLevel,sep=""),coef=vdCoef)
      }
      return(retList)
    }, lsQCCounts=list())
ret1$adP = round(ret1$adP,5)
test_that("Test that the funcMakeContrasts works on 1 factor covariate with 3 levels.",{
  expect_equal(ret1,list(adP=c(1.0,1.0),lsSig=lsSig,lsQCCounts=list()))})

strFormula = "adCur ~ Covariate4 + Covariate5"
strRandomFormula = NULL
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate4"
lsSig[[1]]$orig = "Covariate42"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = covY
lsSig[[1]]$factors = "Covariate4"
lsSig[[1]]$metadata = covX4 #update
vdCoef = c(Covariate42=NA)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(covX4) #update
lsSig[[1]]$allCoefs = vdCoef
lsSig[[2]] = list()
lsSig[[2]]$name = "Covariate5"
lsSig[[2]]$orig = "Covariate52"
lsSig[[2]]$taxon = "adCur"
lsSig[[2]]$data = covY
lsSig[[2]]$factors = "Covariate5"
lsSig[[2]]$metadata = covX5 #update
vdCoef = c(Covariate52=NA)
lsSig[[2]]$value = vdCoef
lsSig[[2]]$std = sd(covX5) #update
lsSig[[2]]$allCoefs = vdCoef
ret1 = funcMakeContrasts(strFormula=strFormula, strRandomFormula=strRandomFormula, frmeTmp=frmeTmp, iTaxon=iTaxon,
    functionContrast=function(x,adCur,dfData)
    {
      retList = list()
      lmodKW = kruskal(adCur,dfData[[x]],group=FALSE,p.adj="holm")
      asLevels = levels(dfData[[x]])
      # The names of the generated comparisons, sometimes the control is first sometimes it is not so
      # We will just check which is in the names and use that
      asComparisons = row.names(lmodKW$comparisons)
      #Get the comparison with the control
      for(sLevel in asLevels[2:length(asLevels)])
      {
        sComparison = intersect(c(paste(asLevels[1],sLevel,sep=" - "),paste(sLevel,asLevels[1],sep=" - ")),asComparisons)
        #Returning NA for the coef in a named vector
        vdCoef = c()
        vdCoef[[paste(x,sLevel,sep="")]]=NA
        retList[[length(retList)+1]]=list(p.value=lmodKW$comparisons[sComparison,"p.value"],SD=NA,name=paste(x,sLevel,sep=""),coef=vdCoef)
      }
      return(retList)
    }, lsQCCounts=list())
ret1$adP = round(ret1$adP,5)
test_that("1. Test that the funcMakeContrasts works on 2 factor covariate with 2 levels.",{
  expect_equal(ret1,list(adP=round(c(0.24434,0.655852),5),lsSig=lsSig,lsQCCounts=list()))})


#Test Model selection


context("Test funcBoostModel")
context("Test funcForwardModel")
context("Test funcBackwardsModel")


#Test Univariates
context("Test funcSpearman")
strFormula = "adCur ~ Covariate1"
adCur = c(.26,  .31,  .25,  .50,  .36,  .40,  .52,  .28,  .38)
x = c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
frmeTmp = data.frame(Covariate1=x, adCur=adCur)
iTaxon = 2
lsQCCounts = list()
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate1"
lsSig[[1]]$orig = "Covariate1"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = adCur
lsSig[[1]]$factors = "Covariate1"
lsSig[[1]]$metadata = x
vdCoef = c(Covariate1=0.6)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(x)
lsSig[[1]]$allCoefs = vdCoef
ret1 = funcSpearman(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=NULL)
ret1$adP = round(ret1$adP,5)
test_that("Test that the spearman test has the correct results for 1 covariate.",{
  expect_equal(ret1,list(adP=round(c(0.09679784),5),lsSig=lsSig,lsQCCounts=list()))
})

strFormula = "adCur ~ Covariate1 + Covariate2"
frmeTmp = data.frame(Covariate1=x, Covariate2=x, adCur=adCur)
iTaxon = 3
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate1"
lsSig[[1]]$orig = "Covariate1"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = adCur
lsSig[[1]]$factors = "Covariate1"
lsSig[[1]]$metadata = x
vdCoef = c(Covariate1=0.6)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(x)
lsSig[[1]]$allCoefs = vdCoef
lsSig[[2]] = list()
lsSig[[2]]$name = "Covariate2"
lsSig[[2]]$orig = "Covariate2"
lsSig[[2]]$taxon = "adCur"
lsSig[[2]]$data = adCur 
lsSig[[2]]$factors = "Covariate2"
lsSig[[2]]$metadata = x
vdCoef = c(Covariate2=0.6)
lsSig[[2]]$value = vdCoef
lsSig[[2]]$std = sd(x)
lsSig[[2]]$allCoefs = vdCoef
lsQCCounts = list()
ret1 = funcSpearman(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=NULL)
ret1$adP = round(ret1$adP,5)
test_that("Test that the spearman test has the correct results for 2 covariates.",{
  expect_equal(ret1,list(adP=round(c(0.09679784,0.09679784),5),lsSig=lsSig,lsQCCounts=list()))
})


context("Test funcWilcoxon")
strFormula = "adCur ~ Covariate1"
x = c(TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE)
frmeTmp = data.frame(Covariate1=x, adCur=adCur)
iTaxon = 2
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate1"
lsSig[[1]]$orig = "Covariate1"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = adCur
lsSig[[1]]$factors = "Covariate1"
lsSig[[1]]$metadata = x
vdCoef = c(Covariate1=13)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(x)
lsSig[[1]]$allCoefs = vdCoef
lsQCCounts = list()
ret1 = funcWilcoxon(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=NULL)
print("ret1")
print(ret1)
print("list(adP=round(c(0.55555556),5),lsSig=lsSig,lsQCCounts=list())")
print(list(adP=round(c(0.55555556),5),lsSig=lsSig,lsQCCounts=list()))
ret1$adP = round(ret1$adP,5)
test_that("Test that the wilcoxon test has the correct results for 1 covariate.",{
  expect_equal(ret1,list(adP=round(c(0.55555556),5),lsSig=lsSig,lsQCCounts=list()))
})


context("Test funcKruskalWallis")
strFormula = "adCur ~ Covariate1"
x = as.factor(c("one","two","three","one","one","three","two","three","two"))
frmeTmp = data.frame(Covariate1=x, adCur=adCur)
iTaxon = 2
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate1"
lsSig[[1]]$orig = "Covariate1three"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = adCur
lsSig[[1]]$factors = "Covariate1"
lsSig[[1]]$metadata = x
vdCoef = c(Covariate1three=NA)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(x)
lsSig[[1]]$allCoefs = vdCoef
lsSig[[2]] = list()
lsSig[[2]]$name = "Covariate1"
lsSig[[2]]$orig = "Covariate1two"
lsSig[[2]]$taxon = "adCur"
lsSig[[2]]$data = adCur
lsSig[[2]]$factors = "Covariate1"
lsSig[[2]]$metadata = x
vdCoef = c(Covariate1two=NA)
lsSig[[2]]$value = vdCoef
lsSig[[2]]$std = sd(x)
lsSig[[2]]$allCoefs = vdCoef
lsQCCounts = list()
ret1 = funcKruskalWallis(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=NULL)
ret1$adP = round(ret1$adP,5)
test_that("Test that the Kruskal Wallis (Nonparameteric anova) has the correct results for 1 covariate.",{
  expect_equal(ret1,list(adP=c(1.0,1.0),lsSig=lsSig,lsQCCounts=list()))
})


context("test funcDoUnivariate")
covX1 = c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
covX2 = c(144.4, 245.9, 141.9, 253.3, 144.7, 244.1, 150.7, 245.2, 160.1)
covX3 = as.factor(c(1,2,3,1,2,3,1,2,3))
covX4 = as.factor(c(1,1,1,1,2,2,2,2,2))
covX5 = as.factor(c(1,2,1,2,1,2,1,2,1))
covX6 = as.factor(c("one","two","three","one","one","three","two","three","two"))
covY = c(.26,  .31,  .25,  .50,  .36,  .40,  .52,  .28,  .38)
frmeTmp = data.frame(Covariate1=covX1, Covariate2=covX2, Covariate3=covX3, Covariate4=covX4, Covariate5=covX5, Covariate6=covX6, adCur= covY)
iTaxon = 7
# 1 cont answer
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate1"
lsSig[[1]]$orig = "Covariate1"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = adCur
lsSig[[1]]$factors = "Covariate1"
lsSig[[1]]$metadata = frmeTmp[["Covariate1"]]
vdCoef = c(Covariate1=0.6)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(frmeTmp[["Covariate1"]])
lsSig[[1]]$allCoefs = vdCoef
ret1 = funcDoUnivariate(strFormula="adCur ~ Covariate1",frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula=NULL)
ret2 = funcDoUnivariate(strFormula=NULL,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula="adCur ~ 1|Covariate1")
ret1$adP = round(ret1$adP,5)
ret2$adP = round(ret2$adP,5)
test_that("2. Test that the funcMakeContrasts works on a continuous variable.",{
  expect_equal(ret1,list(adP=round(c(0.09679784),5),lsSig=lsSig,lsQCCounts=list()))
  expect_equal(ret2,list(adP=round(c(0.09679784),5),lsSig=lsSig,lsQCCounts=list()))
})
lsSig[[2]] = list()
lsSig[[2]]$name = "Covariate2"
lsSig[[2]]$orig = "Covariate2"
lsSig[[2]]$taxon = "adCur"
lsSig[[2]]$data = adCur
lsSig[[2]]$factors = "Covariate2"
lsSig[[2]]$metadata = frmeTmp[["Covariate2"]]
vdCoef = c(Covariate2=0.46666667)
lsSig[[2]]$value = vdCoef
lsSig[[2]]$std = sd(frmeTmp[["Covariate2"]])
lsSig[[2]]$allCoefs = vdCoef
ret1 = funcDoUnivariate(strFormula="adCur ~ Covariate1 + Covariate2",frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula=NULL)
ret2 = funcDoUnivariate(strFormula=NULL,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula="adCur ~ 1|Covariate1 + 1|Covariate2")
ret1$adP = round(ret1$adP,5)
ret2$adP = round(ret2$adP,5)
test_that("Test that the funcMakeContrasts works on 2 continuous variables.",{
  expect_equal(ret1,list(adP=round(c(0.09679784,0.21252205),5),lsSig=lsSig,lsQCCounts=list()))
  expect_equal(ret2,list(adP=round(c(0.09679784,0.21252205),5),lsSig=lsSig,lsQCCounts=list()))
})
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate4"
lsSig[[1]]$orig = "Covariate4"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = adCur
lsSig[[1]]$factors = "Covariate4"
lsSig[[1]]$metadata = frmeTmp[["Covariate4"]]
vdCoef = c(Covariate4=5)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(frmeTmp[["Covariate4"]])
lsSig[[1]]$allCoefs = vdCoef
ret1 = funcDoUnivariate(strFormula="adCur ~ Covariate4",frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula=NULL)
ret2 = funcDoUnivariate(strFormula=NULL,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula="adCur ~ 1|Covariate4")
ret1$adP = round(ret1$adP,5)
ret2$adP = round(ret2$adP,5)
test_that("Test that the funcMakeContrasts works on 1 factor covariate with 2 levels.",{
  expect_equal(ret1,list(adP=round(c(0.2857143),5),lsSig=lsSig,lsQCCounts=list()))
  expect_equal(ret2,list(adP=round(c(0.2857143),5),lsSig=lsSig,lsQCCounts=list()))
})
lsSig[[2]] = list()
lsSig[[2]]$name = "Covariate5"
lsSig[[2]]$orig = "Covariate5"
lsSig[[2]]$taxon = "adCur"
lsSig[[2]]$data = adCur
lsSig[[2]]$factors = "Covariate5"
lsSig[[2]]$metadata = frmeTmp[["Covariate5"]]
vdCoef = c(Covariate5=8)
lsSig[[2]]$value = vdCoef
lsSig[[2]]$std = sd(frmeTmp[["Covariate5"]])
lsSig[[2]]$allCoefs = vdCoef
ret1 = funcDoUnivariate(strFormula="adCur ~ Covariate4 + Covariate5",frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula=NULL)
ret2 = funcDoUnivariate(strFormula=NULL,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula="adCur ~ 1|Covariate4 + 1|Covariate5")
ret3 = funcDoUnivariate(strFormula="adCur ~ Covariate4",frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula="adCur ~ 1|Covariate5")
ret1$adP = round(ret1$adP,5)
ret2$adP = round(ret2$adP,5)
ret3$adP = round(ret3$adP,5)
test_that("2. Test that the funcMakeContrasts works on 2 factor covariate with 2 levels.",{
  expect_equal(ret1,list(adP=round(c(0.2857143,0.73016),5),lsSig=lsSig,lsQCCounts=list()))
  expect_equal(ret2,list(adP=round(c(0.2857143,0.73016),5),lsSig=lsSig,lsQCCounts=list()))
  expect_equal(ret3,list(adP=round(c(0.2857143,0.73016),5),lsSig=lsSig,lsQCCounts=list()))
})
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate4"
lsSig[[1]]$orig = "Covariate4"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = adCur
lsSig[[1]]$factors = "Covariate4"
lsSig[[1]]$metadata = frmeTmp[["Covariate4"]]
vdCoef = c(Covariate4=5)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(frmeTmp[["Covariate4"]])
lsSig[[1]]$allCoefs = vdCoef
lsSig[[2]] = list()
lsSig[[2]]$name = "Covariate1"
lsSig[[2]]$orig = "Covariate1"
lsSig[[2]]$taxon = "adCur"
lsSig[[2]]$data = adCur
lsSig[[2]]$factors = "Covariate1"
lsSig[[2]]$metadata = frmeTmp[["Covariate1"]]
vdCoef = c(Covariate1=0.6)
lsSig[[2]]$value = vdCoef
lsSig[[2]]$std = sd(frmeTmp[["Covariate1"]])
lsSig[[2]]$allCoefs = vdCoef
ret1 = funcDoUnivariate(strFormula="adCur ~ Covariate4 + Covariate1",frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula=NULL)
ret2 = funcDoUnivariate(strFormula=NULL,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula="adCur ~ 1|Covariate4 + 1|Covariate1")
ret3 = funcDoUnivariate(strFormula="adCur ~ Covariate4",frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula="adCur ~ 1|Covariate1")
ret1$adP = round(ret1$adP,5)
ret2$adP = round(ret2$adP,5)
ret3$adP = round(ret3$adP,5)
test_that("Test that the funcMakeContrasts works on 1 factor covariate with 2 levels and a continuous variable.",{
  expect_equal(ret1,list(adP=round(c(0.2857143,0.09679784),5),lsSig=lsSig,lsQCCounts=list()))
  expect_equal(ret2,list(adP=round(c(0.2857143,0.09679784),5),lsSig=lsSig,lsQCCounts=list()))
  expect_equal(ret3,list(adP=round(c(0.2857143,0.09679784),5),lsSig=lsSig,lsQCCounts=list()))
})
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate3"
lsSig[[1]]$orig = "Covariate32"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = adCur
lsSig[[1]]$factors = "Covariate3"
lsSig[[1]]$metadata = frmeTmp[["Covariate3"]]
vdCoef = c(Covariate32=NA)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(frmeTmp[["Covariate3"]])
lsSig[[1]]$allCoefs = vdCoef
lsSig[[2]] = list()
lsSig[[2]]$name = "Covariate3"
lsSig[[2]]$orig = "Covariate33"
lsSig[[2]]$taxon = "adCur"
lsSig[[2]]$data = adCur
lsSig[[2]]$factors = "Covariate3"
lsSig[[2]]$metadata = frmeTmp[["Covariate3"]]
vdCoef = c(Covariate33=NA)
lsSig[[2]]$value = vdCoef
lsSig[[2]]$std = sd(frmeTmp[["Covariate3"]])
lsSig[[2]]$allCoefs = vdCoef
lsSig[[3]] = list()
lsSig[[3]]$name = "Covariate1"
lsSig[[3]]$orig = "Covariate1"
lsSig[[3]]$taxon = "adCur"
lsSig[[3]]$data = adCur
lsSig[[3]]$factors = "Covariate1"
lsSig[[3]]$metadata = frmeTmp[["Covariate1"]]
vdCoef = c(Covariate1=0.6)
lsSig[[3]]$value = vdCoef
lsSig[[3]]$std = sd(frmeTmp[["Covariate1"]])
lsSig[[3]]$allCoefs = vdCoef
lsSig[[4]] = list()
lsSig[[4]]$name = "Covariate2"
lsSig[[4]]$orig = "Covariate2"
lsSig[[4]]$taxon = "adCur"
lsSig[[4]]$data = adCur
lsSig[[4]]$factors = "Covariate2"
lsSig[[4]]$metadata = frmeTmp[["Covariate2"]]
vdCoef = c(Covariate2=0.46666667)
lsSig[[4]]$value = vdCoef
lsSig[[4]]$std = sd(frmeTmp[["Covariate2"]])
lsSig[[4]]$allCoefs = vdCoef
ret1 = funcDoUnivariate(strFormula="adCur ~ Covariate3 + Covariate1 + Covariate2",frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula=NULL)
ret2 = funcDoUnivariate(strFormula=NULL,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula="adCur ~ 1|Covariate3 + 1|Covariate1 + 1|Covariate2")
ret3 = funcDoUnivariate(strFormula="adCur ~ Covariate3 + Covariate1",frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula="adCur ~ 1|Covariate2")
ret1$adP = round(ret1$adP,5)
ret2$adP = round(ret2$adP,5)
ret3$adP = round(ret3$adP,5)
test_that("Test that the funcMakeContrasts works on 1 factor covariate with 3 levels and 2 continuous variables.",{
  expect_equal(ret1,list(adP=round(c(1.0,1.0,0.09679784,0.21252205),5),lsSig=lsSig,lsQCCounts=list()))
  expect_equal(ret2,list(adP=round(c(1.0,1.0,0.09679784,0.21252205),5),lsSig=lsSig,lsQCCounts=list()))
  expect_equal(ret3,list(adP=round(c(1.0,1.0,0.09679784,0.21252205),5),lsSig=lsSig,lsQCCounts=list()))
})
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate4"
lsSig[[1]]$orig = "Covariate4"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = adCur
lsSig[[1]]$factors = "Covariate4"
lsSig[[1]]$metadata = frmeTmp[["Covariate4"]]
vdCoef = c(Covariate4=5)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(frmeTmp[["Covariate4"]])
lsSig[[1]]$allCoefs = vdCoef
lsSig[[2]] = list()
lsSig[[2]]$name = "Covariate2"
lsSig[[2]]$orig = "Covariate2"
lsSig[[2]]$taxon = "adCur"
lsSig[[2]]$data = adCur
lsSig[[2]]$factors = "Covariate2"
lsSig[[2]]$metadata = frmeTmp[["Covariate2"]]
vdCoef = c(Covariate2=0.46666667)
lsSig[[2]]$value = vdCoef
lsSig[[2]]$std = sd(frmeTmp[["Covariate2"]])
lsSig[[2]]$allCoefs = vdCoef
ret1 = funcDoUnivariate(strFormula="adCur ~ Covariate4 + Covariate2",frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula=NULL)
ret2 = funcDoUnivariate(strFormula=NULL,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula="adCur ~ 1|Covariate4 + 1|Covariate2")
ret3 = funcDoUnivariate(strFormula= "adCur ~ Covariate4",frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula="adCur ~ 1|Covariate2")
ret1$adP = round(ret1$adP,5)
ret2$adP = round(ret2$adP,5)
ret3$adP = round(ret3$adP,5)
test_that("3. Test that the funcMakeContrasts works on 2 factor covariate with 2 levels and a continuous variable.",{
  expect_equal(ret1,list(adP=round(c(0.2857143,0.21252205),5),lsSig=lsSig,lsQCCounts=list()))
  expect_equal(ret2,list(adP=round(c(0.2857143,0.21252205),5),lsSig=lsSig,lsQCCounts=list()))
  expect_equal(ret3,list(adP=round(c(0.2857143,0.21252205),5),lsSig=lsSig,lsQCCounts=list()))
})
lsSig = list()
lsSig[[1]] = list()
lsSig[[1]]$name = "Covariate4"
lsSig[[1]]$orig = "Covariate4"
lsSig[[1]]$taxon = "adCur"
lsSig[[1]]$data = adCur
lsSig[[1]]$factors = "Covariate4"
lsSig[[1]]$metadata = frmeTmp[["Covariate4"]]
vdCoef = c(Covariate4=5)
lsSig[[1]]$value = vdCoef
lsSig[[1]]$std = sd(frmeTmp[["Covariate4"]])
lsSig[[1]]$allCoefs = vdCoef
lsSig[[2]] = list()
lsSig[[2]]$name = "Covariate3"
lsSig[[2]]$orig = "Covariate32"
lsSig[[2]]$taxon = "adCur"
lsSig[[2]]$data = adCur
lsSig[[2]]$factors = "Covariate3"
lsSig[[2]]$metadata = frmeTmp[["Covariate3"]]
vdCoef = c(Covariate32=NA)
lsSig[[2]]$value = vdCoef
lsSig[[2]]$std = sd(frmeTmp[["Covariate3"]])
lsSig[[2]]$allCoefs = vdCoef
lsSig[[3]] = list()
lsSig[[3]]$name = "Covariate3"
lsSig[[3]]$orig = "Covariate33"
lsSig[[3]]$taxon = "adCur"
lsSig[[3]]$data = adCur
lsSig[[3]]$factors = "Covariate3"
lsSig[[3]]$metadata = frmeTmp[["Covariate3"]]
vdCoef = c(Covariate33=NA)
lsSig[[3]]$value = vdCoef
lsSig[[3]]$std = sd(frmeTmp[["Covariate3"]])
lsSig[[3]]$allCoefs = vdCoef
lsSig[[4]] = list()
lsSig[[4]]$name = "Covariate2"
lsSig[[4]]$orig = "Covariate2"
lsSig[[4]]$taxon = "adCur"
lsSig[[4]]$data = adCur
lsSig[[4]]$factors = "Covariate2"
lsSig[[4]]$metadata = frmeTmp[["Covariate2"]]
vdCoef = c(Covariate2=0.46666667)
lsSig[[4]]$value = vdCoef
lsSig[[4]]$std = sd(frmeTmp[["Covariate2"]])
lsSig[[4]]$allCoefs = vdCoef
ret1 = funcDoUnivariate(strFormula="adCur ~ Covariate4 + Covariate3 + Covariate2",frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula=NULL)
ret2 = funcDoUnivariate(strFormula=NULL,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=list(),strRandomFormula="adCur ~ 1|Covariate4 +1|Covariate3 + 1|Covariate2") 
ret1$adP = round(ret1$adP,5)
ret2$adP = round(ret2$adP,5)
test_that("Test that the funcMakeContrasts works on 1 factor covariate with 2 levels , 1 factor with 3 levels, and a continuous variable.",{
  expect_equal(ret1,list(adP=round(c(0.2857143,1.0,1.0,0.21252205),5),lsSig=lsSig,lsQCCounts=list()))
  expect_equal(ret2,list(adP=round(c(0.2857143,1.0,1.0,0.21252205),5),lsSig=lsSig,lsQCCounts=list()))
})

#Test multivariates
context("Test funcLasso")


context("Test funcLM")
#This test just makes sure the statistical method is being called correctly for one covariate with the correct return
strFormula = "adCur ~ Covariate1"
strRandomFormula = NULL
x = c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
x2 = c(34.2, 32.5, 22.4, 43, 3.25, 6.4, 7, 87, 9)
xf1 = c(1,1,2,2,1,2,1,1,2)
xf2 = c(1,1,1,1,2,2,2,2,2)
frmeTmp = data.frame(Covariate1=x, Covariate2=x2, FCovariate3=xf1, FCovariate4=xf2, adCur=adCur)
iTaxon = 5
lmRet = lm(as.formula(strFormula), data=frmeTmp, na.action = c_strNA_Action)
test_that("Test that the lm has the correct results for 1 covariate.",{
  expect_equal(funcLM(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
})
#Test for correct call for 2 covariates
strFormula = "adCur ~ Covariate1 + Covariate2"
lmRet = lm(as.formula(strFormula), data=frmeTmp, na.action = c_strNA_Action)
test_that("Test that the lm has the correct results for 2 covariates.",{
  expect_equal(funcLM(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
})
##Test for correct call with 1 random and one fixed covariate
#strFormula = "adCur ~ Covariate1"
#strRandomFormula = "~1|FCovariate3"
#lmRet = glmmPQL(fixed=as.formula(strFormula), random=as.formula(strRandomFormula), family=gaussian(link="identity"), data=frmeTmp)
#test_that("Test that the lm has the correct results for 1 random and one fixed covariate.",{
#  expect_equal(funcLM(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
#})
##Test for correct call with 1 random and 2 fixed covariates
#strFormula = "adCur ~ Covariate1 + Covariate2"
#strRandomFormula = "~1|FCovariate3"
#lmRet = glmmPQL(fixed=as.formula(strFormula), random=as.formula(strRandomFormula), family=gaussian(link="identity"), data=frmeTmp)
#test_that("Test that the lm has the correct results for 1 random and 2 fixed covariates.",{
#  expect_equal(funcLM(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
#})
##Test for correct call with 2 random and 1 fixed covariates
#strFormula = "adCur ~ Covariate1"
#strRandomFormula = "~1|FCovariate4+1|FCovariate3"
#lmRet = glmmPQL(fixed=as.formula(strFormula), random=as.formula(strRandomFormula), family=gaussian(link="identity"), data=frmeTmp)
#test_that("Test that the lm has the correct results for 2 random and 1 fixed covariates.",{
#  expect_equal(funcLM(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
#})


context("Test funcBinomialMult")
strFormula = "adCur ~ Covariate1"
strRandomFormula = NULL
x = c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
x2 = c(34.2, 32.5, 22.4, 43, 3.25, 6.4, 7, 87, 9)
xf1 = c(1,1,2,2,1,2,1,1,2)
xf2 = c(1,1,1,1,2,2,2,2,2)
frmeTmp = data.frame(Covariate1=x, Covariate2=x2, FCovariate3=xf1, FCovariate4=xf2, adCur=adCur)
iTaxon = 5
lmRet = glm(as.formula(strFormula), family=binomial(link=logit), data=frmeTmp, na.action=c_strNA_Action)
test_that("Test that the neg binomial regression has the correct results for 1 covariate.",{
  expect_equal(funcBinomialMult(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
})
#Test for correct call for 2 covariates
strFormula = "adCur ~ Covariate1 + Covariate2"
iTaxon = 5
lmRet = glm(as.formula(strFormula), family=binomial(link=logit), data=frmeTmp, na.action=c_strNA_Action)
test_that("Test that the neg binomial regression has the correct results for 2 covariates.",{
  expect_equal(funcBinomialMult(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
})
##Test for correct call with 1 random and one fixed covariate
#strFormula = "adCur ~ Covariate1"
#strRandomFormula = "~1|FCovariate3"
#lmRet = glmmPQL(fixed=as.formula(strFormula), random=as.formula(strRandomFormula), family=binomial(link=logit), data=frmeTmp)
#test_that("Test that the lm has the correct results for 1 random and one fixed covariate.",{
#  expect_equal(funcBinomialMult(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
#})
##Test for correct call with 1 random and 2 fixed covariates
#strFormula = "adCur ~ Covariate1 + Covariate2"
#strRandomFormula = "~1|FCovariate3"
#lmRet = glmmPQL(fixed=as.formula(strFormula), random=as.formula(strRandomFormula), family=binomial(link=logit), data=frmeTmp)
#test_that("Test that the lm has the correct results for 1 random and 2 fixed covariates.",{
#  expect_equal(funcBinomialMult(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
#})
##Test for correct call with 2 random and 1 fixed covariates
#strFormula = "adCur ~ Covariate1"
#strRandomFormula = "~1|FCovariate4+1|FCovariate3"
#lmRet = glmmPQL(fixed=as.formula(strFormula), random=as.formula(strRandomFormula), family=binomial(link=logit), data=frmeTmp)
#test_that("Test that the lm has the correct results for 2 random and 1 fixed covariates.",{
#  expect_equal(funcBinomialMult(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
#})


context("Test funcQuasiMult")
strFormula = "adCur ~ Covariate1"
strRandomFormula = NULL
x = c(44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1,44.4, 45.9, 41.9, 53.3, 44.7, 44.1, 50.7, 45.2, 60.1)
x2 = c(34.2, 32.5, 22.4, 43, 3.25, 6.4, 7, 87, 9,34.2, 32.5, 22.4, 43, 3.25, 6.4, 7, 87, 9)
xf1 = c(1,1,2,2,1,1,2,2,2,1,1,2,2,1,1,2,2,2)
xf2 = c(1,1,1,1,2,2,2,2,2,1,1,1,1,2,2,2,2,2)
frmeTmp = data.frame(Covariate1=x, Covariate2=x2, FCovariate3=xf1, FCovariate4=xf2, adCur=adCur)
iTaxon = 5
lmRet = glm(as.formula(strFormula), family=quasipoisson, data=frmeTmp, na.action=c_strNA_Action)
test_that("Test that the quasi poisson has the correct results for 1 covariate.",{
  expect_equal(funcQuasiMult(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
})
#Test for correct call for 2 covariates
strFormula = "adCur ~ Covariate1 + Covariate2"
iTaxon = 5
lmRet = glm(as.formula(strFormula), family=quasipoisson, data=frmeTmp, na.action=c_strNA_Action)
test_that("Test that the quasi poisson has the correct results for 2 covariates.",{
  expect_equal(funcQuasiMult(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
})
##Test for correct call with 1 random and one fixed covariate
#strFormula = "adCur ~ Covariate1"
#strRandomFormula = "~1|FCovariate3"
#lmRet = glmmPQL(fixed=as.formula(strFormula), random=as.formula(strRandomFormula), family=quasipoisson, data=frmeTmp)
#test_that("Test that the lm has the correct results for 1 random and one fixed covariate.",{
#  expect_equal(funcQuasiMult(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
#})
##Test for correct call with 1 random and 2 fixed covariates
#strFormula = "adCur ~ Covariate1 + Covariate2"
#strRandomFormula = "~1|FCovariate3"
#lmRet = glmmPQL(fixed=as.formula(strFormula), random=as.formula(strRandomFormula), family=quasipoisson, data=frmeTmp)
#test_that("Test that the lm has the correct results for 1 random and 2 fixed covariates.",{
#  expect_equal(funcQuasiMult(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
#})
##Test for correct call with 2 random and 1 fixed covariates
#strFormula = "adCur ~ Covariate1"
#strRandomFormula = "~1|FCovariate4+1|FCovariate3"
#lmRet = glmmPQL(fixed=as.formula(strFormula), random=as.formula(strRandomFormula), family=quasipoisson, data=frmeTmp)
#test_that("Test that the lm has the correct results for 2 random and 1 fixed covariates.",{
#  expect_equal(funcQuasiMult(strFormula=strFormula,frmeTmp=frmeTmp,iTaxon=iTaxon,lsQCCounts=lsQCCounts,strRandomFormula=strRandomFormula),lmRet)
#})


#Test transforms
context("Test funcNoTransform")
aTest1 = c(NA)
aTest2 = c(NULL)
aTest3 = c(0.5,1.4,2.4,3332.4,0.0,0.0000003)
aTest4 = c(0.1)
test_that("Test that no transform does not change the data.",{
  expect_equal(funcNoTransform(aTest1), aTest1)
  expect_equal(funcNoTransform(aTest2), aTest2)
  expect_equal(funcNoTransform(aTest3), aTest3)
  expect_equal(funcNoTransform(aTest4), aTest4)
})


context("Test funcArcsinSqrt")
aTest1 = c(NA)
aTest2 = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
aTest3 = c(0.000001)
test_that("Test that funcArcsinSqrt performs the transform correctly.",{
  expect_equal(funcArcsinSqrt(NA), as.numeric(NA))
  expect_equal(funcArcsinSqrt(aTest1), asin(sqrt(aTest1)))
  expect_equal(funcArcsinSqrt(aTest2), asin(sqrt(aTest2)))
  expect_equal(funcArcsinSqrt(aTest3), asin(sqrt(aTest3)))
})