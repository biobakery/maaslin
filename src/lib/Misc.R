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

### Modified Code
### This code is from the package agricolae by Felipe de Mendiburu
### Modifications here are minimal and allow one to use the p.values from the post hoc comparisons
### Authors do not claim credit for this solution only needed to modify code to use the output.
kruskal <- function (y, trt, alpha = 0.05, p.adj = c("none", "holm", "hochberg", 
    "bonferroni", "BH", "BY", "fdr"), group = TRUE, main = NULL) 
{
    dfComparisons=NULL
    dfMeans=NULL
    dntStudent=NULL
    dLSD=NULL
    dHMean=NULL
    name.y <- paste(deparse(substitute(y)))
    name.t <- paste(deparse(substitute(trt)))
    p.adj <- match.arg(p.adj)
    junto <- subset(data.frame(y, trt), is.na(y) == FALSE)
    N <- nrow(junto)
    junto[, 1] <- rank(junto[, 1])
    means <- tapply.stat(junto[, 1], junto[, 2], stat = "sum")
    sds <- tapply.stat(junto[, 1], junto[, 2], stat = "sd")
    nn <- tapply.stat(junto[, 1], junto[, 2], stat = "length")
    means <- data.frame(means, replication = nn[, 2])
    names(means)[1:2] <- c(name.t, name.y)
    ntr <- nrow(means)
    nk <- choose(ntr, 2)
    DFerror <- N - ntr
    rs <- 0
    U <- 0
    for (i in 1:ntr) {
        rs <- rs + means[i, 2]^2/means[i, 3]
        U <- U + 1/means[i, 3]
    }
    S <- (sum(junto[, 1]^2) - (N * (N + 1)^2)/4)/(N - 1)
    H <- (rs - (N * (N + 1)^2)/4)/S
#    cat("\nStudy:", main)
#    cat("\nKruskal-Wallis test's\nTies or no Ties\n")
#    cat("\nValue:", H)
#    cat("\ndegrees of freedom:", ntr - 1)
    p.chisq <- 1 - pchisq(H, ntr - 1)
#    cat("\nPvalue chisq  :", p.chisq, "\n\n")
    DFerror <- N - ntr
    Tprob <- qt(1 - alpha/2, DFerror)
    MSerror <- S * ((N - 1 - H)/(N - ntr))
    means[, 2] <- means[, 2]/means[, 3]
#    cat(paste(name.t, ",", sep = ""), " means of the ranks\n\n")
    dfMeans=data.frame(row.names = means[, 1], means[, -1])
    if (p.adj != "none") {
#        cat("\nP value adjustment method:", p.adj)
        a <- 1e-06
        b <- 1
        for (i in 1:100) {
            x <- (b + a)/2
            xr <- rep(x, nk)
            d <- p.adjust(xr, p.adj)[1] - alpha
            ar <- rep(a, nk)
            fa <- p.adjust(ar, p.adj)[1] - alpha
            if (d * fa < 0) 
                b <- x
            if (d * fa > 0) 
                a <- x
        }
        Tprob <- qt(1 - x/2, DFerror)
    }
    nr <- unique(means[, 3])
    if (group) {
        Tprob <- qt(1 - alpha/2, DFerror)
#        cat("\nt-Student:", Tprob)
#        cat("\nAlpha    :", alpha)
        dntStudent=Tprob
        dAlpha=alpha
        if (length(nr) == 1) {
            LSD <- Tprob * sqrt(2 * MSerror/nr)
#            cat("\nLSD      :", LSD, "\n")
            dLSD=LSD
        }
        else {
            nr1 <- 1/mean(1/nn[, 2])
            LSD1 <- Tprob * sqrt(2 * MSerror/nr1)
#            cat("\nLSD      :", LSD1, "\n")
            dLSD =LSD1
#            cat("\nHarmonic Mean of Cell Sizes ", nr1)
            dHMean=nr1
        }
#        cat("\nMeans with the same letter are not significantly different\n")
#        cat("\nGroups, Treatments and mean of the ranks\n")
        output <- order.group(means[, 1], means[, 2], means[, 
            3], MSerror, Tprob, std.err = sqrt(MSerror/means[, 
            3]))
        dfComparisons=order.group(means[, 1], means[, 2], means[, 
            3], MSerror, Tprob, std.err = sqrt(MSerror/means[, 
            3]))
    }
    if (!group) {
        comb <- combn(ntr, 2)
        nn <- ncol(comb)
        dif <- rep(0, nn)
        LCL <- dif
        UCL <- dif
        pvalue <- dif
        sdtdif <- dif
        for (k in 1:nn) {
            i <- comb[1, k]
            j <- comb[2, k]
            if (means[i, 2] < means[j, 2]) {
                comb[1, k] <- j
                comb[2, k] <- i
            }
            dif[k] <- abs(means[i, 2] - means[j, 2])
            sdtdif[k] <- sqrt(S * ((N - 1 - H)/(N - ntr)) * (1/means[i, 
                3] + 1/means[j, 3]))
            pvalue[k] <- 2 * round(1 - pt(dif[k]/sdtdif[k], DFerror), 
                6)
        }
        if (p.adj != "none") 
            pvalue <- round(p.adjust(pvalue, p.adj), 6)
        LCL <- dif - Tprob * sdtdif
        UCL <- dif + Tprob * sdtdif
        sig <- rep(" ", nn)
        for (k in 1:nn) {
            if (pvalue[k] <= 0.001) 
                sig[k] <- "***"
            else if (pvalue[k] <= 0.01) 
                sig[k] <- "**"
            else if (pvalue[k] <= 0.05) 
                sig[k] <- "*"
            else if (pvalue[k] <= 0.1) 
                sig[k] <- "."
        }
        tr.i <- means[comb[1, ], 1]
        tr.j <- means[comb[2, ], 1]
        dfComparisons <- data.frame(Difference = dif, p.value = pvalue, 
            sig, LCL, UCL)
        rownames(dfComparisons) <- paste(tr.i, tr.j, sep = " - ")
#        cat("\nComparison between treatments mean of the ranks\n\n")
#        print(output)
        dfMeans <- data.frame(trt = means[, 1], means = means[, 
            2], M = "", N = means[, 3])
    }
#    invisible(output)
     invisible(list(study=main,test="Kruskal-Wallis test",value=H,df=(ntr - 1),chisq.p.value=p.chisq,p.adj.method=p.adj,ntStudent=dntStudent,alpha=alpha,LSD=dLSD,Harmonic.mean=dHMean,comparisons=dfComparisons,means=dfMeans))
}

### This function is NOT original code but is from the gamlss package.
### It is written here in an effort to over write the gamlss object summary method
### so that I can return information of interest.
estimatesgamlss<-function (object, Qr, p1, coef.p, 
                           est.disp , df.r, 
                           digits = max(3, getOption("digits") - 3),
                           covmat.unscaled , ...)
{
  #covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
  covmat <- covmat.unscaled #in glm is=dispersion * covmat.unscaled, but here is already multiplied by the dispersion
  var.cf <- diag(covmat)
  s.err <- sqrt(var.cf)
  tvalue <- coef.p/s.err
  dn <- c("Estimate", "Std. Error")
  if (!est.disp) 
  {
    pvalue <- 2 * pnorm(-abs(tvalue))
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), c(dn, "z value","Pr(>|z|)"))
  } else if (df.r > 0) {
    pvalue <- 2 * pt(-abs(tvalue), df.r)
    coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef.p), c(dn, "t value","Pr(>|t|)"))
  } else {
    coef.table <- cbind(coef.p, Inf)
    dimnames(coef.table) <- list(names(coef.p), dn)
  }
  return(coef.table)
}

### This function is NOT original code but is from the gamlss package.
### It is written here in an effort to over write the gamlss object summary method
### so that I can return information of interest.
summary.gamlss<- function (object, type = c("vcov", "qr"), save = FALSE, ...) 
{
  return(as.data.frame(estimatesgamlss(object=object,Qr=object$mu.qr, p1=1:(object$mu.df-object$mu.nl.df), 
    coef.p=object$mu.coefficients[object$mu.qr$pivot[1:(object$mu.df-object$mu.nl.df)]], 
    est.disp =TRUE, df.r=(object$noObs - object$mu.df),
    covmat.unscaled=chol2inv(object$mu.qr$qr[1:(object$mu.df-object$mu.nl.df), 1:(object$mu.df-object$mu.nl.df), drop = FALSE]) )) )
}