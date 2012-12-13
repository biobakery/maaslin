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

inlinedocs <- function(
##author<< Curtis Huttenhower <chuttenh@hsph.harvard.edu> and Timothy Tickle <ttickle@hsph.harvard.edu>
##description<< Holds MaAsLin related plotting
) { return( pArgs ) }

funcPDF <- function(
### Function to plot raw data with linear model information.
### Continuous and integer variables are plotted with a line of best fit.
### Other data is plotted as boxplots.
frmeTmp,
lsCur,
### Linear model information
curPValue,
### Pvalue to display
curQValue,
### Qvalue to display
strFilePDF,
### PDF file to create or to which to append
strBaseOut,
### Project directory to place pdf in
strName,
### Name of taxon
fInvert
### Invert the figure so the background is black
){
  if( is.na( strFilePDF ) )
  {
    strFilePDF <- sprintf( "%s-%s.pdf", strBaseOut, strName )
    pdf( strFilePDF, width = 11, useDingbats=FALSE )
  }
  
  #Invert plots
  adColorMin <- c(1, 0, 0)
  adColorMax <- c(0, 1, 0)
  adColorMed <- c(0, 0, 0)
  if( fInvert )
  {
    par( bg = "black", fg = "white", col.axis = "white", col.lab = "white", col.main = "white", col.sub = "white" )
    adColorMin <- c(1, 1, 0)
    adColorMax <- c(0, 1, 1)
    adColorMed <- c(1, 1, 1)
  }

  #Create linear model title data string
  print("ls")
  print(lsCur)
  strTitle <- sprintf( "%s (%.3g sd %.3g, p=%.3g, q=%.3g)", lsCur$orig, lsCur$value, lsCur$std, curPValue, curQValue )
  adMar <- c(5, 4, 4, 2) + 0.1
  dLine <- NA
  strTaxon <- lsCur$taxon
  if( nchar( strTaxon ) > 80 )
  {
    dCEX <- 0.75
    iLen <- nchar( strTaxon )
    if( iLen > 120 )
    {
      dLine <- 2.5
      i <- round( iLen / 2 )
      strTaxon <- paste( substring( strTaxon, 0, i ), substring( strTaxon, i + 1 ), sep = "\n" )
      adMar[2] <- adMar[2] + 1
    }
  } else { dCEX = 1 }

  #Plot 1x2 graphs per page if there are multiple coefficients
  iNumberCoefficients = length(setdiff(names(lsCur$allCoefs),c("(Intercept)")))
  if(iNumberCoefficients>1){ par(mfrow=c(1,2)) }else{ par(mfrow=c(1,1)) }

  # Plot factor data as boxplot if is descrete data
  # Otherwise plot as a line
  adCur <- lsCur$metadata
  adData <- lsCur$data
  adY <- adData
  if( class( lsCur$metadata ) == "factor" )
  {
    if( "NA" %in% levels( lsCur$metadata ) )
    {
      afNA <- adCur == "NA"
      adY <- adY[!afNA]
      adData <- adData[!afNA]
      adCur <- adCur[!afNA]
      adCur <- factor( adCur, levels = setdiff( levels( adCur ), "NA" ) )
    }
    astrNames <- c()
    astrColors <- c()
    dMed <- median( adY[adCur == levels( adCur )[1]], na.rm = TRUE )
    adIQR <- quantile( adY, probs = c(0.25, 0.75), na.rm = TRUE )
    dIQR <- adIQR[2] - adIQR[1]
    if( dIQR <= 0 )
    {
      dIQR <- sd( adY, na.rm = TRUE )
    }
    dIQR <- dIQR / 2

    #Print boxplots/strip charts of raw data. Add model data to it.
    for( strLevel in levels( adCur ) )
    {
      c_iLen <- 32
      strLength <- strLevel
      if( nchar( strLength ) > c_iLen )
      {
        iTmp <- ( c_iLen / 2 ) - 2
        strLength <- paste( substr( strLength, 1, iTmp ), substring( strLength, nchar( strLength ) - iTmp ), sep = "..." )
      }
      astrNames <- c(astrNames, sprintf( "%s (%d)", strLength, sum( adCur == strLevel, na.rm = TRUE ) ))
      astrColors <- c(astrColors, sprintf( "%sAA", funcColor( ( median( adY[adCur == strLevel], na.rm = TRUE ) - dMed ) /
        dIQR, dMax = 3, dMin = -3, adMax = adColorMin, adMin = adColorMax, adMed = adColorMed ) ))
    }
    #Controls boxplot labels
    #(When there are many factor levels some are skipped and not plotted
    #So this must be reduced)
    dBoxPlotLabelCex = dCEX
    if(length(astrNames)>8)
    {
      dBoxPlotLabelCex = dBoxPlotLabelCex * 1.5/(length(astrNames)/8)
    }
    par(cex.axis = dBoxPlotLabelCex)
    boxplot( adY ~ adCur, notch = TRUE, names = astrNames, mar = adMar, col = astrColors,
      main = strTitle, xlab = lsCur$name, ylab = NA, cex.lab = dCEX, outpch = 4, outcex = 0.5 )
    par(cex.axis = dCEX)
    stripchart( adY ~ adCur, add = TRUE, col = astrColors, method = "jitter", vertical = TRUE, pch = 20 )
    title( ylab = strTaxon, cex.lab = dCEX, line = dLine )
  } else {
    #Plot continuous data
    plot( adCur, adY, mar = adMar, main = strTitle, xlab = lsCur$name, pch = 20,
      col = sprintf( "%s99", funcGetColor( ) ), ylab = NA, xaxt = "s" )
    title( ylab = strTaxon, cex.lab = dCEX )
    lmod <- lm( adY ~ adCur )
    dColor <- lmod$coefficients[2] * mean( adCur, na.rm = TRUE ) / mean( adY, na.rm = TRUE )
    strColor <- sprintf( "%sDD", funcColor( dColor, adMax = adColorMin, adMin = adColorMax, adMed = adColorMed ) )
    abline( reg = lmod, col = strColor, lwd = 3 )
  }

  ### Plot the residual plot
#  funcResidualPlot(lsCur=lsCur, frmeTmp= frmeTmp)

  return(strFilePDF)
  ### File to which the pdf was written
}

funcResidualPlot <- function(
### Plot to data after confounding.
### That is, in a linear model with significant coefficient b1 for variable x1,
### that's been sparsified to some subset of terms: y = b0 + b1*x1 + sum(bi*xi)
### Plot x1 on the X axis, and instead of y on the Y axis, instead plot:
### y' = b0 + sum(bi*xi)
lsCur,
frmeTmp
){
  print("lsCur")
  print(lsCur)
  #Now plot residual hat plot
  #Get coefficient names
  lsAllCoefs = setdiff(names(lsCur$allCoefs),c("(Intercept)"))

  #All coefficients except for the one of interest
  lsOtherCoefs = setdiff(lsAllCoefs, c(lsCur$orig))

  #If there are no other coefficients then skip plot
  if(!length(lsOtherCoefs)){return()}

  #Get xi (raw data)
  mtrxRawData = as.matrix(lsCur$metadata)

  lmod = lm(as.formula(paste( "mtrxRawData ~", paste( sprintf( "`%s`", lsOtherCoefs ), collapse = " + " ))),frmeTmp)

  #Plot
  plot(bi ~ mtrxRawData, xlab = lsCur$orig, ylab = paste(lsOtherCoefs,sep="", collapse="+"), main = paste(lsCur$taxon,"~",paste(lsOtherCoefs,sep="+")), pch = 20)
  rug(mtrxRawData, side=1)
  rug(bi, side=2)
}