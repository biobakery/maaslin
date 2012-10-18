

funcPDF <- function( lsCur, curPValue, curQValue, strFilePDF, strBaseOut, strName, fInvert )
{
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
    fGenetics <- length( aiGenetics ) && ( class( adCur ) == "integer" ) &&
      length( intersect( astrFactors, colnames( frmeData )[aiGenetics] ) )
    if( fGenetics )
    {
      astrLabels <- c()
      for( i in 0:2 )
      {
        astrLabels <- c(astrLabels, sprintf( "%d (%d)", i, sum( adCur == i, na.rm = TRUE ) ))
      }
      adCur <- adCur + rnorm( length( adCur ), sd = 0.05 )
    }
    plot( adCur, adY, mar = adMar, main = strTitle, xlab = lsCur$name, pch = 20,
      col = sprintf( "%s99", funcGetColor( ) ), ylab = NA, xaxt = ifelse( fGenetics, "n", "s" ) )
    if( fGenetics )
    {
      axis( 1, at = 0:2, labels = astrLabels )
    }
    title( ylab = strTaxon, cex.lab = dCEX )
    lmod <- lm( adY ~ adCur )
    dColor <- lmod$coefficients[2] * mean( adCur, na.rm = TRUE ) / mean( adY, na.rm = TRUE )
    strColor <- sprintf( "%sDD", funcColor( dColor, adMax = adColorMin, adMin = adColorMax, adMed = adColorMed ) )
    abline( reg = lmod, col = strColor, lwd = 3 )
  }
  return(strFilePDF)
}

#funcPlotDeconfounded <- function()
#{
      #Now plot residual hat plot
      #Get coefficient names
#      lsAllCoefs = setdiff(names(lsCur$allCoefs),c("(Intercept)"))

      #Current b1 coefficient
#      sCurSigData = lsCur$orig

      #Other coefficients names
#      lsOtherCoefs = setdiff(lsAllCoefs, c(sCurSigData))

      #Get xi (raw data)
#      b1 = as.matrix(adData)

      #Get bi coefficients * data
#      bi = as.matrix(adData) %*% t(as.matrix(lsCur$allCoefs[lsOtherCoefs]))
#      bi = bi %*% as.matrix(rep(1,ncol(bi)))

      #Plot
#      plot(bi ~ b1, xlab = sCurSigData, ylab = paste(lsOtherCoefs,sep="", collapse="+"), main = paste(strTaxon,"~",strTitle), pch = 20)
#      rug(b1, side=1)
#      rug(bi, side=2)
#}
