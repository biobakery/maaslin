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
##description<< Global project constants.
) { return( pArgs ) }

#General
c_COMMA = ","
c_DASH = "-"

#For reading IO
c_MATRIX_NAME = "Matrix:"
c_FILE_NAME = "File:"
c_DELIMITER = "Delimiter:"
c_ID_ROW = "Name_Row_Number:"
c_ID_COLUMN = "Name_Column_Number:"
c_ROWS = "Read_Rows:"
c_PCLROWS = "Read_PCL_Rows:"
c_TSVROWS = "Read_TSV_Rows:"
c_COLUMNS = "Read_Columns:"
c_PCLCOLUMNS = "Read_PCL_Columns:"
c_TSVCOLUMNS = "Read_TSV_Columns:"
c_CHARACTER_DATA_TYPE = "DT_Character:"
c_FACTOR_DATA_TYPE = "DT_Factor:"
c_INTEGER_DATA_TYPE = "DT_Integer:"
c_LOGICAL_DATA_TYPE = "DT_Logical:"
c_NUMERIC_DATA_TYPE = "DT_Numeric:"
c_ORDEREDFACTOR_DATA_TYPE = "DT_Ordered_Factor:"

### The name of the data matrix read in using a read.config file
c_strMatrixData	 <- "Abundance"
### The name of the metadata matrix read in using a read.config file
c_strMatrixMetadata <- "Metadata"
# Settings for MFA visualization/ordination
c_iMFA <- 30
c_dHeight <- 9
c_dDefaultScale = 0.5
# The column that is used to determine if information meets a certain significance threshold (dSignificanceLevel) to include in the Summary text file)
c_strKeywordEvaluatedForInclusion <- "Q.value"
#The name of the custom process function
c_strCustomProcessFunction = "processFunction"

#Delimiters
#Feature name delimiter
c_cFeatureDelim = "|"
c_cFeatureDelimRex = "\\|"

#The word used for unclassified
c_strUnclassified = "unclassified"

#Maaslincore settings
#If a metadata does not have more than count of unique values, it is changed to factor data mode.
c_iNonFactorLevelThreshold = 3

#Extensions
c_sDetailFileSuffix = ".txt"
c_sSummaryFileSuffix = ".txt"
c_sLogFileSuffix = "_log"

#Delimiter for output tables
c_cTableDelimiter="\t"

#Testing Related
c_strCorrectAnswers = "answers"
c_strTemporaryFiles = "tmp"
c_strTestingInput = "input"

#Reading matrix defaults
c_strDefaultMatrixDelimiter = "\t"
c_strDefaultMatrixRowID = "1"
c_strDefaultMatrixColID = "1"
c_strDefaultReadRows = "-"
c_strDefaultReadCols = "-"

#Separator used when collapsing factor names
c_sFactorNameSep = ""

#Separator used by the mfa
c_sMFANameSep1 = "_"
c_sMFANameSep2 = "."

#Analysis Module list positioning
c_iSelection = 1
c_iTransform = 2
c_iAnalysis = 3
c_iResults = 4
c_iUnTransform = 5
c_iIsUnivariate = 6

#Count based models
c_vCountBasedModels = c("neg_binomial","quasi")

# Na action in anaylsis, placed here to standardize
c_strNA_Action = "na.omit"
