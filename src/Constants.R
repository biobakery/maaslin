####################################
# Summary: Constants
# Author: Timothy Tickle
# Start Date: 11-02-2011
####################################

#General
c_COMMA = ","
c_DASH = "-"

#For reading IO
c_MATRIX_NAME = "Matrix:"
c_FILE_NAME = "File:"
c_DELIMITER = "Delimiter:"
c_ID_ROW = "Name_Row_Number:"
c_ID_COLUMN = "Name_Column_Number:"
c_ROWS = "Read_Rows"
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
# The column that is used to determine if information meets a certain significance threshold (dSignificanceLevel) to include in the Summary text file)
c_strKeywordEvaluatedForInclusion <- "Q.value"
#The name of the custom process function
c_strCustomProcessFunction = "processFunction"

#Delimiters
#Feature name delimiter
c_cFeatureDelim = "|"
c_cFeatureDelimREx = "\\|"

#The word used for unclassified
c_strUnclassified = "unclassified"

#Outlier related constants
c_dFence <- 0
c_dPOutlier <- 0.05

#Maaslincore settings
#If a metadata does not have more than count of unique values, it is changed to factor data mode.
c_iNonFactorLevelThreshold = 3

#Extensions
c_DETAIL_FILE_SUFFIX = ".txt"

#Delimiter for output tables
c_cTableDelimiter="\t"

#Testing Related
c_strTestingDirectory = "testing"
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