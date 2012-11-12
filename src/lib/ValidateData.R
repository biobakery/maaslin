#######################################################################################
# This file is provided under the Creative Commons Attribution 3.0 license.
#
# You are free to share, copy, distribute, transmit, or adapt this work
# PROVIDED THAT you attribute the work to the authors listed below.
# For more information, please see the following web page:
# http://creativecommons.org/licenses/by/3.0/
#
# This file is a component of the MaAsLin (Multivariate Associations Using Linear Models), 
# authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Timothy Tickle, ttickle@hsph.harvard.edu).
#######################################################################################

inlinedocs <- function(
##author<< Curtis Huttenhower <chuttenh@hsph.harvard.edu> and Timothy Tickle <ttickle@hsph.harvard.edu>
##description<< Minor validation files to check data typing when needed.
) { return( pArgs ) }

funcIsValid <- function(
### Requires a data to not be NA, not be NULL
### Returns True on meeting these requirements, returns false otherwise
### Return boolean Indicator of not being empty (TRUE = not empty)
tempData = NA
### Parameter tempData Is evaluated as not empty
){
  #If the data is not na or null return true
  if(!is.null(tempData))
  {
    if(length(tempData)==1){ return(!is.na(tempData)) }
    return(TRUE)
  }
  return(FALSE)
  ### True (Valid) false (invalid)
}

funcIsValidString <- function(
### Requires a data to not be NA, not be NULL, and to be of type Character
### Returns True on meeting these requirements, returns false otherwise
### Return boolean Indicator of identity as a string
tempData = NA
### Parameter tempData Is evaluated as a string
){
  #If is not a valid data return false
  if(!funcIsValid(tempData))
  {
    return(FALSE)
  }
  #If is a string return true
  if((class(tempData)=="character")&&(length(tempData)==1))
  {
    return(TRUE)
  }
  return(FALSE)
  ### True (Valid) false (invalid)
}

funcIsValidFileName <- function(
### Requires a data to not be NA, not be NULL, and to be a valid string
### which points to an existing file
### Returns True on meeting these requirements, returns false otherwise
### Return boolean Indicator of identity as a file name
tempData = NA,
### Parameter tempData Is evaluated as a file name
fVerbose=FALSE
### Verbose will print the file path when not valid.
){
  #If is not valid string return false
  if(!(funcIsValidString(tempData)))
  {
    if(fVerbose){print(paste("FunctIsValidFileName: InvalidString. Value=",tempData,sep=""))}
    return(FALSE)
  }
  #If is a valid string and points to a file
  if(file.exists(tempData))
  {
    return(TRUE)
  }
  if(fVerbose){print(paste("FunctIsValidFileName: Path does not exist. Value=",tempData,sep=""))}
  return(FALSE)
  ### True (Valid) false (invalid)
}