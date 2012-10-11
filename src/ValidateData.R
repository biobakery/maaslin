##########################
# Summary: Code used to validate the existence and correct typing of variable
# Author: Timothy Tickle
# Date: 09-30-2011
##########################

# Requires a data to not be NA, not be NULL
# Returns True on meeting these requirements, returns false otherwise
# Parameter tempData Is evaluated as not empty
# Return boolean Indicator of not being empty (TRUE = not empty)
funcIsValid <- function(tempData = NA)
{
  #If the data is not na or null return true
  if(!is.null(tempData))
  {
    if(length(tempData)==1){ return(!is.na(tempData)) }
    return(TRUE)
  }
  return(FALSE)
}

# Requires a data to not be NA, not be NULL, and to be of type Character
# Returns True on meeting these requirements, returns false otherwise
# Parameter tempData Is evaluated as a string
# Return boolean Indicator of identity as a string
funcIsValidString <- function(tempData = NA)
{
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
}

# Requires a data to not be NA, not be NULL, and to be a valid string
# which points to an existing file
# Returns True on meeting these requirements, returns false otherwise
# Parameter tempData Is evaluated as a file name
# Return boolean Indicator of identity as a file name
funcIsValidFileName <- function(tempData = NA, fVerbose=FALSE)
{
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
}