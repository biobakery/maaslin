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
