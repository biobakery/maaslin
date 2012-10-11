
#Arcsin sqrt transform data
processFunction = function( frmeData, aiMetadata, aiGenetics, aiData )
{
  #AsinSqrt data
  for(aiDatum in aiData)
  {
    frmeData[,aiDatum] = asin(sqrt(frmeData[,aiDatum]))
  }
  return( list(frmeData = frmeData, aiMetadata = aiMetadata, aiGenetics = aiGenetics, aiData = aiData, lsQCCounts=list()) )
}
