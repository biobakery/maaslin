
#Arcsin sqrt transform data
processFunction = function( frmeData, aiMetadata, aiGenetics, aiData, funcTransform )
{
  #Transform data
  for(aiDatum in aiData)
  {
    frmeData[,aiDatum] = funcTransform(frmeData[,aiDatum])
  }
  return( list(frmeData = frmeData, aiMetadata = aiMetadata, aiGenetics = aiGenetics, aiData = aiData, lsQCCounts=list()) )
}
