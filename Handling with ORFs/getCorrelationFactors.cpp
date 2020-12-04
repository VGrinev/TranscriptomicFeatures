#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
SEXP getCorrelationFactors (SEXP string, SEXP pattern, SEXP blockSize)
{
  //  input data
  Rcpp::CharacterVector strVect(string);
  Rcpp::CharacterVector patVect(pattern);
  Rcpp::NumericVector partSizeVect(blockSize);
  int rowsNum = strVect.size();
  int partSize = partSizeVect[0];
  char firstNucleotide = patVect[0][0];
  char secondNucleotide = patVect[0][1];
  //  output vector
  Rcpp:NumericVector result(rowsNum);
  //  foreach sequence in source data
  for (int i = 0; i < rowsNum; i ++)
  {
    int stringSize = strlen(strVect[i]);
    result[i] = 0;
    //  calculate number of parts -1
    int partsNum = stringSize / partSize;
    //  foreach part (apart from last one)
    int firstMatches = 0;
    int secondMatches = 0;
    for (int currPart = 0; currPart < partsNum; currPart ++)
    {
      //  variables to calculate frequencies of nucleotides
      firstMatches = 0;
      secondMatches = 0;
      for (int currChar = currPart * partSize; currChar < (currPart+1) * partSize; currChar ++)
      {
        if(strVect[i][currChar] == firstNucleotide)
          firstMatches ++;
        if(strVect[i][currChar] == secondNucleotide)
          secondMatches ++;
      }
      //  caculate sum of absolute value of frequencies difference
      result[i] += ((double)abs(firstMatches - secondMatches)) / partSize;
    }
    //  last part processing if it exists
    if (stringSize % partSize != 0)
    {
      firstMatches = 0;
      secondMatches = 0;
      for (int currChar = partsNum * partSize; currChar < stringSize; currChar ++)
      {
        if(strVect[i][currChar] == firstNucleotide)
          firstMatches ++;
        if(strVect[i][currChar] == secondNucleotide)
          secondMatches ++;
      }
      result[i] += ((double)abs(firstMatches - secondMatches)) / (stringSize % partSize);
      result[i] = ((double)result[i]) / (partsNum + 1);
    }
    else
      result[i] = ((double)result[i]) / partsNum;
  }
  return result;
}
