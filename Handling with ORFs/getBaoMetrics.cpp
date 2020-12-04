#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP getBaoMetrics (SEXP sequences, SEXP indexes, SEXP types)
{
  Rcpp::CharacterVector strVect(sequences);
  Rcpp::CharacterVector typesVect(types);
  Rcpp::NumericVector indVect(indexes);
  
  int rowsNum = indVect.size();
  Rcpp::DoubleVector HRR(rowsNum);
  Rcpp::DoubleVector HRY(rowsNum);
  Rcpp::DoubleVector HYR(rowsNum);
  Rcpp::DoubleVector HYY(rowsNum);
  Rcpp::DoubleVector HMM(rowsNum);
  Rcpp::DoubleVector HMK(rowsNum);
  Rcpp::DoubleVector HKM(rowsNum);
  Rcpp::DoubleVector HKK(rowsNum);
  Rcpp::DoubleVector HWW(rowsNum);
  Rcpp::DoubleVector HWS(rowsNum);
  Rcpp::DoubleVector HSW(rowsNum);
  Rcpp::DoubleVector HSS(rowsNum);
  
  for (int i = 0; i < rowsNum; i ++)
  {
    int stringSize = strlen(strVect[i]);
    
    //	first step - transform source sequence into 3 sequences with RY, MK and WS characters
    char** tmpStrings = new char*[3];
    tmpStrings[0] = new char[stringSize];
    tmpStrings[1] = new char[stringSize];
    tmpStrings[2] = new char[stringSize];
    for (int j = 0; j < stringSize; j ++)
    {
      if (strVect[i][j] == 'A')
      {
        tmpStrings[0][j] = 'R';
        tmpStrings[1][j] = 'M';
        tmpStrings[2][j] = 'W';
      }
      if (strVect[i][j] == 'T')
      {
        tmpStrings[0][j] = 'Y';
        tmpStrings[1][j] = 'K';
        tmpStrings[2][j] = 'W';
      }
      if (strVect[i][j] == 'G')
      {
        tmpStrings[0][j] = 'R';
        tmpStrings[1][j] = 'K';
        tmpStrings[2][j] = 'S';
      }
      if (strVect[i][j] == 'C')
      {
        tmpStrings[0][j] = 'Y';
        tmpStrings[1][j] = 'M';
        tmpStrings[2][j] = 'S';
      }
    }
    
    //  second step - transform 3 previous sequences into 12 numeric sequences of local frequences
    Rcpp::CharacterVector patterns = Rcpp::CharacterVector::create
      ("RR","RY","YR","YY","MM","MK","KM","KK","WW","WS","SW","SS");
    
    //  12 because we have 12 Bao metrics
    double** tmpNumbers = new double*[12];
    int* tmpPrevPosition = new int[12];
    for (int k = 0; k < 12; k ++)
    {
      tmpPrevPosition[k] = 0;
      tmpNumbers[k] = new double[stringSize];
    }

    // for all characters
    for (int j = 0; j < stringSize - 1; j ++)
    {
      //  for all patterns
      for (int k = 0; k < 12; k ++)
      {
        //  if current characters matches the current pattern
        if ((tmpStrings[k/4][j] == patterns[k][0]) && (tmpStrings[k/4][j+1] == patterns[k][1]))
        {
          //  calculate local frequency as 1 / (current position - position of previous match)
          //  +1 cause first index of array in C++ is 0, not 1
          tmpNumbers[k][j] = (double) 1 / (j - tmpPrevPosition[k] + 1);
          tmpPrevPosition[k] = j + 1;
        }
        //  if not matched, local frequency is null
        else
          tmpNumbers[k][j] = 0;
      }
    }
    
    //  third step - transform previous sequences into partitial sums 
    //  and normalize on total sum of all partitial sums 
    double** tmpPartSum = new double*[12];
    double* totalSum = new double[12];
    for (int k = 0; k < 12; k ++)
    {
      tmpPartSum[k] = new double[stringSize - 1];
      tmpPartSum[k][0] = tmpNumbers[k][0];
      totalSum[k] = tmpPartSum[k][0];
    }
    for (int j = 1; j < stringSize - 1; j ++)
    {
      for (int k = 0; k < 12; k ++)
      {
        tmpPartSum[k][j] = tmpPartSum[k][j-1] + tmpNumbers[k][j];
        totalSum[k] += tmpPartSum[k][j];
      }
    }
    for (int j = 0; j < stringSize - 1; j ++)
    {
      for (int k = 0; k < 12; k ++)
      {
        //  to not divide by zero
        if (totalSum[k] != 0)
          tmpPartSum[k][j] = (double) tmpPartSum[k][j] / totalSum[k];
      }
    }
    
    //  final step - calculate Shannon's entropies of partitial sum sequence
    double* entropies = new double[12];
    for (int k = 0; k < 12; k ++)
      entropies[k] = 0.0;
    for (int j = 0; j < stringSize - 1; j ++) 
    {
      for (int k = 0; k < 12; k ++)
      {
        //  log(x)*x in limit to 0 is 0, so if tmpPartSum[k][j] == 0 we shoudln't add anything
        if (tmpPartSum[k][j] != 0)
        {
          entropies[k] = entropies[k] - tmpPartSum[k][j] * log(tmpPartSum[k][j]) / log(2.0);
        }
      }
    }
    HRR[i] = entropies[0];
    HRY[i] = entropies[1];
    HYR[i] = entropies[2];
    HYY[i] = entropies[3];
    HMM[i] = entropies[4];
    HMK[i] = entropies[5];
    HKM[i] = entropies[6];
    HKK[i] = entropies[7];
    HWW[i] = entropies[8];
    HWS[i] = entropies[9];
    HSW[i] = entropies[10];
    HSS[i] = entropies[11];
    
    
    for (int k = 0; k < 12; k ++)
    {
      delete tmpPartSum[k];
      delete tmpNumbers[k];
    }
    delete tmpStrings[0];
    delete tmpStrings[1];
    delete tmpStrings[2];
    delete[] tmpStrings;
    delete[] tmpPartSum;
    delete[] tmpNumbers;
    delete[] entropies;
    delete[] totalSum;
    delete[] tmpPrevPosition;
  }
  
  //  output data.frame
  Rcpp::DataFrame result = Rcpp::DataFrame::create
    (Rcpp::Named("index")=indVect,
     Rcpp::Named("types")=typesVect,
     Rcpp::Named("HRR")=HRR,
     Rcpp::Named("HRY")=HRY,
     Rcpp::Named("HYR")=HYR,
     Rcpp::Named("HYY")=HYY,
     Rcpp::Named("HMM")=HMM,
     Rcpp::Named("HMK")=HMK,
     Rcpp::Named("HKM")=HKM,
     Rcpp::Named("HKK")=HKK,
     Rcpp::Named("HWW")=HWW,
     Rcpp::Named("HWS")=HWS,
     Rcpp::Named("HSW")=HSW,
     Rcpp::Named("HSS")=HSS);
  return result;
}