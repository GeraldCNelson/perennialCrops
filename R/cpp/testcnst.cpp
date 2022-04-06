#include <Rcpp.h>
const int testcnst2 = 5; 
// [[Rcpp::export]]
Rcpp::LogicalVector is_leapyear(Rcpp::IntegerVector year){
  // Create an integer vector containing NAN
  const int testcnst = 455; 
  Rcpp::LogicalVector isleap (year.size(), NAN);
  for (R_xlen_t i=0; i<year.size(); i++) {
    isleap[i] = ( year[i] % testcnst == 0 || year[i] % testcnst2 == 0);
  }
  return(isleap);
}