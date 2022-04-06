#include <Rcpp.h>

using namespace Rcpp;  

// // [[Rcpp::export]]
// NumericVector vecpow(const NumericVector base, const NumericVector exp) {
//   NumericVector out(base.size());
//   std::transform(base.begin(), base.end(),
//                  exp.begin(), out.begin(), static_cast<double(*)(double, double)>(::pow));
//   return out;
// }

// [[Rcpp::export]]

NumericVector vecpow(NumericVector base, NumericVector exp) {
  NumericVector out(base.size());
  for (R_xlen_t i=0; i<base.size(); i++) {
    out[i] = pow(base[i], exp[i]);
  }
  return out;
}

/*** R
vecpow(c(0:3), c(4:1))
***/
