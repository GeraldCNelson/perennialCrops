#include <Rcpp.h>

using namespace Rcpp;  

// [[Rcpp::export]]
double dpow(const double lhs, const double rhs) {
  return std::pow(lhs, rhs);
}

// [[Rcpp::export]]
std::vector<double> mypow(const std::vector<double>& base, const std::vector<double>& exp) {
  std::vector<double> res(base.size());
  std::transform(base.begin(), base.end(), exp.begin(), res.begin(), dpow);
  return res;
}

/*** R
mypow(c(0:3), c(4:1))
***/

// // [[Rcpp::export]]
// NumericVector vecpow(const NumericVector base, const NumericVector exp) {
//   NumericVector out(base.size());
//   std::transform(base.begin(), base.end(),
//                  exp.begin(), out.begin(), ::pow);
//   return out;
// }
// 
// /*** R
// vecpow(c(0:3), c(4:1))
// ***/