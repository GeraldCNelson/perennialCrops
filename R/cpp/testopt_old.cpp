#include <Rcpp.h>
const double tolerance = 1e-0;

[[Rcpp::export]]
Rcpp::NumericVector f(Rcpp::NumericVector x) {
  return Rcpp::pow(x-0.33, 2);
}
// 
// Rcpp::List fTg_opt(const double optmin, const double optmax) {
//   Rcpp::Environment base("package:stats");
//   Rcpp::Function optimize_r = base["optimize"];    
//   Rcpp::NumericVector interval = {optmin,optmax};
//   return optimize_r(f, interval, tolerance);
// }

//#include <Rcpp.h>
// 
// // [[Rcpp::export]]
// Rcpp::List fTg_opt(Rcpp::Function f, const double optmin, const double optmax) {
//   Rcpp::Environment base("package:stats");
//   Rcpp::Function optimize_r = base["optimize"];
//   Rcpp::NumericVector interval = {optmin,optmax};
//   Rcpp::List res = optimize_r(f, interval);
//   return res;
// }

/*** R
f <- function (x) (x - .33)^2
xmin <- optimize(f, c(0, 1), tol = 0.0001)
xmin

fTg_opt(f, 0, 1)
*/