#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

std::vector<double> f_gdd(std::vector<double> tmp, double topt_min, double topt_max) {
  std::vector<double> gdd(tmp.size(), NAN);
  
  for (size_t i=0; i<tmp.size(); i++) {
    if (std::isnan(tmp[i])) continue;
    //        Rcpp::Rcout << i << std::endl;
    double max1 = std::min(tmp[i], topt_max) - topt_min;
    //     double max2 = std::max(0.0, max1);
    //    gdd[i] = max2-topt_min;
    gdd[i] = std::max(0.0, max1);
  }
  return gdd;
}
