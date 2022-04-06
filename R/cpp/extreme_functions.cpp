#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

std::vector<double> f_extremeCold(std::vector<double> tmin) {
  std::vector<double> cld(tmin.size(), NAN);
  for (size_t i=0; i<tmin.size(); i++) {
    
    if (std::isnan(tmin[i])) continue;
    
    if (tmin[i] <= - 30) {
      cld[i] =  1;
    } else  {
      cld[i] = 0;
    }
  }
  return cld;
}
