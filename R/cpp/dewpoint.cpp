#include <Rcpp.h>
// using namespace Rcpp;
// Source: https://bmcnoldy.rsmas.miami.edu/Humidity.html
//TD: =243.04*(LN(RH/100)            +((17.625*T)/    (243.04+T)))/     (17.625-LN(RH/100)-       ((17.625*T)/(243.04+T))) 
  
// TtDewp = 243.04 * (log(hurs/100) + ((17.625 * tmp)/(243.04 + tmp))) /(17.625 - log(hurs/100) - ((17.625 * tmp)/(243.04 + tmp))) 

// [[Rcpp::export]]
Rcpp::NumericVector f_tDewp(Rcpp::NumericVector hurs, Rcpp::NumericVector tmp) {
  Rcpp::NumericVector tDew(hurs.size(), Rcpp::NumericVector::get_na());
 hurs = log(hurs/100.0);
 tDew = 243.04 * (hurs + ((17.625 * tmp)/(243.04 + tmp))) /(17.625 - hurs) - ((17.625 * tmp)/(243.04 + tmp));
   return tDew;
}

