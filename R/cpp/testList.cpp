// [[Rcpp::export]]
Rcpp::List fTg_prep(Rcpp::NumericVector tas, Rcpp::NumericVector hurs, Rcpp::NumericVector radiation, Rcpp::NumericVector zenith) {
  Rcpp::NumericVector cza(tas.size(), NAN);
  
  for (R_xlen_t i=0; i<tas.size(); i++) {
    if (zenith[i] <= 0) {
      zenith[i] = 1e-10;
    } else if (radiation[i] > 0 & zenith[i] > 1.57) {
      zenith[i] = 1.57;
    } else if (radiation[i] > 15 & zenith[i] > 1.54) {
      zenith[i] = 1.54;
    } else if (radiation[i] > 900 & zenith[i] > 1.52) { 
      zenith[i] = 1.52;
    } else if (radiation[i] < 10 & zenith[i] == 1.57) { 
      radiation[i] = 0;
    }
    cza[i] = cos(zenith[i]);
  }
 return Rcpp::List::create(Rcpp::Named("cza") = cza, Rcpp::Named("radiation") = radiation);
}