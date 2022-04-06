#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector esat(NumericVector Tk) {
  NumericVector esat_out(Tk.size(), NAN);
  for (size_t i=0; i<Tk.size(); i++) {
    esat_out[i] = 6.1121 * Tk[i];
  }
  return esat_out;
}

// [[Rcpp::export]]
NumericVector h_evap(NumericVector Tk) {
  // Environment env = Environment::global_env();
  // NumericVector f_esat = env["esat"];
  NumericVector h_evap_out(Tk.size(), NAN);
  for (size_t i=0; i<Tk.size(); i++) {
    
    h_evap_out[i] = (313.15 - Tk[i]);
 //   NumericVector temp = esat(Tk[i]);
 //      Rcpp::Rcout << temp << std::endl;
    
    h_evap_out[i] = 313.15 + esat(Tk[i]);
  }
  
  return h_evap_out;
}



// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
h_evap(42)
*/
