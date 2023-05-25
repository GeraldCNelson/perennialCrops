//' Chill portion calculator
//' 
//' \code{chill_func} returns the chill portion count of a temperature vector
//' 
//' @param temp a vector of temperautres \emph{in Celsius}
//' 
//' @return The output is the chill portion calculation. Portions are rounded
//' down to the largest integer.
//' @export

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

double chill_func(DoubleVector temp) {
  /* initial setups */
  DoubleVector tempvec = temp + 273;
  DoubleVector ftmprt = 1.6 * 277 * (tempvec - 277) / tempvec;
  DoubleVector st = exp(ftmprt);
  DoubleVector xi = st / (1 + st);
  
  double bignum1 = -30.5434496955383;
  DoubleVector  log_xs = bignum1 + (12888.8 - 4153.5) / tempvec;
  DoubleVector  xs = exp(log_xs);
 
  DoubleVector  log_ak1 =  42.3892695757818 - 12888.8 / tempvec;
  DoubleVector  ak1 = exp(log_ak1); 
 
  /*initializing the vectors*/
  tempvec.insert(0, 285);
  tempvec.insert(0, 288);
  ftmprt.insert(0, 12.44); 
  ftmprt.insert(0, 16.93);
  st.insert(0, 252887.94);
  st.insert(0, 22471935.51);
  xi.insert(0,1);
  xi.insert(0,1);
  xs.insert(0, 1.11);
  xs.insert(0, 0.81);
  ak1.insert(0, 0.06);
  ak1.insert(0, 0.09);
  
  int n = tempvec.size();
  
  DoubleVector  inter_S(n);
  inter_S[0] = 0.0;
  DoubleVector inter_E(n);
  inter_E[0] = 0.0;
  DoubleVector  delt(n);
  delt[0] = 0.0;
  
  double portions = 0.0;
  inter_E[0] = xs[0] * (xs[0] - inter_S[0]) * exp(-ak1[0]); 
  
  /* chill generating loop */
  for(int i = 1; i < n; ++i) {
    if (inter_E[i-1] < 1) {
      inter_S[i] = inter_E[i-1];
    } else {
      inter_S[i] = inter_E[i-1] * (1 - xi[i - 1]);
    }
    
    inter_E[i] = xs[i] - (xs[i] - inter_S[i]) * exp(-ak1[i]);
    
    if (inter_E[i] < 1){
      delt[i] = 0;
    } else {
      delt[i] = inter_E[i] * xi[i];
    }
    
    portions = portions + delt[i];
    
  }
  
  return floor(portions);
}
