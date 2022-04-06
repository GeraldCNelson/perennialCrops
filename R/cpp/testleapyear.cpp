#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::LogicalVector checkLeapYear(Rcpp::IntegerVector year){
  
  // Create an integer vector containing NA
  Rcpp::LogicalVector isleap (year.size(), NAN);
  for (R_xlen_t i=0; i<year.size(); i++) {
  isleap[i] = ( year[i] % 400 == 0 || year[i] % 4 == 0);
  }
  return(isleap);
}
  
  // // If you use the element of LogicalVector directly in the "if" statement
  // // NA_LOGICAL will be evaluated as TRUE
  // for(int i=0; i<v.size();++i) {
  //   if(v[i]) Rprintf("v[%i] is evaluated as true.\n",i);
  //   else Rprintf("v[%i] is evaluated as false.\n",i);
  // }
  // 
  // // Evaluate the elements of LogicalVector
  // for(int i=0; i<v.size();++i) {
  //   if(v[i]==TRUE) Rprintf("v[%i] is TRUE.\n",i);
  //   else if (v[i]==FALSE) Rprintf("v[%i] is FALSE.\n",i);
  //   else if (v[i]==NA_LOGICAL) Rprintf("v[%i] is NA.\n",i);
  //   else Rcout << "v[" << i << "] is not 1\n";
  // }
  
  // Displays the value of TRUE, FALSE and NA_LOGICAL
  // Rcout << "TRUE " << TRUE << "\n";
  // Rcout << "FALSE " << FALSE << "\n";
  // Rcout << "NA_LOGICAL " << NA_LOGICAL << "\n";
  // 
  // return v;
  