#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

std::vector<double> fr(std::vector<double> Twb_prev, std::vector<double> Tair, std::vector<double> emis_atm, std::vector<double> Tsfc, std::vector<double> radiation, std::vector<double> zenith, std::vector<double> eair, std::vector<double> density, std::vector<double> wind) {
  Environment env = Environment::global_env();
  
  const double emis_wick = env["emis_wick"];
  const double emis_sfc = env["emis_sfc"];
  const double alb_wick = env["alb_wick"];
  const double diam_wick = env["diam_wick"];
  const double len_wick = env["len_wick"];
  const double propDirect = env["propDirect"];
  const double SurfAlbedo = env["SurfAlbedo"];
  const double stefanb = env["stefanb"];
  const double cp = env["cp"];
  const double r_air = env["r_air"];
  const double irad = env["irad"];
  const double ratio = env["ratio"];
  const double Pair = env["Pair"];
  
  double Pr;
  Pr = cp / (cp + (1.25 * r_air));
  
  NumericVector xxx(Twb_prev.size(), NAN);
  
  std::vector<double> fr_out(Twb_prev.size(), NAN);
  std::vector<double> Tref(Twb_prev.size(), NAN);
  std::vector<double> Fatm(Twb_prev.size(), NAN);
  std::vector<double> Sc(Twb_prev.size(), NAN);
  std::vector<double> h(Twb_prev.size(), NAN);
  std::vector<double> ewick(Twb_prev.size(), NAN);
  std::vector<double> evap(Twb_prev.size(), NAN);
  std::vector<double> Twb(Twb_prev.size(), NAN);
  
  // called functions
  NumericVector f_esat = env["esat"];
  NumericVector f_viscosity = env["viscosity"];
  NumericVector f_diffusivity = env["diffusivity"];
  NumericVector f_h_cylinder_in_air = env["h_cylinder_in_air"];
  NumericVector f_h_evap = env["h_evap"];
  
  for (size_t i=0; i<radiation.size(); i++) {
    Tref[i] = 0.5 * (Twb_prev[i] + Tair[i]);
    Fatm[i] = stefanb * emis_wick * (0.5 * (emis_atm[i] * pow(Tair[i], 4) + emis_sfc * pow(Tsfc[i], 4)) - pow(Twb_prev[i], 4)) + (1.0 - alb_wick) * radiation[i] * ((1.0 - propDirect) * (1.0 + 0.25 * diam_wick/len_wick) + ((std::tan(zenith[i])/3.1416) + 0.25 * diam_wick/len_wick) * propDirect + SurfAlbedo);
    Sc[i] = f_viscosity(Tair[i])/(density[i] * f_diffusivity(Tref[i]));
    h[i] = f_h_cylinder_in_air(Twb_prev[i], wind[i]);
    ewick[i] = f_esat(Twb_prev[i]);
    evap[i] = f_h_evap(Twb_prev[i]);
    Twb[i] = Tair[i] - evap[i]/ratio * (ewick[i] - eair[i])/(Pair - ewick[i]) * pow((Pr/Sc[i]), 0.56) + Fatm[i]/h[i] * irad;
    fr_out[i] = abs(Twb[i] - Twb_prev[i]);
  }
  return fr_out;
}

// [[Rcpp::export]]

NumericVector h_cylinder_in_air(std::vector<double> Tk, std::vector<double> speed) { //speed is wind speed
   Environment env = Environment::global_env();
  
  //variables
  std::vector<double> therm_con(Tk.size(), NAN);
  std::vector<double> Re(Tk.size(), NAN);
  std::vector<double> Nu(Tk.size(), NAN);
  std::vector<double> h_cylinder_in_air_val_out(Tk.size(), NAN);
  std::vector<double> density(Tk.size(), NAN);
  
  //functions
  NumericVector f_thermal_cond = env["thermal_cond"];
  NumericVector f_viscosity = env["viscosity"];
  
  // Constants
  const double m_air = env["m_air"];
  const double r_gas = env["r_gas"];
  const double cp = env["cp"];
  const double Pair = env["Pair"];
  const double diam_wick = env["diam_wick"];
  const double min_speed = env["min_speed"];
  double Pr;
  double r_air;
  r_air = r_gas / m_air;
  Pr = cp / (cp + (1.25 * r_air));
  
  for (size_t i=0; i<Tk.size(); i++) {
    
    therm_con[i] = f_thermal_cond(Tk[i]);
    density[i] = Pair * 100 / (r_air * Tk[i]);
    if (speed[i] < min_speed) {speed[i] = min_speed;}
    Re[i] = speed[i] * density[i] * diam_wick / f_viscosity(Tk[i]);
    
    // Nusselt number
    Nu[i] = 0.281 * pow(Re[i], 0.6) * pow(Pr, 0.44);
    
    // Convective heat transfer coefficient in W/(m2 K) for a long cylinder in cross flow
    h_cylinder_in_air_val_out[i] = Nu[i] * therm_con[i] / diam_wick;
  }
  return h_cylinder_in_air_val_out;
}



//      Rcpp::Rcout << i << std::endl;


// [[Rcpp::export]]
std::vector<double> h_evap(std::vector<double> Tk) {
  std::vector<double> h_evap_out(Tk.size(), NAN);
  for (size_t i=0; i<Tk.size(); i++) {
    
    h_evap_out[i] = (313.15 - Tk[i]) / 30.0 * (-71100.0) + 2407300.0;
  }
  return h_evap_out;
}

// [[Rcpp::export]]

std::vector<double> viscosity(std::vector<double> Tk) {
  std::vector<double> viscosity_out(Tk.size(), NAN);
  std::vector<double> omega(Tk.size(), NAN);
  std::vector<double> viscosity(Tk.size(), NAN);
  for (size_t i=0; i<Tk.size(); i++) {
    omega[i] = (Tk[i] / 97 - 2.9) / 0.4 * (-0.034) + 1.048;
    viscosity_out[i] = 0.0000026693 * pow((28.97 * Tk[i]), 0.5) / (pow(3.617, 2) * omega[i]); // here's the original in R (3.617 ^ 2 * omega)
  }
  return viscosity_out;
}

// [[Rcpp::export]]

std::vector<double> emis_atm(std::vector<double> Tk, std::vector<double> RH) {
  Environment env = Environment::global_env();
  
  std::vector<double> emis_atm_out(Tk.size(), NAN);
  std::vector<double> e(Tk.size(), NAN);
  NumericVector f = env["esat"];
  for (size_t i=0; i<Tk.size(); i++) {
    e[i] = RH[i] * f(Tk[i]);
    emis_atm_out[i] = 0.575 * pow(e[i], 0.143);
  }
  return emis_atm_out;
}

// [[Rcpp::export]]

std::vector<double> esat(std::vector<double> Tk) {
  Environment env = Environment::global_env();
  
  std::vector<double> esat_out(Tk.size(), NAN);
  const double kVal = env["kVal"];
  for (size_t i=0; i<Tk.size(); i++) {
    esat_out[i] = 6.1121 * exp(17.502 * (Tk[i] - kVal) / (Tk[i] - 32.18));
    esat_out[i] = 1.004 * esat_out[i];
  }
  return esat_out;
}

// [[Rcpp::export]]

std::vector<double> diffusivity(std::vector<double> Tk) {
  Environment env = Environment::global_env();
  
  std::vector<double> diffusivity_out(Tk.size(), NAN);
  const double pcrit13 = env["pcrit13"];
  const double tcrit512 = env["tcrit512"];
  const double Tcrit12 = env["Tcrit12"];
  const double Mmix = env["Mmix"];
  const double Pair = env["Pair"];
  for (size_t i=0; i<Tk.size(); i++) {
    diffusivity_out[i] = 0.000364 * pow((Tk[i] / Tcrit12),  2.334) * pcrit13 * tcrit512 * Mmix / (Pair / 1013.25) * 0.0001;
  }
  return diffusivity_out;
}

// [[Rcpp::export]]

std::vector<double> thermal_cond(std::vector<double> Tk) {
  Environment env = Environment::global_env();
  NumericVector f_viscosity = env["viscosity"];
  const double cp = env["cp"];
  const double r_air = env["r_air"];
  
  std::vector<double> thermal_cond_out(Tk.size(), NAN);
  
  for (size_t i=0; i<Tk.size(); i++) {
    thermal_cond_out[i] = (cp + 1.25 * r_air) * f_viscosity(Tk[i]);
  }
  return thermal_cond_out;
}

// // [[Rcpp::export]]
// 
// std::vector<double> h_sphere_in_air(std::vector<double> Tk, std::vector<double> speed) { //speed is wind speed
//   Environment env = Environment::global_env();
//   
//   //variables
//   std::vector<double> therm_con(Tk.size(), NAN);
//   std::vector<double> Re(Tk.size(), NAN);
//   std::vector<double> Nu(Tk.size(), NAN);
//   std::vector<double> h_sphere_in_air_val_out(Tk.size(), NAN);
//   std::vector<double> density(Tk.size(), NAN);
//   
//   //functions
//   NumericVector f_thermal_cond = env["thermal_cond"];
//   NumericVector f_viscosity = env["viscosity"];
//   
//   // Constants
//   const double m_air = env["m_air"];
//   const double r_gas = env["r_gas"];
//   const double cp = env["cp"];
//   const double Pair = env["Pair"];
//   const double diam_globe = env["diam_globe"];
//   const double min_speed = env["min_speed"];
//   double Pr;
//   double r_air;
//   r_air = r_gas / m_air;
//   Pr = cp / (cp + (1.25 * r_air));
//   
//   for (size_t i=0; i<Tk.size(); i++) {
//     
//     therm_con[i] = f_thermal_cond(Tk[i]);
//     density[i] = Pair * 100 / (r_air * Tk[i]);
//     if (speed[i] < min_speed) {speed[i] = min_speed;}
//     Re[i] = speed[i] * density[i] * diam_globe / f_viscosity(Tk[i]);
//     
//     // Nusselt number
//     Nu[i] = 2.0 + 0.6 * pow(Re[i], 0.5) * pow(Pr, 0.3333);
//     // Convective heat transfer coefficient in W/(m2 K) for a long cylinder in cross flow
//     h_sphere_in_air_val_out[i] = Nu[i] * therm_con[i] / diam_globe;
//   }
//   return h_sphere_in_air_val_out;
// }


// std::vector<double> h_cylinder_in_air(std::vector<double> Tk, std::vector<double> speed) {
//   Environment env = Environment::global_env();
//   
//   //variables
//   std::vector<double> therm_con(Tk.size(), NAN);
//   std::vector<double> Re(Tk.size(), NAN);
//   std::vector<double> Nu(Tk.size(), NAN);
//   std::vector<double> h_cylinder_in_air_val_out(Tk.size(), NAN);
//   std::vector<double> density(Tk.size(), NAN);
//   //functions
//   NumericVector f_thermal_cond = env["thermal_cond"];
//   NumericVector f_viscosity = env["viscosity"];
//   
//   // Constants
//   const double m_air = env["m_air"];
//   const double r_gas = env["r_gas"];
//   const double cp = env["cp"];
//   const double Pair = env["Pair"];
//   const double diam_wick = env["diam_wick"];
//   const double min_speed = env["min_speed"];
//   double Pr;
//   double r_air;
//   Pr = cp / (cp + (1.25 * r_air));
//   r_air = r_gas / m_air;
//   
//   for (size_t i=0; i<Tk.size(); i++) {
//     
//     therm_con[i] = f_thermal_cond(Tk[i]);
//     density[i] = Pair * 100 / (r_air * Tk[i]);
//     if (speed[i] < min_speed) {speed[i] = min_speed;}
//     Re[i] = speed[i] * density[i] * diam_wick / f_viscosity(Tk[i]);
//     
//     // Nusselt number
//     Nu[i] = 0.281 * pow(Re[i], 0.6) * pow(Pr, 0.44);
//     
//     // Convective heat transfer coefficient in W/(m2 K) for a long cylinder in cross flow
//     h_cylinder_in_air_val_out[i] = Nu[i] * therm_con[i] / diam_wick;
//   }
//   return h_cylinder_in_air_val_out;
// }


//
// NumericVector alb_globe = env["alb_globe"];
// NumericVector alb_wick = env["alb_wick"];
// NumericVector cp = env["cp"];
// NumericVector cp = env["cp"];
// NumericVector diam_globe = env["diam_globe"];
// NumericVector diam_wick = env["diam_wick"];
// NumericVector emis_globe = env["emis_globe"];
// NumericVector emis_sfc = env["emis_sfc"];
// NumericVector emis_wick = env["emis_wick"];
// NumericVector len_wick = env["len_wick"];
// NumericVector m_air = env["m_air"];
// NumericVector m_h2o = env["m_h2o"];
// NumericVector r_gas = env["r_gas"];
// NumericVector stefanb = env["stefanb"];
// NumericVector r_air = env["r_air"];
// NumericVector Pr = env["Pr"];
// NumericVector ratio = env["ratio"];
// NumericVector SurfAlbedo = env["SurfAlbedo"];
// NumericVector irad = env["irad"];
// NumericVector kVal = env["kVal"];
// NumericVector EQTIME1 = env["EQTIME1"];
// NumericVector EQTIME2 = env["EQTIME2"];
// NumericVector EQTIME3 = env["EQTIME3"];
// NumericVector EQTIME4 = env["EQTIME4"];
// NumericVector EQTIME5 = env["EQTIME5"];
// NumericVector EQTIME6 = env["EQTIME6"];
// NumericVector DECL1 = env["DECL1"];
// NumericVector DECL2 = env["DECL2"];
// NumericVector DECL3 = env["DECL3"];
// NumericVector DECL4 = env["DECL4"];
// NumericVector DECL5 = env["DECL5"];
// NumericVector DECL6 = env["DECL6"];
// NumericVector DECL7 = env["DECL7"];
// 
// NumericVector propDirect = env["propDirect"];
// NumericVector Pair = env["Pair"];
//




