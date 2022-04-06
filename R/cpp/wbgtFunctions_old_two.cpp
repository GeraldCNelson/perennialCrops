#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

Rcpp::NumericVector emis_atm(Rcpp::NumericVector Tk, Rcpp::NumericVector RH) {
  Environment env = Environment::global_env();
  
  Rcpp::NumericVector emis_atm_out(Tk.size(), NAN);
  Rcpp::NumericVector e(Tk.size(), NAN);
  Rcpp::NumericVector f = env["esat"];
  for (R_xlen_t i=0; i<Tk.size(); i++) {
    e[i] = RH[i] * f(Tk[i]);
    emis_atm_out[i] = 0.575 * pow(e[i], 0.143);
  }
  return emis_atm_out;
}

// [[Rcpp::export]]

Rcpp::NumericVector esat(Rcpp::NumericVector Tk) {
  Environment env = Environment::global_env();
  
  Rcpp::NumericVector esat_out(Tk.size(), NAN);
  const double kVal = env["kVal"];
  for (R_xlen_t i=0; i<Tk.size(); i++) {
    esat_out[i] = 6.1121 * exp(17.502 * (Tk[i] - kVal) / (Tk[i] - 32.18));
    esat_out[i] = 1.004 * esat_out[i];
  }
  return esat_out;
}

// [[Rcpp::export]]

Rcpp::NumericVector h_sphere_in_air(Rcpp::NumericVector Tk, Rcpp::NumericVector speed) { //speed is wind speed
  Environment env = Environment::global_env();
  
  //variables
  Rcpp::NumericVector therm_con(Tk.size(), NAN);
  Rcpp::NumericVector Re(Tk.size(), NAN);
  Rcpp::NumericVector Nu(Tk.size(), NAN);
  Rcpp::NumericVector h_sphere_in_air_val_out(Tk.size(), NAN);
  Rcpp::NumericVector density(Tk.size(), NAN);
  
  //functions
  Function f_thermal_cond = env["thermal_cond"];
  Rcpp::NumericVector thermal_cond_out = f_thermal_cond(Tk);
  
  Function f_viscosity = env["viscosity"];
  Rcpp::NumericVector viscosity_out = f_viscosity(Tk);
  
  // constants
  const double m_air = env["m_air"];
  const double r_gas = env["r_gas"];
  const double cp = env["cp"];
  const double Pair = env["Pair"];
  const double diam_globe = env["diam_globe"];
  const double min_speed = env["min_speed"];
  double Pr;
  double r_air;
  r_air = r_gas / m_air;
  Pr = cp / (cp + (1.25 * r_air));
  
  for (R_xlen_t i=0; i<Tk.size(); i++) {
    
    therm_con[i] = thermal_cond_out[i];
    density[i] = Pair * 100 / (r_air * Tk[i]);
    if (speed[i] < min_speed) {speed[i] = min_speed;}
    Re[i] = speed[i] * density[i] * diam_globe / viscosity_out[i];
    
    // Nusselt number
    Nu[i] = 2.0 + 0.6 * pow(Re[i], 0.5) * pow(Pr, 0.3333);
    // Convective heat transfer coefficient in W/(m2 K) for a long cylinder in cross flow
    h_sphere_in_air_val_out[i] = Nu[i] * therm_con[i] / diam_globe;
  }
  return h_sphere_in_air_val_out;
}

// [[Rcpp::export]]

Rcpp::NumericVector h_cylinder_in_air(Rcpp::NumericVector Tk, Rcpp::NumericVector speed) { //speed is wind speed
  Environment env = Environment::global_env();
  
  //variables
  Rcpp::NumericVector therm_con(Tk.size(), NAN);
  Rcpp::NumericVector Re(Tk.size(), NAN);
  Rcpp::NumericVector Nu(Tk.size(), NAN);
  Rcpp::NumericVector h_cylinder_in_air_val_out(Tk.size(), NAN);
  Rcpp::NumericVector density(Tk.size(), NAN);
  
  //functions
  Function f_thermal_cond = env["thermal_cond"];
  Rcpp::NumericVector thermal_cond_out = f_thermal_cond(Tk);
  Function f_viscosity = env["viscosity"];
  Rcpp::NumericVector viscosity_out = f_viscosity(Tk);
  
  // constants
  const double m_air = env["m_air"];
  const double r_gas = env["r_gas"];
  const double cp = env["cp"];
  const double Pair = env["Pair"];
  const double diam_wick = env["diam_wick"];
  const double min_speed = env["min_speed"];
  double r_air = r_gas / m_air;
  double Pr = cp / (cp + (1.25 * r_air));
  
  for (R_xlen_t i=0; i<Tk.size(); i++) {
    
    therm_con[i] = thermal_cond_out[i];
    density[i] = Pair * 100 / (r_air * Tk[i]);
    if (speed[i] < min_speed) {speed[i] = min_speed;}
    Re[i] = speed[i] * density[i] * diam_wick / viscosity_out[i];
    
    // Nusselt number
    Nu[i] = 0.281 * pow(Re[i], 0.6) * pow(Pr, 0.44);
    
    // Convective heat transfer coefficient in W/(m2 K) for a long cylinder in cross flow
    h_cylinder_in_air_val_out[i] = Nu[i] * therm_con[i] / diam_wick;
  }
  return h_cylinder_in_air_val_out;
}

// [[Rcpp::export]]

Rcpp::NumericVector viscosity(Rcpp::NumericVector Tk) {
  Rcpp::NumericVector viscosity_out(Tk.size(), NAN);
  Rcpp::NumericVector omega(Tk.size(), NAN);
  Rcpp::NumericVector viscosity(Tk.size(), NAN);
  for (R_xlen_t i=0; i<Tk.size(); i++) {
    omega[i] = (Tk[i] / 97 - 2.9) / 0.4 * (-0.034) + 1.048;
    viscosity_out[i] = 0.0000026693 * pow((28.97 * Tk[i]), 0.5) / (pow(3.617, 2) * omega[i]); // here's the original in R (3.617 ^ 2 * omega)
  }
  return viscosity_out;
}

// [[Rcpp::export]]

Rcpp::NumericVector thermal_cond(Rcpp::NumericVector Tk) {
  Environment env = Environment::global_env();
  Function f_viscosity = env["viscosity"];
  Rcpp::NumericVector viscosity_out = f_viscosity(Tk);
  
  const double cp = env["cp"];
  const double r_air = env["r_air"];
  
  Rcpp::NumericVector thermal_cond_out(Tk.size(), NAN);
  
  for (R_xlen_t i=0; i<Tk.size(); i++) {
    thermal_cond_out[i] = (cp + 1.25 * r_air) * viscosity_out[i];
  }
  return thermal_cond_out;
}

// [[Rcpp::export]]

Rcpp::NumericVector diffusivity(Rcpp::NumericVector Tk) {
  Environment env = Environment::global_env();
  
  Rcpp::NumericVector diffusivity_out(Tk.size(), NAN);
  const double pcrit13 = env["pcrit13"];
  const double tcrit512 = env["tcrit512"];
  const double Tcrit12 = env["Tcrit12"];
  const double Mmix = env["Mmix"];
  const double Pair = env["Pair"];
  for (R_xlen_t i=0; i<Tk.size(); i++) {
    diffusivity_out[i] = 0.000364 * pow((Tk[i] / Tcrit12),  2.334) * pcrit13 * tcrit512 * Mmix / (Pair / 1013.25) * 0.0001;
  }
  return diffusivity_out;
}

// [[Rcpp::export]]

Rcpp::NumericVector fr(Rcpp::NumericVector Twb_prev, Rcpp::NumericVector Tair, Rcpp::NumericVector emis_atm_out, Rcpp::NumericVector Tsfc, Rcpp::NumericVector radiation, Rcpp::NumericVector zenith, Rcpp::NumericVector eair, Rcpp::NumericVector density, Rcpp::NumericVector speed, Rcpp::NumericVector RH) {
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
  
  Rcpp::NumericVector fr_out(Twb_prev.size(), NAN);
  Rcpp::NumericVector Tref(Twb_prev.size(), NAN);
  Rcpp::NumericVector Fatm(Twb_prev.size(), NAN);
  Rcpp::NumericVector Sc(Twb_prev.size(), NAN);
  Rcpp::NumericVector h(Twb_prev.size(), NAN);
  Rcpp::NumericVector ewick(Twb_prev.size(), NAN);
  Rcpp::NumericVector evap(Twb_prev.size(), NAN);
  Rcpp::NumericVector Twb(Twb_prev.size(), NAN);
  
  // called functions
  Function f_esat = env["esat"];
  Rcpp::NumericVector f_esat_out = f_esat(Twb_prev);
  Function f_viscosity = env["viscosity"];
  Rcpp::NumericVector viscosity_tair_out = f_viscosity(Tair);
  Function f_diffusivity = env["diffusivity"];
  Rcpp::NumericVector diffusivity_tref_out = f_esat(Tref);
  Function f_h_cylinder_in_air = env["h_cylinder_in_air"];
  Rcpp::NumericVector h_cylinder_in_air_out = f_h_cylinder_in_air(Twb_prev, speed);
  Function f_h_evap = env["h_evap"];
  Rcpp::NumericVector h_evap_out = f_h_evap(Twb_prev);
  
  // Function f_emis_atm = env["emis_atm"];
  // Rcpp::NumericVector emis_atm_out = f_emis_atm(Tk, RH);
  
  
  for (R_xlen_t i=0; i<radiation.size(); i++) {
    Tref[i] = 0.5 * (Twb_prev[i] + Tair[i]);
    Fatm[i] = stefanb * emis_wick * (0.5 * (emis_atm_out[i] * pow(Tair[i], 4) + emis_sfc * pow(Tsfc[i], 4)) - pow(Twb_prev[i], 4)) + (1.0 - alb_wick) * radiation[i] * ((1.0 - propDirect) * (1.0 + 0.25 * diam_wick/len_wick) + ((std::tan(zenith[i])/3.1416) + 0.25 * diam_wick/len_wick) * propDirect + SurfAlbedo);
    Sc[i] = viscosity_tair_out[i]/(density[i] * diffusivity_tref_out[i]);
    h[i] = h_cylinder_in_air_out(Twb_prev[i], speed[i]);
    ewick[i] = f_esat_out[i];
    evap[i] = h_evap_out[i];
    Twb[i] = Tair[i] - evap[i]/ratio * (ewick[i] - eair[i])/(Pair - ewick[i]) * pow((Pr/Sc[i]), 0.56) + Fatm[i]/h[i] * irad;
    fr_out[i] = abs(Twb[i] - Twb_prev[i]);
  }
  return fr_out;
}

  // [[Rcpp::export]]
  
  Rcpp::NumericVector fr_1(Rcpp::NumericVector Tglobe_prev, Rcpp::NumericVector Tair, Rcpp::NumericVector speed, 
                   Rcpp::NumericVector RH, Rcpp::NumericVector radiation, Rcpp::NumericVector cza) {
    Environment env = Environment::global_env();
    
    // variables
    Rcpp::NumericVector Tsfc(Tglobe_prev.size(), NAN);
    Rcpp::NumericVector Tref(Tglobe_prev.size(), NAN);
    Rcpp::NumericVector h(Tglobe_prev.size(), NAN);
    Rcpp::NumericVector Tglobe(Tglobe_prev.size(), NAN);
    Rcpp::NumericVector fr_1_out(Tglobe_prev.size(), NAN);
    
    const double emis_globe = env["emis_globe"];
    const double emis_sfc = env["emis_sfc"];
    const double stefanb = env["stefanb"];
    const double alb_globe = env["alb_globe"];
    const double propDirect = env["propDirect"];
    const double SurfAlbedo = env["SurfAlbedo"];
    
    // called functions
    Function f_h_sphere_in_air = env["h_sphere_in_air"];
    Rcpp::NumericVector h_sphere_in_air_out =f_h_sphere_in_air(Tref, speed);
    Function f_emis_atm = env["emis_atm"];
    Rcpp::NumericVector emis_atm_out =f_emis_atm(Tair, RH);
    for (R_xlen_t i=0; i<radiation.size(); i++) {
      
  Tsfc = Tair;
  Tref = 0.5 * (Tglobe_prev + Tair);
  h = h_sphere_in_air_out;
  Tglobe[i] = pow(0.5 * (emis_atm_out[i] * pow(Tair[i], 4.0) + emis_sfc * pow(Tsfc[i], 4.0)) - h[i]/(emis_globe * stefanb) * (Tglobe_prev[i] - Tair[i]) + radiation[i]/(2.0 * emis_globe * stefanb) * (1.0 - alb_globe) * (propDirect * (1.0/(2.0 * cza[i]) - 1.0) + 1.0 + SurfAlbedo), 0.25);
  fr_1_out[i] = abs(Tglobe[i] - Tglobe_prev[i]);
}
return fr_1_out;
}

//      Rcpp::Rcout << i << std::endl;


// [[Rcpp::export]]
Rcpp::NumericVector h_evap(Rcpp::NumericVector Tk) {
  Rcpp::NumericVector h_evap_out(Tk.size(), NAN);
  for (R_xlen_t i=0; i<Tk.size(); i++) {
    
    h_evap_out[i] = (313.15 - Tk[i]) / 30.0 * (-71100.0) + 2407300.0;
  }
  return h_evap_out;
}







// Rcpp::NumericVector h_cylinder_in_air(Rcpp::NumericVector Tk, Rcpp::NumericVector speed) {
//  Environment env = Environment::global_env();
//   
//   //variables
//   Rcpp::NumericVector therm_con(Tk.size(), NAN);
//   Rcpp::NumericVector Re(Tk.size(), NAN);
//   Rcpp::NumericVector Nu(Tk.size(), NAN);
//   Rcpp::NumericVector h_cylinder_in_air_val_out(Tk.size(), NAN);
//   Rcpp::NumericVector density(Tk.size(), NAN);
//   //functions
//   Rcpp::NumericVector f_thermal_cond = env["thermal_cond"];
//   Rcpp::NumericVector f_viscosity = env["viscosity"];
//   
//   constants
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
//   for (R_xlen_t i=0; i<Tk.size(); i++) {
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
// Rcpp::NumericVector alb_globe = env["alb_globe"];
// Rcpp::NumericVector alb_wick = env["alb_wick"];
// Rcpp::NumericVector cp = env["cp"];
// Rcpp::NumericVector cp = env["cp"];
// Rcpp::NumericVector diam_globe = env["diam_globe"];
// Rcpp::NumericVector diam_wick = env["diam_wick"];
// Rcpp::NumericVector emis_globe = env["emis_globe"];
// Rcpp::NumericVector emis_sfc = env["emis_sfc"];
// Rcpp::NumericVector emis_wick = env["emis_wick"];
// Rcpp::NumericVector len_wick = env["len_wick"];
// Rcpp::NumericVector m_air = env["m_air"];
// Rcpp::NumericVector m_h2o = env["m_h2o"];
// Rcpp::NumericVector r_gas = env["r_gas"];
// Rcpp::NumericVector stefanb = env["stefanb"];
// Rcpp::NumericVector r_air = env["r_air"];
// Rcpp::NumericVector Pr = env["Pr"];
// Rcpp::NumericVector ratio = env["ratio"];
// Rcpp::NumericVector SurfAlbedo = env["SurfAlbedo"];
// Rcpp::NumericVector irad = env["irad"];
// Rcpp::NumericVector kVal = env["kVal"];
// Rcpp::NumericVector EQTIME1 = env["EQTIME1"];
// Rcpp::NumericVector EQTIME2 = env["EQTIME2"];
// Rcpp::NumericVector EQTIME3 = env["EQTIME3"];
// Rcpp::NumericVector EQTIME4 = env["EQTIME4"];
// Rcpp::NumericVector EQTIME5 = env["EQTIME5"];
// Rcpp::NumericVector EQTIME6 = env["EQTIME6"];
// Rcpp::NumericVector DECL1 = env["DECL1"];
// Rcpp::NumericVector DECL2 = env["DECL2"];
// Rcpp::NumericVector DECL3 = env["DECL3"];
// Rcpp::NumericVector DECL4 = env["DECL4"];
// Rcpp::NumericVector DECL5 = env["DECL5"];
// Rcpp::NumericVector DECL6 = env["DECL6"];
// Rcpp::NumericVector DECL7 = env["DECL7"];
// 
// Rcpp::NumericVector propDirect = env["propDirect"];
// Rcpp::NumericVector Pair = env["Pair"];
//




