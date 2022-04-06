//#include <iostream>
#include <Rcpp.h>

using namespace Rcpp;  

const double alb_globe = 0.05;
const double alb_wick = 0.4;
const double cp = 1003.5;
const double diam_globe = 0.0508;
const double diam_wick = 0.007;
const double emis_globe = 0.95;
const double emis_sfc = 0.999  ;
const double emis_wick = 0.95;
const double len_wick = 0.0254;
const double m_air = 28.97;
const double m_h2o = 18.015;
const double r_gas = 8314.34;
const double stefanb = 0.000000056696;
const double SurfAlbedo = 0.4;
const double alb_sfc = SurfAlbedo;
const double irad = 1.0;
const double kVal = 273.15;
const double EQTIME1 = 229.18;
const double EQTIME2 = 7.5e-05;
const double EQTIME3 = 0.001868;
const double EQTIME4 = 0.032077;
const double EQTIME5 = 0.014615;
const double EQTIME6 = 0.040849;
const double DECL1 = 0.006918;
const double DECL2 = 0.399912;
const double DECL3 = 0.070257;
const double DECL4 = 0.006758;
const double DECL5 = 0.000907;
const double DECL6 = 0.002697;
const double DECL7 = 0.00148;

const double propDirect = 0.8; //Assume a proportion of direct radiation = direct/(diffuse + direct)
const double Pair = 1010; // Atmospheric pressure in hPa

const double pcrit13 =  19.94585; //(36.4 * 218) ^ (1 / 3);
const double tcrit512 = 113.4662; //(132 * 647.3) ^ (5 / 12);
const double Tcrit12 = 292.3074; //(132 * 647.3) ^ 0.5;
const double Mmix =  0.3000463; //(1 / 28.97 + 1 / 18.015) ^ 0.5
const double min_speed = 0.1;

const double r_air = r_gas / m_air;
const double Pr = cp / (cp + (1.25 * r_air));
const double ratio = cp * m_air/m_h2o;
const double pi = 3.14159265358979323846;
const double tolerance = 1e-0;
const double utc_hour = 12.0;
const double TimeOffset = 0.0;

const double CONVERGENCE = 0.02;
const double MAX_ITER =	50;

// // [[Rcpp::export]]
// NumericVector vecpow(const NumericVector base, const NumericVector exp) {
//   NumericVector out(base.size());
//   std::transform(base.begin(), base.end(),
//                  exp.begin(), out.begin(), static_cast<double(*)(double, double)>(::pow));
//   return out;
// }

// [[Rcpp::export]]
NumericVector vecpow(NumericVector base, NumericVector exp) {
  NumericVector out(base.size());
  for (R_xlen_t i=0; i<base.size(); i++) {
    out[i] = pow(base[i], exp[i]);
  }
//  Rcpp::Rcout << "out = " << out << std::endl << "-------------" << std::endl;
  
  return out;
}



// [[Rcpp::export]]
Rcpp::NumericVector degToRad(Rcpp::NumericVector angleDeg) {
  Rcpp::NumericVector degToRad(angleDeg.size(), Rcpp::NumericVector::get_na());
  return(pi * angleDeg / 180.0);
}

// [[Rcpp::export]]
Rcpp::NumericVector radToDeg(Rcpp::NumericVector angleRad) {
  Rcpp::NumericVector degToRad(angleRad.size(), Rcpp::NumericVector::get_na());
  return(180 * angleRad / pi);
}

// [[Rcpp::export]]
Rcpp::IntegerVector get_years(const Rcpp::DateVector& dates) {
  Rcpp::IntegerVector out(dates.size(), Rcpp::NumericVector::get_na());
  for (R_xlen_t i=0; i< dates.size(); i++) {
    Rcpp::Date d = dates[i];
    out[i] = d.getYear();
    // Rcpp::Rcout << "dates[i] = " << dates[i] << std::endl << "-------------" << std::endl;
    // Rcpp::Rcout << "out[i] = " << out[i] << std::endl << "-------------" << std::endl;
  }
  return out;
}

// // [[Rcpp::export]]
// Rcpp::LogicalVector is_leapyear(Rcpp::DateVector& year){
//   Rcpp::LogicalVector isleap (year.size(), Rcpp::NumericVector::get_na());
//   Rcpp::IntegerVector year_num(year.size(), Rcpp::NumericVector::get_na());
//   
//   year_num = get_years(year);
//   for (R_xlen_t i=0; i<year.size(); i++) {
//     isleap[i] = ( year_num[i] % 400 == 0 || year_num[i] % 4 == 0);
//     // Rcpp::Rcout << "year[i] = " << year[i] << std::endl << "-------------" << std::endl;
//     // Rcpp::Rcout << "isleap[i] = " << isleap[i] << std::endl << "-------------" << std::endl;
//   }
//   return(isleap);
// }

// [[Rcpp::export]]
Rcpp::IntegerVector get_days(const Rcpp::DateVector& dates) {
  Rcpp::IntegerVector out(dates.size(), Rcpp::NumericVector::get_na());
  for (R_xlen_t i=0; i< dates.size(); i++) {
    Rcpp::Date d = dates[i];
    out[i] = d.getYearday();
    // Rcpp::Rcout << "dates[i] = " << dates[i] << std::endl << "-------------" << std::endl;
    // Rcpp::Rcout << "out[i] = " << out[i] << std::endl << "-------------" << std::endl;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector calZenith(const Rcpp::DateVector& dates, Rcpp::NumericVector lon, Rcpp::NumericVector lat) { //Calculate zenith angle in degrees
  Rcpp::NumericVector d1(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::IntegerVector year(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::IntegerVector year_num(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector doy(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector dpy(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector RadLon(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector RadLat(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Gamma(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector EquTime(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Decli(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector CosZen(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector SZARad(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector SZA(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector TrueSolarTime(dates.size(), (utc_hour * 60.0) + TimeOffset);
  Rcpp::NumericVector HaDeg(dates.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector HaRad(dates.size(), Rcpp::NumericVector::get_na());
  
  year_num = get_years(dates);
  doy = get_days(dates);
  // Rcpp::Rcout << "dates.size() = " << dates.size() << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "year_num = " << year_num << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "doy = " << doy << std::endl << "-------------" << std::endl;
  
  // check for leap year
  for (R_xlen_t i=0; i<dates.size(); i++) {
    if( year_num[i] % 400 == 0 || year_num[i] % 4 == 0) {
      dpy[i] = 366;
    } else { dpy[i] = 365;
      // Rcpp::Rcout << "dpy[i] = " << dpy[i] << std::endl << "-------------" << std::endl;
      // Rcpp::Rcout << "i = " << i << std::endl << "-------------" << std::endl;
    }
    // Rcpp::Rcout << "year_num = " << year_num << std::endl << "-------------" << std::endl;
    // Rcpp::Rcout << "dpy = " << dpy << std::endl << "-------------" << std::endl;
  }
//  return(dpy);

RadLon = degToRad(lon);
RadLat = degToRad(lat);

// Rcpp::Rcout << "lon = " << lon << std::endl << "-------------" << std::endl;
//  Rcpp::Rcout << " lat = " << lat << std::endl << "-------------" << std::endl;

Gamma = 2 * pi * ((doy - 1.0) + (utc_hour/24.0))/dpy; //Evaluate the fractional year in radians 
// Rcpp::Rcout << "Gamma = " << Gamma << std::endl << "-------------" << std::endl;
// Rcpp::Rcout << "doy = " << doy << std::endl << "-------------" << std::endl;
// Rcpp::Rcout << "dpy = " << doy << std::endl << "-------------" << std::endl;
// Rcpp::Rcout << "utc_hour = " << utc_hour << std::endl << "-------------" << std::endl;
// Rcpp::Rcout << "pi = " << pi << std::endl << "-------------" << std::endl;

EquTime = EQTIME1 * (EQTIME2 + EQTIME3 * cos(Gamma) - EQTIME4 * sin(Gamma) - EQTIME5 * cos(2 * Gamma) - EQTIME6 * sin(2 * Gamma)); //Evaluate the Equation of time in minutes 
Decli = DECL1 - DECL2 * cos(Gamma) + DECL3 * sin(Gamma) - DECL4 * cos(2 * Gamma) + DECL5 * sin(2 * Gamma) - DECL6 * cos(3 * Gamma) + DECL7 * sin(3 * Gamma); //Evaluate the solar declination angle in radians (must be between -23.5 and 23.5 degrees)
// Note: TrueSolarTime defined in variable def above
HaDeg = ((TrueSolarTime/4.0) - 180.0);
HaRad = degToRad(HaDeg);
//  Rcpp::Rcout << "HaDeg = " << HaDeg;
CosZen = (sin(RadLat) * sin(Decli) + cos(RadLat) * cos(Decli) * cos(HaRad));
// Rcpp::Rcout << "CosZen = " << CosZen;
//  Rcpp::Rcout << "RadLat = " << RadLat;
// Rcpp::Rcout << "Decli = " << Decli << std::endl << "-------------" << std::endl;
// Rcpp::Rcout << "Gamma = " << doy;
// Rcpp::Rcout << "doy = " << doy;
// Rcpp::Rcout << "HaDeg = " << HaDeg;
// Rcpp::Rcout << "HaRad = " << HaRad;
// Rcpp::Rcout << "TrueSolarTime = " << TrueSolarTime;
for (R_xlen_t i=0; i<year.size(); i++) {
  if(CosZen[i] > 1.0 ) {
    CosZen[i] = 1.0;
  } else if(CosZen[i] < -1.0){
    CosZen[i] = 1.0;
  }
}
SZARad = acos(CosZen);
SZA = radToDeg(SZARad);
return(SZA);
}

// [[Rcpp::export]]
Rcpp::NumericVector h_evap(Rcpp::NumericVector Tk) { //heat of evaporation, J/(kg K), for temperature in the range 283-313 K.
  Rcpp::NumericVector h_evap_out(Tk.size(), Rcpp::NumericVector::get_na());
  h_evap_out = (313.15 - Tk) / 30.0 * (-71100.0) + 2407300.0;
  return h_evap_out;
}

// [[Rcpp::export]]
Rcpp::NumericVector esat(Rcpp::NumericVector Tk) { //saturation vapor pressure (hPa) over water.  Reference: Buck's (1981) approximation (eqn 3) of Wexler's (1976) formulae over liquid water.
  Rcpp::NumericVector esat_out(Tk.size(), Rcpp::NumericVector::get_na());
  
  esat_out = 1.004 * 6.1121 * exp(17.502 * (Tk - kVal) / (Tk - 32.18));
  return esat_out;
}

// [[Rcpp::export]]
Rcpp::NumericVector emis_atm(Rcpp::NumericVector Tk, Rcpp::NumericVector RH) {
  Rcpp::NumericVector emis_atm_out(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector e(Tk.size(), Rcpp::NumericVector::get_na());
  e = RH * esat(Tk);
  emis_atm_out = 0.575 * pow(e, 0.143);
  return emis_atm_out;
}

// [[Rcpp::export]]
//   // https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html?vA=290&units=K#

Rcpp::NumericVector viscosity(Rcpp::NumericVector Tk) { //viscosity of air, kg/(m s)
  Rcpp::NumericVector viscosity_out(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector omega(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector base_vec(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector temp2(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector test_o(Tk.size(), Rcpp::NumericVector::get_na());
  omega = (((Tk / 97) - 2.9) / 0.4 * (-0.034)) + 1.048;
  const double temp = pow(3.617, 2.0);
  temp2 = 28.97 * Tk;
      // Rcpp::Rcout << "omega = " << omega << std::endl << "-------------" << std::endl;
      // Rcpp::Rcout << "temp = " << temp << std::endl << "-------------" << std::endl;
  viscosity_out = 0.0000026693 * pow(temp2, 0.5) / (temp * omega); 

  return viscosity_out;
}

// [[Rcpp::export]]
Rcpp::NumericVector thermal_cond(Rcpp::NumericVector Tk) { //Thermal conductivity of air, W/(m K). Reference: BSL, page 257.
  Rcpp::NumericVector thermal_cond_out(Tk.size(), Rcpp::NumericVector::get_na());
  thermal_cond_out = (cp + 1.25 * r_air) * viscosity(Tk);
  return thermal_cond_out;
}

// [[Rcpp::export]]
Rcpp::NumericVector h_sphere_in_air(Rcpp::NumericVector Tk, Rcpp::NumericVector speed) { //Convective heat transfer coefficient for flow around a sphere, W/(m2 K). Reference: Bird, Stewart, and Lightfoot (BSL), page 409
  Rcpp::NumericVector therm_con(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Re(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Nu(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector h_sphere_in_air_val_out(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector density(Tk.size(), Rcpp::NumericVector::get_na());
  
  density = Pair * 100 / (r_air * Tk);
  speed[speed < min_speed] = min_speed;
  Re = speed * density * diam_globe / viscosity(Tk); // Reynolds number for sphere
  Nu = 2.0 + 0.6 * pow(Re, 0.5) * pow(Pr, 0.3333); // Nusselt number for sphere
  h_sphere_in_air_val_out = Nu * thermal_cond(Tk) / diam_globe;
  // Rcpp::Rcout << "density = " << density << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "speed = " << speed << std::endl << "-------------" << std::endl;
  //  Rcpp::Rcout << "Re = " << Re << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "Nu = " << Nu << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "h_sphere_in_air_val_out = " << h_sphere_in_air_val_out << std::endl << "-------------" << std::endl;
  // double Tksize = Tk.size();
  // Rcpp::Rcout << "Tksize_sphere = " << Tksize << "; h_sphere_in_air_val_out = " << h_sphere_in_air_val_out << std::endl;
  
  return h_sphere_in_air_val_out;
}

// [[Rcpp::export]]
Rcpp::NumericVector h_cylinder_in_air(Rcpp::NumericVector Tk, Rcpp::NumericVector speed) { //Convective heat transfer coefficient for a long cylinder, W/(m2 K). Bedingfield and Drew, eqn 32
  Rcpp::NumericVector therm_con(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Re(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Nu(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector h_cylinder_in_air_val_out(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector density(Tk.size(), Rcpp::NumericVector::get_na());
  density = Pair * 100 / (r_air * Tk);
  speed[speed < min_speed] = min_speed;
  Rcpp::NumericVector vis = viscosity(Tk);
  Re = speed * density * diam_wick / viscosity(Tk); // Reynolds number,  different than the one in h_sphere_in_air
  Nu = 0.281 * pow(Re, 0.6) * pow(Pr, 0.44);  // Nusselt number, different than the one in h_sphere_in_air
  h_cylinder_in_air_val_out = Nu * thermal_cond(Tk) / diam_wick;
//  Rcpp::Rcout << "vis = " << vis << "; Pr = " << Pr << std::endl;
  return h_cylinder_in_air_val_out;
}

// [[Rcpp::export]]
Rcpp::NumericVector diffusivity(Rcpp::NumericVector Tk) { //diffusivity of water vapor in air, m2/s. Reference: BSL, page 505.
  Rcpp::NumericVector diffusivity_out(Tk.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector holder(Tk.size(), Rcpp::NumericVector::get_na());
  holder = Tk / Tcrit12; 
  diffusivity_out = 0.000364 * pow(holder, 2.334) * pcrit13 * tcrit512 * Mmix / (Pair / 1013.25) * 0.0001;
  return diffusivity_out;
}

// [[Rcpp::export]]
Rcpp::List zenRadCorrect(Rcpp::NumericVector radiation, Rcpp::NumericVector zenith) {
  for (R_xlen_t i=0; i<radiation.size(); i++) {
    if (zenith[i] <= 0.0) {zenith[i] = 1e-10;
    } else if (radiation[i] > 0.0 & zenith[i] > 1.57) {zenith[i] = 1.57;
    } else if (radiation[i] > 15.0 & zenith[i] > 1.54) {zenith[i] = 1.54;
    } else if (radiation[i] > 900.0 & zenith[i] > 1.52) {zenith[i] = 1.52;
    } else if (radiation[i] < 10.0 & zenith[i] == 1.57) {radiation[i] = 0.0;
    }
  }   
  return Rcpp::List::create(Rcpp::Named("zenithMod") = zenith, 
                            Rcpp::Named("radiationMod") = radiation);
}

// // [[Rcpp::export]]
// Rcpp::List fTg_prep(Rcpp::NumericVector tas, Rcpp::NumericVector hurs, Rcpp::NumericVector radiation, Rcpp::NumericVector zenith) {
//   Rcpp::NumericVector cza(tas.size(), Rcpp::NumericVector::get_na());
//   Rcpp::NumericVector RH(tas.size(), Rcpp::NumericVector::get_na());
//   Rcpp::NumericVector Tair(tas.size(), Rcpp::NumericVector::get_na());
//   Rcpp::NumericVector zenithMod(tas.size(), Rcpp::NumericVector::get_na());
//   
//   }
// //  Rcpp::Rcout << "cza = " << cza << std::endl;
//   
//   return Rcpp::List::create(Rcpp::Named("zenithMod") = zenithMod, 
//                              Rcpp::Named("radiationMod") = radiation);
// }

// [[Rcpp::export]]
Rcpp::List fTnwb_prep(Rcpp::NumericVector tas, Rcpp::NumericVector dewp, Rcpp::NumericVector hurs, Rcpp::NumericVector radiation, Rcpp::NumericVector wind, Rcpp::NumericVector zenith) {
  // variables
  Rcpp::NumericVector Tair(tas.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Tdew(tas.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector eair(tas.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector emis_atm_out(tas.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector density(tas.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector RH(tas.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector viscosity_out(tas.size(), Rcpp::NumericVector::get_na());
  
  for (R_xlen_t i=0; i<tas.size(); i++) {
    if (zenith[i] <= 0) {zenith[i] = 1e-10;
    } else if (radiation[i] > 0 & zenith[i] > 1.57) {zenith[i] = 1.57;
    } else if (radiation[i] > 15 & zenith[i] > 1.54) {zenith[i] = 1.54;
    } else if (radiation[i] > 900 & zenith[i] > 1.52) {zenith[i] = 1.52;
    } else if (radiation[i] < 10 & zenith[i] == 1.57) {radiation[i] = 0;
    }
  }
  RH = hurs * 0.01;
  Tdew = dewp + kVal;
  Tair = tas + kVal;
  eair = RH * esat(Tair);
  emis_atm_out = emis_atm(Tair, RH);
  density = Pair * 100/(Tair * r_air);
  viscosity_out = viscosity(Tair);
  
  return Rcpp::List::create(Rcpp::Named("emis_atm_out") = emis_atm_out,
                            Rcpp::Named("zenithMod") = zenith, 
                            Rcpp::Named("radiationMod") = radiation,
                            Rcpp::Named("RH") = RH,
                            Rcpp::Named("density") = density,
                            Rcpp::Named("viscosity_out") = viscosity_out,
                            Rcpp::Named("eair") = eair);
}


// fr_tg is the function to be minimized in fTg.R. Returns globe temperature in degC.
// Tglobe_prev is the value of Tair over which the optimization occurs. The range is Tair-2, Tair+10
// changed tas to Tair, Aug 31, 2021

// [[Rcpp::export]]
Rcpp::NumericVector fr_tg(double Tglobe_prev, Rcpp::NumericVector Tair, Rcpp::NumericVector hurs, Rcpp::NumericVector speed, Rcpp::NumericVector radiation, Rcpp::NumericVector zenith) {
  Rcpp::NumericVector Tglobe(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector h_sphere(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector emis_atm_out(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector RH(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Tsfc(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Tref_globe(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector cza(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Tglobe_base(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector zenithMod(speed.size(), Rcpp::NumericVector::get_na());
  cza = cos(zenithMod);
  
  Tsfc = Tair; // Tsfc is surface temperature; Tair is air temp. Since we don't have separate values for these Liljegren, et al, set them equal. 
  Tref_globe = 0.5 * (Tglobe_prev + Tair);
  h_sphere = h_sphere_in_air(Tref_globe, speed); //Convective heat transfer coefficient for flow around a sphere, W/(m2 K)
  Tglobe = pow((0.5 * (emis_atm_out * pow(Tair, 4) + emis_sfc * pow(Tsfc, 4)) - h_sphere / (emis_globe * stefanb) * (Tglobe_prev - Tair) + radiation / 
    (2 * emis_globe * stefanb) * (1 - alb_globe) * (propDirect * (1 / (2 * cza) - 1) + 1 + alb_sfc)), 0.25);
  return abs(Tglobe - Tglobe_prev);
}

// fr_tnwb is the function to be minimized in fTnwb.R.
// Twb_prev is the value of Tair over which the optimization occurs. The range is Tdew-1 to Tair+1

// [[Rcpp::export]]
Rcpp::NumericVector fr_tnwb(double Twb_prev, Rcpp::NumericVector Tair, Rcpp::NumericVector wind) {
  Rcpp::NumericVector Tref_cylinder(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Fatm(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Sc(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector h_cyl(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector ewick(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector evap(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Twb(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector diffusivity(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector density(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector viscosity_out(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector emis_atm_out(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector zenithMod(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector radiationMod(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector eair(Tair.size(), Rcpp::NumericVector::get_na());
  
  Tref_cylinder = 0.5 * (Twb_prev + Tair);
  diffusivity = 0.000364 * pow((Tref_cylinder / Tcrit12), 2.334) * pcrit13 * tcrit512 * Mmix / (Pair / 1013.25) * 0.0001;
  Sc = viscosity_out/(density * diffusivity);
  h_cyl = h_cylinder_in_air(Twb_prev, wind);
  evap = (313.15 - Twb_prev) / 30.0 * (-71100.0) + 2407300.0; //h_evap moved here
  Fatm =  stefanb * emis_wick * (0.5 * (emis_atm_out * pow(Tair, 4.0) + emis_sfc * pow(Tair, 4.0)) - pow(Twb_prev, 4.0)) + (1.0 - alb_wick) * radiationMod * ((1.0 - propDirect) * 
    (1.0 + 0.25 * diam_wick/len_wick) + ((tan(zenithMod)/3.1416) + 0.25 * diam_wick/len_wick) * propDirect + SurfAlbedo);
  Twb = Tair - evap / ratio * (ewick - eair) / (Pair - ewick) * pow(Pr / Sc, 0.56) + Fatm / h_cyl * irad;
  return abs(Twb - Twb_prev);
}
