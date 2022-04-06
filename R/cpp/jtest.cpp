
// optimization directions from https://cran.r-project.org/web/packages/RcppNumerical/vignettes/introduction.html#numerical-optimization
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
const double utc_hour = 12.0;
const double TimeOffset = 0.0;

// [[Rcpp::export]]
Rcpp::NumericVector esat(Rcpp::NumericVector Tk) { //saturation vapor pressure (hPa) over water.  Reference: Buck's (1981) approximation (eqn 3) of Wexler's (1976) formula over liquid water.
  return 1.004 * 6.1121 * exp(17.502 * (Tk - kVal) / (Tk - 32.18));
}

// [[Rcpp::export]]
Rcpp::NumericVector emis_atm(Rcpp::NumericVector Tk, Rcpp::NumericVector RH) {
  Rcpp::NumericVector e(Tk.size(), Rcpp::NumericVector::get_na());
  e = RH * 1.004 * 6.1121 * exp(17.502 * (Tk - kVal) / (Tk - 32.18)); // replaces call to esat function
  // e = RH * esat(Tk);
  return 0.575 * pow(e, 0.143);
}

// [[Rcpp::export]]
// https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html?vA=290&units=K#
Rcpp::NumericVector viscosity(Rcpp::NumericVector Tk) { //viscosity of air, kg/(m s)
  Rcpp::NumericVector omega(Tk.size(), Rcpp::NumericVector::get_na());
  omega = (((Tk / 97.0) - 2.9) / 0.4 * (-0.034)) + 1.048;
  return 0.0000026693 * pow((28.97 * Tk), 0.5) / (pow(3.617, 2.0) * omega);
}

// [[Rcpp::export]]
Rcpp::NumericVector thermal_cond(Rcpp::NumericVector Tk, Rcpp::NumericVector viscosity) { //Thermal conductivity of air, W/(m K). Reference: BSL, page 257.
  Rcpp::NumericVector thermal_cond_out(Tk.size(), Rcpp::NumericVector::get_na());
  thermal_cond_out = (cp + 1.25 * r_air) * viscosity;
  return thermal_cond_out;
}

// [[Rcpp::export]]
Rcpp::NumericVector h_evap(Rcpp::NumericVector Tk) { //heat of evaporation, J/(kg K), for temperature in the range 283-313 K.
  Rcpp::NumericVector h_evap_out(Tk.size(), Rcpp::NumericVector::get_na());
  h_evap_out = (313.15 - Tk) / 30.0 * (-71100.0) + 2407300.0;
  return h_evap_out;
}

// [[Rcpp::export]]
Rcpp::NumericVector h_cylinder_in_air(double Tk, Rcpp::NumericVector speed, Rcpp::NumericVector viscosity) {
  Rcpp::NumericVector therm_con(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Re(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Nu(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector h_cylinder_in_air_val_out(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector density(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector thermal_con(speed.size(), Rcpp::NumericVector::get_na());
  density = Pair * 100 / (r_air * Tk);
  speed[speed < min_speed] = min_speed;
  Re = speed * density * diam_wick / viscosity; // Reynolds number,  different than the one in h_sphere_in_air
  Nu = 0.281 * pow(Re, 0.6) * pow(Pr, 0.44);  // Nusselt number, different than the one in h_sphere_in_air
  thermal_con = (cp + 1.25 * r_air) * viscosity; //replaces thermal_cond(Tk, viscosity)
  //  Rcpp::Rcout << "Nu = " << Nu << std::endl << "-------------" << std::endl;
  
  return Nu * thermal_con / diam_wick;
}

// [[Rcpp::export]]
Rcpp::NumericVector h_sphere_in_air(Rcpp::NumericVector Tk, Rcpp::NumericVector speed, Rcpp::NumericVector viscosity) { //Convective heat transfer coefficient for flow around a sphere, W/(m2 K). Reference: Bird, Stewart, and Lightfoot (BSL), page 409
  Rcpp::NumericVector therm_con(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Re(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Nu(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector h_sphere_in_air_val_out(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector density(speed.size(), Rcpp::NumericVector::get_na());
  
  density = Pair * 100 / (r_air * Tk);
  speed[speed < min_speed] = min_speed;
  Re = speed * density * diam_globe / viscosity; // Reynolds number for sphere
  Nu = 2.0 + 0.6 * pow(Re, 0.5) * pow(Pr, 0.3333); // Nusselt number for sphere
  // Rcpp::Rcout << "Nu = " << Nu << std::endl << "-------------" << std::endl;
  return Nu * thermal_cond(Tk, viscosity) / diam_globe;
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

// [[Rcpp::export]]
Rcpp::IntegerVector get_days(const Rcpp::DateVector& dates) {
  Rcpp::IntegerVector out(dates.size(), Rcpp::NumericVector::get_na());
  for (R_xlen_t i=0; i< dates.size(); i++) {
    Rcpp::Date d = dates[i];
    out[i] = d.getYearday();
    // Rcpp::Rcout << "dates[i] = " << dates[i] << std::endl << "-------------" << std::endl;
    // Rcpp::Rcout << "d = " << d << std::endl << "-------------" << std::endl;
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
  //   Rcpp::Rcout << "Gamma = " << Gamma << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "doy = " << doy << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "dpy = " << doy << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "utc_hour = " << utc_hour << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "pi = " << pi << std::endl << "-------------" << std::endl;
  
  EquTime = EQTIME1 * (EQTIME2 + EQTIME3 * cos(Gamma) - EQTIME4 * sin(Gamma) - EQTIME5 * cos(2 * Gamma) - 
    EQTIME6 * sin(2 * Gamma)); //Evaluate the Equation of time in minutes 
  Decli = DECL1 - DECL2 * cos(Gamma) + DECL3 * sin(Gamma) - DECL4 * cos(2 * Gamma) + DECL5 * sin(2 * Gamma) - 
    DECL6 * cos(3 * Gamma) + DECL7 * sin(3 * Gamma); //Evaluate the solar declination angle in radians (must be between -23.5 and 23.5 degrees)
  // Note: TrueSolarTime defined in variable def above
  HaDeg = ((TrueSolarTime/4.0) - 180.0);
  HaRad = degToRad(HaDeg);
  //  Rcpp::Rcout << "HaDeg = " << HaDeg;
  CosZen = (sin(RadLat) * sin(Decli) + cos(RadLat) * cos(Decli) * cos(HaRad));
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
  // Rcpp::Rcout << "Gamma = " << Gamma << std::endl;
  // Rcpp::Rcout << "EquTime = " << EquTime << std::endl;
  // Rcpp::Rcout << "doy = " << doy << std::endl;
  // Rcpp::Rcout << "dpy = " << dpy << std::endl;
  // Rcpp::Rcout << "pi = " << pi << std::endl;
  // Rcpp::Rcout << "CosZen = " << CosZen << std::endl;
  // Rcpp::Rcout << "RadLat = " << RadLat << std::endl;
  // Rcpp::Rcout << "Decli = " << Decli << std::endl << "-------------" << std::endl;
  
  SZARad = acos(CosZen);
  SZA = radToDeg(SZARad);
  // Rcpp::Rcout << "SZARad = " << SZARad << std::endl << "-------------" << std::endl;
  return(SZA);
}

// [[Rcpp::export]]
Rcpp::List fwbgt_prep(Rcpp::NumericVector tas, Rcpp::NumericVector dewp, Rcpp::NumericVector hurs, Rcpp::NumericVector radiation, 
                      Rcpp::NumericVector speed, Rcpp::NumericVector zenith) {
  // variables
  Rcpp::NumericVector Tair(tas.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Tdew(tas.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector eair(tas.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector e(tas.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector emis_atm_out(tas.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector density(tas.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector RH(tas.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector viscosity_out(tas.size(), Rcpp::NumericVector::get_na());
  
  for (R_xlen_t i=0; i<tas.size(); i++) {
    if (zenith[i] <= 0) {zenith[i] = 1e-10;
    } else if ((radiation[i] > 0) & (zenith[i] > 1.57)) {zenith[i] = 1.57;  // 90°
    } else if ((radiation[i] > 15) & (zenith[i] > 1.54)) {zenith[i] = 1.54; // 88°
    } else if ((radiation[i] > 900) & (zenith[i] > 1.52)) {zenith[i] = 1.52;// 87°
    } else if ((radiation[i] < 10) & (zenith[i] == 1.57)) {radiation[i] = 0;
    }
  }
  // Rcpp::Rcout << "fwbgt_prep zenith adj = " << zenith << std::endl << "-------------" << std::endl;
  
  RH = hurs * 0.01;
  Tdew = dewp + kVal;
  Tair = tas + kVal;
  //  eair = RH * esat(Tair);
  eair = RH * 1.004 * 6.1121 * exp(17.502 * (Tair - kVal) / (Tair - 32.18)); // replaces call to esat function
  
  // emis_atm_out = emis_atm(Tair, RH);
  
  e = RH * 1.004 * 6.1121 * exp(17.502 * (Tair - kVal) / (Tair - 32.18)); // replaces call to esat function
  // e = RH * esat(Tk);
  emis_atm_out = 0.575 * pow(e, 0.143); // replaces call to emis_atm function
  
  density = Pair * 100/(Tair * r_air);
  viscosity_out = viscosity(Tair);
  
  // Rcpp::Rcout << "emis_atm_out = " << emis_atm_out << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "zenith = " << zenith << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "radiation = " << radiation << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "RH = " << RH << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "density = " << density << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "viscosity_out = " << viscosity_out << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "eair = " << eair << std::endl << "-------------" << std::endl;
  
  return Rcpp::List::create(Rcpp::Named("emis_atm_out") = emis_atm_out,
                            Rcpp::Named("zenith") = zenith, 
                            Rcpp::Named("radiation") = radiation,
                            Rcpp::Named("RH") = RH,
                            Rcpp::Named("density") = density,
                            Rcpp::Named("viscosity_out") = viscosity_out,
                            Rcpp::Named("eair") = eair);
}

// fr_tnwb is the function to be minimized in fTnwb.R.
// Twb_prev is the value of Tair over which the optimization occurs. The range is Tdew-1 to Tair+1
// [[Rcpp::export]]
Rcpp::NumericVector fr_tnwb(double Twb_prev, Rcpp::NumericVector Tair, Rcpp::NumericVector speed, Rcpp::NumericVector radiation, Rcpp::NumericVector zenith, 
                            Rcpp::NumericVector viscosity_out, Rcpp::NumericVector emis_atm_out, Rcpp::NumericVector eair, Rcpp::NumericVector density) {
  Rcpp::NumericVector Tref_cylinder(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Fatm(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Sc(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector h_cyl(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector ewick(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector evap(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Twb(Tair.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector diffusivity(Tair.size(), Rcpp::NumericVector::get_na());
  
  Tref_cylinder = 0.5 * (Twb_prev + Tair);
  diffusivity = 0.000364 * pow((Tref_cylinder / Tcrit12), 2.334) * pcrit13 * tcrit512 * Mmix / (Pair / 1013.25) * 0.0001; // replaces diffusivityR(Tref_cylinder, Pair)
  Sc = viscosity_out/(density * diffusivity);
  h_cyl = h_cylinder_in_air(Twb_prev, speed, viscosity_out);
  // Rcpp::Rcout << "Tair.size() = " << Tair.size() << ", speed = " << speed << ", viscosity_out = " << viscosity_out << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "Sc = " << Sc << "Twb_prev = " << Twb_prev  << std::endl << "-------------" << std::endl;
  
  ewick = 1.004 * 6.1121 * exp(17.502 * (Twb_prev - kVal) / (Twb_prev)); //esat formula used here directly
  evap = (313.15 - Twb_prev) / 30.0 * (-71100.0) + 2407300.0; //h_evap moved here
  Fatm =  stefanb * emis_wick * (0.5 * (emis_atm_out * pow(Tair, 4.0) + emis_sfc * pow(Tair, 4.0)) - pow(Twb_prev, 4.0)) + (1.0 - alb_wick) * radiation * ((1.0 - propDirect) * 
    (1.0 + 0.25 * diam_wick/len_wick) + ((tan(zenith)/3.1416) + 0.25 * diam_wick/len_wick) * propDirect + SurfAlbedo);
  Twb = Tair - evap / ratio * (ewick - eair) / (Pair - ewick) * pow(Pr / Sc, 0.56) + Fatm / h_cyl * irad;
  //  Rcpp::Rcout << "fr_tnwb consts " << " " << emis_wick <<  " " << stefanb <<  " alb_wick " << alb_wick <<  " " << propDirect <<   " zenith " << zenith << std::endl;
  
  // Rcpp::Rcout << "Tnwb = " << Twb << " Tnwb_prev = " << Twb_prev << " Fatm = " << Fatm << std::endl;
  // Rcpp::Rcout << "Twb_prev = " << Twb_prev << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "Fatm = " << Fatm << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "Sc = " << Sc << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "Fatm = " << Fatm << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "h_cyl = " << h_cyl << std::endl << "-------------" << std::endl;
  // Rcpp::Rcout << "irad = " << irad << std::endl << "-------------" << std::endl;
  
  return abs(Twb - Twb_prev);
}

// fr_tg is the function to be minimized in fTg.R. Returns globe temperature in degC.
// Tglobe_prev is the value of Tair over which the optimization occurs. The range is Tair-2, Tair+10
// [[Rcpp::export]]
Rcpp::NumericVector fr_tg(double Tglobe_prev, Rcpp::NumericVector Tair, Rcpp::NumericVector hurs, Rcpp::NumericVector speed, 
                          Rcpp::NumericVector radiation, Rcpp::NumericVector zenith, Rcpp::NumericVector viscosity_out, Rcpp::NumericVector emis_atm_out) {
  Rcpp::NumericVector Tglobe(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector h_sphere(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector RH(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Tsfc(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Tref_globe(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector cza(speed.size(), Rcpp::NumericVector::get_na());
  Rcpp::NumericVector Tglobe_base(speed.size(), Rcpp::NumericVector::get_na());
  cza = cos(zenith);
  //  Rcpp::Rcout << "fr_tg zenith = " << zenith << std::endl << "-------------" << std::endl;
  
  Tsfc = Tair; // Tsfc is surface temperature; Tair is air temp. Since we don't have separate values for these Liljegren, et al, set them equal. 
  Tref_globe = 0.5 * (Tglobe_prev + Tair);
  //Rcpp::Rcout << "Tref_globe = " << Tref_globe << "T globe_prev = " << Tglobe_prev << " Tair = " << Tair << std::endl << "-------------" << std::endl;
  
  h_sphere = h_sphere_in_air(Tref_globe, speed, viscosity_out); //Convective heat transfer coefficient for flow around a sphere, W/(m2 K)
  Tglobe = pow((0.5 * (emis_atm_out * pow(Tair, 4) + emis_sfc * pow(Tsfc, 4)) - h_sphere / (emis_globe * stefanb) * (Tglobe_prev - Tair) + radiation / 
    (2 * emis_globe * stefanb) * (1 - alb_globe) * (propDirect * (1 / (2 * cza) - 1) + 1 + alb_sfc)), 0.25);
  //  Rcpp::Rcout << "Tglobe = " << Tglobe << " Tglobe_prev = " << Tglobe_prev << " emis_atm_out = " << emis_atm_out << " radiation = " << radiation << std::endl;
  // Rcpp::Rcout << "emis_atm_out = " << emis_atm_out << std::endl << "-------------" << std::endl;
  
  // Rcpp::Rcout << "fr_tg consts " << " " << emis_globe <<  " " << stefanb<<  " " << alb_globe <<  " " << propDirect <<  " cza " << cza <<   " zenith " << zenith << std::endl;
  return abs(Tglobe - Tglobe_prev);
}


/*** R

library(terra)
library(HeatStress)
#terraOptions(memfrac = 2,  ncopies = 1, progress = 10, tempdir =  "data/ISIMIP", verbose = FALSE)
#constants
kVal = 273.15
tolerance <- 1e-4

data("data_obs") 
head(data_obs)
dates_test <- data_obs$Dates
tas_test <- data_obs$tasmean
dewp_test <- data_obs$dewp
speed_test <- data_obs$wind
radiation_test <- data_obs$solar

# generate hurs from dewpoint test values
hurs_test <- 100*(exp((17.625*dewp_test)/(243.04+dewp_test))/exp((17.625*tas_test)/(243.04+tas_test))) 
dewp_test_compare <- 243.04*(log(hurs_test/100)+((17.625*tas_test)/(243.04+tas_test)))/(17.625-log(hurs_test/100)-((17.625*tas_test)/(243.04+tas_test))) # should be the same as dewp_test
#compareOutputs(dewp_test_compare, dewp_test)

r <- rast(nrows=5, ncols=5, xmin = -180, xmax = 180, ymin = -90, ymax = 90, nlyrs = 92) 
n <- ncell(r)
tas_r <- setValues(r, rep(data_obs$tasmean, each=n))
dewp_r <- setValues(r, rep(data_obs$dewp, each=n))
hurs_r <- setValues(r, rep(data_obs$hurs, each=n))
wind_r <- setValues(r, rep(data_obs$wind, each=n))
solar_r <- setValues(r, rep(data_obs$solar, each=n))

lat_r <- init(r, "y")
lon_r <- init(r, "x")
lat_test <- values(lat_r)[,1]
lon_test <- values(lon_r)[,1]
lat_test <- c(lat_test, lat_test, lat_test, lat_test[1:17])
lon_test <- c(lon_test, lon_test, lon_test, lon_test[1:17])
hurs_test <- 100*(exp((17.625*dewp_test)/(243.04+dewp_test))/exp((17.625*tas_test)/(243.04+tas_test))) 

# optimization, fTnwb -----
opt_fTnwb <- function(tas, dewp, hurs, speed, radiation, zenith, viscosity_out, emis_atm_out, eair, density, lon, lat, dates, tolerance) {
  #      browser() 
  Tair = tas + kVal
  r <- rep(NA, length(Tair))
  Tdew <- dewp + kVal
  #print(str(tas))
  not_na <- which(!(is.na(tas) | is.na(hurs) | is.na(speed) | is.na(radiation)))
  for (i in not_na) {
    temp <- stats::optimize(f = fr_tnwb, interval = range(Tdew[i] - 1.0, Tair[i] + 1.0), tol = tolerance, Tair[i],  speed[i], radiation[i], zenith[i], viscosity_out[i], emis_atm_out[i], eair[i], density[i])
    r[i] <- temp$minimum - kVal
    outt <- round(temp$minimum )
    #    print(paste0("range fTnwb start: ", round((Tdew[i] - 1.0), 2),  ", range fTnwb end: ", round((Tair[i] + 1.0), 2), ", minimum fTnwb: ", outt))
  }
  return(r)
}

# optimization, fTg -----
opt_fTg <- function(tas = tas, hurs = hurs, speed = speed, radiation = radiation, zenith = zenith, viscosity_out = viscosity_out, emis_atm_out, lon, lat, dates, tolerance) {
  Tair = tas + kVal
  r <- rep(NA, length(Tair))
  not_na <- which(!(is.na(tas) | is.na(hurs) | is.na(speed) | is.na(radiation)))
  for (i in not_na) {
    temp <- stats::optimize(f = fr_tg, interval = range(Tair[i] - 2, Tair[i] + 10), tol = tolerance,  Tair[i], hurs[i], speed[i], radiation[i], zenith[i], viscosity_out[i], emis_atm_out[i])
    r[i] <- temp$minimum - kVal
    outt <- temp$minimum 
  }
  return(r)
}

# A new wbgt function -----
f_wbgt.Liljegren_cpp <- function(tas, dewp, hurs, radiation, speed, lon, lat, dates, tolerance) { # raw radiation read in here
  zenithDeg <- calZenith(dates, lon, lat)
  ZenithAngle <- degToRad(zenithDeg)
  zenith = ZenithAngle # needed for debugging
  
  test <- fwbgt_prep(tas, dewp, hurs, radiation, speed, zenith)
  emis_atm_out <- unlist(test[1])
  zenith_adj <- unlist(test[2])
  radiation <- unlist(test[3]) # adjusted radiation here
  RH <- unlist(test[4])
  density <- unlist(test[5])
  viscosity_out <- unlist(test[6])
  eair <- unlist(test[7])
  Tnwb <- opt_fTnwb(tas = tas,  dewp = dewp,  hurs = hurs, speed = speed, radiation = radiation, zenith = zenith_adj, viscosity_out = viscosity_out, emis_atm_out = emis_atm_out, eair = eair, density = density, lon = lon, lat = lat, 
                    dates = dates, tolerance = tolerance)
  Tg <- opt_fTg(tas = tas, hurs = hurs, speed = speed, radiation = radiation, zenith = zenith_adj, viscosity_out = viscosity_out, 
                emis_atm_out = emis_atm_out, lon = lon, lat = lat, dates = dates, tolerance = tolerance)
  wbgt <- list(data = 0.7 * Tnwb + 0.2 * Tg + 0.1 * tas, Tnwb = Tnwb, Tg = Tg)
}

# try with one day of data
i = 92
out_cpp <- f_wbgt.Liljegren_cpp(tas = tas_test[i], hurs = hurs_test[i], dewp = dewp_test[i], speed = speed_test[i], radiation = radiation_test[i], lon = lon_test[i], lat = lat_test[i], dates = dates_test[i], tolerance = tolerance)
str(out_cpp)

# use wbgt.Liljegren from heat stress
wbgt_outdoors <- wbgt.Liljegren(tas = tas_test[i], dewp = dewp_test[i], wind = speed_test[i],radiation = radiation_test[i], lon = lon_test[i], lat = lat_test[i], dates = dates_test[i], tolerance)
str(wbgt_outdoors)

# try with larger data set, 92 obs
out_cpp_larger <- f_wbgt.Liljegren_cpp(tas = tas_test, hurs = hurs_test, dewp = dewp_test, speed = speed_test,radiation = radiation_test, lon = lon_test, lat = lat_test, dates = dates_test, tolerance = tolerance)
str(out_cpp_larger)

# try lapp with the 92 day rasters
ff <- function(...) {
  f_wbgt.Liljegren_cpp(...)$Tg
}
combined_r <- sds(tas_r, dewp_r, hurs_r,  wind_r, solar_r, lon_r, lat_r)
out_yr <- lapp(combined_r, fun = ff, dates=dates_test, tolerance=tolerance)




*/
