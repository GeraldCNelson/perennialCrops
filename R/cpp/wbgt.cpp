
// optimization directions from https://cran.r-project.org/web/packages/RcppNumerical/vignettes/introduction.html#numerical-optimization
#include <Rcpp.h>

const double alb_globe = 0.05;
const double alb_wick = 0.4;
const double diam_globe = 0.0508;
const double diam_wick = 0.007;
const double emis_globe = 0.95;
const double emis_sfc = 0.999  ;
const double emis_wick = 0.95;
const double len_wick = 0.0254;
const double m_h2o = 18.015;
const double stefanb = 0.000000056696;
const double SurfAlbedo = 0.4;
const double propDirect = 0.8; //Assume a proportion of direct radiation = direct/(diffuse + direct)
const double Pair = 1010; // Atmospheric pressure in hPa

const double cp = 1003.5;
const double m_air = 28.97;
const double r_gas = 8314.34;
const double r_air = r_gas / m_air;

const double Pr = cp / (cp + (1.25 * r_air));
const double pi = 3.14159265358979323846;
const double min_speed = 0.1;

const double kVal = 273.15;



inline double esat(const double& Tk) { 
  //saturation vapor pressure (hPa) over water.  
  // Reference: Buck's (1981) approximation (eqn 3) of Wexler's (1976) formula over liquid water.
  return 1.004 * 6.1121 * exp(17.502 * (Tk - kVal) / (Tk - 32.18));
}


inline double emis_atm(const double& Tk, const double& RH) {
  // replaces call to esat function e = RH * esat(Tk);
    double e = RH * 1.004 * 6.1121 * exp(17.502 * (Tk - kVal) / (Tk - 32.18)); 
    return 0.575 * pow(e, 0.143);
}

// https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html?vA=290&units=K#
inline double viscosity(const double&  Tk) { 
  //viscosity of air, kg/(m s)
  double visc = (((Tk / 97.0) - 2.9) / 0.4 * (-0.034)) + 1.048;
  return 0.0000026693 * pow((28.97 * Tk), 0.5) / (pow(3.617, 2.0) * visc);
}


inline double thermal_cond(const double& viscosity) {
  //Thermal conductivity of air, W/(m K). Reference: BSL, page 257.
  return viscosity * (cp + 1.25 * r_air);
}

inline double h_evap(const double& Tk) {
  return (313.15 - Tk) / 30.0 * (-71100.0) + 2407300.0;
}


inline double h_cylinder_in_air(const double& Tk, const double& speed, const double& viscosity) {
   double density = Pair * 100 / (r_air * Tk);
   double okspeed = speed < min_speed ?  min_speed : speed;
   // Reynolds number,  different than the one in h_sphere_in_air
   double Re = okspeed * density * diam_wick / viscosity; 
   // Nusselt number, different than the one in h_sphere_in_air     
   double Nu = 0.281 * pow(Re, 0.6) * pow(Pr, 0.44); 
   //replaces thermal_cond(Tk, viscosity)  
   double thermal_con = (cp + 1.25 * r_air) * viscosity;
   return Nu * thermal_con / diam_wick;
}

inline double h_sphere_in_air(const double& Tk, const double& speed, const double& viscosity) { 
  //Convective heat transfer coefficient for flow around a sphere, W/(m2 K). 
  // Reference: Bird, Stewart, and Lightfoot (BSL), page 409
   double tcond = thermal_cond(viscosity);
   double density = Pair * 100 / (r_air * Tk);
   double okspeed = speed < min_speed ?  min_speed : speed;
   // Reynolds number for sphere
   double Re = okspeed * density * diam_globe / viscosity;
   // Nusselt number for sphere
   double Nu = 2.0 + 0.6 * pow(Re, 0.5) * pow(Pr, 0.3333); 
   return tcond * Nu / diam_globe;     
}

const double toRad = pi/180;
const double toDeg = 180/pi;

inline void degToRad(double& angle) {
  angle *= toRad;
}

inline void radToDeg(double& angle) {
  angle *= toDeg;
}



double calZenith(const int doy, const int year_num, double lon, double lat) { 
//  const double EQTIME1 = 229.18;
//  const double EQTIME2 = 7.5e-05;
//  const double EQTIME3 = 0.001868;
//  const double EQTIME4 = 0.032077;
//  const double EQTIME5 = 0.014615;
//  const double EQTIME6 = 0.040849;
  const double DECL1 = 0.006918;
  const double DECL2 = 0.399912;
  const double DECL3 = 0.070257;
  const double DECL4 = 0.006758;
  const double DECL5 = 0.000907;
  const double DECL6 = 0.002697;
  const double DECL7 = 0.00148;
	const double utc_hour = 12.0;
	const double TimeOffset = 0.0;
	const double TrueSolarTime = (utc_hour * 60.0) + TimeOffset;

  //Calculate zenith angle in degrees
  int dpy = 365;
  if( year_num % 400 == 0 || year_num % 4 == 0) {
      dpy = 366;
  }  

  degToRad(lat);
  
  //Evaluate the fractional year in radians 
  double Gamma = 2 * pi * ((doy - 1.0) + (utc_hour/24.0))/dpy; 
  //Evaluate the Equation of time in minutes 
  // not used
  // double EquTime = EQTIME1 * (EQTIME2 + EQTIME3 * cos(Gamma) - EQTIME4 * sin(Gamma) - EQTIME5 * cos(2 * Gamma) - EQTIME6 * sin(2 * Gamma)); 
  //Evaluate the solar declination angle in radians (must be between -23.5 and 23.5 degrees)
  double Decli = DECL1 - DECL2 * cos(Gamma) + DECL3 * sin(Gamma) - DECL4 * cos(2 * Gamma) + DECL5 * sin(2 * Gamma) - 
    DECL6 * cos(3 * Gamma) + DECL7 * sin(3 * Gamma); 
  double Ha = ((TrueSolarTime/4.0) - 180.0);
  degToRad(Ha);
  double CosZen = (sin(lat) * sin(Decli) + cos(lat) * cos(Decli) * cos(Ha));
  CosZen = CosZen > 1 ? 1 : (CosZen < -1 ? 1 : CosZen);
  double SZA = acos(CosZen);
  radToDeg(SZA);
  return(SZA);
}


// [[Rcpp::export]]
Rcpp::List fwbgt_prep(const Rcpp::DateVector& dates, const double& lon, const double& lat, 
                      const Rcpp::NumericVector& tas, const Rcpp::NumericVector& hurs, 
                      const Rcpp::NumericVector& dewp, const Rcpp::NumericVector& speed, 
                      Rcpp::NumericVector radiation) {
  // variables
  Rcpp::NumericVector zenith(tas.size());
  Rcpp::NumericVector density(tas.size());
  Rcpp::NumericVector emis_atm_out(tas.size());
  Rcpp::NumericVector eair(tas.size());
  Rcpp::NumericVector viscosity_out(tas.size());
  Rcpp::NumericVector RH(tas.size());

  size_t n = tas.size();
  
  for (size_t i=0; i<n; i++) {
    Rcpp::Date d = dates[i];
    
    zenith[i] = calZenith(d.getYearday(), d.getYear(), lon, lat);
  
    if (zenith[i] <= 0) { zenith[i] = 1e-10;
    } else if ((radiation[i] > 0) & (zenith[i] > 1.57)) {zenith[i] = 1.57;  // 90°
    } else if ((radiation[i] > 15) & (zenith[i] > 1.54)) {zenith[i] = 1.54; // 88°
    } else if ((radiation[i] > 900) & (zenith[i] > 1.52)) {zenith[i] = 1.52;// 87°
    } else if ((radiation[i] < 10) & (zenith[i] == 1.57)) {radiation[i] = 0;
    }
  
    RH[i] = hurs[i] * 0.01;
    double KTair = tas[i] + kVal;
    eair[i] = RH[i] * esat(KTair); 
    // replaces call to emis_atm function
    emis_atm_out[i] = 0.575 * pow(eair[i], 0.143);
    density[i] = Pair * 100/(KTair * r_air);
    viscosity_out[i] = viscosity(KTair);
    Rcpp::Rcout << "emis_atm_out = " << emis_atm_out << ", eair[i] = " << eair[i] << ", RH[i] = " << RH[i] << ", emis_atm_out[i] = " << emis_atm_out[i] << std::endl;
  }
  
 //   ", emis_sfc = " << emis_sfc << ", Fatm = " << Fatm  << ", zenith[i]) = " << zenith[i] << ", ewick = " << ewick  << ", h_cyl = " << h_cyl << ", h_cyl = " << h_cyl  << std::endl;
  
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
Rcpp::NumericVector fr_tnwb(double Twb_prev, Rcpp::NumericVector Tair, 
                            Rcpp::NumericVector speed, Rcpp::NumericVector radiation, 
                            Rcpp::NumericVector zenith, Rcpp::NumericVector viscosity_out, 
                            Rcpp::NumericVector emis_atm_out, Rcpp::NumericVector eair, 
                            Rcpp::NumericVector density) {

	const double pcrit13 =  19.94585; //(36.4 * 218) ^ (1 / 3);
	const double tcrit512 = 113.4662; //(132 * 647.3) ^ (5 / 12);
	const double Tcrit12 = 292.3074; //(132 * 647.3) ^ 0.5;
	const double Mmix =  0.3000463; //(1 / 28.97 + 1 / 18.015) ^ 0.5
	const double ratio = cp * m_air/m_h2o;
	const double irad = 1.0;

  size_t n = Tair.size();
  Rcpp::NumericVector out(n);
  for (size_t i=0; i<n; i++) {
     double Tref_cylinder = 0.5 * (Twb_prev + Tair[i]);
     double diffusivity = 0.000364 * pow((Tref_cylinder / Tcrit12), 2.334) * pcrit13 * tcrit512 * Mmix / (Pair / 1013.25) * 0.0001; // replaces diffusivityR(Tref_cylinder, Pair)
     double Sc = viscosity_out[i] / (density[i] * diffusivity);

     double h_cyl = h_cylinder_in_air(Twb_prev, speed[i], viscosity_out[i]);

     double ewick = 1.004 * 6.1121 * exp(17.502 * (Twb_prev - kVal) / (Twb_prev)); //esat formula used here directly
     double evap = (313.15 - Twb_prev) / 30.0 * (-71100.0) + 2407300.0; //h_evap moved here
     double Fatm =  stefanb * emis_wick * (0.5 * (emis_atm_out[i] * pow(Tair[i], 4.0) + emis_sfc *
                    pow(Tair[i], 4.0)) - pow(Twb_prev, 4.0)) + (1.0 - alb_wick) * radiation[i] * 
                    ((1.0 - propDirect) * (1.0 + 0.25 * diam_wick/len_wick) +
                    ((tan(zenith[i])/3.1416) + 0.25 * diam_wick/len_wick) * propDirect + SurfAlbedo);
  
     double Twb = Tair[i] - evap / ratio * (ewick - eair[i]) / (Pair - ewick) * 
                  pow(Pr / Sc, 0.56) + Fatm / h_cyl * irad;
     
     //RH: this seems odd. Should Twb_prev be out[i-1]?
     out[i] = (Twb - Twb_prev);
   }   
  return out;
}

// fr_tg is the function to be minimized in fTg.R. Returns globe temperature in degC.
// Tglobe_prev is the value of Tair over which the optimization occurs. The range is Tair-2, Tair+10
// [[Rcpp::export]]
Rcpp::NumericVector fr_tg(double Tglobe_prev, Rcpp::NumericVector Tair, 
                          Rcpp::NumericVector hurs, Rcpp::NumericVector speed, 
                          Rcpp::NumericVector radiation, Rcpp::NumericVector zenith, 
                          Rcpp::NumericVector viscosity_out, Rcpp::NumericVector emis_atm_out) {

  size_t n = Tair.size();
  Rcpp::NumericVector out(n);
  for (size_t i=0; i<n; i++) {
     double cza = cos(zenith[i]);
    // Tsfc is surface temperature; Tair is air temp.
    // Since we don't have separate values for these Liljegren, et al, set them equal. 
    // double Tsfc = Tair[i];
     double Tref_globe = 0.5 * (Tglobe_prev + Tair[i]);
     //Convective heat transfer coefficient for flow around a sphere, W/(m2 K)
     double h_sphere = h_sphere_in_air(Tref_globe, speed[i], viscosity_out[i]);

      double Tglobe = pow((0.5 * (emis_atm_out[i] * pow(Tair[i], 4) + emis_sfc * 
                      pow(Tair[i], 4)) - h_sphere / (emis_globe * stefanb) * 
                      (Tglobe_prev - Tair[i]) + radiation[i] / (2 * emis_globe * stefanb) 
                      * (1 - alb_globe) * (propDirect * (1 / (2 * cza) - 1) + 1 + SurfAlbedo)), 0.25);
     out[i] = abs(Tglobe - Tglobe_prev);
  }
  return out;
}

// 
// 
// /*** R
// #remotes::install_github("anacv/HeatStress")
// library(HeatStress)
// library(Rcpp)
// #constants
// kVal = 273.15
// tolerance <- 1e-4
// 
// data("data_obs") 
// d <- data_obs
// d[apply(is.na(data_obs), 1, any), ] <- NA
// test <- wbgt.Liljegren(tas=data_obs$tasmean, dewp=data_obs$dewp, 
//                                 wind=data_obs$wind, radiation=data_obs$solar, dates= data_obs$Dates, lon=-5.66, lat=40.96)
// 
// 
// # generate hurs from dewpoint test values
// #hurs_test <- 100*(exp((17.625*dewp_test)/(243.04+dewp_test))/exp((17.625*tas_test)/(243.04+tas_test))) 
// 
// # optimization, fTnwb -----
// opt_fTnwb <- function(tas, hurs, dewp, speed, radiation, zenith, viscosity_out, 
//                       emis_atm_out, eair, density, tolerance) {
//   # browser() 
//   Tair = tas + kVal
//   r <- rep(NA, length(Tair))
//   Tdew <- dewp + kVal 
//   notna <- which(!is.na(tas))
//   for (i in notna) {
//       rng <- range(Tdew[i] - 1.0, Tair[i] + 1.0)
//       temp <- stats::optimize(f = fr_tnwb, interval = rng, tol = tolerance, 
//                      Tair[i],  speed[i], radiation[i], zenith[i], viscosity_out[i],
//                      emis_atm_out[i], eair[i], density[i])
//       r[i] <- temp$minimum - kVal
//   }
//   return(r)
// }
// 
// # optimization, fTg -----
// opt_fTg <- function(tas = tas, hurs = hurs, speed = speed, radiation = radiation, zenith = zenith, viscosity_out = viscosity_out, emis_atm_out, tolerance) {
//   Tair = tas + kVal
//   r <- rep(NA, length(Tair))
//   not_na <- which(!is.na(tas))
//   for (i in not_na) {
//     rng <- range(Tair[i] - 2, Tair[i] + 10)
//     temp <- stats::optimize(f = fr_tg, interval = rng, tol = tolerance,  
//                       Tair[i], hurs[i], speed[i], radiation[i], zenith[i], 
//                       viscosity_out[i], emis_atm_out[i])
//     r[i] <- temp$minimum - kVal
//   }
//   return(r)
// }
// 
// 
// # A new wbgt function -----
// f_wbgt.Liljegren_cpp <- function(tas, hurs, dewp, radiation, speed, lon, lat, dates, tolerance) { # raw radiation read in here
//   prep <- fwbgt_prep(dates, lon, lat, tas, hurs, dewp, radiation, speed)
//   Tnwb <- opt_fTnwb(tas = tas, hurs = hurs, dewp = dewp,  speed = speed, radiation = radiation, 
//                     zenith = prep$zenith, viscosity_out = prep$viscosity_out, 
//                     emis_atm_out = prep$emis_atm_out, eair = prep$eair, density = prep$density, 
//                     tolerance = tolerance)
//   Tg <- opt_fTg(tas = tas, hurs = hurs, speed = speed, radiation = radiation, zenith = prep$zenith, 
//               viscosity_out = prep$viscosity_out, emis_atm_out = prep$emis_atm_out, tolerance = tolerance)
//   cbind(data = 0.7 * Tnwb + 0.2 * Tg + 0.1 * tas, Tnwb = Tnwb, Tg = Tg)
// }
// 
// 
// i = 1:nrow(d)
// x <- f_wbgt.Liljegren_cpp(tas = d$tasmean[i], hurs = d$hurs[i], 
//                   dewp = d$dewp[i], speed = d$wind[i], radiation = d$solar[i], 
//                      lon = -5.66, lat = 40.96, dates = d$Dates[i], 
//                      tolerance = tolerance)
// tail(x)
// par(mfrow=c(2,2))
// for (i in 1:3) {
//   plot(test[[i]], x[,i], xlab="HeatStress", ylab="cpp")
//   abline(0,1)
// }
// */


/* 
#terra stuff --- not ready for that
library(terra)
 
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
lat_r <- setValues(r, rep(lat_test, each=n))
lon_r <- setValues(r, rep(lon_test, each=n))
hurs_test <- 100*(exp((17.625*dewp_test)/(243.04+dewp_test))/exp((17.625*tas_test)/(243.04+tas_test))) 

# try lapp with the 92 day rasters
#fwbgt <- function (vtas, vhurs, vdewp, vwind, vrad, vlon, vlat, dates) {
#   browser()
#  not_na <- unname(which(!is.na(vtas[1,]) & !is.na(vdewp[1,]) & !is.na(vhurs[1,]) & !is.na(vwind[1,])))
#  for (i in not_na) {
#     x <- f_wbgt.Liljegren_cpp(tas = vtas[,i], hurs = vhurs[,i], dewp = vdewp[,i], speed = vwind[,i], radiation = vrad[,i], vlon[,i], vlat[,i], dates, tolerance)
#     print(paste0("x: ", x))
#  }
#  return(x)
#}
#
#combined_r <- sds(tas_r, hurs_r, dewp_r,  wind_r, solar_r, lon_r, lat_r)
#out_92 <- lapp(combined_r, fun = fwbgt, dates_test)

*/
