
Tair = tas_test + kVal;
Tsfc <- Tair
Tref = 0.5 * (Tair + Tair);
h = h_sphere_in_air(Tref, speed_test);
RH = hurs_test * 0.01;
emis_atm_out = emis_atm(Tair, RH) 
zenith_test <- calZenith(dates_test, lon_test, lat_test)

cza = cos(zenith_test)

Tglobe_prev  <- mean(Tair) #


Tglobe <- (0.5 * (emis_atm(Tair, RH) * Tair ^ 4 + emis.sfc * Tsfc ^ 4) - h / (emis.globe * stefanb) * (Tglobe_prev - Tair) + radiation / (2 * emis.globe * stefanb) * (1 - alb.globe) * (propDirect * (1 / (2 * cza) - 1) + 1 + alb.sfc)) ^ 0.25


i = 10

Tglobe <- (0.5 * (emis_atm(Tair[i], RH[i]) * Tair[i] ^ 4 + emis.sfc * Tsfc[i] ^ 4) - h / (emis.globe * stefanb) * (Tglobe_prev - Tair[i]) + radiation_test[i] / (2 * emis.globe * stefanb) * (1 - alb.globe) * (propDirect * (1 / (2 * cza) - 1) + 1 + alb.sfc)) ^ 0.25



Tglobe <- (0.5 * (emis_atm(Tair, RH) * Tair ^ 4 + emis.sfc * Tsfc ^ 4) - h / (emis.globe * stefanb) * (Tglobe_prev - Tair) + radiation / (2 * emis.globe * stefanb) * (1 - alb.globe) * (propDirect * (1 / (2 * cza) - 1) + 1 + alb.sfc)) ^ 0.25

emis_atm_out = emis_atm(Tair, RH) ;
Rcpp::NumericVector TglobeDelta = Tglobe_prev - Tair;
Rcpp::NumericVector propString = propDirect * (1.0 / (2.0 * cza) - 1.0) + 1.0 + alb.sfc; 
Rcpp::NumericVector endValues = h / (emisXstefanb) * (TglobeDelta) + radiation / (2 * emisXstefanb) * (1.0 - alb.globe) * (propString);
Tglobe_base = 0.5 * ((emis_atm_out + emis.sfc) * pow(Tair, 4)) - endValues;