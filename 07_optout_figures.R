# Run this file from the scripts folder to reproduce the Supplement examples.
# Edit eta_values, or run an individual script after setting metatournament.eta.
# Numerical settings: optout.npop is population integration, optout.nalpha is
# the plotted talent raster, optout.nscan is the boundary candidate scan.
# Requires Rcpp and Rtools on Windows; model calculations are in R and C++.
eta_values<-c(.5,1)
scripts<-c('07_known_talents_strategy_map_optout_fixed_E_HDR.R',
  '07_known_talent_optout_utility_contours_fixed_E.R',
  '07_known_talents_strategy_map_optout_HDR_responsive_E.R',
  '07_known_talent_optout_utility_contours_responsive_E.R',
  '07_asymmetric_tournaments_optout_HDR.R',
  '07_asymmetric_tournaments_optout_utility_contours.R')
for(eta in eta_values) {
  options(metatournament.eta=eta)
  for(file in scripts) {
    message('eta=',eta,': ',file)
    source(file,local=new.env(parent=globalenv()))
  }
}
