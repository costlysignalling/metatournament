# Separate numerical engine: original example calculations remain unchanged.
# Global returns are evaluated as E/M*w directly, even in rare pools.
source('10_convincing_parameter_levels.R')
campaign_efficiency_compile<-local({e<-NULL;function(){if(is.null(e)){
 e<<-new.env();Rcpp::sourceCpp('10_campaign_efficiency_core.cpp',env=e,
 cacheDir=getOption('metatournament.efficiency.cpp.cache',file.path(tempdir(),'10_campaign_efficiency_cpp')))};e}})
campaign_efficiency_cell<-campaign_level_cell
.campaign_efficiency_environment<-new.env(parent=environment(campaign_level_cell))
.campaign_efficiency_environment$convincing_compile<-campaign_efficiency_compile
environment(campaign_efficiency_cell)<-.campaign_efficiency_environment
