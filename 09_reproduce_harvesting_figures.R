# Run from the project scripts folder. Saved population data are reused by default.
# Requires Rcpp, sf and isoband; compiling requires Rtools on Windows.
# install.packages(c('Rcpp','sf','isoband'))
source('09_strategy_examples.R')
source('09_prepare_harvesting_phispace.R')
source('09_harvesting_phispace_figures.R')
source('09_fixed_performance_harvesting_figure.R')

reproduce_harvesting_figures<-function(recalculate=FALSE){
 if(recalculate||!file.exists('09_harvesting_results/finite_examples.rds'))prepare_finite_pool_examples()
 if(recalculate||!file.exists('09_harvesting_results/examples.rds'))prepare_harvest_examples()
 prepare_harvesting_phispace(overwrite=recalculate)
 refresh_harvesting_presence()
 if(recalculate&&file.exists('09_harvesting_results/plot_fields.rds'))unlink('09_harvesting_results/plot_fields.rds')
 finite<-readRDS('09_harvesting_results/finite_examples.rds')
 draw_strategy_examples(finite,'Figure_finite_pool_without_boredom',finite=TRUE)
 draw_all_harvesting_examples()
 draw_all_harvesting_phispace()
 # The explanatory figure's function is defined in its own script.
 plot_fixed_performance_harvesting()
 invisible(NULL)
}
if(sys.nframe()==0L)reproduce_harvesting_figures()
