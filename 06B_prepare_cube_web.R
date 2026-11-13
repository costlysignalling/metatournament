# Run after rebuilding the two manuscript figures to refresh the web bundle.
# No model calculation is repeated here. The app folder is independently
# deployable after this preparation; it does not read parent-folder files.
dir.create('06B_cube_web/engine',recursive=TRUE,showWarnings=FALSE)
code <- readLines('06B_predict_many_tournaments.R',warn=FALSE)
code <- sub("Rcpp::sourceCpp('06B_many_tournaments_core.cpp',",
  "Rcpp::sourceCpp('engine/06B_many_tournaments_core.cpp',env=environment(many_compile),",
  code,fixed=TRUE)
writeLines(code,'06B_cube_web/engine/06B_predict_many_tournaments.R',useBytes=TRUE)
for(name in c('06B_many_tournaments_core.cpp','06B_three_tournament_geometry.R'))
  stopifnot(file.copy(name,file.path('06B_cube_web/engine',name),overwrite=TRUE))
source('06B_three_tournament_geometry.R')
z <- readRDS('06B_many_results/three_tournament_population.rds')
r <- tail(z$checks,1)[[1]]$responsive
g <- readRDS('06B_many_results/allocation_cube_meshes.rds')
stopifnot(r$converged)
saveRDS(list(beta=z$beta,phi=z$phi,n=r$n,n_cube=length(g$axis),
  residual=r$residual,max_gap=r$max_allocation_gap,E=r$E,
  data=many_mesh_payload(g$meshes)),'06B_cube_web/default_cube.rds')
cat('Refreshed the self-contained 06B_cube_web folder.\n')
