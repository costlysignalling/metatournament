# Run from the scripts folder to export both additional HDR strategy figures.
# Edit the A-D examples in 05_asymmetric_supplement_examples.R.
source('05_asymmetric_tournaments_figure_HDR.R')
source('05_asymmetric_supplement_examples.R')
cache_file <- '05_asymmetric_supplement_population.rds'
signature <- list(scenarios=asymmetric_supplement_examples,B=50,nquad=64L,nxpop=801L,
  engine=unname(tools::md5sum('05_energy_allocation_population_phispace.cpp')),
  solver=unname(tools::md5sum('05_explorer_population.R')))
saved <- if(file.exists(cache_file))readRDS(cache_file) else NULL
pool <- list()
if(!is.null(saved) && identical(saved$signature[names(signature)[-1]],signature[-1]))
  pool <- unlist(saved$population,recursive=FALSE,use.names=FALSE)
if(is.null(saved) || !identical(saved$signature,signature)) {
  saved <- list(signature=signature,population=list())
}
ids <- commandArgs(trailingOnly=TRUE)
if(!length(ids))ids <- names(asymmetric_supplement_examples)
stopifnot(all(ids %in% names(asymmetric_supplement_examples)))
for(id in ids) {
  scenarios <- asymmetric_supplement_examples[[id]]
  if(is.null(saved$population[[id]])) {
    # Display limits do not affect population means. Reuse matching parameters
    # when columns are moved or the same example is shown at different zooms.
    fields <- c('beta1','beta2','phi1','phi2')
    prep <- lapply(scenarios,function(sc) {
      matches <- which(vapply(pool,function(p)identical(p$pars[fields],sc[fields]),TRUE))
      if(length(matches)) p <- pool[[matches[1]]] else
        p <- prepare_asymmetric_examples(list(sc),50,64L,801L)[[1]]
      p$pars <- sc
      pool[[length(pool)+1L]] <<- p
      p
    })
    attr(prep,'signature') <- list(scenarios=scenarios,B=50,nquad=64L,nxpop=801L)
    saved$population[[id]] <- prep
    saveRDS(saved,cache_file)
  }
  result <- plot_asymmetric_hdr(scenarios=scenarios,B=50,
    population=saved$population[[id]],png=TRUE,pdf=TRUE,
    file=paste0('Figure_asymmetric_tournaments_',id,'.png'))
  message('Saved ',id)
}
