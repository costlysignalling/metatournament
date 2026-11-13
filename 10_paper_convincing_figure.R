source('10_local_global_figures.R')
source('10_convincing_effort_limits.R')

# Reuse the plotted Supplement maps and their boundary geometry for A–C.
# D uses its own global baseline and recipient response at the requested parameters.
prepare_paper_convincing_figure <- function(png=TRUE,pdf=TRUE) {
  selections <- data.frame(family=c('local','global','global'),
                           set=c('symmetric_beta04','symmetric_beta04','symmetric_beta08'))
  data <- list()
  geometry <- list(fields=list(),equal=list(),alpha=list())
  limits <- readRDS('10_local_global_results/effort_limits.rds')
  benchmarks <- vector('list',4)
  for (i in seq_len(nrow(selections))) {
    family <- selections$family[i]
    set <- selections$set[i]
    local <- family=='local'
    data_file <- if(local) file.path('10_convincing_results',paste0(set,'.rds')) else
      file.path('10_comparison_results',paste0('popularizing_',set,'.rds'))
    geom_file <- if(local) file.path('10_convincing_results',paste0(set,'_geometry.rds')) else
      file.path('10_comparison_results',paste0('geometry_popularizing_',set,'.rds'))
    saved <- readRDS(geom_file)
    data[[i]] <- readRDS(data_file)[[2]]
    stopifnot(all(data[[i]]$p$phi==1))
    data[[i]]$map <- readRDS(file.path('10_local_global_results','high_resolution',
                                      paste0(family,'_',set,'_2.rds')))
    geometry$fields[[i]] <- saved$fields[[2]]
    geometry$equal[[i]] <- saved$equal[[2]]
    geometry$alpha[[i]] <- saved$alpha[[2]]
    benchmarks[[i]] <- limits[[paste(family,set,sep='_')]][[2]]
  }

  cache <- '10_local_global_results/paper_asymmetric_panel.rds'
  if(file.exists(cache)) {
    panel <- readRDS(cache)
  } else {
    p <- convincing_prepare(beta=c(.8,.4),phi=c(1.5,.66),N=1e60,
                            n=8192L,npressure=161L,nx=80L,adaptive=TRUE)
    message('Calculating the asymmetric example on the 361 by 361 talent grid')
    a <- seq(0,2.5,length.out=361)
    allocation <- convincing_optimize(p,as.matrix(expand.grid(a,a)))
    d <- list(p=p,map=list(alpha=a,allocation=allocation,recipient_abs=numeric()))
    panel <- list(data=d,fields=list(NULL,convincing_boundary_fields(d,2)),
                  equal=list(NULL,convincing_equal_curves(d,2)),
                  benchmarks=campaign_specialist_limits(p))
    saveRDS(panel,cache)
  }
  data[[4]] <- panel$data
  geometry$fields[[4]] <- panel$fields
  geometry$equal[[4]] <- panel$equal
  geometry$alpha[[4]] <- panel$data$map$alpha
  benchmarks[[4]] <- panel$benchmarks
  for(d in data) {
    o <- d$map$allocation
    stopifnot(all(is.finite(o)),max(abs(rowSums(o[,1:2])+abs(o[,3])-d$p$B))<1e-6)
  }
  draw_local_global_set(png=png,pdf=pdf,file='Figure_paper_convincing_examples',
    data=data,geometry=geometry,
    effort_geometry=list(readRDS('10_local_global_results/effort_geometry_local_symmetric_beta04_2.rds'),NULL,NULL,NULL),
    benchmarks=benchmarks,effort_limit=1,response_limit=1,rows=2:3,
    column_labels=c('N = 30','Global','Global','Global'))
}

if(sys.nframe()==0L)prepare_paper_convincing_figure()
