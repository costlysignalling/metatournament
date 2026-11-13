# Shared examples for the strategy and utility figures. Edit this list once,
# then source either figure script. Columns follow this order in both figures.
asymmetric_examples <- list(
  list(beta1=0.6,beta2=0.6,phi1=1,phi2=0.66),
  list(beta1=0.8,beta2=0.6,phi1=1.5,phi2=1),
  list(beta1=0.8,beta2=0.6,phi1=1.5,phi2=1.5),
  list(beta1=0.6,beta2=0.4,phi1=0.66,phi2=1)
)

validate_asymmetric_examples <- function(scenarios,B,alim,nalpha,nx,nroot,nquad,nxpop,png,pdf) {
  if(!is.list(scenarios) || !length(scenarios)) stop('Supply a nonempty list of examples.')
  for(sc in scenarios) {
    keys <- c('beta1','beta2','phi1','phi2')
    if(!is.list(sc) || !all(keys %in% names(sc)) ||
       any(!vapply(sc[keys],function(x) is.numeric(x)&&length(x)==1L&&is.finite(x),TRUE)))
      stop('Each example needs finite scalar beta1, beta2, phi1 and phi2.')
    if(min(sc$beta1,sc$beta2)<=0 || max(sc$beta1,sc$beta2)>=1 || min(sc$phi1,sc$phi2)<=0)
      stop('Require 0 < beta1,beta2 < 1 and phi1,phi2 > 0.')
    if(!is.null(sc$alim) && (length(sc$alim)!=2 || any(!is.finite(sc$alim)) ||
       sc$alim[1]<0 || diff(sc$alim)<=0))
      stop('Each optional scenario alim must be an increasing nonnegative interval.')
  }
  if(length(B)!=1 || !is.finite(B) || B<=0) stop('B must be positive.')
  if(length(alim)!=2 || any(!is.finite(alim)) || alim[1]<0 || diff(alim)<=0)
    stop('alim must be an increasing nonnegative interval.')
  resolutions <- c(nalpha,nx,nroot,nquad,nxpop)
  if(length(resolutions)!=5 || any(!is.finite(resolutions)) || any(resolutions<3) ||
     any(resolutions!=as.integer(resolutions))) stop('Resolutions must be integers >= 3.')
  for(flag in list(png,pdf)) if(!is.logical(flag)||length(flag)!=1L||is.na(flag))
    stop('png and pdf must be TRUE or FALSE.')
}

prepare_asymmetric_examples <- function(scenarios,B,nquad,nxpop,population=NULL) {
  signature <- list(scenarios=scenarios,B=B,nquad=max(64L,nquad),nxpop=nxpop)
  if(!is.null(population)) {
    if(!identical(attr(population,'signature'),signature))
      stop('Reused population must have matching scenarios, B, nquad and nxpop.')
    return(population)
  }
  prep <- lapply(scenarios,function(sc) {
    fixed <- random_allocation_reference(sc$beta1,sc$beta2,sc$phi1,sc$phi2,B)
    pop <- explorer_population(sc$beta1,sc$beta2,sc$phi1,sc$phi2,B,fixed$s,
      nquad=max(64L,nquad),nx=nxpop)
    list(pars=sc,fixed=fixed,response=pop$first,eq=unlist(pop$responsive))
  })
  attr(prep,'signature') <- signature
  prep
}

asymmetric_figure_output <- function(draw,ncolumns,file,png,pdf) {
  # Same physical panel size as the original four-column figure.
  width <- (3580/300)*(ncolumns+.46)/4.46
  height <- 1740/300
  files <- character()
  if(png) {
    grDevices::png(file,width=width,height=height,units='in',res=600,
      type=if(capabilities('cairo')) 'cairo' else getOption('bitmapType'))
    tryCatch(draw(),finally=grDevices::dev.off())
    files <- c(files,png=file)
    if(basename(file)=='Figure_asymmetric_tournaments.png') {
      alias <- file.path(dirname(file),'Figure_asymmetric_tournaments_strategy_map_HDR.png')
      if(!file.copy(file,alias,overwrite=TRUE)) stop('Could not save the HDR filename alias.')
      files <- c(files,png_hdr=alias)
    }
  }
  if(pdf) {
    pdf_file <- paste0(tools::file_path_sans_ext(file),'.pdf')
    if(capabilities('cairo')) grDevices::cairo_pdf(pdf_file,width=width,height=height,family='Arial')
    else grDevices::pdf(pdf_file,width=width,height=height,useDingbats=FALSE)
    tryCatch(draw(),finally=grDevices::dev.off())
    files <- c(files,pdf=pdf_file)
  }
  if(interactive() || grDevices::dev.cur()!=1L) draw()
  files
}
