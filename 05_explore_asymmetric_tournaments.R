# Edit these settings, then source this file from the scripts folder.
# By default, draw in R's plotting window, without a legend or file exports.
# Set png=TRUE and/or pdf=TRUE to export; legend=TRUE adds the optional legend.
# Requires Rcpp, a C++ compiler, and the supplied 05_* helper files.
# Fixed: random allocations independent of talent, E1=E2=B/2, separate M1,M2.
# Responsive: a self-consistent population of talent-informed best responses.
# Investment grids resolve interior maxima; increase nx/nxpop to check accuracy.
# Population integration splits at strategy switches rather than sampling a fixed talent grid.
# nquad controls population quadrature (minimum 64); try 96 or 128 to check convergence.
# The first call compiles the C++ helper; later calls reuse its cached library.
settings <- list(beta1=.6,beta2=.4,phi1=.66,phi2=1,B=50,
                 mode="both",alim=c(0,2.5),nalpha=301L,nx=801L,
                 nxpop=801L,nquad=31L,nroot=501L,
                 legend=FALSE,png=FALSE,pdf=FALSE)

source("05_random_allocation_reference.R")
source("05_asymmetric_model.R")
source("05_explorer_population.R")

explore_asymmetric <- function(beta1,beta2,phi1,phi2,B=50,
    mode=c("both","fixed","responsive"),alim=c(0,2.5),
    nalpha=301L,nx=801L,nxpop=801L,nquad=31L,nroot=501L,
    file="Figure_asymmetric_tournaments_explorer.png",
    legend=FALSE,png=FALSE,pdf=FALSE) {
  mode <- match.arg(mode)
  for(flag in list(legend=legend,png=png,pdf=pdf))
    if(!is.logical(flag) || length(flag)!=1L || is.na(flag))
      stop("legend, png and pdf must each be TRUE or FALSE.")
  if(any(!is.finite(c(beta1,beta2))) || min(beta1,beta2)<=0 || max(beta1,beta2)>=1)
    stop("Training exponents beta1 and beta2 must lie strictly between 0 and 1.")
  if(any(!is.finite(c(phi1,phi2,B))) || min(phi1,phi2,B)<=0)
    stop("phi1, phi2 and B must be finite and positive.")
  if(length(alim)!=2 || any(!is.finite(alim)) || alim[1]<0 || diff(alim)<=0)
    stop("alim must be an increasing nonnegative interval.")
  if(any(c(nalpha,nx,nxpop,nquad,nroot)<3)) stop("Resolutions must be at least 3.")
  nx <- as.integer(nx); nxpop <- as.integer(nxpop)
  # Include B/2 exactly when testing whether the equal split is optimal.
  if(nx%%2==0) nx <- nx+1L
  if(nxpop%%2==0) nxpop <- nxpop+1L
  model <- new_asymmetric_model(nquad)
  fixed <- random_allocation_reference(beta1,beta2,phi1,phi2,B)
  population <- explorer_population(beta1,beta2,phi1,phi2,B,fixed$s,
    responsive=mode!="fixed",nquad=max(64L,nquad),nx=nxpop)
  states <- list(fixed=fixed)
  if(mode!="fixed") states$responsive <- population$responsive
  if(mode=="responsive") states$fixed <- NULL
  investments <- lapply(names(states),function(name) {
    st <- states[[name]]
    if(name=="responsive") return(st)
    population$first
  })
  names(investments) <- names(states)
  alpha <- seq(alim[1],alim[2],length.out=nalpha)
  maps <- lapply(states,function(st)
    model$strategy_map(alpha,B,beta1,beta2,phi1,phi2,st$s,nx,nroot))
  hdr <- model$talent_hdr_curves()
  cols <- colorRampPalette(c("#315A9B","white","#B54343"))(801)
  draw <- function() {
    layout(matrix(seq_len(length(states)+as.integer(legend)),nrow=1),
           widths=c(rep(1,length(states)),if(legend) .26 else numeric()))
    par(mar=c(4,4,5,1.6),oma=c(0,0,0,0),
        mgp=c(2.5,.7,0),pty="s",cex=1)
    for(i in seq_along(states)) {
      sm <- maps[[i]]; st <- states[[i]]
      image(alpha,alpha,sm$z,col=cols,zlim=c(0,1),useRaster=TRUE,
            xlim=alim,ylim=alim,xaxs="i",yaxs="i",
            xlab=expression(alpha[1]~"(talent in tournament 1)"),
            ylab=if(i==1) expression(alpha[2]~"(talent in tournament 2)") else "")
      # Fill the global maxxer regions as polygons of their computed
      # boundary curves. The muted allocation field remains a raster, but
      # the bright specialization edges no longer follow pixel-cell edges.
      clip <- function(z) pmax(alim[1],pmin(alim[2],z))
      polygon(c(alim[1],clip(sm$b2),alim[1]),
              c(alim[1],alpha,alim[2]),col="#08F",border=NA)
      polygon(c(alim[1],alpha,alim[2]),
              c(alim[1],clip(sm$b1),alim[1]),col="#E06",border=NA)
      lines(sm$b2,alpha,col=adjustcolor("#08F",.9),lwd=1.05)
      lines(alpha,sm$b1,col=adjustcolor("#E06",.9),lwd=1.05)
      for(z in hdr) lines(z[1,],z[2,],col="#EEE",lwd=.8)
      abline(0,1,lty=2,lwd=.8)
      hc <- model$half_curve(alpha,B,beta1,beta2,phi1,phi2,st$s,alim,nx)
      ii <- which(is.finite(hc$x)&hc$global)
      if(length(ii)) for(j in split(ii,cumsum(c(1,diff(ii)!=1))))
        if(length(j)>1) lines(hc$x[j],hc$y[j],col="#FF0",lwd=3)
      box()
      bar <- investments[[i]]
      x0 <- alim[2]+.025*diff(alim); x1 <- alim[2]+.075*diff(alim)
      split <- alim[1]+(bar$E1/B)*diff(alim)
      rect(x0,alim[1],x1,split,col="#E06",border=NA,xpd=NA)
      rect(x0,split,x1,alim[2],col="#08F",border=NA,xpd=NA)
      title(main=if(names(states)[i]=="fixed") expression("Fixed "*E[1]==E[2])
        else expression("Responsive "*E[t]),line=3.7,cex.main=1.15)
      tx1 <- alim[1]+.30*diff(alim); tx2 <- alim[1]+.70*diff(alim)
      mtext(bquote(beta[1]==.(formatC(beta1,format="f",digits=2))),
            side=3,line=1.7,at=tx1,cex=1.05,col="#E06")
      mtext(bquote(beta[2]==.(formatC(beta2,format="f",digits=2))),
            side=3,line=1.7,at=tx2,cex=1.05,col="#08F")
      mtext(bquote(varphi[1]==.(formatC(phi1,format="f",digits=2))),
            side=3,line=.3,at=tx1,cex=1.05,col="#E06")
      mtext(bquote(varphi[2]==.(formatC(phi2,format="f",digits=2))),
            side=3,line=.3,at=tx2,cex=1.05,col="#08F")
      text((x0+x1)/2,alim[1],expression(E[1]),pos=3,cex=.65,col="white",xpd=NA)
      text((x0+x1)/2,alim[2],expression(E[2]),pos=1,cex=.65,col="white",xpd=NA)
    }
    if(legend) {
      # Same narrow bar, coloured endpoint/equal-split strips, ticks and
      # vertical label as Figure_asymmetric_tournaments.
      par(mar=c(4,.55,5,1.6),pty="m",xpd=NA,cex=1)
      plot.new()
      plot.window(xlim=c(-.08,1),ylim=c(0,1),xaxs="i",yaxs="i")
      zz <- seq(0,1,length.out=801)
      for(k in 1:800)
        rect(.34,zz[k],.64,zz[k+1],col=cols[k],border=NA)
      rect(.34,0,.64,1,border="black",lwd=.8)
      wd <- .006
      rect(.34,-wd,.64,wd,col="#08F",border=1,lwd=.8)
      rect(.34,.5-wd,.64,.5+wd,col="#FF0",border=1,lwd=.8)
      rect(.34,1-wd,.64,1+wd,col="#E06",border=1,lwd=.8)
      ticks <- seq(0,1,.25)
      segments(.64,ticks,.75,ticks,lwd=.8)
      text(.88,ticks,sprintf("%.2f",ticks),adj=c(0,.5),cex=.9)
      text(.10,.5,"Relative time invested in tournament 1",srt=90,cex=0.85)
      text(.49,1.05,"1maxxers",cex=0.85,col="#E06")
      text(.49,-.05,"2maxxers",cex=0.85,col="#08F")
    }
  }
  width <- 3.8*length(states)+if(legend) 1.0 else 0
  height <- 4.3
  files <- character()
  if(png) {
    grDevices::png(file,width=width,height=height,units="in",res=600,
      type=if(capabilities('cairo')) 'cairo' else getOption('bitmapType'))
    tryCatch(draw(),finally=grDevices::dev.off())
    files <- c(files,png=file)
  }
  if(pdf) {
    pdf_file <- paste0(tools::file_path_sans_ext(file),'.pdf')
    if(capabilities('cairo')) {
      grDevices::cairo_pdf(pdf_file,width=width,height=height,family='Arial')
    } else {
      grDevices::pdf(pdf_file,width=width,height=height,useDingbats=FALSE)
    }
    tryCatch(draw(),finally=grDevices::dev.off())
    files <- c(files,pdf=pdf_file)
  }
  # In R/RStudio this opens or redraws the plotting window. In a headless
  # Rscript run, draw only if the caller supplied a device, so the null
  # device cannot silently create Rplots.pdf when exports are disabled.
  if(interactive() || grDevices::dev.cur()!=1L) draw()
  print(do.call(rbind,lapply(states,function(st) unlist(st[c("E1","E2","M1","M2","s")]))))
  cat("Mean investments after the plotted best responses:\n")
  print(do.call(rbind,lapply(investments,function(st) unlist(st[c("E1","E2")]))))
  # Point query: compare every grid-detected local maximum and both endpoints;
  # refine each interior candidate before comparing utilities. Increase nx if
  # unusually narrow competing peaks are suspected. Ties return one maximizer.
  optimal_investment <- function(alpha1,alpha2,state=names(states)[1]) {
    if(!state %in% names(states)) stop("Choose a computed state.")
    if(length(alpha1)!=1 || length(alpha2)!=1 ||
       any(!is.finite(c(alpha1,alpha2))) || min(alpha1,alpha2)<0)
      stop("Supply two finite nonnegative talents.")
    st <- states[[state]]
    u <- function(x) st$s*model$w(x,alpha1,beta1,phi1)+
      model$w(B-x,alpha2,beta2,phi2)
    xx <- seq(0,B,length.out=nx); yy <- u(xx)
    k <- which(yy[2:(nx-1)]>=yy[1:(nx-2)] & yy[2:(nx-1)]>=yy[3:nx])+1L
    if(alpha1==0 && alpha2==0) return(c(x1=0,x2=B,utility=0))
    candidates <- c(0,B,vapply(k,function(i)
      optimize(u,c(xx[i-1],xx[i+1]),maximum=TRUE,tol=1e-10)$maximum,1))
    x <- candidates[which.max(u(candidates))]
    c(x1=x,x2=B-x,utility=(st$E2/st$M2)*u(x))
  }
  invisible(list(population=states,investments=investments,alpha=alpha,maps=maps,
                 optimal_investment=optimal_investment,file=if(png) file else NULL,files=files))
}

#result <- do.call(explore_asymmetric,settings)

#explore_asymmetric(beta1=0.8,beta2=0.4,phi1=1,phi2=1,legend=T,png=T)

# for(i in c(0.2,0.4,0.6,0.8,0.9)){
#   explore_asymmetric(beta1=i,beta2=0.5,phi1=0.8,phi2=0.8)
# }
# 
# for(i in c(0.2,0.4,0.6,0.8,0.9)){
#   explore_asymmetric(beta1=i,beta2=0.5,phi1=1.2,phi2=1.2)
# }

# explore_asymmetric(beta1=0.8,beta2=0.4,phi1=0.66,phi2=0.66) #probably will not be used



# explore_asymmetric(beta1=0.92,beta2=0.4,phi1=1.5,phi2=1)
# explore_asymmetric(beta1=0.8,beta2=0.4,phi1=1.5,phi2=1)
# explore_asymmetric(beta1=0.8,beta2=0.5,phi1=1.5,phi2=1)
# explore_asymmetric(beta1=0.9,beta2=0.9,phi1=2.1,phi2=0.5)

# explore_asymmetric(beta1=0.9,beta2=0.5,phi1=2,phi2=1.5)
# explore_asymmetric(beta1=0.5,beta2=0.9,phi1=1.4,phi2=0.3)
# explore_asymmetric(beta1=0.5,beta2=0.9,phi1=1.3,phi2=2.2)
# explore_asymmetric(beta1=0.5,beta2=0.1,phi1=1.8,phi2=2,alim=c(0,0.5))
# explore_asymmetric(beta1=0.5,beta2=0.1,phi1=1.8,phi2=2,alim=c(0,10000))