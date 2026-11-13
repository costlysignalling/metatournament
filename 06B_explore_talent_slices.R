# Interactive counterpart of Figure S16. Run from the scripts folder.
# source('06B_explore_talent_slices.R')
# many_slices(beta=c(.9,.6,.6),phi=c(.8,1.5,1),B=50)
# Four or more tournaments: alpha_1 and alpha_2 are plotted, alpha_3 varies
# between columns; other_talents fixes alpha_4,... (default: all equal to 1).
# Population feedback always uses ALL tournaments and the full population.
# A returned many_example() result can be reused: many_slices(result=example).
# Nothing is written to disk unless png=TRUE or pdf=TRUE.

if(!exists('many_allocate',mode='function'))
  source('06B_predict_many_tournaments.R')

many_slices <- function(beta=c(.9,.6,.6),phi=c(.8,1.5,1),B=50,
  cuts=qlnorm(c(.1,.3,.5,.7,.9),sdlog=sqrt(log((1+sqrt(2))/2))),
  other_talents=NULL,talent_limit=2.5,n_grid=61L,n_population=4096L,
  result=NULL,png=FALSE,pdf=FALSE,file='Figure_many_talent_slices_preview') {
  if(!is.null(result)) {
    if((!missing(beta) && !isTRUE(all.equal(beta,result$beta))) ||
       (!missing(phi) && !isTRUE(all.equal(phi,result$phi))) ||
       (!missing(B) && !isTRUE(all.equal(B,result$B))))
      stop('A reused result must have the same beta, phi and B.')
    beta <- result$beta;phi <- result$phi;B <- result$B
  }
  T <- length(beta)
  if(T<3 || length(phi)!=T || any(!is.finite(c(beta,phi,B,cuts,talent_limit,n_grid))) ||
     any(beta<=0 | beta>=1) || any(phi<=0) || length(B)!=1 || B<=0 ||
     length(cuts)<1 || any(cuts<0) || length(talent_limit)!=1 || talent_limit<=0 ||
     length(n_grid)!=1 || n_grid<11 || n_grid!=as.integer(n_grid))
    stop('Use at least three tournaments, matching beta/phi, positive B, nonnegative cuts and n_grid >= 11.')
  if(is.null(other_talents))other_talents <- rep(1,T-3)
  if(length(other_talents)!=T-3 || any(!is.finite(other_talents)) || any(other_talents<0))
    stop('other_talents must contain nonnegative talents for tournaments 4,...,T.')
  if(is.null(result)) {
    population <- many_lognormal_population(n_population,T)
    fixed <- many_random_reference(beta,phi,B)
    message('Solving the responsive population...')
    responsive <- many_responsive(population,beta,phi,B,start=fixed$scaler)
  } else {fixed <- result$fixed;responsive <- result$responsive}
  if(!isTRUE(responsive$converged))
    stop('Responsive population did not converge; inspect many_responsive() diagnostics.')
  scalers <- list(fixed$scaler,responsive$scaler)
  specialist_colours <- c('#EE0066','#0088FF','#40A56A',
    if(T>3)hcl.colors(T-3,'Dark 3') else character())
  muted <- c('#B54343','#315A9B','#46765E',
    if(T>3)hcl.colors(T-3,'Dark 3') else character())
  boundaries <- c('#A70047','#005BBB','#22663F',
    if(T>3)hcl.colors(T-3,'Dark 3') else character())
  alpha <- seq(0,talent_limit,length.out=n_grid)
  plane <- as.matrix(expand.grid(alpha1=alpha,alpha2=alpha))
  paint <- function(X) {
    shares <- X/B
    rgbmat <- shares%*%t(col2rgb(muted)/255)
    saturation <- apply(shares,1,max)-apply(shares,1,min)
    rgbmat <- 1-(1-rgbmat)*saturation
    col <- rgb(pmin(1,pmax(0,rgbmat[,1])),pmin(1,pmax(0,rgbmat[,2])),
               pmin(1,pmax(0,rgbmat[,3])))
    # Solid specialist colours are drawn only as vector polygons below.
    # Keeping the raster muted avoids bright, jagged pixel fringes.
    col
  }
  d <- function(a,x,t,scaler) scaler[t]*phi[t]*a*beta[t]*(1+x)^(beta[t]-1)*
    (1+a*expm1(beta[t]*log1p(x)))^(phi[t]-1)
  # In concave slices, the entry conditions give exact specialist polygons.
  # Nonconcave slices locate specialist boundaries by global best responses;
  # local marginal conditions alone are not sufficient there.
  polygons <- function(talents,scaler) {
    entry <- scaler*phi*beta;out <- vector('list',T)
    fixed_entry <- entry[3:T]*talents[3:T]
    for(t in 3:T) {
      final <- d(talents[t],B,t,scaler)
      rivals <- setdiff(3:T,t)
      if(length(rivals) && any(final<entry[rivals]*talents[rivals]))next
      xmax <- min(talent_limit,final/entry[1]);ymax <- min(talent_limit,final/entry[2])
      if(xmax>0 && ymax>0)out[[t]] <- cbind(c(0,xmax,xmax,0),c(0,0,ymax,ymax))
    }
    for(t in 1:2) {
      cutoff <- max(fixed_entry)
      if(d(talent_limit,B,t,scaler)<=cutoff)next
      start <- if(cutoff==0)0 else uniroot(function(a)d(a,B,t,scaler)-cutoff,
        c(0,talent_limit),tol=1e-12)$root
      own <- seq(start,talent_limit,length.out=601)
      rival <- pmin(talent_limit,d(own,B,t,scaler)/entry[3-t])
      P <- cbind(c(start,own,talent_limit),c(0,rival,0))
      out[[t]] <- if(t==2)P[,2:1] else P
    }
    out
  }
  numerical_polygons <- function(A,X,scaler) {
    out <- vector('list',T)
    # For fixed alpha1, specialization in tournament 2 increases with alpha2;
    # specialization in any other tournament decreases with alpha2. This
    # follows from positive talent-investment increasing differences.
    for(t in seq_len(T)) {
      bottom <- X[seq_len(n_grid),t]>=B*(1-1e-8)
      top <- X[(n_grid-1L)*n_grid+seq_len(n_grid),t]>=B*(1-1e-8)
      rising <- t==2
      threshold <- if(rising)ifelse(bottom,0,talent_limit) else ifelse(top,talent_limit,0)
      pending <- which(bottom!=top)
      if(length(pending)) {
        lo <- rep(0,length(pending));hi <- rep(talent_limit,length(pending))
        for(iter in seq_len(23)) {
          mid <- (lo+hi)/2;test <- A[pending,,drop=FALSE];test[,2] <- mid
          fit <- many_allocate(test,beta,phi,B,scaler)
          if(!all(fit$within_tolerance))stop('Boundary calculation exceeds global utility tolerance.')
          inside <- fit$allocation[,t]>=B*(1-1e-8)
          move_low <- if(rising)!inside else inside
          lo[move_low] <- mid[move_low];hi[!move_low] <- mid[!move_low]
        }
        threshold[pending] <- (lo+hi)/2
      }
      if((rising && any(threshold<talent_limit)) || (!rising && any(threshold>0)))
        out[[t]] <- cbind(c(alpha,rev(alpha)),c(threshold,rep(if(rising)talent_limit else 0,n_grid)))
    }
    out
  }
  slices <- lapply(1:2,function(row)lapply(seq_along(cuts),function(j) {
    message(sprintf('%s: talent-3 slice %d/%d',if(row==1)'Fixed' else 'Responsive',j,length(cuts)))
    A <- cbind(plane,cuts[j])
    if(T>3)A <- cbind(A,matrix(rep(other_talents,each=nrow(A)),nrow=nrow(A)))
    fit <- many_allocate(A,beta,phi,B,scalers[[row]])
    if(!all(fit$within_tolerance))stop('Slice allocations exceed the global utility tolerance.')
    margin <- sapply(seq_len(T),function(t) {
      rival <- sapply(setdiff(seq_len(T),t),function(k)d(A[,k],0,k,scalers[[row]]))
      d(A[,t],B,t,scalers[[row]])-apply(rival,1,max)
    })
    list(allocation=fit$allocation,colours=paint(fit$allocation),
         concave=all(fit$concave),margins=margin,
         polygons=if(all(fit$concave))polygons(A[1,],scalers[[row]]) else
           numerical_polygons(A,fit$allocation,scalers[[row]]),
         max_relative_gap=max(fit$relative_gap_bound))
  }))
  sdlog <- sqrt(log((1+sqrt(2))/2));th <- seq(0,2*pi,length.out=601)
  hdr <- lapply(qchisq(c(.5,.9,.99),df=2,ncp=2*sdlog^2),function(q)
    exp(-sdlog^2+sdlog*sqrt(q)*rbind(cos(th),sin(th))))
  contour_zero <- function(field,col,lwd,level=0) {
    finite <- field[is.finite(field)]
    if(length(finite) && min(finite)<level && max(finite)>level)
      contour(alpha,alpha,matrix(field,n_grid,n_grid),levels=level,
        add=TRUE,drawlabels=FALSE,col=col,lwd=lwd)
  }
  draw <- function() {
    old <- par(no.readonly=TRUE);on.exit(par(old))
    size <- dev.size('in');W <- size[1];H <- size[2];K <- length(cuts)
    left <- .96;right <- .22;gap <- .33;top <- .8;bottom <- .48;rowgap <- .68
    panel <- min((W-left-right-gap*(K-1))/K,(H-top-bottom-rowgap)/2)
    if(panel<=.2)stop('Enlarge the plotting window to display these panels.')
    xstart <- left+(W-left-right-gap*(K-1)-K*panel)/2
    ystart <- bottom+(H-top-bottom-rowgap-2*panel)/2
    centres <- numeric(2)
    for(row in 1:2)for(j in seq_along(cuts)) {
      x <- xstart+(j-1)*(panel+gap);y <- ystart+(2-row)*(panel+rowgap)
      centres[row] <- (y+panel/2)/H
      par(fig=c(x/W,(x+panel)/W,y/H,(y+panel)/H),mar=rep(0,4),
          new=!(row==1 && j==1),mgp=c(1.5,.45,0),tcl=-.2,xpd=FALSE,cex=1)
      s <- slices[[row]][[j]];X <- s$allocation
      plot(NA,xlim=range(alpha),ylim=range(alpha),xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='')
      rasterImage(as.raster(t(matrix(s$colours,n_grid,n_grid)[,n_grid:1])),
        0,0,talent_limit,talent_limit,interpolate=TRUE)
      for(t in seq_len(T))if(!is.null(s$polygons[[t]]))
        polygon(s$polygons[[t]],col=specialist_colours[t],border=NA)
      for(curve in hdr)lines(curve[1,],curve[2,],col='#EEEEEE',lwd=.8)
      abline(0,1,lty=3,lwd=.65,col='grey40')
      # Use the same vector paths for the fill and its outline.
      for(t in seq_len(T))if(!is.null(s$polygons[[t]]))
        polygon(s$polygons[[t]],col=NA,border=boundaries[t],lwd=1.2)
      fields <- lapply(combn(T,2,simplify=FALSE),function(pair) {
        other <- setdiff(seq_len(T),pair)
        field <- X[,pair[1]]-X[,pair[2]]
        field[X[,pair[1]]+X[,pair[2]]<1e-4*B |
          pmax(X[,pair[1]],X[,pair[2]])<apply(X[,other,drop=FALSE],1,max)-.02] <- NA
        field
      })
      # Draw all outlines first, then all yellow strokes, including junctions.
      for(field in fields)contour_zero(field,'black',2)
      for(field in fields)contour_zero(field,'#FFFF00',1.3)
      ticks <- pretty(c(0,talent_limit),n=5);ticks <- ticks[ticks>=0 & ticks<=talent_limit]
      axis(1,at=ticks,cex.axis=.85)
      axis(2,at=ticks,labels=if(j==1)TRUE else FALSE,las=1,cex.axis=.85);box()
      if(row==2)mtext(expression(alpha[1]),1,line=1.45,xpd=NA)
      if(j==1)mtext(expression(alpha[2]),2,line=2.05,xpd=NA)
      if(row==1)mtext(bquote(alpha[3]==.(formatC(cuts[j],format='f',digits=2))),
        3,line=.35,adj=0,cex=1.1,xpd=NA)
      # The plotting region, not merely its enclosing device, must be square.
      stopifnot(abs(diff(par('pin')))<1e-6)
    }
    par(fig=c(0,1,0,1),new=TRUE,mar=rep(0,4),xpd=NA);plot.new()
    plot.window(c(0,1),c(0,1),xaxs='i',yaxs='i')
    text(.12/W,centres[1],expression(Fixed~equal~E[t]~reference),srt=90,cex=.95)
    text(.12/W,centres[2],expression(Responsive~E[t]),srt=90,cex=.95)
    labels <- c(paste('Tournament',seq_len(T)),'Equal leading investments')
    cex <- .95;symbol <- .027;space <- .01;between <- .025
    widths <- strwidth(labels,cex=cex);total <- sum(widths)+(T+1)*(symbol+space)+T*between
    if(total>.97){cex <- cex*.97/total;widths <- strwidth(labels,cex=cex)}
    total <- sum(widths)+(T+1)*(symbol+space)+T*between
    x <- (1-total)/2;y <- 1-.20/H
    for(t in seq_along(labels)) {
      if(t<=T)points(x+symbol/2,y,pch=15,col=specialist_colours[t],cex=1.4)
      else {segments(x,y,x+symbol,y,col='black',lwd=2);segments(x,y,x+symbol,y,col='#FFFF00',lwd=1.3)}
      text(x+symbol+space,y,labels[t],adj=c(0,.5),cex=cex)
      x <- x+symbol+space+widths[t]+between
    }
  }

  draw()
  if(png) {
    grDevices::png(paste0(file,'.png'),width=2.2*length(cuts),height=5.5,units='in',res=600,type='cairo')
    tryCatch(draw(),finally=dev.off())
  }
  if(pdf) {
    grDevices::pdf(paste0(file,'.pdf'),width=2.2*length(cuts),height=5.5)
    tryCatch(draw(),finally=dev.off())
  }
  invisible(list(beta=beta,phi=phi,B=B,cuts=cuts,other_talents=other_talents,
    alpha=alpha,slices=slices,fixed=fixed,responsive=responsive,draw=draw))
}


# Three tournaments, as in S16
many_slices(beta=c(.6,.4,.8), phi=c(.66,1,1.2), B=50)
many_slices(beta=c(.9,.6,.6), phi=c(.8,1.5,1), B=50)

# Your four-tournament example, holding talent 4 at 1.5
many_slices(beta=c(.6,.8,.4,.8), phi=c(.66,1,1.2,1.4), B=50, other_talents=1.5)

