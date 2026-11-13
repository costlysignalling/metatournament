# Geometry of the 50:50 marginal-return curve and the equal-talent diagonal.
# Run from the scripts folder. Calculations are cached separately from plotting.
# The supplied input RDS shares Script 06's deterministic population equilibria.
# Linear cases agree with the exact one-dimensional population integrals.
# Unlike Script 06's old fixed-scaler reference, this fixed block recalculates
# M1 and M2 from uniform random allocations with E1=E2=B/2.
for(pkg in c('sf','isoband')) if(!requireNamespace(pkg,quietly=TRUE))
  stop('Install package ',pkg,' to draw the vector region boundaries.')
source('05_random_allocation_reference.R')
source('05_equal_talent_geometry.R')
source('05_interpolate_equal_talent_beta.R')

B <- 50
# Move the COMPLETE bottom gradient-legend row (bars, titles, ticks, labels).
# Increase this one number to move it lower; 0 restores the previous position.
gradient_legend_drop <- .06
n_plot <- 501L
force_recompute <- identical(Sys.getenv('RECOMPUTE_5050_GEOMETRY'),'1')
input_file <- '05_equal_talent_geometry_inputs.rds'
cache_file <- '05_equal_talent_5050_geometry_cache.rds'
inputs <- readRDS(input_file)
shared_population <- '06_coexistence_beta_inputs.rds'
if(file.exists(shared_population) &&
   !identical(inputs$source$md5,unname(tools::md5sum(shared_population)))) {
  source('05_prepare_equal_talent_geometry_inputs.R')
  inputs <- readRDS(input_file)
  force_recompute <- TRUE
}
stopifnot(inputs$B==B,identical(inputs$phis,c(.66,1,1.5)))
phis <- inputs$phis; betas <- inputs$betas
plot_betas <- seq(min(betas),max(betas),length.out=n_plot)
input_md5 <- unname(tools::md5sum(input_file))

bilinear <- function(z,x,y,xs=betas,ys=betas) {
  i <- pmax(1L,pmin(length(xs)-1L,findInterval(x,xs)))
  j <- pmax(1L,pmin(length(ys)-1L,findInterval(y,ys)))
  tx <- (x-xs[i])/(xs[i+1]-xs[i]); ty <- (y-ys[j])/(ys[j+1]-ys[j])
  (1-tx)*(1-ty)*z[cbind(i,j)]+tx*(1-ty)*z[cbind(i+1,j)]+
    (1-tx)*ty*z[cbind(i,j+1)]+tx*ty*z[cbind(i+1,j+1)]
}

colours <- c(colorRampPalette(c('#315A9B','white','#B54343'))(801),
             '#E06','#08F','white')
interpolation_md5 <- unname(tools::md5sum('05_interpolate_equal_talent_beta.R'))
if(file.exists(cache_file)) {
  old_cache <- readRDS(cache_file)
  if(!identical(old_cache$version,2L) ||
     !identical(old_cache$interpolation_md5,interpolation_md5)) force_recompute <- TRUE
}
if(file.exists(cache_file) && !force_recompute) {
  cache <- readRDS(cache_file)
  if(!identical(cache$input_md5,input_md5) || cache$n_plot!=n_plot || cache$version!=2L)
    stop('Geometry cache differs from the settings. Set force_recompute=TRUE.')
  panels <- cache$panels
} else {
  moments <- sapply(phis,function(phi) vapply(plot_betas,function(beta)
    random_allocation_reference(beta,beta,phi,phi,B)$M1,1))
  panels <- list()
  fine_grid <- expand.grid(beta1=plot_betas,beta2=plot_betas)
  coarse_grid <- expand.grid(beta1=betas,beta2=betas)
  for(mode in c('fixed','responsive')) for(row in 1:3) for(column in 1:3) {
    message(mode,': phi1=',phis[column],', phi2=',phis[row])
    if(mode=='fixed') {
      ratio <- outer(moments[,column],moments[,row],function(a,b)b/a)
      coarse_m1 <- approx(plot_betas,moments[,column],xout=betas)$y
      coarse_m2 <- approx(plot_betas,moments[,row],xout=betas)$y
      coarse_ratio <- outer(coarse_m1,coarse_m2,function(a,b)b/a)
    } else {
      coarse_ratio <- inputs$ratios[[column,row]]
      if(column==row) {
        direct <- if(file.exists(shared_population))
          readRDS(shared_population)$panels[[column,row]]$direct else NULL
        low <- equal_talent_beta_low_field(coarse_ratio,betas,plot_betas,B,direct)
        low_coarse <- equal_talent_beta_low_field(coarse_ratio,betas,betas,B,direct)
        ratio <- exp(low-outer(log(plot_betas),log(plot_betas),'-')-
          outer(plot_betas,plot_betas,'-')*log1p(B/2))
        coarse_ratio <- exp(low_coarse-outer(log(betas),log(betas),'-')-
          outer(betas,betas,'-')*log1p(B/2))
      } else ratio <- matrix(exp(bilinear(log(coarse_ratio),fine_grid$beta1,fine_grid$beta2)),n_plot)
    }
    if(column==row) { diag(ratio)<-1;diag(coarse_ratio)<-1 }
    coarse <- equal_talent_geometry(coarse_grid$beta1,coarse_grid$beta2,
                                  phis[column],phis[row],as.vector(coarse_ratio),B)
    coarse$global_equal_split <- check_equal_split_global(coarse,phis[column],phis[row],B)
    fine <- equal_talent_geometry(fine_grid$beta1,fine_grid$beta2,
                                phis[column],phis[row],as.vector(ratio),B)
    colour_id <- rep(804L,nrow(fine))
    colour_id[fine$type==1] <- 802L
    colour_id[fine$type==2] <- 803L
    w <- 1/(1+fine$alpha_cross)
    k <- which(fine$type==3);colour_id[k] <- 401L+as.integer(round(400*w[k]))
    k <- which(fine$type==4);colour_id[k] <- 401L-as.integer(round(400*w[k]))
    key <- paste(mode,row,column,sep='_')
    panels[[key]] <- list(mode=mode,row=row,column=column,
      phi1=phis[column],phi2=phis[row],colours=matrix(colour_id,n_plot),
      always_above=matrix(fine$type==1,n_plot),always_below=matrix(fine$type==2,n_plot),
      coarse=coarse,bad=matrix(!is.na(coarse$global_equal_split)&!coarse$global_equal_split,length(betas)))
  }
  saveRDS(list(version=2L,input_md5=input_md5,interpolation_md5=interpolation_md5,
               n_plot=n_plot,panels=panels),cache_file)
  tab <- do.call(rbind,lapply(panels,function(panel) {
    z <- panel$coarse;z$mode<-panel$mode;z$phi1<-panel$phi1;z$phi2<-panel$phi2;z
  }))
  con <- gzfile('equal_talent_5050_crossings.csv.gz','wt')
  write.csv(tab,con,row.names=FALSE);close(con)
}

band_geometry <- function(mask,x=plot_betas) {
  if(!any(mask)) return(NULL)
  bands <- isoband::isobands(x,x,t(1*mask),.5,1.5)
  if(!length(bands[[1]]$x)) return(NULL)
  sf::st_sfc(isoband::iso_to_sfg(bands)[[1]])
}
draw_solid_region <- function(mask,colour) {
  geom <- band_geometry(mask)
  if(!is.null(geom)) plot(geom,add=TRUE,col=colour,border=NA)
}
draw_hatching_one <- function(mask,crossing_mask,colour) {
  geom <- band_geometry(mask,betas)
  if(is.null(geom)) return(invisible(NULL))
  crossing <- band_geometry(crossing_mask)
  if(is.null(crossing)) return(invisible(NULL))
  geom <- suppressWarnings(sf::st_intersection(geom,crossing))
  offsets <- seq(-1,1,by=.055)
  hatch <- sf::st_sfc(sf::st_multilinestring(lapply(offsets,function(d)
    rbind(c(-1,-1+d),c(2,2+d)))))
  clipped <- suppressWarnings(sf::st_intersection(hatch,geom))
  if(length(clipped)) plot(clipped,add=TRUE,col=adjustcolor(colour,.8),lwd=.85)
}
# Colour the hatch by the better specialist at the stationary crossing.
draw_hatching <- function(panel) {
  z <- panel$coarse
  u1 <- z$ratio*expm1(panel$phi1*log1p(z$alpha_cross*expm1(z$beta1*log1p(B))))
  u2 <- expm1(panel$phi2*log1p(z$alpha_cross*expm1(z$beta2*log1p(B))))
  winner1 <- matrix(!is.na(u1)&!is.na(u2)&u1>=u2,length(betas))
  draw_hatching_one(panel$bad & winner1,panel$colours<=801L,'#E06')
  draw_hatching_one(panel$bad & !winner1,panel$colours<=801L,'#08F')
}

# The total of two independent reference lognormal talents exceeds this
# cutoff for 50% of actors. Integrate the convolution; no simulation needed.
talent_sdlog <- sqrt(log((1+sqrt(2))/2))
talent_sum_cdf <- function(total) integrate(function(a)
  dlnorm(a,0,talent_sdlog)*plnorm(total-a,0,talent_sdlog),
  0,total,rel.tol=1e-11)$value
talent_sum_cutoff <- uniroot(function(total) talent_sum_cdf(total)-.5,
                             c(2,5),tol=1e-10)$root

# Trace crossings at the 50th percentile of total reference talent.
# Evaluate the marginal-return equation on the numerical input grid, then
# interpolate its zero contour rather than contouring clipped display colours.
draw_talent_sum_contour <- function(panel) {
  grid <- expand.grid(beta1=betas,beta2=betas)
  ratio <- panel$coarse$ratio
  L <- 1+B/2
  alpha <- talent_sum_cutoff/2
  z <- log(ratio)+log(panel$phi1*grid$beta1/(panel$phi2*grid$beta2))+
    (grid$beta1-grid$beta2)*log(L)+
    (panel$phi1-1)*log1p(alpha*expm1(grid$beta1*log(L)))-
    (panel$phi2-1)*log1p(alpha*expm1(grid$beta2*log(L)))
  z <- matrix(z,length(betas))
  if(panel$phi1==panel$phi2) {
    # Equal tournament parameters make the whole diagonal coincident. Divide
    # out this trivial zero so it is not mistaken for a talent-sum contour.
    delta <- outer(betas,betas,'-')
    z <- z/delta
    for(i in seq_along(betas)) {
      neighbours <- c(max(1L,i-1L),min(length(betas),i+1L))
      neighbours <- neighbours[neighbours!=i]
      z[i,i] <- mean(z[i,neighbours])
    }
  }
  fine <- expand.grid(beta1=plot_betas,beta2=plot_betas)
  z <- matrix(bilinear(z,fine$beta1,fine$beta2),n_plot)
  segments <- contourLines(plot_betas,plot_betas,z,levels=0)
  region <- band_geometry(panel$colours<=801L)
  if(!length(segments) || is.null(region)) return(invisible(NULL))
  lines <- sf::st_sfc(lapply(segments,function(seg)
    sf::st_linestring(cbind(seg$x,seg$y))))
  clipped <- suppressWarnings(sf::st_intersection(lines,region))
  if(length(clipped)) {
    plot(clipped,add=TRUE,col='white',lwd=1.1)
    plot(clipped,add=TRUE,col='#333333',lwd=1.1,lty=2)
  }
}

source('05_parameter_space_example_markers.R')
examples <- parameter_space_example_sets('beta',phis,range(betas))
example_cols <- parameter_space_marker_colours

draw_panel <- function(panel) {
  par(mar=c(1.5,1.5,.5,.5),mgp=c(1.35,.36,0),tcl=-.24,
      xaxs='i',yaxs='i',pty='s',xpd=FALSE)
  image(plot_betas,plot_betas,panel$colours,col=colours,
        breaks=seq(.5,804.5,1),axes=FALSE,xlab='',ylab='',useRaster=TRUE)
  draw_solid_region(panel$always_above,'#E06')
  draw_solid_region(panel$always_below,'#08F')
  draw_hatching(panel)
  draw_talent_sum_contour(panel)
  for(set in names(examples)) {
    z <- examples[[set]]
    k <- which(z[,3]==panel$phi1 & z[,4]==panel$phi2)
    if(length(k)) {
      points(z[k,1],z[k,2],pch=4,cex=1.3,lwd=3.5,col='#333333')
      points(z[k,1],z[k,2],pch=4,cex=1.3,lwd=2.5,col=example_cols[set])
    }
  }
  if(panel$row==3) {
    axis(1,at=c(.2,.4,.6,.8),cex.axis=1.02,lwd=.8,lwd.ticks=.8)
    mtext(expression(beta[1]),side=1,line=1.9,cex=1.18*.9,xpd=NA)
  }
  if(panel$column==1) {
    axis(2,las=0,at=c(.2,.4,.6,.8),cex.axis=1.02,lwd=.8,lwd.ticks=.8)
    mtext(expression(beta[2]),side=2,line=1.4,cex=1.18*.9,las=1)
    usr <- par('usr')
    text(usr[1]-.36,mean(usr[3:4]),
         bquote(varphi[2]==.(formatC(panel$phi2,format='f',digits=2))),
         srt=90,xpd=NA,cex=1.18*1.5)
  }
  if(panel$row==1) mtext(bquote(varphi[1]==.(formatC(panel$phi1,format='f',digits=2))),
                        side=3,line=.22,cex=1.18)
  box(lwd=.8)
}
draw_legend <- function() {
  par(mar=c(.03,.03,.03,.03),xpd=NA,pty='m')
  plot.new()
  width <- par('pin')[1]
  plot.window(xlim=c(0,width),ylim=c(0,1),xaxs='i',yaxs='i')
  # Measure the complete text-plus-symbol items in inches. The gaps between
  # items stay constant regardless of how long their labels are.
  legend_row <- function(labels,kinds,cols,y,cex=1.5) {
    symbol_width <- .15; text_gap <- .08; item_gap <- .24
    measure <- function() strwidth(labels,cex=cex,units='inches')
    text_width <- measure()
    available <- width-.25-length(labels)*(symbol_width+text_gap)-
      (length(labels)-1)*item_gap
    if(sum(text_width)>available) {
      cex <- cex*available/sum(text_width);text_width <- measure()
    }
    x <- (width-sum(text_width)-length(labels)*(symbol_width+text_gap)-
          (length(labels)-1)*item_gap)/2
    for(i in seq_along(labels)) {
      if(kinds[i]=='box') rect(x,y-.047,x+symbol_width,y+.047,
                               col=cols[i],border=NA)
      if(kinds[i]=='hatch') rect(x,y-.047,x+symbol_width,y+.047,
                      density=18,angle=45,col=cols[i],border=NA)
      if(kinds[i]=='cross') {
        points(x+symbol_width/2,y,pch=4,cex=1.3,lwd=3.5,col='black')
        points(x+symbol_width/2,y,pch=4,cex=1.3,lwd=2.5,col=cols[i])
      }
      if(kinds[i]=='line') {
        segments(x,y,x+symbol_width,y,col='white',lwd=1.1)
        segments(x,y,x+symbol_width,y,col=cols[i],lwd=1.1,lty=2)
      }
      text(x+symbol_width+text_gap,y,labels[i],adj=c(0,.5),cex=cex)
      x <- x+symbol_width+text_gap+text_width[i]+item_gap
    }
  }
  legend_row(c('All equal talents lean 1','All equal talents lean 2',
               '50:50 Crosspoint covered by the other maxxing region'),
             c('box','box','hatch'),c('#E06','#08F','#333333'),.60)
  draw_parameter_example_legend(examples,y=.80,
    line_label='50th percentile of talent sum')
  gradient_y <- function(y) y-gradient_legend_drop
  # Red precedes blue, matching the tournament order of the solid colours.
  bars <- list(list(left=.12*width,right=.45*width,side=1,
                    label='Larger sums start leaning 1'),
               list(left=.55*width,right=.88*width,side=-1,
                    label='Larger sums start leaning 2'))
  for(bar in bars) {
    xx <- seq(bar$left,bar$right,length.out=401)
    idx <- 401+bar$side*round(400*seq(1,0,length.out=401))
    for(i in 1:400) rect(xx[i],gradient_y(.21),xx[i+1],gradient_y(.35),col=colours[idx[i]],border=NA)
    text(mean(c(bar$left,bar$right)),gradient_y(.44),bar$label,cex=1.5)
    sums <- c(0,1,2,4,10,Inf)
    at <- bar$left+(bar$right-bar$left)*(1-2/(2+sums))
    segments(at,gradient_y(.21),at,gradient_y(.18),lwd=.6)
    text(at,gradient_y(.12),c('0','1','2','4','10',expression(infinity)),cex=1.3)
  }
  text(.025*width,gradient_y(.28),'Talent sum',adj=0,cex=1.5)
  text(.89*width,gradient_y(.28),'Later crossing\npaler color',pos=4,cex=1.5)
}

draw_equal_talent_geometry <- function() {
  layout(rbind(c(1,2,3,0,10,11,12),c(4,5,6,0,13,14,15),
               c(7,8,9,0,16,17,18),rep(19,7)),
         widths=c(1,1,1,.32,1,1,1),heights=c(1,1,1,.95))
  par(oma=c(.45,3.35,6.55,.45))
  for(mode in c('fixed','responsive')) for(r in 1:3) for(c in 1:3)
    draw_panel(panels[[paste(mode,r,c,sep='_')]])
  draw_legend()
  mtext(expression('Fixed '*E[1]==E[2]),side=3,outer=TRUE,at=.225,line=3.95,cex=1.55,font=2)
  mtext('What actors would do if everyone else allocated time randomly',side=3,outer=TRUE,at=.225,line=2.25,cex=1.1)
  mtext(expression('Responsive '*E[t]),side=3,outer=TRUE,at=.775,line=3.95,cex=1.55,font=2)
  mtext('What actors end up doing in Nash equilibrium',side=3,outer=TRUE,at=.775,line=2.25,cex=1.1)
}
# Figure 06 uses 3500 x (3500*1900/3100) pixels at 300 dpi.
# Its outer margins consume 7 lines at .132 inches per line after layout().
# Add .39 row-height units solely for the taller legend; the panel rows,
# margins, headings, fonts and block gap retain Figure 06's dimensions.
figure_width <- 3500/300
reference_height <- figure_width*1900/3100
reference_row_height <- (reference_height-7*.132)/3.56
figure_height <- reference_height+(.95-.56)*reference_row_height
png('Figure_equal_talent_5050_crossings.png',width=figure_width,height=figure_height,units='in',res=600)
tryCatch(draw_equal_talent_geometry(),finally=dev.off())
cairo_pdf('Figure_equal_talent_5050_crossings.pdf',width=figure_width,height=figure_height,family='Arial')
tryCatch(draw_equal_talent_geometry(),finally=dev.off())
