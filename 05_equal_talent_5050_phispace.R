# Phi-plane equivalent of Figure_equal_talent_5050_crossings.
# Axes: phi1, phi2 in [0.10,3.10]; panel betas: 0.10,0.50,0.90.
# Dot overlays identify two crossings; hues and saturation describe the LAST.
# Coloured hatching tests global optimality at that last stationary crossing.
# Uses the validated Script 06 continuous-population equilibrium ratios.
# Run from the scripts folder. First build needs Rcpp/Rtools; cached redraws do not.
for(pkg in c('sf','isoband')) if(!requireNamespace(pkg,quietly=TRUE)) stop('Install ',pkg)
source('05_random_allocation_reference.R')
B <- 50
n_plot <- 501L
gradient_legend_drop <- .06
input_file <- '06_coexistence_phi_inputs.rds'
cache_file <- '05_equal_talent_5050_phispace_geometry_cache.rds'
force_recompute <- identical(Sys.getenv('RECOMPUTE_5050_PHISPACE'),'1')
inputs <- readRDS(input_file)
phis <- inputs$axis;panel_betas <- inputs$panel_vals
stopifnot(inputs$B==B,identical(panel_betas,c(.1,.5,.9)))
plot_phis <- seq(min(phis),max(phis),length.out=n_plot)
colours <- c(colorRampPalette(c('#315A9B','white','#B54343'))(801),'#E06','#08F','white')
signature <- list(input=unname(tools::md5sum(input_file)),n_plot=n_plot,version=3L,
  geometry=unname(tools::md5sum('05_equal_talent_phi_geometry.cpp')),
  objective=unname(tools::md5sum('05_energy_allocation_population_phispace.cpp')))
bilinear <- function(z,x,y,xs=phis,ys=phis) {
 i<-pmax(1L,pmin(length(xs)-1L,findInterval(x,xs)))
 j<-pmax(1L,pmin(length(ys)-1L,findInterval(y,ys)))
 tx<-(x-xs[i])/(xs[i+1]-xs[i]);ty<-(y-ys[j])/(ys[j+1]-ys[j])
 (1-tx)*(1-ty)*z[cbind(i,j)]+tx*(1-ty)*z[cbind(i+1,j)]+
 (1-tx)*ty*z[cbind(i,j+1)]+tx*ty*z[cbind(i+1,j+1)]
}
# Reconstruct interpolated values from directly evaluated population nodes.
# Interpolating raw log(rho) on earlier meshes biases small-phi ratios;
# interpolate log(rho/rho_zero) instead, preserving every direct evaluation.
responsive_low_field <- function(panel,b1,b2) {
 g<-expand.grid(phi1=phis,phi2=phis)
 exact<-matrix(log(as.vector(panel$rho))+log(g$phi1*b1/(g$phi2*b2))+
                 (b1-b2)*log(1+B/2),length(phis))
 base<-seq(1L,length(phis),by=4L)
 stopifnot(all(panel$direct[base,base]))
 if(b1==.1 && b2==.1) {
   # This panel is exceptionally smooth: cubic interpolation of its direct
   # base mesh agrees with all additional direct checks to < 8e-8 in log ratio.
   xs<-phis[base];mat<-exact[base,base]
   tmp<-vapply(seq_along(xs),function(j)splinefun(xs,mat[,j])(phis),numeric(length(phis)))
   return(t(vapply(seq_len(nrow(tmp)),function(i)splinefun(xs,tmp[i,])(phis),numeric(length(phis)))))
 }
 ids<-base;field<-exact[ids,ids]
 for(step in c(2L,1L)) {
  next_ids<-seq(1L,length(phis),by=step);xy<-expand.grid(x=phis[next_ids],y=phis[next_ids])
  pred<-matrix(bilinear(field,xy$x,xy$y,phis[ids],phis[ids]),length(next_ids))
  direct<-panel$direct[next_ids,next_ids];truth<-exact[next_ids,next_ids]
  pred[direct]<-truth[direct];field<-pred;ids<-next_ids
 }
 field
}

talent_sdlog <- sqrt(log((1+sqrt(2))/2))
talent_sum_cdf <- function(total) integrate(function(a)
 dlnorm(a,0,talent_sdlog)*plnorm(total-a,0,talent_sdlog),0,total,rel.tol=1e-11)$value
talent_sum_cutoff <- uniroot(function(total)talent_sum_cdf(total)-.5,c(2,5),tol=1e-10)$root
if(file.exists(cache_file) && !force_recompute && identical(readRDS(cache_file)$signature,signature)) {
 panels <- readRDS(cache_file)$panels
} else {
 source('05_equal_talent_phi_geometry.R')
 moments <- sapply(panel_betas,function(b)vapply(phis,function(p)
   random_allocation_reference(b,b,p,p,B)$M1,1))
 fine_grid <- expand.grid(phi1=plot_phis,phi2=plot_phis)
 coarse_grid <- expand.grid(phi1=phis,phi2=phis)
 panels <- list()
 for(mode in c('fixed','responsive'))for(row in 1:3)for(column in 1:3){
  b1<-panel_betas[column];b2<-panel_betas[row]
  message(mode,': beta1=',b1,', beta2=',b2)
  coarse_offset <- log(coarse_grid$phi1*b1/(coarse_grid$phi2*b2))+(b1-b2)*log(1+B/2)
  fine_offset <- log(fine_grid$phi1*b1/(fine_grid$phi2*b2))+(b1-b2)*log(1+B/2)
  if(mode=='fixed') {
   cr <- outer(moments[,column],moments[,row],function(a,b)b/a)
   low_field <- matrix(log(as.vector(cr))+coarse_offset,length(phis))
  } else {
   low_field <- responsive_low_field(inputs$panels[[column,row]],b1,b2)
   cr <- matrix(exp(as.vector(low_field)-coarse_offset),length(phis))
  }
  ratio <- matrix(exp(bilinear(low_field,fine_grid$phi1,fine_grid$phi2)-fine_offset),n_plot)
  if(column==row){diag(ratio)<-1;diag(cr)<-1}
  coarse <- equal_talent_phi_geometry_cpp(coarse_grid$phi1,coarse_grid$phi2,as.vector(cr),b1,b2,B,TRUE)
  fine <- equal_talent_phi_geometry_cpp(fine_grid$phi1,fine_grid$phi2,as.vector(ratio),b1,b2,B,FALSE)
  colour_id <- rep(804L,nrow(fine));colour_id[fine$type==1]<-802L;colour_id[fine$type==2]<-803L
  w <- 1/(1+fine$alpha_cross)
  k<-which(fine$type==3);colour_id[k]<-401L+as.integer(round(400*w[k]))
  k<-which(fine$type==4);colour_id[k]<-401L-as.integer(round(400*w[k]))
  a<-talent_sum_cutoff/2;L<-1+B/2
  score <- log(as.vector(ratio))+log(fine_grid$phi1*b1/(fine_grid$phi2*b2))+
   (b1-b2)*log(L)+(fine_grid$phi1-1)*log1p(a*expm1(b1*log(L)))-
   (fine_grid$phi2-1)*log1p(a*expm1(b2*log(L)))
  panels[[paste(mode,row,column,sep='_')]] <- list(mode=mode,row=row,column=column,
   beta1=b1,beta2=b2,colours=matrix(colour_id,n_plot),
   always_above=matrix(fine$type==1,n_plot),always_below=matrix(fine$type==2,n_plot),
   two_crossings=matrix(fine$n_crossings==2L,n_plot),median_score=matrix(score,n_plot),
   coarse=coarse,bad=matrix(!is.na(coarse$global_equal_split)&!coarse$global_equal_split,length(phis)))
 }
 saveRDS(list(signature=signature,panels=panels),cache_file)
 tab <- do.call(rbind,lapply(panels,function(p){z<-p$coarse;z$mode<-p$mode;z$beta1<-p$beta1;z$beta2<-p$beta2;z$sum_first<-2*z$alpha_first;z$sum_cross<-2*z$alpha_cross;z}))
 con<-gzfile('equal_talent_5050_crossings_phispace.csv.gz','wt');write.csv(tab,con,row.names=FALSE);close(con)
}

band_geometry <- function(mask,x=plot_phis) {
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
  geom <- band_geometry(mask,phis)
  if(is.null(geom)) return(invisible(NULL))
  crossing <- band_geometry(crossing_mask)
  if(is.null(crossing)) return(invisible(NULL))
  geom <- suppressWarnings(sf::st_intersection(geom,crossing))
  offsets <- seq(-6,6,by=.055*3/.98)
  hatch <- sf::st_sfc(sf::st_multilinestring(lapply(offsets,function(d)
    rbind(c(-3,-3+d),c(7,7+d)))))
  clipped <- suppressWarnings(sf::st_intersection(hatch,geom))
  if(length(clipped)) plot(clipped,add=TRUE,col=adjustcolor(colour,.8),lwd=.85)
}
# Colour the hatch by the better specialist at the stationary crossing.
draw_hatching <- function(panel) {
  z <- panel$coarse
  u1 <- z$ratio*expm1(z$phi1*log1p(z$alpha_cross*expm1(panel$beta1*log1p(B))))
  u2 <- expm1(z$phi2*log1p(z$alpha_cross*expm1(panel$beta2*log1p(B))))
  winner1 <- matrix(!is.na(u1)&!is.na(u2)&u1>=u2,length(phis))
  draw_hatching_one(panel$bad & winner1,panel$colours<=801L,'#E06')
  draw_hatching_one(panel$bad & !winner1,panel$colours<=801L,'#08F')
}

# Trace crossings at the 50th percentile of total reference talent.
# Evaluate the marginal-return equation on the numerical input grid, then
# interpolate its zero contour rather than contouring clipped display colours.
draw_talent_sum_contour <- function(panel) {
  z <- panel$median_score
  if(panel$beta1==panel$beta2) {
    # Remove the coincident phi1=phi2 diagonal before finding nontrivial roots.
    delta <- outer(plot_phis,plot_phis,'-')
    z <- z/delta
    for(i in seq_along(plot_phis)) {
      k <- c(max(1L,i-1L),min(length(plot_phis),i+1L));k<-k[k!=i]
      z[i,i] <- mean(z[i,k])
    }
  }
  segments <- contourLines(plot_phis,plot_phis,z,levels=0)
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
examples <- parameter_space_example_sets('phi',panel_betas,range(phis))
example_cols <- parameter_space_marker_colours

draw_panel <- function(panel) {
  par(mar=c(1.5,1.5,.5,.5),mgp=c(1.35,.36,0),tcl=-.24,
      xaxs='i',yaxs='i',pty='s',xpd=FALSE)
  image(plot_phis,plot_phis,panel$colours,col=colours,
        breaks=seq(.5,804.5,1),axes=FALSE,xlab='',ylab='',useRaster=TRUE)
  draw_solid_region(panel$always_above,'#E06')
  draw_solid_region(panel$always_below,'#08F')
  if(any(panel$two_crossings)) {
    ids <- seq(1,n_plot,by=18L)
    xy <- expand.grid(i=ids,j=ids)
    xy <- xy[panel$two_crossings[cbind(xy$i,xy$j)],,drop=FALSE]
    if(nrow(xy)) points(plot_phis[xy$i],plot_phis[xy$j],pch=20,cex=.35,col='#333333')
  }
  draw_hatching(panel)
  draw_talent_sum_contour(panel)
  for(set in names(examples)) {
    z <- examples[[set]]
    k <- which(z[,3]==panel$beta1 & z[,4]==panel$beta2)
    if(length(k)) {
      points(z[k,1],z[k,2],pch=4,cex=1.3,lwd=3.5,col='#333333')
      points(z[k,1],z[k,2],pch=4,cex=1.3,lwd=2.5,col=example_cols[set])
    }
  }
  if(panel$row==3) {
    axis(1,at=c(.5,1,1.5,2,2.5,3),cex.axis=1.02,lwd=.8,lwd.ticks=.8,gap.axis=.35)
    mtext(expression(varphi[1]),side=1,line=1.9,cex=1.18*.9,xpd=NA)
  }
  if(panel$column==1) {
    axis(2,las=0,at=c(.5,1,1.5,2,2.5,3),cex.axis=1.02,lwd=.8,lwd.ticks=.8,gap.axis=.35)
    mtext(expression(varphi[2]),side=2,line=1.4,cex=1.18*.9,las=1)
    usr <- par('usr')
    text(usr[1]-.36*3/.98,mean(usr[3:4]),
         bquote(beta[2]==.(formatC(panel$beta2,format='f',digits=2))),
         srt=90,xpd=NA,cex=1.18*1.5)
  }
  if(panel$row==1) mtext(bquote(beta[1]==.(formatC(panel$beta1,format='f',digits=2))),
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
      if(kinds[i]=='dots') {
        # Physical square, filled with a regular dot pattern.
        half_height <- symbol_width/(2*par('pin')[2])
        rect(x,y-half_height,x+symbol_width,y+half_height,
             col='white',border=cols[i],lwd=.65)
        dots <- expand.grid(x=x+symbol_width*c(.2,.5,.8),
                            y=y+half_height*c(-.6,0,.6))
        points(dots$x,dots$y,pch=20,cex=.4,col=cols[i])
      }
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
               '50:50 not optimal at last crossing','Two crossings'),
             c('box','box','hatch','dots'),c('#E06','#08F','#333333','#333333'),.60)
  draw_parameter_example_legend(examples,y=.80,
    line_label='50th percentile of talent sum')
  gradient_y <- function(y) y-gradient_legend_drop
  # Red precedes blue, matching the tournament order of the solid colours.
  bars <- list(list(left=.12*width,right=.45*width,side=1,
                    label='After last crossing: leaning 1'),
               list(left=.55*width,right=.88*width,side=-1,
                    label='After last crossing: leaning 2'))
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
png('Figure_equal_talent_5050_crossings_phispace.png',width=figure_width,height=figure_height,units='in',res=600)
tryCatch(draw_equal_talent_geometry(),finally=dev.off())
cairo_pdf('Figure_equal_talent_5050_crossings_phispace.pdf',width=figure_width,height=figure_height,family='Arial')
tryCatch(draw_equal_talent_geometry(),finally=dev.off())
