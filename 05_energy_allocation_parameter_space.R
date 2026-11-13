# Population energy after the first response and at Nash equilibrium.
# Replot from the supplied RDS; run 05_prepare_energy_allocation_inputs.R
# only to rebuild the numerical data (that step requires Rcpp and Rtools).
# The reference is random allocations with E1=E2=B/2 and separate M1,M2.
# Fixed-panel colours are the induced first-response energies, not the inputs.
# A-D classify the change from first response to equilibrium, using the
# magnitude and sign of E1-E2. The SAME unfilled polygons appear in both blocks.
for(pkg in c('sf','isoband'))if(!requireNamespace(pkg,quietly=TRUE))stop('Install ',pkg)
source('05_energy_allocation_quadrature.R')
input_file<-if(file.exists('05_energy_allocation_refined_inputs.rds'))
 '05_energy_allocation_refined_inputs.rds' else '05_energy_allocation_inputs.rds'
data<-readRDS(input_file)
B<-data$B;betas<-data$betas;phis<-data$phis
n_plot<-1001L
# Increase this one number to move the entire bottom energy scale lower.
gradient_legend_drop<-.06
show_examples<-TRUE
region_circle_cex<-2.7
region_letter_cex<-1.05
png_resolution<-600L
# Preserve the preceding figure's label centres while refining its boundaries.
label_positions<-if(file.exists('05_energy_allocation_label_positions.rds'))
 readRDS('05_energy_allocation_label_positions.rds') else NULL
plot_betas<-seq(min(betas),max(betas),length.out=n_plot)
energy_colours<-colorRampPalette(c('#08F','white','#E06'))(801)
bilinear<-function(z,x,y,xs=betas,ys=betas){
 i<-pmax(1L,pmin(length(xs)-1L,findInterval(x,xs)))
 j<-pmax(1L,pmin(length(ys)-1L,findInterval(y,ys)))
 tx<-(x-xs[i])/(xs[i+1]-xs[i]);ty<-(y-ys[j])/(ys[j+1]-ys[j])
 (1-tx)*(1-ty)*z[cbind(i,j)]+tx*(1-ty)*z[cbind(i+1,j)]+
 (1-tx)*ty*z[cbind(i,j+1)]+tx*ty*z[cbind(i+1,j+1)]
}
fine_grid<-expand.grid(x=plot_betas,y=plot_betas)
interpolate_field<-function(z)matrix(bilinear(z,fine_grid$x,fine_grid$y),n_plot)
region_geometry<-function(field){
 if(!any(field>1e-10))return(NULL)
 iso<-isoband::isobands(plot_betas,plot_betas,t(field),1e-10,max(field)+1)
 if(!length(iso[[1]]$x))return(NULL)
 sf::st_sfc(isoband::iso_to_sfg(iso)[[1]])
}
# Approximate pole of inaccessibility, as in the established coexistence figure.
visual_centre<-function(poly,avoid=matrix(numeric(),ncol=2),n=25,refine=3){
 bb<-sf::st_bbox(poly);xmin<-bb[1];ymin<-bb[2];xmax<-bb[3];ymax<-bb[4]
 best<-sf::st_coordinates(sf::st_point_on_surface(poly))[1,1:2]
 for(k in seq_len(refine)){
  gr<-expand.grid(x=seq(xmin,xmax,length.out=n),y=seq(ymin,ymax,length.out=n))
  pts<-sf::st_as_sf(gr,coords=c('x','y'));inside<-lengths(sf::st_within(pts,poly))>0
  if(any(inside)){
   pin<-pts[inside,,drop=FALSE];d<-as.numeric(sf::st_distance(pin,sf::st_boundary(poly)))
   if(nrow(avoid)){
    xy<-sf::st_coordinates(pin);clearance<-rep(Inf,nrow(xy))
    for(mark in seq_len(nrow(avoid)))clearance<-pmin(clearance,
      sqrt((xy[,1]-avoid[mark,1])^2+(xy[,2]-avoid[mark,2])^2))
    if(any(clearance>=.075))d[clearance<.075]<- -Inf
   }
   best<-sf::st_coordinates(pin[which.max(d),])[1,1:2]
  }
  dx<-(xmax-xmin)/(n-1);dy<-(ymax-ymin)/(n-1)
  xmin<-max(bb[1],best[1]-2*dx);xmax<-min(bb[3],best[1]+2*dx)
  ymin<-max(bb[2],best[2]-2*dy);ymax<-min(bb[4],best[2]+2*dy)
 }
 best
}
source('05_parameter_space_example_markers.R')
examples <- parameter_space_example_sets('beta',phis,range(betas))
example_cols <- parameter_space_marker_colours


cache_file<-'05_energy_allocation_plot_cache.rds'
signature<-list(input=unname(tools::md5sum(input_file)),n_plot=n_plot,version=7L,examples=if(show_examples)examples else NULL,
 label_positions=if(!is.null(label_positions))unname(tools::md5sum('05_energy_allocation_label_positions.rds'))else NULL)
cache<-if(file.exists(cache_file))readRDS(cache_file)else NULL
if(!is.null(cache)&&identical(cache$signature,signature))geometries<-cache$geometries else{
 geometries<-vector('list',9);dim(geometries)<-c(3,3)
 for(i in 1:3)for(j in 1:3){
  p<-data$panels[[i,j]];first<-interpolate_field(p$first);eq<-interpolate_field(p$equilibrium)
  category<-matrix(energy_feedback_class(as.vector(first),as.vector(eq),B),n_plot)
  # Continuous inequalities give vector zero-contours, avoiding raster-step
  # outlines. Positive field values define the corresponding category.
  d0<-first-B/2;d1<-eq-B/2
  class_fields<-list(
   pmax(pmin(d0,d1-d0),pmin(-d0,d0-d1)),
   pmax(pmin(d1,d0-d1),pmin(-d1,d1-d0)),
   pmax(pmin(-d1,d0+d1),pmin(d1,-d0-d1)),
   pmax(pmin(d0,-d0-d1),pmin(-d0,d0+d1)))
  regions<-lapply(class_fields,region_geometry)
  avoid<-matrix(numeric(),ncol=2)
  if(show_examples)for(ex in examples){
   keep<-which(ex[,3]==p$phi1 & ex[,4]==p$phi2)
   if(length(keep))avoid<-rbind(avoid,ex[keep,1:2,drop=FALSE])
  }
  labels<-if(!is.null(label_positions))label_positions[[i+(j-1)*3]] else list()
  if(is.null(label_positions))for(k in 1:4)if(!is.null(regions[[k]])){
   parts<-suppressWarnings(sf::st_cast(regions[[k]],'POLYGON'))
   for(part in seq_along(parts)){
    polygon<-parts[part]
    if(as.numeric(sf::st_area(polygon))<.003)next
    centre<-visual_centre(polygon,avoid)
    distance<-as.numeric(sf::st_distance(sf::st_sfc(sf::st_point(centre)),sf::st_boundary(polygon)))
    if(distance>.018)labels[[length(labels)+1L]]<-c(centre,k)
   }
  }
  # Retain existing centres; add A on both leader sides only off the phi diagonal.
  if(p$phi1!=p$phi2)for(leader in c(-1,1)){
   existing<-Filter(function(v)v[3]==1,labels)
   if(any(vapply(existing,function(v)leader*bilinear(first-B/2,v[1],v[2],plot_betas,plot_betas)>0,logical(1))))next
   side<-region_geometry(pmin(leader*d0,leader*(d1-d0)))
   if(is.null(side))next
   parts<-suppressWarnings(sf::st_cast(side,'POLYGON'))
   areas<-as.numeric(sf::st_area(parts));polygon<-parts[which.max(areas)]
   if(max(areas)<.003)next
   centre<-visual_centre(polygon,avoid)
   distance<-as.numeric(sf::st_distance(sf::st_sfc(sf::st_point(centre)),sf::st_boundary(polygon)))
   if(distance>.018)labels[[length(labels)+1L]]<-c(centre,1)
  }
  geometries[[i,j]]<-list(first=first,equilibrium=eq,regions=regions,labels=labels)
 }
 saveRDS(list(signature=signature,geometries=geometries),cache_file)
}
draw_energy_panel<-function(panel,geometry,r,c,mode){
 par(mar=c(1.5,1.5,.5,.5),mgp=c(1.35,.36,0),tcl=-.24,
     xaxs='i',yaxs='i',pty='s',xpd=FALSE)
 z<-if(mode=='fixed')geometry$first else geometry$equilibrium
 image(plot_betas,plot_betas,z,col=energy_colours,breaks=seq(0,B,length.out=802),
       axes=FALSE,xlab='',ylab='',useRaster=TRUE)
 # Dotted diagonal underneath region curves and labels.
 abline(0,1,lwd=.7,lty=3,col='#333333')
 # Equal-phi reflection preserves A-D, so solid diagonal edges are redundant.
 diagonal<-sf::st_sfc(sf::st_linestring(matrix(c(0,0,1,1),ncol=2,byrow=TRUE)))
 for(region in geometry$regions)if(!is.null(region)){
  boundary<-sf::st_boundary(region)
  if(panel$phi1==panel$phi2)
   boundary<-suppressWarnings(sf::st_difference(boundary,sf::st_buffer(diagonal,1e-7)))
  plot(boundary,add=TRUE,col='#222222',lwd=.9)
 }
 markers<-matrix(numeric(),ncol=2)
 if(show_examples)for(ex in examples){
  keep<-which(ex[,3]==panel$phi1 & ex[,4]==panel$phi2)
  if(length(keep))markers<-rbind(markers,ex[keep,1:2,drop=FALSE])
 }
 for(label in geometry$labels){
  # Keep a newly selected example cross clear of an existing region letter.
  if(nrow(markers) && any(sqrt((markers[,1]-label[1])^2+
                              (markers[,2]-label[2])^2)<.075)){
   parts<-suppressWarnings(sf::st_cast(geometry$regions[[label[3]]],'POLYGON'))
   anchor<-sf::st_sfc(sf::st_point(as.numeric(label[1:2])))
   nearest<-which.min(as.numeric(sf::st_distance(parts,anchor)))
   label[1:2]<-visual_centre(parts[nearest],markers)
  }
  points(label[1],label[2],pch=21,bg='black',col='black',cex=region_circle_cex,lwd=.9)
  text(label[1],label[2],LETTERS[label[3]],cex=region_letter_cex,font=2,col='white')
 }
 if(show_examples)for(set in names(examples)){
  ex<-examples[[set]];k<-which(ex[,3]==panel$phi1 & ex[,4]==panel$phi2)
  if(length(k)){
   points(ex[k,1],ex[k,2],pch=4,cex=1.3,lwd=3.5,col='#333333')
   points(ex[k,1],ex[k,2],pch=4,cex=1.3,lwd=2.5,col=example_cols[set])
  }
 }
 if(r==3){axis(1,at=c(.2,.4,.6,.8),cex.axis=1.02,lwd=.8,lwd.ticks=.8)
  mtext(expression(beta[1]),side=1,line=1.9,cex=1.18*.9,xpd=NA)}
 if(c==1){axis(2,las=0,at=c(.2,.4,.6,.8),cex.axis=1.02,lwd=.8,lwd.ticks=.8)
  mtext(expression(beta[2]),side=2,line=1.4,cex=1.18*.9,las=1)
  usr<-par('usr');text(usr[1]-.36,mean(usr[3:4]),
   bquote(varphi[2]==.(formatC(panel$phi2,format='f',digits=2))),srt=90,xpd=NA,cex=1.18*1.5)}
 if(r==1)mtext(bquote(varphi[1]==.(formatC(panel$phi1,format='f',digits=2))),side=3,line=.22,cex=1.18)
 box(lwd=.8)
}
draw_legend<-function(){
 par(mar=c(.03,.03,.03,.03),xpd=NA,pty='m');plot.new();width<-par('pin')[1]
 plot.window(xlim=c(0,width),ylim=c(0,1),xaxs='i',yaxs='i')
 row<-function(labels,symbols,y,cex=1.5){
  sw<-.16;gap<-.09;item_gap<-.30
  tw<-strwidth(labels,cex=cex,units='inches')
  available<-width-.3-length(labels)*(sw+gap)-(length(labels)-1)*item_gap
  if(sum(tw)>available){cex<-cex*available/sum(tw);tw<-strwidth(labels,cex=cex,units='inches')}
  x<-(width-sum(tw)-length(labels)*(sw+gap)-(length(labels)-1)*item_gap)/2
  for(i in seq_along(labels)){
   if(symbols[i]%in%LETTERS[1:4]){
    points(x+sw/2,y,pch=21,bg='black',col='black',cex=region_circle_cex,lwd=.9)
    text(x+sw/2,y,symbols[i],font=2,cex=region_letter_cex,col='white')
   }else{
    points(x+sw/2,y,pch=4,col='#333333',lwd=3.5,cex=1.3)
    points(x+sw/2,y,pch=4,col=example_cols[as.integer(symbols[i])],lwd=2.5,cex=1.3)
   }
   text(x+sw+gap,y,labels[i],adj=c(0,.5),cex=cex)
   x<-x+sw+gap+tw[i]+item_gap
  }
 }
 row(c('Same leading tournament; tilt grows','Same leading tournament; tilt shrinks',
       'New leader; smaller tilt','New leader; larger tilt'),c('A','B','C','D'),.56)
 if(show_examples)draw_parameter_example_legend(examples,y=.80)
 gy<-function(y)y-gradient_legend_drop
 left<-.18*width;right<-.82*width;xx<-seq(left,right,length.out=802)
 for(i in 1:801)rect(xx[i],gy(.21),xx[i+1],gy(.30),col=energy_colours[802-i],border=NA)
 energy<-c(0,B/4,B/2,3*B/4,B);at<-left+(right-left)*energy/B
 segments(at,gy(.21),at,gy(.18),lwd=.6);segments(at,gy(.30),at,gy(.33),lwd=.6)
 text(at,gy(.12),format(B-energy,trim=TRUE),cex=1.3)
 text(at,gy(.40),format(energy,trim=TRUE),cex=1.3)
 text(.13*width,gy(.12),expression(E[1]),cex=1.5)
 text(.13*width,gy(.40),expression(E[2]),cex=1.5)
 text(.84*width,gy(.255),bquote(E[1]+E[2]==.(B)),pos=4,cex=1.4)
}
draw_energy_allocation <- function() {
  layout(rbind(c(1,2,3,0,10,11,12),c(4,5,6,0,13,14,15),
               c(7,8,9,0,16,17,18),rep(19,7)),
         widths=c(1,1,1,.32,1,1,1),heights=c(1,1,1,.95))
  par(oma=c(.45,3.35,6.55,.45))
  for(mode in c('fixed','responsive')) for(r in 1:3) for(c in 1:3)
    draw_energy_panel(data$panels[[c,r]],geometries[[c,r]],r,c,mode)
  draw_legend()
  mtext(expression('Fixed '*E[1]==E[2]),side=3,outer=TRUE,at=.225,line=3.95,cex=1.55,font=2)
  mtext('Invested energy after the first talent-informed response',side=3,outer=TRUE,at=.225,line=2.25,cex=1.1)
  mtext(expression('Responsive '*E[t]),side=3,outer=TRUE,at=.775,line=3.95,cex=1.55,font=2)
  mtext('Invested energy at Nash equilibrium',side=3,outer=TRUE,at=.775,line=2.25,cex=1.1)
}
# Figure 06 uses 3500 x (3500*1900/3100) pixels at 300 dpi.
# Its outer margins consume 7 lines at .132 inches per line after layout().
# Add .39 row-height units solely for the taller legend; the panel rows,
# margins, headings, fonts and block gap retain Figure 06's dimensions.
figure_width <- 3500/300
reference_height <- figure_width*1900/3100
reference_row_height <- (reference_height-7*.132)/3.56
figure_height <- reference_height+(.95-.56)*reference_row_height
png('Figure_energy_allocation_fixed_vs_responsive.png',width=figure_width,height=figure_height,units='in',res=png_resolution)
tryCatch(draw_energy_allocation(),finally=dev.off())
cairo_pdf('Figure_energy_allocation_fixed_vs_responsive.pdf',width=figure_width,height=figure_height,family='Arial')
tryCatch(draw_energy_allocation(),finally=dev.off())

if(interactive())draw_energy_allocation()
