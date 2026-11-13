# Parameter figures for first local/global interventions, not campaign Nash maps.
# Filled bands are vector polygons; no raster blur or arbitrary smoothing is used.
assemble_campaign_fields<-function(require_refined=TRUE){
 if(file.exists('10_campaign_precision_results/complete.rds')){
  if(!exists('assemble_campaign_precision_fields',mode='function'))source('10_campaign_precision_maps.R')
  return(assemble_campaign_precision_fields(.1,save_output=TRUE))
 }

 axis<-seq(.1,.9,length.out=13);phis<-c(.66,1,1.5)
 fields<-C1<-C2<-setNames(lapply(1:2,function(k)array(NA_real_,c(13,13,3,3))),c('local','global'))
 panels<-setNames(lapply(1:2,function(k)matrix(vector('list',9),3,3)),names(fields))
 for(family in names(fields))for(c in 1:3)for(r in 1:3){
  for(i in 1:13)for(j in 1:13){
   cc<-c;rr<-r;ii<-i;jj<-j;swap<-FALSE
   if(c<r){cc<-r;rr<-c;ii<-j;jj<-i;swap<-!swap}
   if(cc==rr&&jj<ii){tmp<-ii;ii<-jj;jj<-tmp;swap<-!swap}
   file<-sprintf('10_local_global_results/parameter_cells/%s_p%.2f_%.2f_b%.6f_%.6f.rds',family,phis[cc],phis[rr],axis[ii],axis[jj])
   if(!file.exists(file))stop('Missing parameter cell: ',file)
   z<-readRDS(file);if(is.list(z$result))stop('Unresolved parameter cell: ',file)
   v<-unname(z$result[c('campaign1','campaign2')]);if(swap)v<-rev(v)
   if(c==r&&i==j)v[]<-mean(v)
   C1[[family]][i,j,c,r]<-v[1];C2[[family]][i,j,c,r]<-v[2]
   fields[[family]][i,j,c,r]<-sum(v)
  }
  mid<-(head(axis,-1)+tail(axis,-1))/2;mid1<-mid2<-matrix(NA_real_,12,12)
  for(i in 10:12)for(j in 10:12){
   cc<-c;rr<-r;ii<-i;jj<-j;swap<-FALSE
   if(c<r){cc<-r;rr<-c;ii<-j;jj<-i;swap<-!swap}
   if(cc==rr&&jj<ii){tmp<-ii;ii<-jj;jj<-tmp;swap<-!swap}
   file<-sprintf('10_local_global_results/parameter_center_cells/%s_p%.2f_%.2f_b%.6f_%.6f.rds',family,phis[cc],phis[rr],mid[ii],mid[jj])
   if(!file.exists(file)){if(require_refined)stop('Additional parameter calculation missing: ',file);next}
   z<-readRDS(file);v<-unname(z$result[c('campaign1','campaign2')]);if(swap)v<-rev(v)
   if(c==r&&i==j)v[]<-mean(v)
   mid1[i,j]<-v[1];mid2[i,j]<-v[2]
  }
  panels[[family]][[c,r]]<-list(axis=axis,C1=C1[[family]][,,c,r],C2=C2[[family]][,,c,r],center_C1=mid1,center_C2=mid2,actual_grid=13L)
 }
 refined<-'10_local_global_results/linear_beta25.rds'
 if(file.exists(refined)){z<-readRDS(refined);stopifnot(isTRUE(z$complete));panels$global[[2,2]]<-list(axis=z$axis,C1=z$C1,C2=z$C2,actual_grid=length(z$axis))
 }else if(require_refined)stop('Denser linear-spoils grid is not prepared.')
 d<-list(axis=axis,phis=phis,fields=fields,campaign1=C1,campaign2=C2,panels=panels,
  range=range(unlist(fields)),B=50,efficiency=.1,local_N=30)
 actual<-unlist(lapply(panels,function(mm)lapply(mm,function(p)c(as.vector(p$C1+p$C2),as.vector(p$center_C1+p$center_C2)))))
 d$range<-range(actual,na.rm=TRUE);d$display_max<-ceiling(max(actual,na.rm=TRUE)/.05)*.05
 saveRDS(d,'10_local_global_results/parameter_levels.rds');d
}
campaign_bilinear<-function(axis,z,fine,centers=NULL){
 v<-vapply(seq_along(axis),function(j)approx(axis,z[,j],xout=fine)$y,numeric(length(fine)))
 out<-t(vapply(seq_along(fine),function(i)approx(axis,v[i,],xout=fine)$y,numeric(length(fine))))
 if(!is.null(centers))for(i in 1:(length(axis)-1))for(j in 1:(length(axis)-1)){
  if(!is.finite(centers[i,j]))next
  ix<-which(fine>=axis[i]&fine<=axis[i+1]);iy<-which(fine>=axis[j]&fine<=axis[j+1]);if(!length(ix)||!length(iy))next
  u<-outer((fine[ix]-axis[i])/(axis[i+1]-axis[i]),rep(1,length(iy)))
  w<-outer(rep(1,length(ix)),(fine[iy]-axis[j])/(axis[j+1]-axis[j]));m<-out[ix,iy,drop=FALSE];cc<-centers[i,j]
  bottom<-w<=pmin(u,1-u);right<-u>=pmax(w,1-w);top<-w>=pmax(u,1-u);left<-u<=pmin(w,1-w)
  val<-(1-u-w)*z[i,j]+(u-w)*z[i+1,j]+2*w*cc;m[bottom]<-val[bottom]
  val<-(u-w)*z[i+1,j]+(u+w-1)*z[i+1,j+1]+2*(1-u)*cc;m[right]<-val[right]
  val<-(w-u)*z[i,j+1]+(u+w-1)*z[i+1,j+1]+2*(1-w)*cc;m[top]<-val[top]
  val<-(1-u-w)*z[i,j]+(w-u)*z[i,j+1]+2*u*cc;m[left]<-val[left]
  out[ix,iy]<-m
 }
 out
}
campaign_direction_rgb<-function(c1,c2,max_total){
 total<-c1+c2;direction<-ifelse(total>1e-12,(c1-c2)/total,0)
 hue<-colorRamp(c('#004E85','#482769','#8F004D'))((pmax(-1,pmin(1,as.vector(direction)))+1)/2)/255
 strength<-pmin(1,pmax(0,as.vector(total)/max_total))
 col<-rgb(1-strength+strength*hue[,1],1-strength+strength*hue[,2],1-strength+strength*hue[,3])
 matrix(col,nrow(c1),ncol(c1))
}
# Align the fine mesh with every calculated vertex and cell centre.
campaign_rendering_axis<-function(d,minimum=600L){
 denominator<-unique(unlist(lapply(d$panels,function(block)lapply(block,function(p)
  (length(p$axis)-1L)*(if(is.null(p$center_C1))1L else 2L)))))
 gcd<-function(a,b){while(b!=0){z<-b;b<-a%%b;a<-z};a}
 multiple<-Reduce(function(a,b)a*b/gcd(a,b),denominator)
 seq(.1,.9,length.out=1L+multiple*ceiling(minimum/multiple))
}
# Figure S12 (05_energy_allocation_phispace.R) uses a dark outline
# under each coloured cross. Only examples on the plotted phi slices appear here.
campaign_example_markers <- function() {
  symmetric <- function(beta, phi) cbind(beta,beta,phi,phi)
  asymmetric <- rbind(c(.6,.6,1,.66),c(.8,.6,1.5,1),
                      c(.8,.6,1.5,1.5),c(.6,.4,.66,1))
  # Figure 13 takes precedence at the shared coordinates.
  list(
    local=list('Figure 13'=symmetric(.4,1),
               'Figure S46'=symmetric(.4,c(.66,1.5)),
               'Figure S47'=symmetric(.8,c(.66,1,1.5)),
               'Figure S48'=asymmetric),
    global=list('Figure 13'=rbind(c(.4,.4,1,1),c(.8,.8,1,1),c(.8,.4,1.5,.66)),
                'Figure S49'=symmetric(.4,c(.66,1.5)),
                'Figure S50'=symmetric(.8,c(.66,1.5)),
                'Figure S51'=asymmetric)
  )
}

campaign_example_cross <- function(x,y,colour) {
  points(x,y,pch=4,cex=1.3,lwd=3.5,col='#333333')
  points(x,y,pch=4,cex=1.3,lwd=2.5,col=colour)
}

plot_campaign_levels<-function(png=TRUE,pdf=TRUE,palette='viridis',
 interpolation=c('bands','none','continuous'),show_nodes=FALSE,
 file='Figure_convincing_local_global_parameter_levels',direction=FALSE,band_width=.005,
 data=readRDS('10_local_global_results/parameter_levels.rds')){
 interpolation<-match.arg(interpolation);d<-data;phis<-d$phis
 examples<-if(isTRUE(all.equal(d$efficiency,.10)))campaign_example_markers()else NULL
 example_colours<-c(
  "Figure 13"="#FFFFFF",
  "Figure S46"="#F2BB3B", "Figure S47"="#41CC9A", "Figure S48"="#DAB6EF",
  "Figure S49"="#FF755C", "Figure S50"="#65DFFF", "Figure S51"="#CDEB52"
 )
 fine<-campaign_rendering_axis(d)
 breaks<-seq(0,d$display_max+band_width/2,by=band_width)
 cols<-hcl.colors(length(breaks)-1,palette)
 draw<-function(){
  layout(rbind(c(1,2,3,0,10,11,12),c(4,5,6,0,13,14,15),c(7,8,9,0,16,17,18),rep(19,7)),
   widths=c(1,1,1,.32,1,1,1),heights=c(1,1,1,.8))
  par(oma=c(.45,3.35,6.55,.45))
  for(mode in 1:2)for(r in 1:3)for(c in 1:3){
   p<-d$panels[[mode]][[c,r]];grid<-if(interpolation=='none')p$axis else fine
   c1<-if(interpolation=='none')p$C1 else campaign_bilinear(p$axis,p$C1,fine,p$center_C1)
   c2<-if(interpolation=='none')p$C2 else campaign_bilinear(p$axis,p$C2,fine,p$center_C2)
   z<-c1+c2
   par(mar=c(1.5,1.5,.5,.5),mgp=c(1.35,.36,0),tcl=-.24,xaxs='i',yaxs='i',pty='s',xpd=FALSE)
   plot(NA,xlim=range(grid),ylim=range(grid),axes=FALSE,xlab='',ylab='',asp=1)
   if(direction){cm<-campaign_direction_rgb(c1,c2,d$display_max)
    rasterImage(as.raster(t(cm)[nrow(cm):1,,drop=FALSE]),min(grid),min(grid),max(grid),max(grid),interpolate=FALSE)
   }else if(interpolation=='bands'){
    bands<-isoband::isobands(grid,grid,t(z),head(breaks,-1),tail(breaks,-1))
    shape<-sf::st_sfc(isoband::iso_to_sfg(bands));keep<-!sf::st_is_empty(shape)
    if(any(keep))plot(shape[keep],col=cols[keep],border=NA,add=TRUE)
   }else image(grid,grid,z,breaks=breaks,col=cols,add=TRUE,axes=FALSE,useRaster=TRUE)
   abline(0,1,lty=3,lwd=.7,col='#FFFFFF70')
   if(show_nodes){points(expand.grid(p$axis,p$axis),pch=16,cex=.25,col='#FFFFFF88')
    if(!is.null(p$center_C1)){ij<-which(is.finite(p$center_C1),arr.ind=TRUE);mid<-(head(p$axis,-1)+tail(p$axis,-1))/2
     if(nrow(ij))points(mid[ij[,1]],mid[ij[,2]],pch=3,cex=.45,col='white')}}
   if(r==3){axis(1,at=c(.2,.4,.6,.8),cex.axis=1.02);mtext(expression(beta[1]),1,line=1.9,cex=1.062)}
   if(c==1){axis(2,at=c(.2,.4,.6,.8),cex.axis=1.02);mtext(expression(beta[2]),2,line=1.4,cex=1.062,las=1)
    text(par('usr')[1]-.27,mean(par('usr')[3:4]),bquote(varphi[2]==.(sprintf('%.2f',phis[r]))),srt=90,xpd=NA,cex=1.77)}
   if(r==1)mtext(bquote(varphi[1]==.(sprintf('%.2f',phis[c]))),3,line=.22,cex=1.18)
   for(i in seq_along(examples[[mode]])){
    marks<-examples[[mode]][[i]]
    here<-marks[,3]==phis[c] & marks[,4]==phis[r]
    if(any(here))campaign_example_cross(marks[here,1],marks[here,2],example_colours[names(examples[[mode]])[i]])
   }
   box(lwd=.8)
  }
  par(mar=c(.1,.1,.1,.1),pty='m',xpd=NA);plot.new();plot.window(xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i')
  if(!direction){
   for(k in seq_along(cols))rect(.24+.52*breaks[k]/d$display_max,.30,.24+.52*breaks[k+1]/d$display_max,.55,col=cols[k],border=NA)
   rect(.24,.30,.76,.55,lwd=.8);at<-seq(0,d$display_max,length.out=6);xx<-.24+.52*at/d$display_max
   segments(xx,.30,xx,.24);text(xx,.14,sprintf('%.0f',100*at),cex=1.1)
   text(.5,.75,'Mean potential-actor convincing time (% of total budget)',cex=1.45)
  }else{
   dd<-seq(-1,1,length.out=301);tt<-seq(0,d$display_max,length.out=121)
   c1<-outer((1+dd)/2,tt);c2<-outer((1-dd)/2,tt);cm<-campaign_direction_rgb(c1,c2,d$display_max)
   rasterImage(as.raster(t(cm)[ncol(cm):1,,drop=FALSE]),.27,.19,.73,.80,interpolate=FALSE);rect(.27,.19,.73,.80)
   xx<-c(.27,.5,.73);segments(xx,.19,xx,.15)
   text(xx,.04,c('Only tournament 2','Equal campaigns','Only tournament 1'),cex=1.1)
   at<-seq(0,d$display_max,length.out=4);yy<-.19+.61*at/d$display_max
   segments(.73,yy,.742,yy);text(.75,yy,sprintf('%.1f%%',100*at),adj=0,cex=1.05)
   text(.21,.495,'Total convincing',srt=90,cex=1.25)
  }
  # Keep the existing colour legends in place. The example keys use the
  # free space to their left and right.
  example_rows<-seq(.75,.25,length.out=length(examples$local))
  for(mode in seq_along(examples)){
   labels<-names(examples[[mode]])
   keep<-if(mode==1L)seq_along(labels)else which(labels!='Figure 13')
   x<-c(.06,.84)[mode]
   y<-example_rows[seq_along(keep)]
   campaign_example_cross(rep(x,length(y)),y,example_colours[labels[keep]])
   text(x+.017,y,labels[keep],adj=c(0,.5),cex=1.25)
  }
  mtext('Local tournaments',3,outer=TRUE,at=.225,line=3.95,cex=1.55)
  mtext('Personal influence; N = 30',3,outer=TRUE,at=.225,line=2.25,cex=1.1)
  mtext('Global tournaments',3,outer=TRUE,at=.775,line=3.95,cex=1.55)
  mtext('Single global intervention',3,outer=TRUE,at=.775,line=2.25,cex=1.1)
  mtext(bquote(epsilon==.(sprintf('%.2f',d$efficiency))~~B==.(d$B)),3,outer=TRUE,at=.5,line=4.3,cex=1.55)
 }
 if(png){grDevices::png(paste0(file,'.png'),width=3500/300,height=7.6,units='in',res=600);tryCatch(draw(),finally=dev.off())}
 if(pdf){cairo_pdf(paste0(file,'.pdf'),width=3500/300,height=7.6,family='Arial');tryCatch(draw(),finally=dev.off())}
 if(interactive())draw();invisible(d)
}
plot_campaign_direction<-function(png=TRUE,pdf=TRUE,file='Figure_convincing_campaign_direction',...)
 plot_campaign_levels(png=png,pdf=pdf,file=file,direction=TRUE,...)
if(sys.nframe()==0L)plot_campaign_levels()
