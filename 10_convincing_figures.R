# source('10_convincing_figures.R'); draw_convincing_set('symmetric_beta04')
# Replotting reads saved calculations. PNG/PDF are optional, FALSE by default.
source('10_convincing_model.R')
convincing_hdr<-function(){s<-boredom_slog();th<-seq(0,2*pi,length.out=1001)
 lapply(qchisq(c(.5,.9,.99),2,ncp=2*s*s),function(q)exp(-s*s+s*sqrt(q)*rbind(cos(th),sin(th))))}
convincing_polygon<-function(a,f,col){
 iso<-isoband::isobands(a,a,t(f),0,max(1,max(f)+1))
 if(length(iso[[1]]$x))plot(sf::st_sfc(isoband::iso_to_sfg(iso)[[1]]),add=TRUE,col=col,border=NA)
}
convincing_strategy_legend<-function(cols){
 # Approved S14 legend: same colours, wording, ticks and endpoint labels.
 par(mar=c(2.5,.55,.5,1.60),pty='m',xpd=NA);plot.new()
 plot.window(xlim=c(-.08,1),ylim=c(0,1),xaxs='i',yaxs='i')
 zz<-seq(0,1,length.out=801)
 for(k in 1:800)rect(.34,zz[k],.64,zz[k+1],col=cols[k],border=NA)
 rect(.34,0,.64,1,border='black',lwd=.8);wd<-.006
 rect(.34,-wd,.64,wd,col='#08F',lwd=.8,border=1)
 rect(.34,.5-wd,.64,.5+wd,col='#FF0',lwd=.8,border=1)
 rect(.34,1-wd,.64,1+wd,col='#E06',lwd=.8,border=1)
 ticks<-seq(0,1,.25);segments(.64,ticks,.75,ticks,lwd=.8)
 text(.88,ticks,sprintf('%.2f',ticks),adj=c(0,.5),cex=1.2)
 text(.10,.5,'Relative productive time in tournament 1',srt=90,cex=1.4)
 text(.49,1.035,'1maxxers',cex=1.5,col='#E06')
 text(.49,-.035,'2maxxers',cex=1.5,col='#08F')
}
convincing_scale<-function(cols,limits,label){
 par(mar=c(2.75,.55,.5,1.6),pty='m',xpd=NA);plot.new()
 plot.window(xlim=c(-.08,1),ylim=limits,xaxs='i',yaxs='i')
 z<-seq(limits[1],limits[2],length.out=length(cols)+1)
 for(k in seq_along(cols))rect(.34,z[k],.64,z[k+1],col=cols[k],border=NA)
 rect(.34,limits[1],.64,limits[2],lwd=.8)
 ticks<-seq(limits[1],limits[2],length.out=5);segments(.64,ticks,.75,ticks,lwd=.8)
 text(.88,ticks,sprintf('%.2f',ticks),adj=c(0,.5),cex=1.2)
 text(.10,mean(limits),label,srt=90,cex=1.25)
}
convincing_equal_curves<-function(d,row){
 p<-d$p;a<-d$map$alpha;o<-d$map$allocation
 background<-convincing_background_function(p)
 q<-if(row==1)o[,5]/p$B else o[,1]/pmax(1e-8,o[,1]+o[,2])
 n<-length(a);mask<-matrix(q,n)>=.5
 ij<-which(mask[-1,,drop=FALSE]!=mask[-n,,drop=FALSE],arr.ind=TRUE)
 roots<-vector('list',n)
 if(nrow(ij)){
  lo<-a[ij[,1]];hi<-a[ij[,1]+1];yy<-a[ij[,2]];side<-mask[cbind(ij[,1]+1,ij[,2])]
  for(k in 1:12){mid<-(lo+hi)/2
   if(row==1)above<-finite_pool_alloc(p$baseline,'responsive',cbind(mid,yy))>=p$B/2 else{
    vv<-convincing_optimize(p,cbind(mid,yy));above<-vv[,1]>=vv[,2]}
   movehi<-above==side;hi[movehi]<-mid[movehi];lo[!movehi]<-mid[!movehi]
  }
  for(j in seq_len(n))roots[[j]]<-(lo+hi)[ij[,2]==j]/2
 }
 curves<-list();nr<-max(lengths(roots))
 if(nr)for(k in seq_len(nr))curves[[k]]<-list(x=vapply(roots,function(v)if(length(v)>=k)v[k]else NA_real_,0.),y=a)
 # A jump from one global maximum to another is not an optimal 50:50 split.
 # Check the equal-allocation candidate against the unconstrained global search.
 lapply(curves,function(cl){
  keep<-which(is.finite(cl$x));if(!length(keep))return(cl)
  xy<-cbind(cl$x[keep],cl$y[keep]);opt<-convincing_optimize(p,xy)
  eqvalue<-function(alpha,c){
   m<-background(c)
   x<-(p$B-abs(c))/2;w1<-boredom_w(x,alpha[1],p$beta[1],p$phi[1]);w2<-boredom_w(x,alpha[2],p$beta[2],p$phi[2])
   den1<-(p$N-1)*m[,3]+w1;den2<-(p$N-1)*m[,4]+w2
   ifelse(den1>0,((p$N-1)*m[,1]+x+pmax(c,0))*w1/pmax(den1,1e-300),0)+
    ifelse(den2>0,((p$N-1)*m[,2]+x+pmax(-c,0))*w2/pmax(den2,1e-300),0)
  }
  best<-vapply(seq_len(nrow(xy)),function(i){
   if(row==1)return(eqvalue(xy[i,],0))
   cc<-p$effort
   if(isTRUE(p$aligned))cc<-cc[cc*(2*opt[i,5]-p$B)>= -1e-7]
   v<-eqvalue(xy[i,],cc);ans<-max(v)
   kk<-which(v[2:(length(v)-1)]>=v[1:(length(v)-2)]&v[2:(length(v)-1)]>=v[3:length(v)])+1
   for(k in kk)ans<-max(ans,optimize(function(c)eqvalue(xy[i,],c),cc[c(k-1,k+1)],maximum=TRUE,tol=1e-7)$objective)
   ans
  },0.)
  valid<-best>=opt[,if(row==1)6 else 4]-1e-5*(1+p$B)
  cl$x[keep[!valid]]<-NA_real_;cl
 })
}
convincing_boundary_fields<-function(d,row){
 p<-d$p;a<-d$map$alpha;o<-d$map$allocation;n<-length(a)
 q<-matrix(if(row==1)o[,5]/p$B else o[,1]/pmax(1e-8,o[,1]+o[,2]),n)
 allocate<-if(row==1)function(actor)finite_pool_alloc(p$baseline,'responsive',actor)/p$B else function(actor){
  v<-convincing_optimize(p,actor);v[,1]/pmax(1e-8,v[,1]+v[,2])}
 fields<-lapply(c(1e-5,1-1e-5),function(target){
  mask<-q>=target;ij<-which(mask[-1,,drop=FALSE]!=mask[-n,,drop=FALSE],arr.ind=TRUE)
  roots<-vector('list',n)
  if(nrow(ij)){
   lo<-a[ij[,1]];hi<-a[ij[,1]+1];y<-a[ij[,2]];side<-mask[cbind(ij[,1]+1,ij[,2])]
   for(k in 1:10){mid<-(lo+hi)/2;above<-allocate(cbind(mid,y))>=target;movehi<-above==side;hi[movehi]<-mid[movehi];lo[!movehi]<-mid[!movehi]}
   for(j in seq_len(n))roots[[j]]<-(lo+hi)[ij[,2]==j]/2
  }
  field<-matrix(0,n,n)
  for(j in seq_len(n))field[,j]<-(if(length(roots[[j]]))apply(abs(outer(a,roots[[j]],'-')),1,min)else rep(diff(range(a)),n))*ifelse(mask[,j],1,-1)
  field
 })
 # Ordered scalar fields define disjoint specialist polygons on the same mesh.
 # In particular, two independent contour reconstructions cannot overlap at
 # a direct jump between the two pure strategies.
 list(pmax(fields[[1]],fields[[2]]),pmin(fields[[1]],fields[[2]]))
}
draw_convincing_set<-function(set='asymmetric',png=FALSE,pdf=FALSE,
 file=paste0('Figure_convincing_',set),data=readRDS(file.path('10_convincing_results',paste0(set,'.rds'))),
 sidebars=identical(set,'asymmetric'),geometry=NULL){
 stopifnot(requireNamespace('isoband',quietly=TRUE),requireNamespace('sf',quietly=TRUE))
 nc<-length(data);cols<-colorRampPalette(c('#315A9B','white','#B54343'))(801)
 efforts<-colorRampPalette(c('#08F','white','#E06'))(801);greys<-gray(seq(1,.06,length.out=801))
 cmax<-max(.05,ceiling(max(vapply(data,function(d)max(abs(d$map$allocation[,3]))/d$p$B,0.))*20)/20)
 rmax<-max(.05,ceiling(max(vapply(data,function(d)max(d$map$recipient_abs)/d$p$B,0.))*20)/20)
 if(is.null(geometry)){
  equal<-lapply(data,function(d)lapply(1:2,function(r)convincing_equal_curves(d,r)))
  fields<-lapply(data,function(d)lapply(1:2,function(r)convincing_boundary_fields(d,r)))
 }else{
  stopifnot(identical(geometry$alpha,lapply(data,function(d)d$map$alpha)))
  equal<-geometry$equal;fields<-geometry$fields
 }
 if((png||pdf)&&identical(file,paste0('Figure_convincing_',set))&&set%in%c('symmetric_beta04','symmetric_beta08','asymmetric'))
  saveRDS(list(fields=fields,equal=equal,alpha=lapply(data,function(d)d$map$alpha)),file.path('10_convincing_results',paste0(set,'_geometry.rds')))
 draw<-function(){
  lay<-matrix(0,4,nc+1);for(r in 1:4)lay[r,1:nc]<-(r-1)*nc+seq_len(nc)
  lay[1:2,nc+1]<-4*nc+1;lay[3,nc+1]<-4*nc+2;lay[4,nc+1]<-4*nc+3
  layout(lay,widths=c(rep(1,nc),.46));par(oma=c(2,3.8,3,0))
  labels<-c('Without convincing','With convincing','Convincing effort','Recipient response')
  for(r in 1:4)for(j in seq_len(nc)){
   d<-data[[j]];p<-d$p;a<-d$map$alpha;o<-d$map$allocation;n<-length(a)
   z<-switch(r,o[,5]/p$B,o[,1]/pmax(1e-8,o[,1]+o[,2]),o[,3]/p$B,d$map$recipient_abs/p$B)
   z<-matrix(z,n,n)
   par(mar=c(2.75,2.75,.5,.2),mgp=c(1.75,.48,0),tcl=-.23,pty='s',xpd=FALSE)
   plot(NA,xlim=range(a),ylim=range(a),xaxs='i',yaxs='i',axes=FALSE,xlab='',ylab='',asp=1)
   pal<-if(r<=2)cols else if(r==3)efforts else greys
   ix<-if(r<=2)1+round(800*z)else if(r==3)1+round(400*(1+z/cmax))else 1+round(800*z/rmax)
   cm<-matrix(pal[pmax(1,pmin(801,ix))],n,n)
   rasterImage(as.raster(t(cm)[n:1,,drop=FALSE]),min(a),min(a),max(a),max(a),interpolate=TRUE)
   if(r<=2){
    # Vector polygons cover pure specialists; no pixel outlines.
    convincing_polygon(a,-fields[[j]][[r]][[1]],'#08F');convincing_polygon(a,fields[[j]][[r]][[2]],'#E06')
   }
   abline(0,1,lty=2,lwd=.65)
   for(h in convincing_hdr())lines(h[1,],h[2,],col='#EEE',lwd=.8)
   if(r<=2){
    for(f in fields[[j]][[r]])for(cl in contourLines(a,a,f,levels=0))lines(cl$x,cl$y,lwd=.75,col='black')
    eql<-equal[[j]][[r]]
    for(cl in eql)lines(cl$x,cl$y,col='black',lwd=3.5)
    for(cl in eql)lines(cl$x,cl$y,col='#FF0',lwd=2.3)
   }
   box()
   if(sidebars&&r<=2){
    budget<-if(r==1)c(p$N*p$baseline$responsive$moment[1],p$N*(p$B-p$baseline$responsive$moment[1]),0)else p$group_budget
    ends<-c(0,cumsum(budget)/(p$N*p$B))*max(a)
    for(k in 1:3)rect(2.54,ends[k],2.69,ends[k+1],col=c('#E06','#08F','black')[k],border=NA,xpd=NA)
    rect(2.54,0,2.69,max(a),lwd=.7,xpd=NA)
   }
   if(r==4){axis(1,at=seq(0,2.5,.5),cex.axis=1.2);mtext(expression(alpha[1]~~'(talent in tournament 1)'),1,line=2.4,cex=.95)}
   if(j==1){axis(2,at=seq(0,2.5,.5),las=1,cex.axis=1.2);mtext(expression(alpha[2]~~'(talent in tournament 2)'),2,line=2.4,cex=.95)}
   if(r==1){
    mtext(LETTERS[j],3,line=1.8,at=0,adj=0,cex=1.25,font=2,xpd=NA)
    for(t in 1:2){at<-2.5*c(.48,.88)[t];cc<-c('#E06','#08F')[t]
     mtext(bquote(beta[.(t)]==.(sprintf('%.2f',p$beta[t]))),3,line=1.55,at=at,cex=1.12,col=cc)
     mtext(bquote(varphi[.(t)]==.(sprintf('%.2f',p$phi[t]))),3,line=.15,at=at,cex=1.12,col=cc)
    }
   }
  }
  for(r in 1:4)mtext(labels[r],2,outer=TRUE,at=1-(r-.5)/4,line=1.8,cex=1.25)
  convincing_strategy_legend(cols)
  convincing_scale(efforts,c(-cmax,cmax),'Signed convincing time / budget')
  convincing_scale(greys,c(0,rmax),'Mean absolute allocation change / budget')
 }
 width<-(3580/300)*(nc+.46)/4.46;height<-1740/300*2
 if(png){grDevices::png(paste0(file,'.png'),width=width,height=height,units='in',res=600,type='cairo');tryCatch(draw(),finally=dev.off())}
 if(pdf){cairo_pdf(paste0(file,'.pdf'),width=width,height=height,family='Arial');tryCatch(draw(),finally=dev.off())}
 if(interactive()||dev.cur()!=1)draw();invisible(data)
}
