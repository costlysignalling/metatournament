# One-row equivalents of Supplement S11, S12 and S13.
# Tournament 1 has fixed performance; beta2=.10,.50,.90 labels the columns.
# Regenerate input files with 09_prepare_harvesting_phispace.R.
source('09_harvesting_model.R')
for(pkg in c('sf','isoband'))if(!requireNamespace(pkg,quietly=TRUE))stop('Install ',pkg)
h_interp<-function(z,axis,fine){
 v<-vapply(seq_along(axis),function(j)splinefun(axis,z[,j],method='monoH.FC')(fine),numeric(length(fine)))
 t(vapply(seq_along(fine),function(i)splinefun(axis,v[i,],method='monoH.FC')(fine),numeric(length(fine))))
}
h_bilinear<-function(z,axis,x,y){
 i<-pmax(1,pmin(length(axis)-1,findInterval(x,axis)));j<-pmax(1,pmin(length(axis)-1,findInterval(y,axis)))
 u<-(x-axis[i])/(axis[i+1]-axis[i]);v<-(y-axis[j])/(axis[j+1]-axis[j])
 (1-u)*(1-v)*z[cbind(i,j)]+u*(1-v)*z[cbind(i+1,j)]+(1-u)*v*z[cbind(i,j+1)]+u*v*z[cbind(i+1,j+1)]
}
h_mask_geometry<-function(mask,axis){
 if(!any(mask))return(NULL)
 b<-isoband::isobands(axis,axis,t(1*mask),.5,1.5)
 if(!length(b[[1]]$x))return(NULL)
 sf::st_sfc(isoband::iso_to_sfg(b)[[1]])
}
h_partition<-function(fields,axis,levels=rep(0,4)){
 families<-Map(function(f,l)contourLines(axis,axis,f,levels=l),fields,levels)
 lim<-range(axis);square<-sf::st_linestring(rbind(c(lim[1],lim[1]),c(lim[2],lim[1]),c(lim[2],lim[2]),c(lim[1],lim[2]),c(lim[1],lim[1])))
 net<-list(square);for(fam in families)for(cl in fam)net[[length(net)+1]]<-sf::st_linestring(cbind(cl$x,cl$y))
 poly<-sf::st_collection_extract(sf::st_polygonize(sf::st_node(sf::st_union(sf::st_sfc(net)))),'POLYGON')
 pts<-sf::st_coordinates(sf::st_point_on_surface(poly));v<-matrix(vapply(fields,function(z)h_bilinear(z,axis,pts[,1],pts[,2]),numeric(nrow(pts))),nrow=nrow(pts),ncol=4)
 ids<-rowSums(sweep(v,2,levels,'>')*rep(2^(0:3),each=nrow(v)))
 dissolved<-lapply(sort(unique(ids)),function(id)list(id=id,geometry=sf::st_union(poly[ids==id])))
 list(regions=dissolved,lines=families)
}
h_centre<-function(poly){
 bb<-sf::st_bbox(poly);best<-sf::st_coordinates(sf::st_point_on_surface(poly))[1,1:2];box<-unname(bb)
 for(k in 1:3){
  g<-expand.grid(x=seq(box[1],box[3],length.out=21),y=seq(box[2],box[4],length.out=21));pts<-sf::st_as_sf(g,coords=c('x','y'));ok<-lengths(sf::st_within(pts,poly))>0
  if(any(ok)){p<-pts[ok,];dist<-as.numeric(sf::st_distance(p,sf::st_boundary(poly)));best<-sf::st_coordinates(p[which.max(dist),])[1,1:2]}
  dx<-(box[3]-box[1])/20;dy<-(box[4]-box[2])/20;box<-c(max(bb[1],best[1]-2*dx),max(bb[2],best[2]-2*dy),min(bb[3],best[1]+2*dx),min(bb[4],best[2]+2*dy))
 }
 best
}
h_layout<-function(){layout(rbind(c(1,2,3,0,4,5,6),rep(7,7)),widths=c(1,1,1,.32,1,1,1),heights=c(1,.95));par(oma=c(.45,3.35,6.55,.45))}
h_panel<-function(beta2,c,diagonal=TRUE){
 par(mar=c(1.5,1.5,.5,.5),mgp=c(1.35,.36,0),tcl=-.24,xaxs='i',yaxs='i',pty='s',xpd=FALSE)
 plot(NA,xlim=c(.1,3.1),ylim=c(.1,3.1),axes=FALSE,xlab='',ylab='',asp=1)
}
h_axes<-function(beta2,column,diagonal=TRUE){
 if(diagonal)abline(0,1,lty=3,lwd=.7,col='#333333')
 axis(1,at=seq(.5,3,.5),cex.axis=1.02,gap.axis=.35);mtext(expression(varphi[1]),1,line=1.9,cex=1.062)
 if(column==1){axis(2,at=seq(.5,3,.5),cex.axis=1.02,gap.axis=.35);mtext(expression(varphi[2]),2,line=1.4,cex=1.062,las=1)}
 mtext(bquote(beta[2]==.(sprintf('%.2f',beta2))),3,line=.22,cex=1.18);box(lwd=.8)
}
h_titles<-function(){
 mtext(expression('Fixed '*E[1]==E[2]),3,outer=TRUE,at=.225,line=3.95,cex=1.55,font=2)
 mtext('First response to random allocations',3,outer=TRUE,at=.225,line=2.25,cex=1.1)
 mtext(expression('Responsive '*E[t]),3,outer=TRUE,at=.775,line=3.95,cex=1.55,font=2)
 mtext('Self-consistent population',3,outer=TRUE,at=.775,line=2.25,cex=1.1)
}
h_legend_start<-function(){par(mar=c(.03,.03,.03,.03),xpd=NA,pty='m');plot.new();w<-par('pin')[1];plot.window(xlim=c(0,w),ylim=c(0,1),xaxs='i',yaxs='i');w}
h_legend_row<-function(labels,colors,y,kinds=rep('box',length(labels)),cex=1.4){
 width<-par('usr')[2];sw<-.16;gap<-.09;between<-.25;tw<-strwidth(labels,cex=cex,units='inches');available<-width-.25-length(labels)*(sw+gap)-(length(labels)-1)*between
 if(sum(tw)>available){cex<-cex*available/sum(tw);tw<-strwidth(labels,cex=cex,units='inches')}
 x<-(width-sum(tw)-length(labels)*(sw+gap)-(length(labels)-1)*between)/2
 for(i in seq_along(labels)){
  if(kinds[i]=='box')rect(x,y-.045,x+sw,y+.045,col=colors[i],border=NA)
  if(kinds[i]=='hatch')rect(x,y-.045,x+sw,y+.045,density=18,angle=45,col=colors[i],border=NA)
  if(kinds[i]=='dot')points(x+sw/2,y,pch=20,cex=.8,col='black')
  if(kinds[i]=='circle')points(x+sw/2,y,pch=21,bg=colors[i],col='black',cex=1.3,lwd=.9)
  if(kinds[i]%in%c('solid','dash'))segments(x,y,x+sw,y,col=colors[i],lty=if(kinds[i]=='dash')2 else 1,lwd=2.1)
  if(kinds[i]%in%LETTERS[1:4]){points(x+sw/2,y,pch=21,bg='black',col='black',cex=2.7);text(x+sw/2,y,kinds[i],col='white',font=2,cex=1.05)}
  text(x+sw+gap,y,labels[i],adj=c(0,.5),cex=cex);x<-x+sw+gap+tw[i]+between
 }
}
h_export<-function(draw,stem){
 png(paste0(stem,'.png'),width=3500/300,height=4.35,units='in',res=600);tryCatch(draw(),finally=dev.off())
 cairo_pdf(paste0(stem,'.pdf'),width=3500/300,height=4.35,family='Arial');tryCatch(draw(),finally=dev.off())
}
prepare_harvest_plot_fields<-function(n=351L){
 file<-'09_harvesting_results/plot_fields.rds'
 if(file.exists(file)){saved<-readRDS(file);if(identical(saved$interpolation,'shape-preserving-signed-v4'))return(saved)}
 s<-boredom_slog();median_sum<-uniroot(function(z)integrate(function(a)dlnorm(a,0,s)*plnorm(z-a,0,s),0,z)$value-.5,c(1,4))$root
 fine<-seq(.1,3.1,length.out=n);g<-expand.grid(phi1=fine,phi2=fine);out<-list(axis=fine,median_sum=median_sum,panels=list(),interpolation='shape-preserving-signed-v4')
 for(b in c(.1,.5,.9)){
  p<-readRDS(paste0('09_harvesting_results/phi_beta_',b,'.rds'));res<-list(beta2=b)
  for(mode in c('first','responsive')){
   ratio<-exp(h_interp(log(p[[paste0('ratio_',mode)]]),p$axis,fine))
   cr<-harvest_core()$harvesting_crossings_cpp(g$phi1,g$phi2,as.vector(ratio),b,50,median_sum/2)
   presence<-lapply(p[[paste0('presence_',mode)]],h_interp,axis=p$axis,fine=fine)
   energy<-h_interp(p[[paste0('E_',mode)]],p$axis,fine)
   res[[mode]]<-list(cross=lapply(1:7,function(k)matrix(cr[,k],n)),presence=presence,energy=pmax(matrix(0,n,n),pmin(50,energy)))
  }
  out$panels[[as.character(b)]]<-res
 }
 saveRDS(out,file);out
}
draw_harvest_crossings<-function(d){
 axis<-d$axis;palette<-colorRampPalette(c('#315A9B','white','#B54343'))(801)
 draw<-function(){h_layout()
  for(mode in c('first','responsive'))for(c in 1:3){p<-d$panels[[c]];f<-p[[mode]]$cross;h_panel(p$beta2,c,FALSE)
   count<-f[[1]];last<-f[[2]];sign<-f[[3]];idx<-401+sign*round(400/(1+exp(pmin(700,ifelse(is.na(last),0,last)))))
   colors<-matrix(palette[pmax(1,pmin(801,idx))],length(axis));colors[count==0&sign>0]<-'#E06';colors[count==0&sign<0]<-'#08F'
   rasterImage(as.raster(t(colors)[ncol(colors):1,,drop=FALSE]),.1,.1,3.1,3.1,interpolate=FALSE)
   for(winner in c(-1,1)){
    geom<-h_mask_geometry(f[[4]]>.5&f[[5]]==winner&count>0,axis)
    if(!is.null(geom)){
     strokes<-sf::st_sfc(sf::st_multilinestring(lapply(seq(-5,5,by=.17),function(off)rbind(c(-3,-3+off),c(7,7+off)))))
     seg<-suppressWarnings(sf::st_intersection(strokes,geom));if(length(seg))plot(seg,add=TRUE,col=if(winner==1)'#E06'else '#08F',lwd=.8)
    }
   }
   ids<-seq(1,length(axis),by=10);ij<-expand.grid(i=ids,j=ids);keep<-count[cbind(ij$i,ij$j)]>=2;if(any(keep))points(axis[ij$i[keep]],axis[ij$j[keep]],pch=20,cex=.3,col='#333333')
   for(cl in contourLines(axis,axis,f[[6]],levels=0)){
    good<-h_bilinear(f[[7]],axis,cl$x,cl$y)>.5
    cl$x[!good]<-NA;lines(cl$x,cl$y,col='white',lwd=1.1);lines(cl$x,cl$y,col='black',lwd=1.1,lty=2)
   }
   h_axes(p$beta2,c,FALSE)
  }
  width<-h_legend_start()
  h_legend_row(c('No crossing: marginal tilt to 1','No crossing: marginal tilt to 2','50:50 not optimal at last crossing','Two crossings'),c('#E06','#08F','black','black'),.84,c('box','box','hatch','dot'))
  h_legend_row('50th percentile of talent sum','black',.64,'dash',cex=1.3)
  drop<-.06
  for(side in 1:2){left<-c(.12,.55)[side]*width;right<-c(.45,.88)[side]*width;xx<-seq(left,right,length.out=401);ind<-401+c(1,-1)[side]*round(400*seq(1,0,length.out=401))
   for(k in 1:400)rect(xx[k],.21-drop,xx[k+1],.35-drop,col=palette[ind[k]],border=NA)
   text(mean(c(left,right)),.44-drop,paste('After last crossing: marginal tilt to',side),cex=1.35)
   sums<-c(0,1,2,4,10,Inf);at<-left+(right-left)*(1-2/(2+sums));segments(at,.21-drop,at,.18-drop,lwd=.6);text(at,.12-drop,c('0','1','2','4','10',expression(infinity)),cex=1.2)
  }
  text(.025*width,.28-drop,'Talent sum',adj=0,cex=1.4);text(.89*width,.28-drop,'Later crossing\npaler color',pos=4,cex=1.4)
  h_titles()
 }
 h_export(draw,'Figure_harvesting_equal_talent_phispace')
}
draw_harvest_energy<-function(d){
 axis<-d$axis;palette<-colorRampPalette(c('#08F','white','#E06'))(801)
 geometry<-lapply(d$panels,function(p){a<-p$first$energy-25;b<-p$responsive$energy-25
  fields<-list(pmax(pmin(a,b-a),pmin(-a,a-b)),pmax(pmin(b,a-b),pmin(-b,b-a)),pmax(pmin(-b,a+b),pmin(b,-a-b)),pmax(pmin(a,-a-b),pmin(-a,a+b)))
  h_partition(fields,axis)
 })
 draw<-function(){h_layout()
  for(mode in c('first','responsive'))for(c in 1:3){p<-d$panels[[c]];h_panel(p$beta2,c)
   image(axis,axis,p[[mode]]$energy/50,col=palette,zlim=c(0,1),add=TRUE,useRaster=TRUE,axes=FALSE)
   h_axes(p$beta2,c)
   for(reg in geometry[[c]]$regions){if(reg$id==0)next
    plot(sf::st_boundary(reg$geometry),add=TRUE,col='#222222',lwd=.9)
    parts<-suppressWarnings(sf::st_cast(reg$geometry,'POLYGON'))
    for(k in seq_along(parts))if(as.numeric(sf::st_area(parts[k]))>.055){pt<-h_centre(parts[k]);clearance<-as.numeric(sf::st_distance(sf::st_sfc(sf::st_point(pt)),sf::st_boundary(parts[k])));if(clearance<.045||min(pt-.1,3.1-pt)<.06)next;id<-which(bitwAnd(reg$id,2^(0:3))>0)[1];points(pt[1],pt[2],pch=21,bg='black',col='black',cex=2.7);text(pt[1],pt[2],LETTERS[id],col='white',font=2,cex=1.05)}
   }
  }
  w<-h_legend_start();h_legend_row(c('Same leading tournament; tilt grows','Same leading tournament; tilt shrinks','New leader; smaller tilt','New leader; larger tilt'),rep('black',4),.82,LETTERS[1:4])
  xx<-seq(.24*w,.76*w,length.out=801);for(k in 1:800)rect(xx[k],.25,xx[k+1],.44,col=rev(palette)[k],border=NA)
  rect(.24*w,.25,.76*w,.44,lwd=.8);text(.24*w,.57,'Tournament 1',col='#E06',adj=0,cex=1.5);text(.76*w,.57,'Tournament 2',col='#08F',adj=1,cex=1.5)
  at<-.24*w+.52*w*c(0,.5,1);segments(at,.25,at,.19);text(at,.10,c('1','0.5','0'),cex=1.2);text(.8*w,.3,expression(E[1]/B),cex=1.4,adj=0)
  h_titles()
 }
 h_export(draw,'Figure_harvesting_invested_energy_phispace')
}
draw_harvest_coexistence<-function(d){
 axis<-d$axis;fill<-c('#F2F2F2','#AFC7EE','#7FAEE6','#A9D2FF','#E6B0B0','#D8C7C5','#D9D9D9','#CFE2F8','#F2A9C2','#CBB6E7','#E4C0D6','#C2C7F1','#F6C2CF','#E8C1D6','#F0D7D7','#CFE8C7')
 catcols<-c('#08F','#6F8FCF','#D27A7A','#E06');linecols<-c('#005BBB','#3E7CB1','#C06C84','#D0006F');lt<-c(1,2,2,1)
 geometry<-lapply(d$panels,function(p)lapply(c('first','responsive'),function(m)h_partition(p[[m]]$presence,axis,c(0,1e-4,1e-4,0))))
 draw<-function(){h_layout()
  for(im in 1:2)for(c in 1:3){p<-d$panels[[c]];g<-geometry[[c]][[im]];h_panel(p$beta2,c)
   for(reg in g$regions)plot(reg$geometry,add=TRUE,col=fill[reg$id+1],border=NA)
   h_axes(p$beta2,c)
   for(k in 1:4)for(cl in g$lines[[k]])lines(cl$x,cl$y,col=linecols[k],lty=lt[k],lwd=1.55)
   for(reg in g$regions){parts<-suppressWarnings(sf::st_cast(reg$geometry,'POLYGON'))
    for(k in seq_along(parts))if(as.numeric(sf::st_area(parts[k]))>.056){pt<-h_centre(parts[k]);dx<-.0645;pos<-rbind(c(pt[1],pt[2]+dx),c(pt[1]-dx,pt[2]),c(pt[1]+dx,pt[2]),c(pt[1],pt[2]-dx))
     for(bit in 1:4)if(bitwAnd(reg$id,2^(bit-1))>0)points(pos[bit,1],pos[bit,2],pch=21,bg=catcols[bit],col='black',cex=.92,lwd=.85)
    }
   }
  }
  h_legend_start();h_legend_row(c('2maxxers','2-leaning diversifiers','1-leaning diversifiers','1maxxers'),catcols,.70,rep('circle',4),cex=1.5)
  h_legend_row(c('2maxxer boundary','2-leaning boundary','1-leaning boundary','1maxxer boundary'),linecols,.30,c('solid','dash','dash','solid'),cex=1.5)
  h_titles()
 }
 h_export(draw,'Figure_harvesting_strategy_coexistence_phispace')
}
draw_all_harvesting_phispace<-function(){d<-prepare_harvest_plot_fields();draw_harvest_crossings(d);draw_harvest_energy(d);draw_harvest_coexistence(d);invisible(d)}
if(sys.nframe()==0L)draw_all_harvesting_phispace()

