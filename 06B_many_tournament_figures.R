# Five talent slices and a genuine three-dimensional allocation cube.
# Run 06B_verify_many_tournaments.R and 06B_refine_many_population.R first.
# Layout, specialist colours, pale density contours and yellow allocation
# boundaries follow the two-tournament figures. PNG output is 600 dpi.
for(e in parse('06B_predict_many_tournaments.R'))
  if(is.call(e) && identical(e[[1]],as.name('<-')))eval(e)
z <- readRDS('06B_many_results/three_tournament_population.rds')
latest <- tail(z$checks,1)[[1]]
scalers <- list(z$fixed$scaler,latest$responsive$scaler)
beta <- z$beta;phi <- z$phi;B <- z$B
colours <- c('#EE0066','#0088FF','#40A56A')
# Load definitions only, preserving editable demonstration calls without
# executing them while rebuilding manuscript figures.
for(e in parse('06B_explore_talent_slices.R'))
  if(is.call(e) && identical(e[[1]],as.name('<-')))eval(e)
result <- list(beta=beta,phi=phi,B=B,fixed=z$fixed,responsive=latest$responsive)
pdf(NULL)
slice_result <- many_slices(result=result,n_grid=181L,png=TRUE,pdf=TRUE,
  file='Figure_many_asymmetric_talent_slices')
dev.off()
saveRDS(list(alpha=slice_result$alpha,cuts=slice_result$cuts,slices=slice_result$slices),
  '06B_many_results/talent_slices.rds')

# Marching tetrahedra: transparent specialist boundary surfaces and surfaces
# on which two jointly leading activities receive equal positive investment.
# These are surfaces, not generally planes. The cube has the SAME 0--2.5 axes.
source('06B_three_tournament_geometry.R')
geometry <- many_cube_meshes(beta,phi,B,scalers,n_cube=35L,progress=message)
axis3 <- geometry$axis;meshes <- geometry$meshes
saveRDS(list(axis=axis3,meshes=meshes),'06B_many_results/allocation_cube_meshes.rds')

project3 <- function(P) {
  # Near-front view: alpha2 runs right, alpha1 down-left, alpha3 upward.
  # These two angles alone control the static camera.
  az <- -1.77;el <- .43
  R <- matrix(c(cos(az),-sin(az),0,sin(az)*sin(el),cos(az)*sin(el),cos(el),sin(az)*cos(el),cos(az)*cos(el),-sin(el)),3,3,byrow=TRUE)
  P%*%t(R)
}
cube_vertices <- as.matrix(expand.grid(c(0,2.5),c(0,2.5),c(0,2.5)))
cube_edges <- which(as.matrix(dist(cube_vertices))==2.5,arr.ind=TRUE);cube_edges <- cube_edges[cube_edges[,1]<cube_edges[,2],]
draw_cubes <- function() {
  old <- par(no.readonly=TRUE);on.exit(par(old));par(mfrow=c(1,2),mar=c(.5,.5,2.5,.5),oma=c(1.9,0,.2,0))
  for(row in 1:2) {
    v <- project3(cube_vertices)
    plot(NA,xlim=range(v[,1])+c(-.3,.35),ylim=range(v[,2])+c(-.3,.35),asp=1,axes=FALSE,xlab='',ylab='')
    parts <- meshes[[row]]
    triangles <- lapply(do.call(c,lapply(parts,function(p)p$triangles)),project3)
    fills <- unlist(lapply(parts,function(p)rep(adjustcolor(p$colour,alpha.f=p$opacity),length(p$triangles))))
    depth <- vapply(triangles,function(Q)mean(Q[,3]),numeric(1))
    o <- order(depth)
    xx <- unlist(lapply(triangles[o],function(Q)c(Q[,1],NA_real_)))
    yy <- unlist(lapply(triangles[o],function(Q)c(Q[,2],NA_real_)))
    polygon(xx,yy,col=fills[o],border=NA)
    for(e in seq_len(nrow(cube_edges)))lines(v[cube_edges[e,],1:2],col='grey35',lwd=.65)
    par(xpd=NA) # Keep mathematical axis labels fully visible.
    for(t in 1:3) {
      P <- matrix(0,1,3);P[t] <- 2.72;Q <- project3(P)
      text(Q[1,1],Q[1,2],bquote(alpha[.(t)]),cex=1.05)
      for(a in c(0,1,2)) {
        P <- matrix(0,1,3);P[t] <- a;Q <- project3(P)
        if(a>0 || t==1)text(Q[1,1],Q[1,2]-.08,a,cex=.7)
      }
    }
    par(xpd=FALSE)
    mtext(if(row==1)expression(Fixed~equal~E[t]~reference) else expression(Responsive~E[t]),3,line=.3,cex=1.15)
  }
  par(oma=rep(0,4),fig=c(0,1,0,1),new=TRUE,mar=rep(0,4));plot.new()
  legend('bottom',horiz=TRUE,bty='n',inset=.01,pch=15,col=c(colours,'#D9AE00'),
    legend=c('Specialist 1 boundary','Specialist 2 boundary','Specialist 3 boundary','Equal leading investments'),cex=.9)
}
png('Figure_many_asymmetric_allocation_cube.png',width=9,height=4.6,units='in',res=600,type='cairo');draw_cubes();dev.off()
pdf('Figure_many_asymmetric_allocation_cube.pdf',width=9,height=4.6);draw_cubes();dev.off()

draw_cubes()

# Export the same triangles for a dependency-free interactive HTML viewer.
payload <- lapply(meshes,function(parts)lapply(parts,function(p)list(colour=p$colour,
  opacity=p$opacity,triangles=lapply(p$triangles,function(P)unname(split(t(P),rep(1:3,each=3)))))))
jsonlite::write_json(payload,'06B_many_results/allocation_cube_meshes.json',auto_unbox=TRUE,digits=7)
cat('Both figures and interactive-viewer mesh data written.\n')
