# Shared three-tournament mesh calculation for the manuscript and web app.
# Population scalers are supplied by the same many_responsive() solver.
many_cube_meshes <- function(beta,phi,B,scalers,n_cube=35L,talent_limit=2.5,
                             progress=function(detail)invisible(NULL)) {
  stopifnot(length(beta)==3,length(phi)==3,length(scalers)==2,n_cube>=9)
  # Extra resolution near the origin, where optimal shares vary fastest.
  axis3 <- talent_limit*seq(0,1,length.out=n_cube)^1.5
  cubeA <- as.matrix(expand.grid(alpha1=axis3,alpha2=axis3,alpha3=axis3))
  colours <- c('#EE0066','#0088FF','#40A56A')
tetra <- rbind(c(1,2,4,8),c(1,4,3,8),c(1,3,7,8),c(1,7,5,8),c(1,5,6,8),c(1,6,2,8))
edges <- rbind(c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4))
iso_mesh <- function(value,X=NULL,pair=NULL) {
  n <- length(axis3);polys <- list()
  for(k in 1:(n-1))for(j in 1:(n-1))for(i in 1:(n-1)) {
    start <- i+(j-1)*n+(k-1)*n*n
    ids <- start+c(0,1,n,n+1,n*n,n*n+1,n*n+n,n*n+n+1)
    v <- value[ids];v[abs(v)<1e-12] <- 1e-12
    if(min(v)>0 || max(v)<0)next
    for(l in seq_len(nrow(tetra))) {
      id <- ids[tetra[l,]];vv <- v[tetra[l,]]
      if(min(vv)>0 || max(vv)<0)next
      points <- matrix(numeric(),0,3);alloc <- matrix(numeric(),0,3)
      for(e in seq_len(6)) {
        a <- edges[e,1];b <- edges[e,2]
        if(vv[a]*vv[b]>=0)next
        fraction <- vv[a]/(vv[a]-vv[b])
        points <- rbind(points,cubeA[id[a],]+fraction*(cubeA[id[b],]-cubeA[id[a],]))
        if(!is.null(X))alloc <- rbind(alloc,X[id[a],]+fraction*(X[id[b],]-X[id[a],]))
      }
      if(nrow(points)<3)next
      if(nrow(points)==4) {
        centre <- colMeans(points);u <- points[1,]-centre;u <- u/sqrt(sum(u*u))
        normal <- c((points[2,2]-points[1,2])*(points[3,3]-points[1,3])-(points[2,3]-points[1,3])*(points[3,2]-points[1,2]),
          (points[2,3]-points[1,3])*(points[3,1]-points[1,1])-(points[2,1]-points[1,1])*(points[3,3]-points[1,3]),
          (points[2,1]-points[1,1])*(points[3,2]-points[1,2])-(points[2,2]-points[1,2])*(points[3,1]-points[1,1]))
        if(sum(normal^2)<1e-24)next
        normal <- normal/sqrt(sum(normal^2));v2 <- c(normal[2]*u[3]-normal[3]*u[2],normal[3]*u[1]-normal[1]*u[3],normal[1]*u[2]-normal[2]*u[1])
        P <- sweep(points,2,centre);o <- order(atan2(P%*%v2,P%*%u));points <- points[o,]
        if(!is.null(X))alloc <- alloc[o,]
      }
      if(!is.null(pair)) {
        # Clip EXACTLY at the shared three-way tie within each tetrahedron,
        # instead of keeping/discarding whole triangles near that junction.
        other <- setdiff(1:3,pair)
        margin <- rowMeans(alloc[,pair,drop=FALSE])-alloc[,other]
        clipped <- matrix(numeric(),0,3)
        for(a in seq_len(nrow(points))) {
          b <- if(a==nrow(points))1L else a+1L
          if(margin[a]>=0)clipped <- rbind(clipped,points[a,])
          if((margin[a]>=0)!=(margin[b]>=0)) {
            f <- margin[a]/(margin[a]-margin[b])
            clipped <- rbind(clipped,points[a,]+f*(points[b,]-points[a,]))
          }
        }
        points <- clipped
      }
      if(nrow(points)>=3)for(a in 2:(nrow(points)-1))
        polys[[length(polys)+1]] <- points[c(1,a,a+1),]
    }
  }
  polys
}
  equal_mesh <- function(pair,scaler) {
    # Parameterize a genuine equal-leading optimum by its shared investment
    # u and marginal return lambda. With all three active, u ranges from B/3
    # to B/2. The remaining patch has u=B/2 and the third activity excluded.
    # The marginal is strictly increasing in talent, so its inverse is unique.
    C <- scaler/exp(mean(log(scaler)));other <- setdiff(1:3,pair)
    marginal <- function(a,x,t)C[t]*phi[t]*a*beta[t]*(1+x)^(beta[t]-1)*
      (1+a*expm1(beta[t]*log1p(x)))^(phi[t]-1)
    talent <- function(lambda,x,t) {
      lo <- rep(0,length(lambda));hi <- rep(talent_limit,length(lambda))
      for(iter in 1:34){mid<-(lo+hi)/2;lower<-marginal(mid,x,t)<lambda;
        lo[lower]<-mid[lower];hi[!lower]<-mid[!lower]}
      a <- (lo+hi)/2;a[lambda==0]<-0;a
    }
    patch <- function(interior) {
      grid <- expand.grid(v=seq(0,1,length.out=n_cube)^1.5,
        u=if(interior)seq(B/3,B/2,length.out=n_cube) else seq(0,1,length.out=n_cube))
      u <- if(interior)grid$u else rep(B/2,nrow(grid))
      x3 <- if(interior)B-2*u else rep(0,length(u))
      upper <- pmin(marginal(talent_limit,u,pair[1]),marginal(talent_limit,u,pair[2]))
      if(interior)upper <- pmin(upper,marginal(talent_limit,x3,other))
      lambda <- grid$v*upper
      A <- matrix(0,nrow(grid),3);X <- A
      for(t in pair){A[,t]<-talent(lambda,u,t);X[,t]<-u}
      if(interior){A[,other]<-talent(lambda,x3,other);X[,other]<-x3}
      else A[,other]<-grid$u*pmin(talent_limit,lambda/(C[other]*phi[other]*beta[other]))
      fit <- many_allocate(A,beta,phi,B,C)
      if(!all(fit$within_tolerance))stop('Equal-investment surface accuracy check failed.')
      value <- rep(0,nrow(A))
      for(t in 1:3)value <- value+C[t]*expm1(phi[t]*log1p(A[,t]*expm1(beta[t]*log1p(X[,t]))))
      valid <- fit$utility-value<=1e-8*(1+abs(fit$utility))
      triangles <- list()
      for(j in 1:(n_cube-1))for(i in 1:(n_cube-1)) {
        k <- i+(j-1)*n_cube
        for(ids in list(c(k,k+1,k+n_cube+1),c(k,k+n_cube+1,k+n_cube)))
          if(all(valid[ids])) {
            P <- A[ids,,drop=FALSE]
            if(sum((P[1,]-P[2,])^2)>1e-20 && sum((P[1,]-P[3,])^2)>1e-20)
              triangles[[length(triangles)+1]] <- P
          }
      }
      triangles
    }
    c(patch(TRUE),patch(FALSE))
  }
  meshes <- lapply(1:2,function(row){
    progress(paste(if(row==1)'Fixed' else 'Responsive','cube boundaries'))
    specialist <- lapply(1:3,function(t) {
      # For fixed rival talents, increasing own talent can only favour full
      # specialization (positive talent-investment increasing differences).
      progress(paste('Tracing specialist',t,'in row',row))
      rival <- setdiff(1:3,t)
      pairs <- as.matrix(expand.grid(axis3,axis3))
      test <- matrix(0,nrow(pairs),3);test[,rival] <- pairs
      low <- many_allocate(test,beta,phi,B,scalers[[row]])
      test[,t] <- talent_limit
      high <- many_allocate(test,beta,phi,B,scalers[[row]])
      if(!all(low$within_tolerance) || !all(high$within_tolerance))stop('Boundary accuracy check failed.')
      inside_low <- low$allocation[,t]>=B*(1-1e-8)
      inside_high <- high$allocation[,t]>=B*(1-1e-8)
      threshold <- ifelse(inside_low,-talent_limit/(n_cube-1),talent_limit+talent_limit/(n_cube-1))
      pending <- which(inside_low!=inside_high)
      if(length(pending)) {
        lo <- rep(0,length(pending));hi <- rep(talent_limit,length(pending))
        for(iter in 1:23) {
          mid <- (lo+hi)/2;A <- test[pending,,drop=FALSE];A[,t] <- mid
          z <- many_allocate(A,beta,phi,B,scalers[[row]])
          if(!all(z$within_tolerance))stop('Boundary accuracy check failed.')
          active <- z$allocation[,t]>=B*(1-1e-8)
          hi[active] <- mid[active];lo[!active] <- mid[!active]
        }
        threshold[pending] <- (lo+hi)/2
      }
      # expand.grid varies its first coordinate fastest.
      indices <- match(cubeA[,rival[1]],axis3)+(match(cubeA[,rival[2]],axis3)-1)*n_cube
      cubeA[,t]-threshold[indices]
    })
    parts <- lapply(1:3,function(t)list(triangles=iso_mesh(specialist[[t]]),colour=colours[t],opacity=.5))
    for(pair in list(c(1,2),c(1,3),c(2,3)))parts[[length(parts)+1]] <- list(
      triangles=equal_mesh(pair,scalers[[row]]),colour='#D9AE00',opacity=.26)
    parts
  })
  list(axis=axis3,meshes=meshes)
}

many_mesh_payload <- function(meshes) {
  lapply(meshes,function(parts)lapply(parts,function(p)list(colour=p$colour,
    opacity=p$opacity,triangles=lapply(p$triangles,function(P)unname(split(t(P),rep(1:3,each=3)))))))
}
