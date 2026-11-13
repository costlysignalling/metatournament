# Stable interpolation for the equal-phi beta-plane panels.
# The population ratio contains a steep beta2/beta1 factor near the origin.
# Interpolate the smooth log marginal-return ratio after cancelling that factor,
# rather than interpolating log(ratio) and cancelling it afterwards.
# The 101-node input mesh was refined from a 51-node population mesh; use that
# original mesh to avoid re-interpolating its raw-log-ratio midpoint estimates.
equal_talent_beta_low_field <- function(ratio,beta_axis,target=beta_axis,B=50,direct=NULL) {
  ids <- seq(1L,length(beta_axis),by=2L)
  xs <- beta_axis[ids]
  z <- log(ratio[ids,ids])+outer(log(xs),log(xs),'-')+
    outer(xs,xs,'-')*log1p(B/2)
  # Exchanging tournament labels reverses the sign exactly.
  z <- (z-t(z))/2
  diag(z) <- 0
  tensor <- function(z,xs,target) {
    tmp <- vapply(seq_along(xs),function(j)
      splinefun(xs,z[,j],method='fmm')(target),numeric(length(target)))
    t(vapply(seq_len(nrow(tmp)),function(i)
      splinefun(xs,tmp[i,],method='fmm')(target),numeric(length(target))))
  }
  field <- tensor(z,xs,beta_axis)
  if(!is.null(direct)) {
    # Preserve additional directly evaluated population nodes on the fine mesh.
    truth <- log(ratio)+outer(log(beta_axis),log(beta_axis),'-')+
      outer(beta_axis,beta_axis,'-')*log1p(B/2)
    field[direct] <- truth[direct]
  }
  field <- (field-t(field))/2
  diag(field) <- 0
  out <- tensor(field,beta_axis,target)
  out <- (out-t(out))/2
  diag(out) <- 0
  out
}
