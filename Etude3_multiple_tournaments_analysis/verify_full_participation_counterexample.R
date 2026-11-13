# Individual-level counterexample; not a population probability calculation.
# At beta=1/2 and phi=3, with equal talents a and B=1, write
# u=sqrt(1+x), v=sqrt(2-x). The derivative of two-tournament utility is
# (3*a/2)*(u-v)*(a^2-(1-a)^2/(u*v)).
# Since u*v >= sqrt(2), a^2 > (1-a)^2/sqrt(2) makes its sign equal
# to that of x-1/2. The midpoint is the minimum and the endpoints are
# the global maxima. No numerical search for the two-tournament optimum
# is needed. If equal three-way allocation beats that endpoint value,
# every global optimum with three tournaments must use all three.
a <- .457
B <- 1
w <- function(x) (1+a*(sqrt(1+x)-1))^3-1
stopifnot(a^2 > (1-a)^2/sqrt(2))
u2 <- w(B)
u3_equal <- 3*w(B/3)
stopifnot(u3_equal > u2)
print(c(two_tournament_optimum=u2,three_way_equal_allocation=u3_equal,
        improvement=u3_equal-u2),digits=12)
