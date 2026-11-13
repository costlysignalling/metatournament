# Figure 02: specialization, generalism, and utility
# Six panels (A-F). Panel C uses beta = 0.1 and phi = 6.84.

set.seed(1)

# Parameters and simulated population --------------------------------------
n <- 4000
budget <- 50
beta <- 0.4
median_talent <- 1
sd_talent <- 0.5
phis <- c(2.34, 2.36, 2.39)

logmean <- log(median_talent)
logsd <- sqrt(log((1 + sqrt(1 + 4 * (sd_talent / median_talent)^2)) / 2))

r <- runif(n)
time <- cbind(budget * r, budget * (1 - r))
talent <- matrix(rlnorm(2 * n, logmean, logsd), ncol = 2)

training <- function(x, beta_value = beta) (1 + x)^beta_value - 1
weight <- function(P, phi) (1 + P)^phi - 1

perf <- talent * training(time)
E <- colSums(time)

utility <- function(phi) {
  spoils <- sapply(1:2, function(t) {
    w <- weight(perf[, t], phi)
    E[t] * w / sum(w)
  })
  rowSums(spoils)
}

utility1 <- utility(1)
entropy <- ifelse(r == 0, 0, -r * log(r)) +
           ifelse(r == 1, 0, -(1 - r) * log(1 - r))

# Expected utility for phi = 1 ---------------------------------------------
mean_training <- ((1 + budget)^(beta + 1) - 1) /
                 (budget * (beta + 1)) - 1

expected1 <- function(x) {
  budget / (2 * mean_training) * (training(x) + training(budget - x))
}

binary_entropy <- function(x) {
  ifelse(x == 0 | x == 1, 0,
         -x * log(x) - (1 - x) * log(1 - x))
}

allocation_from_entropy <- function(h) {
  if (h <= 0) return(0)
  if (h >= log(2)) return(0.5)
  uniroot(function(x) binary_entropy(x) - h, c(0, 0.5))$root
}

# Expected utility for arbitrary phi (Gauss-Hermite quadrature) -------------
gauss_hermite <- function(n = 80) {
  J <- matrix(0, n, n)
  off <- sqrt((1:(n - 1)) / 2)
  J[cbind(1:(n - 1), 2:n)] <- off
  J[cbind(2:n, 1:(n - 1))] <- off
  eig <- eigen(J, symmetric = TRUE)
  order_idx <- order(eig$values)
  list(z = eig$values[order_idx],
       w = eig$vectors[1, order_idx]^2)
}

gh <- gauss_hermite()
a <- exp(logmean + sqrt(2) * logsd * gh$z)

mean_weight <- function(x, phi, beta_value = beta) {
  vapply(x, function(xx) sum(gh$w * weight(a * training(xx, beta_value), phi)), numeric(1))
}

mean_weight_prime <- function(x, phi, beta_value = beta) {
  vapply(x, function(xx) {
    A <- 1 + a * training(xx, beta_value)
    sum(gh$w * phi * a * beta_value * (1 + xx)^(beta_value - 1) * A^(phi - 1))
  }, numeric(1))
}

mean_weight_second <- function(x, phi, beta_value = beta) {
  vapply(x, function(xx) {
    A <- 1 + a * training(xx, beta_value)
    bracket <- (beta_value - 1) * (1 - a) + a * (beta_value * phi - 1) * (1 + xx)^beta_value
    sum(gh$w * phi * a * beta_value * (1 + xx)^(beta_value - 2) * A^(phi - 2) * bracket)
  }, numeric(1))
}

expected_phi <- function(x, phi, beta_value = beta) {
  normalization <- integrate(function(z) mean_weight(budget * z, phi, beta_value),
                             0, 1, rel.tol = 1e-10)$value
  budget / (2 * normalization) *
    (mean_weight(x, phi, beta_value) + mean_weight(budget - x, phi, beta_value))
}

# Find the two local minima, when they exist.
phi_end <- uniroot(function(phi) {
  mean_weight_prime(0, phi) - mean_weight_prime(budget, phi)
}, c(1, 4))$root

phi_mid <- uniroot(function(phi) {
  mean_weight_second(budget / 2, phi)
}, c(1, 4))$root

local_minima <- function(phi) {
  if (phi <= phi_end || phi >= phi_mid) return(numeric(0))
  f <- function(x) mean_weight_prime(x, phi) - mean_weight_prime(budget - x, phi)
  xmin <- uniroot(f, c(1e-10, budget / 2 - 1e-8))$root
  c(xmin, budget - xmin)
}

# Plot ---------------------------------------------------------------------
label_panel <- function(letter) {
  usr <- par("usr")
  mtext(letter, side = 3, at = usr[1] - 0.22 * diff(usr[1:2]),
        line = 0.5, font = 2, cex = 1.1)
}

point_col <- "#00808020"
curve_col <- 2
curve_lwd <- 2
xgrid <- seq(0, budget, length.out = 1001)

w<-2
png("Figure02.png", width = w*3, height = w*2, units = "in", res = 600)
layout(matrix(1:6, nrow = 2, byrow = TRUE))
par(mar = c(3.5, 3.5, 2, 0.5), mgp = c(2, 0.65, 0))

# A: expected utility against allocation entropy
plot(entropy, utility1, pch = 19, col = point_col,
     ylim = c(0, max(utility1)),
     xlab = "Time allocation entropy", ylab = "Total utility")
title(expression(beta == 0.4 * ", " ~ varphi == 1), adj = 0)

text(range(entropy), 4, c("specialists", "generalists"),
     pos = c(4, 2), cex = 0.8, offset = 0)
hgrid <- seq(0, log(2), length.out = 501)
rgrid <- vapply(hgrid, allocation_from_entropy, numeric(1))
lines(hgrid, expected1(budget * rgrid), lwd = curve_lwd, col = curve_col)
label_panel("A")

# B: total utility against training allocation for beta = 0.4, phi = 1
plot(time[, 1], utility1, pch = 19, col = point_col,
     ylim = c(0, max(utility1)), xlab = "Time invested in tournament 1",
     ylab = "")
title(expression(beta == 0.4 * "," ~ varphi == 1), adj = 0)
text(c(0, budget), 4, c("spec 2", "spec 1"),
     pos = c(4, 2), cex = 0.8)
text(budget / 2, 4, "generalists", cex = 0.8)
lines(xgrid, expected1(xgrid), lwd = curve_lwd, col = curve_col)
label_panel("B")

# C: intermediate allocations are optimal (middle scenario of Figure 02B)
# The random draws are identical to those in the supplied 02B script:
# both scripts set.seed(1), then draw r, then two columns of lognormal talents.
beta_c <- 0.1
phi_c <- 6.84
perf_c <- talent * training(time, beta_c)
utility_c <- {
  spoils <- sapply(1:2, function(t) {
    weights <- weight(perf_c[, t], phi_c)
    E[t] * weights / sum(weights)
  })
  rowSums(spoils)
}

# Off-center stationary points can be maxima or minima in this regime.
offcenter_extrema <- function(phi, beta_value) {
  f <- function(x) mean_weight_prime(x, phi, beta_value) -
                   mean_weight_prime(budget - x, phi, beta_value)
  xs <- seq(1e-8, budget / 2 - 1e-6, length.out = 1001)
  ys <- f(xs)
  k <- which(ys[-length(ys)] * ys[-1] < 0)
  if (!length(k)) return(numeric(0))
  x <- uniroot(f, c(xs[k[1]], xs[k[1] + 1]))$root
  c(x, budget - x)
}

plot(time[, 1], utility_c, pch = 19, col = point_col,
     xlim = c(0, budget), ylim = c(49.4, 50.8),
     xlab = "Time invested in tournament 1", ylab = "",
     xaxs = "i", yaxs = "i")
lines(xgrid, expected_phi(xgrid, phi_c, beta_c),
      lwd = curve_lwd, col = curve_col)
ex <- offcenter_extrema(phi_c, beta_c)
if (length(ex)) {
  abline(v = ex, lty = 2)
  points(ex, expected_phi(ex, phi_c, beta_c), pch = 19, col = curve_col)
}
abline(v = budget / 2, lty = 3)
title(bquote(beta == .(beta_c) * "," ~ varphi == .(phi_c)), adj = 0)
title("(zoom)", adj = 1, font.main = 1)
label_panel("C")

# D-F: varying the performance-to-spoils exponent
for (j in seq_along(phis)) {
  phi <- phis[j]
  plot(time[, 1], utility(phi), pch = 19, col = point_col,
       xlim = c(0, budget), ylim = c(49, 51),
       xlab = "Time invested in tournament 1",
       ylab = if (j == 1) "Total utility" else "",
       xaxs = "i", yaxs = "i")
  lines(xgrid, expected_phi(xgrid, phi), lwd = curve_lwd, col = curve_col)
  mins <- local_minima(phi)
  if (length(mins)) {
    abline(v = mins, lty = 2)
    points(mins, expected_phi(mins, phi), pch = 19, col = curve_col)
  }
  abline(v = budget / 2, lty = 3)
  title(bquote(beta == 0.4 * "," ~ varphi == .(phi)), adj = 0)
  label_panel(LETTERS[j + 3])
}

dev.off()
