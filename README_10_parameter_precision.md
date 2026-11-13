# Reproducing the convincing parameter maps

This technical note describes the original working files. In a source-only
checkout, generate the RDS and CSV inputs with the relevant stage of
[reproduce.R](reproduce.R) before using a redraw command.

For a redraw using the saved calculations, run from the scripts folder:

```r
source("10_convincing_efficiency_figures.R")
convincing_efficiency_figures()
```

For a new calculation of the regular 25-by-25 parameter grid:

```r
source("10_compute_convincing_efficiencies.R")
compute_convincing_efficiencies(export_figures=TRUE)
```

The calculation is checkpointed in `10_campaign_precision_results/cells`.
Successful cells are reused. All three convincing efficiencies share the
recipient response profile for each tournament pair. A completed calculation
is marked by `10_campaign_precision_results/complete.rds`. Plotting then reads
the compact `parameter_levels.rds` files and does not repeat the model work.

Recipient expectations use stratified deterministic talent integration of
order 32768, with adaptive Gaussian integration of rare participation tails.
The integrated energy and claim moments are represented by monotone cubic
curves in campaign pressure. This avoids amplifying small integration errors
through local differences between neighbouring pressure values. Individual
campaign decisions are then optimized and averaged over the talent density.
The linear-spoils global case uses its independent one-dimensional adaptive
talent-ratio reduction.

The numerical checks compare the method with full adaptive recipient
integration and denser focal integration at sensitive parameters. The displayed
colour bands are 0.5 percentage points wide. Parameter interpolation is
bilinear and introduces no new extrema; no blur is applied to the final maps.
In the checked difficult local cases, mean effort agrees with full adaptive
integration to about 0.05 percentage points. This is convergence evidence at
those parameters, not a uniform error bound. Very small contours on a nearly
flat plateau should therefore not be interpreted as additional regimes.

## Curvature certificates used to accelerate individual decisions

For a claim `w(x)`, its log marginal curvature is

`w''/w' = [(beta-1)+(phi-1)*alpha*beta*R(x)]/(1+x)`,

where `R(x)=(1+x)^beta/[1+alpha*((1+x)^beta-1)]` is monotone. Its endpoint
values therefore bound the numerator over the entire budget interval.

In a finite pool, with competitor energy `E` and claims `M`, the tournament
payoff is `(E+x)*w/(M+w)`. Its second derivative has the sign of

`2+(E+x)*(w''/w'-2*w'/(M+w))`.

The C++ code uses sufficient bounds on this expression to certify concavity
or convexity over the whole interval. Two concave components permit a unique
marginal-return solution, found with safeguarded Newton steps. Two convex
components require only endpoint comparison. If neither certificate applies,
the original global candidate search is retained. Cached claim values speed
up the population baseline iterations without changing their equations.

This is the same first-intervention model used in the Supplement. These maps
start from the responsive allocation equilibrium without convincing.
