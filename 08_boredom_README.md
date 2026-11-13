# Local locked-time and boredom figures

This technical note describes the original working files. In a source-only
checkout, generate the RDS and CSV inputs with the relevant stage of
[reproduce.R](reproduce.R) before using a redraw command.

These scripts implement Supplement 8. They use a smooth lognormal competitor
density with mass `N-1`, together with the focal actor's own effect on both pools
and spoils totals. This is a specified expected-density model, rather than an
exact average of random finite-classroom equilibria.

Defaults are `N=30`, `B=50`, compulsory proportion `0.5`, and independent
lognormal talents with median 1 and SD 0.5. Boredom redirects pool contribution
to tournament 2 without improving performance there.

Run from the scripts folder. To redraw the supplied results without compiling:

```r
source('08_boredom_examples.R')
plot_boredom_examples()
source('08_boredom_energy_parameter_space.R')
```

The example function accepts `png=FALSE, pdf=FALSE` to display only in R.
For the energy figure, set `options(metatournament.png=FALSE,
metatournament.pdf=FALSE)` before sourcing it. PNG output uses 600 dpi.

To rebuild all data (requires Rcpp and Rtools on Windows):

```r
source('08_prepare_boredom_figures.R')
prepare_boredom_figures()
# Example change; recompilation is not needed when changing parameters:
# prepare_boredom_figures(compulsory=.6, N=15)
```

To calculate one population reference:

```r
source('08_boredom_model.R')
p <- boredom_population(beta=c(.8,.6), phi=c(1.5,1),
     compulsory=.5, N=30, density=boredom_qmc_density(8192))
p$responsive$converged
p$responsive$residual
p$first$mapped
p$responsive$mapped
```

The four `mapped` values are mean productive investment in tournament 1,
mean spoils weight in tournament 1, mean spoils weight in tournament 2,
and mean boredom. Mean productive time in tournament 2 is `B - mapped[1] -
mapped[4]`. Its mean pool contribution is `B - mapped[1]`.

The global one-variable search includes both endpoints, the compulsory-time
kink and interior maxima on both sides. Responsive density calculations solve
three consistency equations, using damped updates and a Newton correction.
The returned residual measures consistency of mean pool contribution and both
mean spoils weights; it does not measure integration error.

The supplied example maps use 32,768 deterministic, exchange-symmetric Halton
points in lognormal quantiles. The parameter grid uses 2,048 points initially,
with higher-resolution retries where necessary. No sampling of individual
classrooms is involved. Independent checks at increased integration and
allocation resolution are saved with the results. To assess precision after
changing parameters, increase both resolutions and compare the returned means.
`08_refine_boredom_density_grid.R` increases the integration resolution around
the sensitive corner transitions before the final plot is drawn.

To read the supplied checks, or recalculate them:

```r
source('08_verify_boredom.R')
verify_boredom()
# verify_boredom(recalculate=TRUE)
```

Both figures show productive time and boredom, which sum to the time budget.
The Aâ€“D regions compare the imbalance of productive tournament investments.
Boredom remains part of tournament 2's pool, unlike the self-directed time
in Supplement 7. The lost symmetry between tournaments is intentional:
compulsory attendance applies only to tournament 1.

The example maps root their boundary locations directly, rather than tracing
pixel masks. For `c=0.5`, equal productive investments can fill an area, whose
edges are yellow; dashed black bounds full boredom and pink full investment
in tournament 1. Sidebars show the integrated time components.

Files are named automatically with compulsory proportion and local group mass.
`08_boredom_core.cpp` is the model implementation; the RDS files contain everything
needed for redraws. All numerical calculations use R and C++.

The example layout follows Figure S26, including its right-hand legend and
font sizes. The full equal-productive-investment region is yellow. The additional examples use 32,768 integration points. Additional
examples can be rebuilt with `source('08_prepare_more_boredom.R');
prepare_more_boredom()` and drawn with
`plot_boredom_examples(variant='more_boredom')`.

The C++ backend remains the default. If Windows blocks the compiled library,
the model automatically uses `08_boredom_R_fallback.R`, a vectorized R
implementation of the same optimization and consistency equations. Select it
explicitly with `options(metatournament.boredom.backend='R')`. It is slower
for large integrations but does not require changing security settings.

`08_refine_boredom_boundaries.R` evaluates equality of the two productive-time
imbalances directly at higher integration resolution. These roots refine the
classification boundaries in the two sensitive panels; population time colours
are retained from the numerical means.

## Interactive boredom explorer

Run `source('08_explore_boredom.R')` from this folder. In interactive R this
draws the single default example. Change `settings` near the top of the file,
or call the function directly:

```r
result <- explore_boredom(beta=c(.8,.8), phi=c(1.5,2),
                          compulsory=.5, N=30, B=50)
explore_boredom(beta1=.6,beta2=.8,phi1=1,phi2=1.25,
               compulsory=.65,N=10,legend=TRUE)
```

Use `mode='fixed'`, `'responsive'` or `'both'`. The default draws in R and
exports no files; enable `png=TRUE` and/or `pdf=TRUE` to save output (600 dpi
PNG). A custom `file='my_example'` sets the output stem. `legend=FALSE` by
default, matching the Script 05 explorer. The fixed condition uses the
compulsory-attendance reference, rather than equal tournament pools.

`result$summary` reports mean productive investments, bored time, pool
contributions and percentages of actors with some or full compulsory-time
boredom. Pool 2 includes the redirected bored-time contribution.
`result$population$responsive$residual` reports consistency error.

`plot(result,legend=TRUE)` redraws immediately. The result contains data only
and can be saved with `saveRDS(result,'my_example.rds')`, then redrawn after
sourcing the explorer with `plot(readRDS('my_example.rds'))`. Repeated calls
with unchanged population parameters reuse the last population calculation.

Use `npop=32768` and larger `nx`/`nxpop` for a precision check; `nalpha`
controls talent-map resolution and `alim=c(0,2.5)` the visible talent range.
The responsive calculation retries at denser integration resolutions when
needed and stops if consistency remains unresolved. The two-dimensional
50:50 region is yellow when `compulsory=.5`; otherwise globally optimal
equal-productive-investment curves are checked before plotting them.

## Main-text boredom examples and separate boredom scale

Run `source("08_prepare_paper_boredom.R")` and then `prepare_paper_boredom()`
to calculate the four main-text examples. `draw_paper_boredom()` redraws both
the allocation maps and their Supplement utility contours from the saved RDS.
The four columns use the specified beta/phi combinations and N=30,30,30,10.
Population integration uses 32768 deterministic points; each equilibrium is
checked before saving. The R/C++ model sources above are shared with the explorer.

Run `source("08_boredom_level_parameter_space.R")` and `plot_boredom_level()`
to redraw the separate greyscale plot of population-average bored time.
It reads `level_refined_c_0.5_N_30.rds`, using the same tournament parameters
and population-reference model as the combined invested-time figure. To rebuild,
run `source("08_refine_boredom_level.R"); refine_boredom_level()` first.
This needs `08_boredom_conditional.cpp` as well as the shared model sources.
The small bored-time mean is integrated conditionally on talent 2: roots find
the complete-disengagement boundary and the onset of boredom in talent 1.
The complete-disengagement contribution is integrated using the lognormal CDF;
Gauss–Legendre quadrature integrates the partial-boredom strip and Gauss–Hermite
quadrature integrates talent 2. This avoids estimating a very small region from
only a few ordinary QMC talent points. Responsive backgrounds retain the checked
QMC calculation. `conditional_precision_checks.rds` records independent
quadrature- and background-resolution checks.
Its common scale uses the actual calculated range across both conditions,
in percent of total B, rather than the full compulsory-time allowance.
