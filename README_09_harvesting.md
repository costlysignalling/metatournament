# Fixed performance and finite competitor pools

This technical note describes the original working files. In a source-only
checkout, generate the RDS and CSV inputs with the relevant stage of
[reproduce.R](reproduce.R) before using a redraw command.

Run R from this scripts folder. Model calculations use Rcpp/C++; plots use `sf` and `isoband`. Windows compilation requires Rtools matching the R version.

```r
# install.packages(c("Rcpp", "sf", "isoband"))
source("09_reproduce_harvesting_figures.R")
reproduce_harvesting_figures()                 # reuse supplied population data
reproduce_harvesting_figures(recalculate=TRUE) # recompute all populations first
```

Outputs are 600-dpi PNGs and vector PDFs. The parameter-plane images use a 41-by-41 grid of calculated populations and shape-preserving interpolation for drawing; the saved input grid is separate from the drawing cache. Delete `09_harvesting_results/plot_fields.rds` to rebuild only the drawing cache.

## Individual allocations

```r
source("09_harvesting_model.R")
p <- harvest_population(beta2=.4, phi1=.66, phi2=1.5, B=50)
actors <- rbind(c(1,1), c(.5,2), c(2,.5))
harvest_alloc(p, "first", actors)      # time in tournament 1
harvest_alloc(p, "responsive", actors)
```

Tournament 1 always has fixed performance. Set `both_fixed=TRUE` for fixed performance in tournament 2 too; `beta2` is then irrelevant. Otherwise tournament 2 uses the original training curve. Both talent distributions are independent lognormals with median 1 and SD .5. The fixed reference uses uniform random time allocations with equal mean investments and separately calculated spoils totals. The responsive reference reproduces both investments and spoils totals.

The harvesting solver constructs the upper concave envelope of the training claim. Its contact points give global individual optima, including discontinuous switches. Conditional integration uses exact lognormal probabilities for specialist regions, Gaussian quadrature over the other talent, and split Gauss–Legendre integration over the interior region. Responsive normalization solves a scalar consistency equation for the ratio of the two population scalers. A checked root is returned; uniqueness of nonlinear responsive solutions is not asserted.

## Finite pools without boredom

```r
source("08_finite_pool_model.R")
p <- finite_pool_population(N=30, beta=c(.8,.4), phi=c(1,1.5), B=50)
finite_pool_alloc(p, "first", rbind(c(1,1),c(.5,2)))
finite_pool_alloc(p, "responsive", rbind(c(1,1),c(.5,2)))
```

This uses an expected competitor density of mass `N-1`, with the focal actor's pool and spoils contributions included exactly. It is not an exact average of equilibria in random finite groups. The solver compares local maxima and endpoints over fully flexible time allocations. Its computational midpoint has no compulsory-attendance interpretation. The supplied examples use 32,768 deterministic, exchange-symmetric Halton talent nodes.

## Files and checks

- `09_fixed_performance_harvesting_figure.R`: continuous and interspersed rectangles.
- `09_harvesting_model.R`, `09_harvesting_core.cpp`: harvesting populations and individual choices.
- `08_finite_pool_model.R`, `08_finite_pool_core.cpp`: finite-pool comparison without boredom.
- `09_strategy_examples.R`: all four-column strategy figures.
- `09_prepare_harvesting_phispace.R`: parameter-grid populations and signed coexistence fields.
- `09_harvesting_phispace_figures.R`: equal-talent geometry, invested energy, and coexistence figures.
- `09_harvesting_results/`: cached populations, drawing fields, and example precision results.
- `09_verify_harvesting_comparisons.R`: independent global-optimization and higher-order quadrature checks.
- `09_verify_harvesting_boundaries.R`: independent envelope-contact checks and higher-order checks at parameter-grid corners.

The model scripts reuse the existing Gaussian, lognormal and deterministic integration helpers in `08_boredom_model.R`. Keep that file in the folder; the harvesting calculation itself does not include boredom.

The example checks doubled Gaussian/Legendre orders from 48/16 to 96/32; population time-share changes were below 0.0000002. At the checked parameter-grid corners and centres, the largest change was about 0.000076. Individual allocations were compared with an independent grid-plus-refinement global search. The signed coexistence criteria were checked against envelope contact intervals at 9,455 talent/parameter combinations. For the finite-pool examples, doubling talent-density resolution changed the best-response population time share by at most 0.000044 at the saved responsive backgrounds. These checks assess numerical resolution; they do not establish uniqueness of a responsive equilibrium.


## Toggle comparison with training at beta 1

Open `Compare_harvesting_P_alpha_vs_beta1.html` in a browser. It is a self-contained offline document with an independent P = α / β = 1 toggle for the fixed-performance examples in the Paper and Supplement. Its figure numbers are read from the current manuscript sources. Both images use the same dimensions, plotting code, approved S14 legend and parameter layout.

In the β = 1 alternative, each tournament originally defined by fixed performance is replaced by the original training claim with β = 1. Both tournaments change in columns A–C of the Paper’s fixed-performance figure; only tournament 1 changes in its column D and in the other two comparison figures. All φ values, budgets and talent distributions stay the same. Both the random-allocation reference spoils totals and the responsive population are recalculated for the alternative. This compares different models, rather than changing a label on the same calculations.

The claims being compared are x[(1+α)^φ−1] and (1+αx)^φ−1. They coincide at φ = 1; at other φ values the strategy maps can differ.

To reproduce the comparison from this folder:

```r
source("09_training_comparison.R")
draw_training_comparison()
```

Then run `09_build_training_comparison.py` with Python to embed the six PNGs into the offline HTML. That builder uses only Python's standard library. Saved comparison populations are in `09_harvesting_results/training_comparison.rds`; `09_verify_training_comparison.R` checks higher-resolution population integration and independent global individual optimization. To recalculate populations before drawing, run `prepare_training_comparison(overwrite=TRUE)`.


## Explanatory rectangle and transformed talent

`09_harvesting_explanations.R` draws the compact Paper figure and provides exact transformed-talent CDF/quantile functions, using the original independent lognormal talent convention (median 1, SD .5).

```r
source("09_harvesting_explanations.R")
plot_harvesting_paper_explanation(phi=1.25)
q_harvesting_talent(.5, phi=2) # transformed median = 3
p_harvesting_talent(3, phi=2) # probability = .5
```

The transformed talent is Z=(1+alpha)^phi-1. It gives an equivalent model with beta=1 and linear spoils. It is different from replacing fixed performance by beta=1 training while retaining the original phi, as in the toggle comparison. The Supplement also derives the convexity of a cumulative claim formed by integrating the spoils rate during training; that alternative is normalized by its population-average cumulative claim.


## Original and transformed talent densities

```r
source("09_transformed_talent_density.R")
plot_transformed_talent_density(phis=c(1.25,2))
d_harvesting_talent(c(.5,1,2), phi=1.25)
```

The density is calculated exactly by the change-of-variable formula, using the inverse alpha=(1+z)^(1/phi)-1 and its Jacobian. No density estimation or simulation is used. Both panels retain the original value scale of the corresponding talents; they share a density scale and use separate horizontal limits reaching at least the transformed 99.5th percentile. The verification helper `09_verify_transformed_talent_density.R` checks unit probability mass, consistency with quantiles and equality with the original lognormal density at phi=1.
