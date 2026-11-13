# Reproducing the asymmetric example figures

Run from this scripts folder. Edit the single list in `05_asymmetric_examples.R`
to choose the examples and their column order. Both figures use this list.
The fixed row uses uniform random allocations independent of talents, with
E1 = E2 = B/2 and separately calculated M1 and M2. Its bars show the first
best response; the responsive row uses the self-consistent population.

```r
source("05_asymmetric_tournaments_figure_HDR.R")
source("05_asymmetric_tournaments_utility_contours.R")

strategy <- plot_asymmetric_hdr()
utility <- plot_asymmetric_utility(population = strategy$population)

# Optional exports (PNG: 600 dpi):
strategy <- plot_asymmetric_hdr(png = TRUE, pdf = TRUE)
utility <- plot_asymmetric_utility(population = strategy$population,
                                  png = TRUE, pdf = TRUE)
```

Alternatively, source `05_draw_asymmetric_examples.R` to draw both figures.
Its two switches control exports. Each individual figure function draws into
R's plotting window by default and creates no files unless requested. In a
headless Rscript session it draws only to an explicitly opened device or export.
Use `strategy$draw()` or `utility$draw()` to redraw without recalculating the
population. Any nonempty list of examples is accepted; the panel sizes are
preserved as the number of columns changes.

For a different set, pass the same `scenarios = my_examples` to both calls.
Each entry is a list with `beta1`, `beta2`, `phi1`, `phi2`. For example:

```r
my_examples <- list(
  list(beta1 = .8, beta2 = .7, phi1 = 1.5, phi2 = 1.5),
  list(beta1 = .9, beta2 = .1, phi1 = 2.8, phi2 = 2)
)
strategy <- plot_asymmetric_hdr(scenarios = my_examples)
utility <- plot_asymmetric_utility(scenarios = my_examples,
                                  population = strategy$population)
```

`population` reuse requires matching examples, B and population resolutions.
`nquad` controls the outer population quadrature (minimum 64); try 96 or 128
for a convergence check. `nxpop` controls the best-response allocation search,
whereas `nalpha`, `nx` and `nroot` control the displayed maps. The solver uses
`05_explorer_population.R`, `05_energy_allocation_quadrature.R` and
`05_energy_allocation_population_phispace.cpp`, plus Rcpp and Rtools for the
first compilation. Compiled code is cached for later calls.

In the utility figure, solid yellow marks the globally optimal 50:50 split.
The thin dashed yellow continuation is the equal-marginal-return curve where
50:50 is not globally optimal, including inside specialization regions. It
must not be interpreted as an optimal-allocation boundary there. Utility
normalization and population mean utility use the same deterministic moments
as the population solver. Utility maxima are refined between grid points.

The utility refinement is supplied in `05_asymmetric_utility_refinement.R`;
it compiles the shared C++ objective and checks every grid-detected local
maximum. Computed utility panels and strategy maps are cached within the
returned drawing functions, so PDF export and redraws reuse those maps.
