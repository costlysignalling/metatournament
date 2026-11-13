# Opting out of tournaments

This technical note describes the original working files. In a source-only
checkout, generate the RDS and CSV inputs with the relevant stage of
[reproduce.R](reproduce.R) before using a redraw command.

Run from this scripts folder in R. The shared model requires Rcpp and a C++ compiler (Rtools on Windows). C++ is compiled into R's local temporary folder to avoid loading a DLL from Dropbox.

To reproduce all twelve Supplement figures, including PNG and PDF versions:

```r
source("07_optout_figures.R")
```

For one figure:

```r
options(metatournament.eta=0.5)
source("07_known_talents_strategy_map_optout_fixed_E_HDR.R")
options(metatournament.eta=1.0)
source("07_asymmetric_tournaments_optout_HDR.R")
```

The value of eta appears automatically in the figure and output filenames. The PNG resolution is 600 dpi. Disable an export with `options(metatournament.png=FALSE)` or `options(metatournament.pdf=FALSE)`. In interactive R, figures also appear in the plotting window; `options(metatournament.plot=FALSE)` disables this additional drawing.

The fixed population reference is `X1 ~ Uniform(0,B)`, `X2=B-X1`, independently of talents. Thus `E1=E2=B/2`, but `M1` and `M2`, and hence their ratios with E, are calculated separately. Bars show the first response, not the initial random allocations. Responsive calculations solve for both absolute E/M ratios. At eta=1, the inactive population is represented by zero effective tournament payoffs, with all time self-directed; E/M itself is undefined there.

The six canonical figure scripts are:

- `07_known_talents_strategy_map_optout_fixed_E_HDR.R`
- `07_known_talent_optout_utility_contours_fixed_E.R`
- `07_known_talents_strategy_map_optout_HDR_responsive_E.R`
- `07_known_talent_optout_utility_contours_responsive_E.R`
- `07_asymmetric_tournaments_optout_HDR.R`
- `07_asymmetric_tournaments_optout_utility_contours.R`

The former four filenames containing `kappa` forward to these current versions. The asymmetric scripts read the same `05_asymmetric_examples.R` list as the main-text figures. Plotting helpers are in `00_helper_functions.R`. Model calculations are shared in `07_optout_model.R` and `07_optout_core.cpp`.

Numerical controls can be set before sourcing a script:

```r
options(metatournament.optout.npop=8192L,  # deterministic population integration
        metatournament.optout.nalpha=251L, # displayed talent raster
        metatournament.optout.nscan=129L)  # stationary-point boundary search
```

The continuous solver includes every stationary candidate in the slack-budget interior by splitting each marginal-return function at its possible inflection point. On the full-budget edge it brackets stationary points using a regular scan supplemented by points close to both endpoints, refines roots, and compares utilities with all other candidates. Increase `nscan` to check unusual parameter combinations. Identical population solutions and talent maps are reused within an R session. The root tolerance is a numerical consistency check; population integration still needs its own resolution check.

Run `source("07_verify_optout.R")` for independent reference, dense-grid allocation, search-resolution, resource-accounting and endpoint checks.

## Main-text examples and invested energy

`source("07_paper_optout_examples.R")` draws the selected A--D talent-space columns, with eta in each column header. `07_pilot_asymmetric_tournaments_optout.R` defines `pilot_asymmetric_optout()` for a single exploratory panel without automatically drawing examples.

The supplied `07_optout_energy_results/energy_eta_*.rds` files allow the three energy figures to be redrawn without recomputing population means:

```r
options(metatournament.eta=0.5)
source("07_optout_energy_parameter_space.R")
```

Use eta=0.3, eta=0.5 or eta=0.7 for the three current figures. To rebuild the numerical data, source `07_prepare_optout_energy.R` and call `prepare_optout_energy()`. Its arguments control the beta grid, talent-population size, allocation grid, residual tolerance, threads and budget. The current eta=0.3, 0.5 and 0.7 figures use 61 beta values per axis, 4096 exchange-symmetric QMC talent pairs and 301 investment points; their base responsive residual tolerance is 0.0001. A quadratic investment grid resolves small entries near zero; exact full-budget complementary allocations are also evaluated. Local quadratic refinement improves continuous allocations before population means are accumulated. This grid engine is checked against the separate continuous solver, and against a larger talent population. All plotted cells must have a resolved responsive residual. The file `07_energy_accuracy_check.rds` records the resolution checks.

In these plots, A--D compare the absolute energy imbalance E1-E2, including changes in total tournament participation. The base pink--white--blue palette is mixed with black according to the self-directed share. This preserves white for balanced tournament investments when no energy is self-directed.

`07_optout_endpoint_integration.R` replaces certified endpoint-only regions with direct one-dimensional normal-density integration. It checks beta_t*phi_t >= 1 and (E_t/M_t)*phi_t*beta_t <= eta before using the 0-or-B reduction; the remaining regions use the general allocation solver. This prevents narrow classification artefacts caused by counting discrete sample winners near equal tournament parameters. It also gives an independent check on near-total withdrawal at eta=0.8.

Cells missing the tighter tolerance are retried individually with larger QMC populations by `07_repair_optout_energy.R`. The retries use at least 401 investment points. Their population size is recorded in `refinement_population`; the main grid still uses 4096 pairs. Run `source("07_verify_optout_endpoint.R")` for independent global-optimum and exact-probability checks of the direct-integration reduction.


The eta=0.5 figure now has additional high-resolution calculations around the B/C line and the small high-phi corner regions. To rebuild this refinement, source `07_refine_optout_energy_boundaries.R` and call `refine_optout_energy_boundaries()`. It uses 32768 talent pairs at the selected cells, increasing the population if needed. The independent checks use up to 65536 pairs and 601 investment points. Run `source("07_verify_optout_energy_boundaries.R"); verify_optout_energy_boundaries()` to check the saved diagnostics, or pass `recalculate=TRUE` to repeat them. Corner regions persist as category B; they are retained in the plot.

The eta=0.3 and eta=0.7 variants use the same layout, colour scale and high-resolution boundary refinement as eta=0.5. Rebuild them with `prepare_optout_energy(eta_values=c(.3,.7))`, followed by `refine_optout_energy_boundaries(eta=.3)` and `refine_optout_energy_boundaries(eta=.7)`. Boundary diagnostics are saved with eta in the filename so variants do not overwrite one another. The separate independent validation file for eta=0.5 remains unchanged.

Endpoint-only conditions are checked separately for the first response and the responsive equilibrium. The responsive population can therefore be directly integrated even when the fixed reference still requires the general allocation solver. `endpoint_exact` marks responsive direct integration; `endpoint_first_exact` marks first-response direct integration.
