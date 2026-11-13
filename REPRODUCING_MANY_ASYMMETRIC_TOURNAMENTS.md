# Multiple asymmetric tournaments

This technical note describes the original working files. In a source-only
checkout, generate the RDS and CSV inputs with the relevant stage of
[reproduce.R](reproduce.R) before using a redraw command.

The additions are **Supplement 6.9 Multiple asymmetric tournaments** and
**Supplement 6.10 Empirical tests when some tournaments are unobserved**.

Scripts beginning with `06B_` cover the generalization in Supplement 6.9;
`06C_` scripts cover the observed-pair application in Supplement 6.10.
They share the supplied `06B_many_results` data folder.

## Predict one actor and draw a pie chart

Run R from this scripts folder. Install `Rcpp` if needed; compiling its supplied
C++ helper requires Rtools on Windows. The helper is freshly compiled in R's
local temporary folder once per session, rather than loading a DLL from a
synchronized project cache. If Windows still blocks the newly compiled DLL,
the application-control policy needs review by the computer administrator. 

```r
source('06B_predict_many_tournaments.R')
example <- many_example(beta=c(.6,.8,.4), phi=c(.66,1,1.2),
                        B=50, actor=c(1,1,1))
```

For your own inputs:

```r
beta <- c(.6,.4,.4)
phi <- c(.66,1,1.2)
B <- 50
actor <- c(1.2,.8,1.5)
meanlog <- c(0,0,0)
sdlog <- rep(sqrt(log((1+sqrt(2))/2)),3)
population <- many_lognormal_population(16384,3,meanlog,diag(sdlog^2))
reference <- many_random_reference(beta,phi,B,meanlog,sdlog)
first <- many_allocate(actor,beta,phi,B,reference$scaler)
responsive <- many_responsive(population,beta,phi,B,start=reference$scaler)
stopifnot(responsive$converged)
equilibrium <- many_allocate(actor,beta,phi,B,responsive$scaler)
first$allocation
equilibrium$allocation
many_pie(first,'First response')
many_pie(equilibrium,'Responsive population')
```

`Sigma` in `many_lognormal_population()` specifies the covariance of **log
talents**. Joint dependence matters. For other distributions, replace
`population` with a matrix: one actor per row, one tournament per column.
`many_responsive()` accepts optional sampling weights and either one common
budget or a vector of actor-specific budgets. For a different initial
allocation rule, obtain the reference scalers with
`many_means(population,initial_X,beta,phi,weights)$scaler`.
The helper `many_random_reference()` assumes a constant budget, lognormal
marginal talents, and allocations uniform over the simplex independently of
talents. It does **not** assume equal M across tournaments.

The utility model uses `0 < beta < 1`, `phi > 0`, and nonnegative talents and
budgets. There is no fixed maximum number of tournaments in the allocation
engine. The included Halton population generator has 20 dimensions; supply
your own population matrix beyond that.

## Explore talent slices like Figure S16

```r
source('06B_explore_talent_slices.R')
many_slices(beta=c(.9,.6,.6),phi=c(.8,1.5,1),B=50)
many_slices(beta=c(.6,.8,.4,.8),phi=c(.66,1,1.2,1.4),B=50,other_talents=1.5)
```

The axes are talents 1 and 2; columns fix talent 3 at its 10th, 30th, 50th,
70th and 90th reference percentiles by default. `cuts=c(.5,1,2)` supplies
specific talent-3 values. For four or more tournaments, `other_talents`
fixes talents 4 onward (all 1 by default). Both population references still
include every tournament and the full population talent distribution.

The default `n_grid=61` is for quick exploration; `n_grid=181` gives finer
boundaries. Concave slices use analytical specialist polygons. Nonconcave
slices trace the specialist cutoffs by repeated global allocation checks;
the same vector paths supply fills and outlines. The display grid sets
resolution along a boundary; bisection refines its crossing location. Yellow curves indicate equal positive investment in jointly
leading tournaments, not necessarily half the total budget.

Plots appear in R. Set `png=TRUE` or `pdf=TRUE` to save under the preview
filename (`file` changes it); manuscript figures are not overwritten.
`preview <- many_slices(...)` returns the numerical allocations in
`preview$slices[[row]][[column]]$allocation`. If `example` already contains a
`many_example()` result, `many_slices(result=example)` reuses its population
solution; changing only cuts or the display grid does not require solving
that population again. A changed beta, phi or budget needs a new solution.

## Accuracy and equilibrium interpretation

- Actor output includes `utility_gap_bound`, `relative_gap_bound`, and
  `within_tolerance`. A bound certifies utility accuracy, not necessarily
  closeness of allocations near a tie. Nonconcave cases use global budget-grid
  dynamic programming. Raise `max_grid` or explicitly choose a different
  `tolerance` when needed. The code never equates a local maximum with the
  global solution. Concave cases use a tangent bound and marginal equalization.
  Nonconcave grid solutions are refined continuously and compared with
  refined specialist starts. Both the grid interpolation bound and a separable
  Lagrange-multiplier bound remain available; local refinement is not used as
  a substitute for the global accuracy check. Concave grid convolutions use
  an exact faster search over monotone maximizing indices.
- `many_responsive()` refuses to use actor solutions outside its requested
  utility tolerance. It reports its equilibrium residual and convergence flag.
  Broad finite differences first handle population strategy switches; smaller
  differences resolve the local root. A damped consistency update is available
  when a Newton step stalls. The convergence tolerance is unchanged.
  A failed solve is not an equilibrium. Inactive E=M=0 tournaments need a
  separate entry analysis. Alternative starting scalers can reveal multiple
  roots; root finding is not a proof of existence or uniqueness.
- Population integration is deterministic Halton quadrature over joint
  talents. Increase population size and compare the results; a tiny consistency
  residual by itself says nothing about integration accuracy. For the plotted
  example we checked 4,096, 16,384 and 65,536 talent vectors. The last increase
  changed mean investments by at most 0.019 budget units (B=50).
- Check the fixed reference separately by increasing `nquad` in
  `many_random_reference()`. Orders 64 and 128 agree for the plotted example;
  extreme parameters can require higher orders. Increasing the responsive
  population size does not check this separate quadrature.
- The nonconcave parameter screens use a larger permitted actor utility gap
  (0.001). One screening scenario did not converge and is explicitly retained
  as a failed solve. No claims or plots treat it as equilibrium.

## Reproduce the figures and checks

Run these R scripts in order:

```r
source('06B_verify_many_tournaments.R')
source('06B_refine_many_population.R')
source('06B_many_tournament_figures.R')
source('06C_check_conditional_equilibrium.R')
source('06B_check_solver_edges.R')
source('06B_validate_final.R')
```

Additional feedback examples:

```r
source('06B_check_feedback_examples.R')
source('06B_check_reduced_imbalance.R')
```

Outputs include 600-dpi PNGs, vector PDFs, the underlying population means and
checks in `06B_many_results`, and cube mesh data. The C++ helper is
`06B_many_tournaments_core.cpp`. Source files are the complete calculation;
no other project scripts are needed for these new results.
The figure script additionally uses `jsonlite` to export the HTML viewer's
mesh data. It is already installed in the project environment.

The cube uses marching tetrahedra with exact linear clipping of the surfaces
at three-way ties within each tetrahedron. Specialist surfaces use analytical
marginal-entry scores, since this entire displayed cube is concave. Increase
`axis3` resolution in the figure script for a finer mesh. Bright solid slice fills are reserved for complete specialists; diversifiers
use muted colours. Yellow slice curves
and gold cube surfaces mark equal **jointly leading** investments. They do not
mean that each activity gets half the full budget. The third tournament colour
is green; the first two retain their established pink and blue.

To rebuild the standalone rotatable HTML cube, run
`python 06B_build_cube_viewer.py` after the R figure script. This Python script
uses only the standard library to package the already calculated meshes; it
does not optimize allocations or calculate population means. The HTML works
offline, with no external libraries or network requests. R can invoke it too:

```r
# system2('python', '06B_build_cube_viewer.py')
```

## Empirical data with unobserved outside tournaments

`06C_observed_pair_prediction.R` provides the conditional predictor and a global
utility-regret test. Its commented example documents the required CSV columns.
For an independently known fixed reference supply `reference_scaler`. For a
responsive prediction leave it unset: the routine solves for the relative
scaler from the joint observed talents and pair budgets, holding those budgets
fixed, without using the observed target splits. 

The synthetic recovery check solves that equation from scratch for 16,384
actors. It reproduces their full three-tournament within-pair allocations to
about 1.2e-7 budget units. 

## Adjustable three-tournament cube

From the scripts folder:

```r
shiny::runApp('06B_cube_web')
```

The app uses the same checked R/C++ solver and the same boundary geometry
as the manuscript figures. Change the three beta and phi sliders and press
**Recalculate cubes**. Budget is fixed at 50; the talent distribution is
independent lognormal with median 1 and SD 0.5. Preview uses 4,096 population
vectors and a 21-by-21-by-21 display grid; the finer option uses 16,384 and
a 35-by-35-by-35 grid. The bundled manuscript example uses 65,536 population
vectors. These settings control integration and surface resolution separately.

**Download HTML snapshot** saves the currently displayed, rotatable cubes.
Changing parameters needs the running Shiny app; a downloaded snapshot needs
only a browser. Calculation failures retain the last successful figure and
report the failure rather than silently substituting an unchecked result.

The folder `06B_cube_web` is self-contained for deployment to a Shiny host.
It requires R packages `shiny`, `jsonlite` and `Rcpp`, plus a C++ compiler
(Rtools on Windows). The app compiles its bundled engine into R's temporary
folder; no DLL from another computer is distributed. After rebuilding the
manuscript figures, run `source('06B_prepare_cube_web.R')` to refresh its
bundled engine and default cube. See the app folder's README for hosting.
