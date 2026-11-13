# Reproducing the Script 06 coexistence figures

This technical note describes the original working files. In a source-only
checkout, generate the RDS and CSV inputs with the relevant stage of
[reproduce.R](reproduce.R) before using a redraw command.

The two maintained plotting scripts are:

- `06_strategy_coexistence_parameter_space_fixed_and_responsive_examples.R`
- `06_strategy_coexistence_phi_space_by_beta_fixed_and_responsive.R`

Their original layouts, colours, boundary styles, symbols, example markers and
legends are retained. PNG exports are 600 dpi with the same physical panel sizes;
matching PDFs are also produced. Both scripts show the plot interactively when sourced in R.

Fixed E1=E2 means that focal actors optimize against a population allocating
uniformly at random, independently of its talents. E1=E2=B/2, while M1 and M2
are calculated separately for the two training/returns functions. Thus the
utility scaler ratio is M2/M1, generally not 1. Colours and coexistence refer
to the resulting first talent-informed responses.

Responsive E_t uses the self-consistent equilibrium of the continuous lognormal
population, with both E_t and M_t responding to allocations. The population
integrals use the deterministic marginal-cutoff / strategy-switch and specialist-
departure splitting routines from the new energy figures, rather than a finite
QMC population. The population equilibrium uses the full talent distribution;
the 99% HDR restricts only the subsequent strategy-coexistence test.

The HDR and four presence-score definitions are unchanged: 2max, 2-leaning,
1-leaning and 1max. Specialists use global endpoint comparisons, and leaning
categories use the global crossover across B/2. The zero tolerances are retained
(0 for specialists and 0.0001 log-talent overlap for leaning diversifiers).
The half-interval search now also refines peaks immediately next to interval
endpoints; otherwise a near-50:50 optimum can be missed. Exact curvature and
monotone-entry shortcuts accelerate the same global decision problem.

When entry jumps directly to B/2 or farther, the leaning-diversifier interval
is absent. Its original width score is then flat at zero, which makes an
interpolated tiny-positive contour follow mesh edges. On such slices the new
engine extends the score negatively using log(b_global/b_left), where b_left is
the cheapest entry talent through allocations up to B/2. This preserves the
positive-width presence criterion and all positive slice scores, while providing
a signed absent side for contour extraction. This is a mathematical absence
margin, not a spline or an image filter.

The supplied `06_coexistence_beta_inputs.rds` and `06_coexistence_phi_inputs.rds`
permit redraws without compiling C++. To rebuild, delete the relevant input RDS
or run `prepare_coexistence('beta', force=TRUE)` / `prepare_coexistence('phi', force=TRUE)`
after sourcing `06_prepare_coexistence_deterministic.R`. Rcpp, Rtools on Windows,
and the existing energy-input datasets are needed for rebuilding. The common
preparation script calls `06_coexistence_hdr_scores.cpp`,
`05_energy_allocation_population_phispace.cpp`, `05_random_allocation_reference.R`,
and `05_energy_allocation_quadrature.R`. Row checkpoints are stored separately
for each plane in `coexistence_deterministic_cache`.

The beta plane is refined from 51 to 101 grid values, with phi panels 0.66, 1.00,
1.50. The phi plane uses beta panels 0.10, 0.50, 0.90 and is refined from 31
through 61 to 121 axis values over 0.10 to 3.10. Direct calculations are required
at every detected presence-boundary cell's corners; background interpolation is
flagged by `direct` in each panel. Diversifier threshold crossings are then bisected directly along parameter-grid
edges using the actual HDR score, separately calculated random-reference moments,
and (in the responsive case) an interpolated equilibrium as a predictor. Both
endpoints of every final bracket are checked with freshly solved population
equilibria; correction and full-model bisection are used until the actual
threshold is bracketed. Nine bisections bound the along-edge location error by half a grid step
/512. The resulting line network is polygonized, and each polygon's class is
recalculated directly at an interior point. This removes contour-interpolation
steps without filtering or visually smoothing the image. The original plotting
styles and region-symbol placement are retained.

`helper_06_verify_deterministic_scores.R` compares the accelerated scores with
a separately compiled global-search reference (with its endpoint defect fixed),
checks analytic linear-spoils cases, and increases the internal search resolution.
The reference source is `helper_06_reference_hdr_scores.cpp`, and results are saved
in `06_coexistence_score_verification.rds`.

The older adaptive/QMC/fallback Script 06 variants remain historical versions;
the two maintained scripts above use the corrected fixed-E condition and the
shared deterministic engine. Manuscript theory will be revised separately.

helper_06_verify_coexistence_grids.R checks saved-field finiteness, exact tournament-exchange symmetry and direct boundary coverage, and exports compressed score CSVs. helper_06_verify_coexistence_diagonals.R independently recalculates directly flagged diagonal scores.

The supplied `06_coexistence_beta_geometry.rds` and
`06_coexistence_phi_geometry.rds` hold these refined vector boundaries and polygon
classes. `06_refine_coexistence_curves.R` rebuilds them using
`06_coexistence_curve_evaluations.cpp` (which uses the HDR score implementation).
Per-curve bisection checkpoints are kept in `coexistence_curve_checkpoints`.

For the difficult mixed-beta phi panel (beta1=0.50, beta2=0.90), the lower-phi
activity is used as the outer integration coordinate. This avoids the poor
convergence found in the opposite direction at very asymmetric phi values.
Self-symmetric panels use the same exchange convention in their field generation
and curve evaluation. `helper_06_verify_orientation_convergence.R` records checks
at quadrature orders 48, 96, 128 and 192; the later order changes in its difficult
sample cases are below 0.005% in the equilibrium scaler ratio. These are numerical
convergence checks, not a claim of arbitrary-precision evaluation.

The same corrected population integration is propagated to the companion phi-energy figure by 05_correct_energy_phispace_population.R and its cached 05_energy_allocation_phispace_population_patch.rds. Both figures agree at shared directly calculated nodes.

Validation of the delivered geometry found 0 polygon-classification differences at 32200 directly calculated beta-space nodes and 51862 phi-space nodes. The verified vector caches contain 1305 beta-space and 2225 phi-space threshold brackets. helper_06_verify_polygon_classes.R and helper_06_verify_curve_geometry.R reproduce these checks.
