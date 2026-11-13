# Reproducing the energy-allocation phi-space figure

This technical note describes the original working files. In a source-only
checkout, generate the RDS and CSV inputs with the relevant stage of
[reproduce.R](reproduce.R) before using a redraw command.

Run `source("05_energy_allocation_phispace.R")` with this folder as the working directory.
The figure uses continuous phi1 and phi2 axes from 0.10 to 3.10, and a 3 by 3
array of beta1,beta2 combinations from 0.10, 0.50 and 0.90 in each population-reference block.
All other model parameters, the colour scale, A-D definitions, layout and export dimensions
match `05_energy_allocation_parameter_space.R`. PNG exports use 600 dpi.

The supplied `05_energy_allocation_phispace_inputs.rds` and
`05_energy_allocation_phispace_refined_inputs.rds` let readers redraw without a compiler.
To recompute, run `05_prepare_energy_allocation_phispace_inputs.R` followed by
`05_refine_energy_allocation_phispace_boundaries.R` (Rcpp and Rtools required).
It uses `05_energy_allocation_population_phispace.cpp`,
`05_energy_allocation_quadrature.R`, and `05_random_allocation_reference.R`.
Row checkpoints in `energy_allocation_phispace_checkpoints` and
`energy_allocation_phispace_refined_checkpoints` make both calculations resumable.
Change `phis`, `panel_betas`, or the quadrature settings in the preparation script to change
the numerical calculation; its checkpoint signature detects such changes automatically.
The default base grid has 31 phi values (0.10 steps), including 1.00 exactly.
The adaptive refinement uses a 61-value lattice (0.05 steps); every cell crossed
by an A-D boundary has directly recalculated corners. Other interpolated
background nodes are flagged by `direct` in the RDS and `direct_high_accuracy` in the refined CSV.

The first-response block starts from random allocations, E1=E2=B/2,
and calculates M1 and M2 separately. Its colours show the mean investments AFTER
actors optimize against that reference. The equilibrium block solves the same
population consistency equation as the beta-space figure. As phi now ranges farther,
the phi-space engine expands its root bracket when necessary; the utility model is unchanged.
Every retained node must pass the equilibrium residual check.

The plotted vector boundaries use continuous inequality contours between calculated nodes.
Interpolation creates no new population calculations; use a finer `phis` grid for resolution checks.

`helper_05_verify_energy_phispace.R` compares representative cases against doubled
quadrature orders and a finer allocation search, and checks tournament-exchange symmetry.
Its saved results are in `05_energy_allocation_phispace_accuracy_check.rds`.

Run helper_05_verify_energy_phispace_grid.R to check the saved grid's domain, exchange symmetry and direct-node coverage at every detected region boundary.

Numerical update: the main plotting script now applies 05_energy_allocation_phispace_population_patch.rds, constructed by 05_correct_energy_phispace_population.R from the corrected Script 06 population reference. It repairs the mixed beta=0.50,0.90 panel (and its mirror) using lower-phi outer integration. The original RDS inputs remain starting meshes; the plot assembles the corrected panel automatically. helper_05_verify_energy_phispace_patch.R checks direct boundary coverage and agreement with the coexistence data.
