# Reproducing the participation figures

This technical note describes the original working files. In a source-only
checkout, generate the RDS and CSV inputs with the relevant stage of
[reproduce.R](reproduce.R) before using a redraw command.

The numerical calculations use **C++ called from R through Rcpp**. Python and reticulate are not needed. All required code and saved data are supplied in this scripts folder and its `Etude3_multiple_tournaments_analysis` subfolder.

Set R's working directory to this scripts folder. Redraw both figures from the saved summaries with:

```r
source("reproduce_participation_figures.R")
reproduce_participation_figures()
```

To recompute every probability from the mathematical model, then regenerate the summaries, convergence tables and figures:

```r
# install.packages("Rcpp")  # only if needed
source("reproduce_participation_figures.R")
reproduce_participation_figures(recompute=TRUE)
```

Full recomputation requires R, Rcpp and a working C++17 compiler (Rtools for the installed R version on Windows; an appropriate compiler on macOS/Linux). It can take hours. It overwrites raw results, summaries, convergence tables and figures; retain a copy first if comparing versions. Use `figures="fixed"` or `figures="proportional"` to select only one figure. Sourcing the helper only defines the function; it does not start computation.

The settings reproduce the supplied figures: independent lognormal talents with median 1 and original-scale SD 0.5; beta 0.4 and 0.8; phi 0.66, 1, 1.25 and 1.5; T from 2 to 10. Fixed budget is B=50; proportional budget is B=25T. The engine uses four deterministically shifted Halton rules, 262144 points each, with nested 65536-point results for convergence comparisons. The source code implements analytic/concave solutions where applicable and global comparisons in nonconcave cases.

## Included files

- `04_known_talents_multiple_tournaments_T10.R`: fixed-budget plotting script.
- `04_known_talents_multiple_tournaments_proportional_budget.R`: proportional-budget plotting script.
- `reproduce_participation_figures.R`: portable R entry point for full calculation and plotting.
- `Etude3_multiple_tournaments_analysis/analytic_participation.R`: independent one-dimensional analytic benchmarks for phi=1, T=2 or 3.
- `Etude3_multiple_tournaments_analysis/T10_phi125/`: fixed-budget C++ engine, `compute_T10.R`, `summarize_T10.R`, `verify_T10.R`, raw and summary CSVs, convergence table and README.
- `Etude3_multiple_tournaments_analysis/proportional_budget/`: proportional-budget C++ engine, `compute_proportional_budget.R`, `summarize_proportional_budget.R`, `verify_proportional_budget.R`, raw and summary CSVs, convergence table and README.
- The analysis folder also contains the earlier T<=7 engine, calculation and validation scripts and results used in developing the figures.

The original `compute_*.R` scripts update an existing calculation incrementally using the supplied baseline CSVs. The new entry point instead recomputes **all** probabilities; baseline files are retained only for the summarizers' consistency checks. Each `verify_*.R` script provides separate checks against dynamic programming; run it from its own calculation folder. The proportional summarizer additionally compares T=3, phi=1 results against the analytic benchmark. These convergence and benchmark checks are numerical evidence, not rigorous bounds on integration error.
