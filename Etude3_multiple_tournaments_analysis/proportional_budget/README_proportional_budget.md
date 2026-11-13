# Proportional budget participation figure

This calculation makes the second participation figure in Supplement 5.6.3 (Figure S6). It uses B_T=25T, with beta=0.4,0.8; phi=0.66,1,1.25,1.5; and T=2,...,10. Thus B_2=50 and B_10=250. Its plotting conventions match the fixed-budget T10 figure.

Run `04_known_talents_multiple_tournaments_proportional_budget.R` from the main scripts folder to redraw the figure from cached probabilities. Adjust fonts and layout there without repeating the integration.

To recompute probabilities, run `compute_proportional_budget.R` and `summarize_proportional_budget.R` from this folder or the main scripts folder. Rcpp and a working C++ toolchain are required and are installed on this computer. The separate source `known_talents_participation.cpp` uses the same global strategy classification as the fixed-budget calculation, at budgets B=25T. The T=2 probabilities are retained unchanged from the fixed-budget calculation because its B is already 50.

All other probabilities use four deterministically shifted Halton rules with 262144 points each, plus nested 65536-point comparisons. The largest change in an averaged probability when increasing resolution fourfold was 0.07372 percentage points. The largest range between fine-resolution shifted rules was 0.20409 percentage points. These are convergence diagnostics, not rigorous error bounds. The phi=1,T=3,B=75 probabilities agreed with independent one-dimensional normal integrals within 0.0031 percentage points.

The global allocation classification was independently checked against simplex dynamic programming in thirteen high-, low- and mixed-talent cases, including ten tournaments at B=250. `verify_proportional_budget.R` reproduces those checks. The nonconcave branch logic is validated for the displayed parameter combinations and B=25T. Changing the budget constant requires changing the numerical engine and rerunning these checks, not only editing the plotting script's label.
