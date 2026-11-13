# Expanded participation figure

This is a separate version of Figure S5. It adds phi=1.25 and extends available tournament counts from T=2 through T=10. The fixed total budget is B=50; beta=0.4,0.8; phi=0.66,1,1.25,1.5. Every series uses filled dots, with colour indicating T. The legend is slightly higher than in the original figure. The existing Supplement and its figure are preserved.

The plotting script is `04_known_talents_multiple_tournaments_T10.R` in the parent scripts folder. It reads this folder's `participation_many_summary.csv`. Change figure dimensions, fonts, colours or legend placement there without recomputing probabilities.

`compute_T10.R` calculates the additions using the separate Rcpp source `known_talents_participation.cpp`, then combines them with `participation_T7_baseline.csv`. The retained baseline uses four deterministic shifts with 262144 points each, plus nested 65536-point checks. `summarize_T10.R` produces the plotted averages and resolution comparison. Both scripts can run from this folder or the main scripts folder. Rcpp and Rtools are installed on this computer.

The original T<=7 values are retained unchanged. Only phi=1.25 and T=8,9,10 are newly integrated. The same lognormal talent distribution, deterministic quadrature rule and global optimization approach are used throughout. Numerical zeros mean that no corresponding strategy region was resolved at this integration resolution, not that its probability has been proved to be exactly zero.

The maximum change between the nested 65536-point and 262144-point rules, averaged across four shifts, was 0.0949 percentage points. The largest range between individual fine-resolution shifted rules was 0.2511 percentage points. These are observed convergence diagnostics, not rigorous bounds on integration error.

For beta=0.8,phi=1.25, the curvature of each tournament's utility has the sign of alpha-1. A strongest talent at least 1 makes specialization globally optimal by convexity and superadditivity; when every talent is below 1, the problem is strictly concave. For beta=0.4,phi=1.25, the utility is concave unless alpha>6, and the largest possible inflection point is below 0.578. The strongest talent receives at least B/T>=5, so it is on its concave branch. The existing global branch comparison therefore extends to this case. At beta=0.4,phi=1.5, the largest inflection remains below 1.756, so the same argument remains valid through T=10. The nonconcave solver is specialized to these parameter combinations at B=50 and T<=10.

`verify_T10.R` checks the added cases against an independent dynamic program, including ten-tournament examples with high, low and mixed talents. The thirteen checks passed. Integration resolution diagnostics are recorded in `participation_many_convergence.csv`; the completion check also verifies that all strategy probabilities sum to one and that all retained baseline values are unchanged.
