# Equal-talent 50:50 crossings in phi space

This technical note describes the original working files. In a source-only
checkout, generate the RDS and CSV inputs with the relevant stage of
[reproduce.R](reproduce.R) before using a redraw command.

Run `05_equal_talent_5050_phispace.R` from the scripts folder. It writes
`Figure_equal_talent_5050_crossings_phispace.png` (600 dpi) and the matching PDF.
It uses the same panel dimensions, colours and legend layout as the beta-plane
figure. Phi1 and phi2 span 0.10 to 3.10; the three panel beta values are
0.10, 0.50 and 0.90. The two blocks are the random-allocation reference
(E1=E2=B/2 with separate M1,M2) and the responsive population equilibrium.
The 50th-percentile talent-sum contour retains the original definition.

The phi plane can contain TWO positive crossings. Small black dots mark these
regions. In all crossing regions, the colour, saturation and gradient talent
sum refer to the LAST crossing. Thus the single-crossing case has the same
interpretation as the original figure. In a dotted region the preference first
reverses, then reverses back; use the CSV for both crossing talent sums.
Coloured diagonal hatching marks a stationary last crossing where equal
allocation is not globally optimal. It uses the colour of the better specialist.
It does not assert that the first crossing is globally optimal; both checks are
included in the CSV.

The crossing calculation is an analytical reduction followed by scalar root
finding, rather than talent sampling. For equal talent a, write Gt=(1+B/2)^beta_t-1.
The log marginal-return ratio is

    log(R(a)) = c + (phi1-1)*log(1+a*G1) - (phi2-1)*log(1+a*G2),
    c = log(rho) + log(phi1*beta1/(phi2*beta2)) + (beta1-beta2)*log(1+B/2).

Its derivative has the sign of

    (phi1-1)*G1 - (phi2-1)*G2 + (phi1-phi2)*G1*G2*a.

The denominator is positive. There is at most one positive turning point and
therefore at most two positive crossings. The C++ helper brackets each monotone
piece separately in log talent. It also handles coincident curves and zero or
one crossing. The global-optimum check compares the equal split against the
full investment interval, including endpoints and refined interior candidates.

Population ratios come from `06_coexistence_phi_inputs.rds`, the validated
continuous-population data used by the newer coexistence figures. This file has
both direct and previously interpolated nodes. To avoid importing interpolation
ripples, the new script rebuilds the interpolation from DIRECT nodes in the
log ratio relative to the exact zero-talent cutoff. For beta1=beta2=0.10 the
smooth direct base mesh permits cubic interpolation; its error against all
additional direct nodes is below 8e-8 in log ratio for the supplied data.
Other panels preserve every direct evaluation and interpolate the same residual
hierarchically. These are numerical population approximations; the analytical
crossing reduction does not make the population equilibrium itself analytical.

Files needed to rebuild the geometry:

- `05_equal_talent_5050_phispace.R`
- `05_equal_talent_phi_geometry.R` and `.cpp`
- `05_energy_allocation_population_phispace.cpp`
- `05_random_allocation_reference.R`
- `06_coexistence_phi_inputs.rds`

The first build needs Rcpp and Rtools; plotting also needs sf and isoband.
The supplied `05_equal_talent_5050_phispace_geometry_cache.rds` permits redraws
without compilation. Input and geometry hashes protect against stale caches.
Set `RECOMPUTE_5050_PHISPACE=1` in the environment to force rebuilding.
`gradient_legend_drop` moves the bottom legend row as in the original script.

`equal_talent_5050_crossings_phispace.csv.gz` contains all coarse-grid values,
including both crossing talents, crossing count and both global-optimum checks.
The verification scripts compare against the existing one-crossing helper,
independent dense root searches, dense investment searches, tournament exchange
symmetry and directly evaluated population nodes.
