# Metatournament reproduction code

Run the commands below from the repository root. Generated plots, numerical
checkpoints, and compiled libraries are created locally.

## Set up R

Use R with a working C++ toolchain for Rcpp. On Windows, use the Rtools
release that matches your R version. On macOS or Linux, install your system's
R-compatible compiler and the system libraries required by sf if installing
it from source. Install the R packages in your normal user library:

    install.packages(c("Rcpp", "sf", "isoband", "shiny", "jsonlite",
                       "viridisLite", "viridis"))

The two offline viewer builders use Python 3 and its standard library.

## Regenerate the analysis

Start with a fresh checkout for a from-scratch calculation. To see whether the
basic R plotting path works before starting the larger calculations:

    Rscript reproduce.R foundations

The full entry point is:

    Rscript reproduce.R

Individual stages can be requested by name. Run dependent stages in the
default order recorded in reproduce.R. For instance, the energy calculations
provide the population inputs for coexistence and the convincing grid.
Participation feeds the combined multiple-tournament figure. The full run can
take hours. Per-step logs go to reproduction_logs. Some expensive stages
save checkpoints and resume from successfully written cells. Keep a fresh
checkout for an independent rerun.

Set the environment variable METATOURNAMENT_CAMPAIGN_WORKERS to a positive
integer to parallelize the campaign grid. Each worker runs a separate R
process and compiles its own engine.

The figures.txt file lists the required manuscript images. A full run checks
that each one was produced. The numerical assertions and convergence checks
live in the analysis scripts.

The clean rebuild produced all 70 manuscript figures. Of these, 57 were
byte-identical to the working exports. The remaining comparisons showed small
rendering differences, a layout difference in one boredom example, and
numerical sensitivity near energy-region boundaries. The beta-space energy
classification agreed at every grid point. Across the three opt-out energy
figures, 56 of 100,467 panel grid points changed classification near a
boundary. The rebuilt opt-out grids conserve the total budget and all cells
meet the solver's residual tolerance of 0.0001.

## Experiment with the model

From the project root in R, source the experiment you want to run.

| Experiment | Entry point |
| --- | --- |
| Two asymmetric tournaments | 05_explore_asymmetric_tournaments.R |
| Individual decisions in many tournaments | 06B_predict_many_tournaments.R |
| Talent slices | 06B_explore_talent_slices.R |
| Predictions from an observed pair of tournaments | 06C_observed_pair_prediction.R |
| Opt-out decisions | 07_pilot_asymmetric_tournaments_optout.R |
| Boredom | 08_explore_boredom.R |
| Convincing | 10_explore_convincing.R |

The chapter-specific reproduction notes explain the individual generators and
resolution controls. When a note refers to a saved numerical output, generate
it with its preparation script before running a redraw.

## Interactive figures

After the many-tournament stage, run python 06B_build_cube_viewer.py to
produce the offline rotatable cube. For the harvesting comparison, run
the following in R after the harvesting stage:

    source("09_training_comparison.R")
    draw_training_comparison()

Then run
python 09_build_training_comparison.py.

Open those generated HTML files locally, or publish them as static files via
GitHub Pages. A project Pages site can serve files from a publishing branch
or folder under https://costlysignalling.github.io/metatournament/. The two
builders produce self-contained pages.

For a cube that recalculates when parameters change, run the many-tournament
stage to build 06B_cube_web/default_cube.rds, then call
shiny::runApp("06B_cube_web") from R. Deploy that app on a Shiny-capable host
to share a public live link.
