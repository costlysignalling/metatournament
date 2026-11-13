# Three tournament allocation cube

This folder is a self-contained Shiny application. It uses the same R/C++
allocation solver and mesh geometry as Supplement Figures S16 and S17.

From the parent scripts folder, run:

```r
shiny::runApp("06B_cube_web")
```

Or, from inside this folder:

```r
shiny::runApp(".")
```

The required packages are `shiny`, `jsonlite` and `Rcpp`. Compiling the
bundled helper requires a C++ compiler, including Rtools on Windows.
The helper is built locally in R's temporary folder; precompiled DLLs are
not distributed.

Change any of the six sliders and select **Recalculate cubes**. Both
population references are recalculated. The budget stays at 50; talents
are independent lognormal variables with median 1 and SD 0.5.

The default manuscript cube is precomputed with 65,536 population vectors.
New calculations use 4,096 vectors and a 21³ surface grid in Preview mode,
or 16,384 vectors and a 35³ grid in the finer mode. The full unbounded
talent distribution is used for population feedback. The displayed cube
is only the 0–2.5 viewing region. The calculation details report the
consistency residual and actor-utility error bound.

Specialist boundaries are traced with the global allocation solver.
Gold surfaces show equal positive investments in two jointly leading
tournaments, rather than necessarily half the budget in either one.
Dragging rotates both cubes; visibility and opacity controls do not
recompute the model.

**Download HTML snapshot** saves the currently displayed rotatable figure.
A snapshot opens directly in a browser. Adjustable parameters require the
running app and its R backend.

## Hosting

Deploy this entire directory to an R Shiny host, retaining `engine/` and
`www/`. There are no parent-folder dependencies or external JavaScript
CDN dependencies. R packages and a compiler must be available on the host.
Keep the server's normal access controls when publishing an app.

For example, after configuring your own hosting account:

```r
rsconnect::deployApp("06B_cube_web") # from the parent scripts folder
```

Hosting instructions: [Posit Shiny deployment documentation](https://docs.posit.co/connect/user/shiny/index.html).

This project has only been prepared for deployment; no public hosting
account or deployment is created by these files.

After changing the manuscript's model engine or rebuilding the manuscript
cube, run `source("06B_prepare_cube_web.R")` from the parent scripts folder
to refresh the bundled sources and default data.
