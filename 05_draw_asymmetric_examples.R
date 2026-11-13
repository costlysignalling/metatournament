# Source from the scripts folder to draw both figures with identical examples.
# Edit 05_asymmetric_examples.R to choose the columns.
# Change these export switches if desired; PNG exports use 600 dpi.
export_png <- T
export_pdf <- T
source('05_asymmetric_tournaments_figure_HDR.R')
source('05_asymmetric_tournaments_utility_contours.R')
strategy <- plot_asymmetric_hdr(png=export_png,pdf=export_pdf)
utility <- plot_asymmetric_utility(population=strategy$population,
                                  png=export_png,pdf=export_pdf)
# The R plotting history contains both figures. Redraw either on demand:
# strategy$draw()
# utility$draw()
