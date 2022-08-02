rm(list=ls())

library(here)
library(choros)

data_dir <- file.path(here(), "data")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

load(file.path(data_dir, "tunney_2018", "tunney", "tunney_diagnostic_plot.Rda"))

# generate plot -----------------------------------------------------------

ggsave(filename=file.path(figures_dir, "figure_S1.pdf"),
       plot=tunney_diagnostic_plot, device="pdf", width=6.5, units="in")
