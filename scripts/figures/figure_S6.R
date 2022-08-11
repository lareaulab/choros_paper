rm(list=ls())

library(here)
library(choros)

data_dir <- file.path(here(), "data")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

load(file.path(data_dir, "meydan_2020", "disome", "disome_diagnostic_plot.Rda"))

# generate plot -----------------------------------------------------------

disome_diagnostic_plot <- disome_diagnostic_plot &
  theme_classic(base_size=6)

ggsave(filename=file.path(figures_dir, "figure_S7.pdf"),
       plot=disome_diagnostic_plot, device="pdf", width=9, height=6, units="in")
