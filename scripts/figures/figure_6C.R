rm(list=ls())

library(here)
library(choros)

data_dir <- file.path(here(), "data", "meydan_2020", "monosome")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

load(file.path(data_dir, "monosome_codon_corr.Rda"))

# make plot ---------------------------------------------------------------

plot_max <- max(unlist(monosome_codon_corr))

figure_6C <- plot_bias(monosome_codon_corr[[2]]) +
  facet_grid(~"Corrected counts") + theme_classic(base_size=6) +
  theme(legend.position="none") + coord_cartesian(ylim=c(0, plot_max)) + 
  xlab("codon position")

ggsave(filename=file.path(figures_dir, "figure_6C.pdf"),
       plot=figure_6C, device="pdf", width=2, height=2, units="in")
