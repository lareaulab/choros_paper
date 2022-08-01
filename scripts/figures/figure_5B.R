rm(list=ls())

library(here)
library(choros)

data_dir <- file.path(here(), "data")
figures_dir <- file.path(here(), "figures")

# load data ---------------------------------------------------------------

ottr_dir <- "~/footprint-bias/expts/ferguson_2021/ottrUMI_Asite_f5_3_f3_3"

load(file.path(ottr_dir, "ottrUMI_f5_3_training_codon_corr.Rda"))
load(file.path(ottr_dir, "test_gc", "ottrUMI_f5_3_gc_corr.Rda"))

# make plot ---------------------------------------------------------------

ottr_max <- max(c(ottrUMI_f5_3_training_codon_corr[, "uncorrected"],
                  gc_codon_corr))

figure_5B <- plot_bias(ottrUMI_f5_3_training_codon_corr[, "uncorrected"]) +
  facet_grid(~"Raw counts") + theme_classic(base_size=6) +
  theme(legend.position="none") + coord_cartesian(ylim=c(0, ottr_max)) +
  xlab("codon position")

ggsave(filename=file.path(figures_dir, "figure_5B.pdf"),
       plot=figure_5B, device="pdf", width=2, height=2, units="in")
