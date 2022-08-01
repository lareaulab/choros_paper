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

figure_5C <- plot_bias(gc_codon_corr) +
  facet_grid(~"Corrected counts") + theme_classic(base_size=6) +
  theme(legend.position="none") +
  xlab("codon position")

ggsave(filename=file.path(figures_dir, "figure_5C.pdf"),
       plot=figure_5C, device="pdf", width=2, height=2, units="in")
