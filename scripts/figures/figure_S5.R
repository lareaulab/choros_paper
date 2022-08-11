rm(list=ls())

library(here)
library(choros)
library(ggplot2)

data_dir <- file.path(here(), "data", "lecanda_2016")
figures_dir <- file.path(here(), "figures")

expts <- c("fixedLinker_fixedPrimer", "randomLinker_randomPrimer")
names(expts) <- c("fixed linker\nfixed primer", "random linker\nrandomprimer")

scer_lengths_fname <- file.path(here(), "reference_data", "scer.transcripts.20cds20.lengths.txt")
scer_lengths <- load_lengths(scer_lengths_fname)

# load data ---------------------------------------------------------------

load(file.path(data_dir, "corr_by_transcript.Rda"))

# generate plot -----------------------------------------------------------

figure_S5 <- ggplot(corr_by_transcript, aes(x=mean_raw+0.001, y=delta_corr, col=which_set)) +
  theme_classic(base_size=8) + geom_hline(yintercept=0, color="grey25") +
  # geom_smooth(method="loess", formula=y~x, se=F, size=0.5) +
  scale_x_log10() + # coord_cartesian(xlim=c(0.5, 1200)) +
  geom_point(size=0.75, alpha=0.5, stroke=0) + labs(col="") +
  xlab("mean coverage") + ylab(expression(Delta*"(correlation)")) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  scale_color_manual(values=setNames(c("orange", "purple", "grey45"),
                                     c("training", "test", "other"))) +
  geom_vline(xintercept=50, col="blue", linetype="dashed", alpha=0.5) + 
  annotate(geom="text", x=55, y=-0.28, label="mean = 50", col="blue", 
           alpha=0.5, size=2, hjust=0)

ggsave(filename=file.path(figures_dir, "figure_S5.pdf"),
       plot=figure_S5, device="pdf", width=6.5, height=3)
