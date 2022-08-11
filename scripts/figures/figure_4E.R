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

figure_4E <- ggplot(subset(corr_by_transcript, good & mean_raw > 50),
                    aes(x=raw_corr, y=corrected_corr, color=which_set)) +
  geom_abline(slope=1, intercept=0, col="blue", alpha=0.5) +
  geom_point(stroke=0) + theme_classic(base_size=8) + labs(col="") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  xlab(expression(rho*"(raw counts)")) + ylab(expression(rho*"(corrected counts)")) +
  scale_color_manual(values=alpha(setNames(c("orange", "purple", "grey45"),
                                     c("training", "test", "other")), 0.5))

ggsave(filename=file.path(figures_dir, "figure_4E.pdf"),
       plot=figure_4E, device="pdf", width=3, height=2)
