rm(list=ls())

library(here)
library(choros)
library(ggplot2)
library(patchwork)

data_dir <- file.path(here(), "data", "lecanda_2016")
figures_dir <- file.path(here(), "figures")

expts <- c("fixedLinker_fixedPrimer", "randomLinker_randomPrimer")
names(expts) <- c("fixed linker\nfixed primer", "random linker\nrandomprimer")

scer_lengths_fname <- file.path(here(), "reference_data", "scer.transcripts.20cds20.lengths.txt")
scer_lengths <- load_lengths(scer_lengths_fname)

# load data ---------------------------------------------------------------

for(expt in expts) {
  load(file.path(data_dir, expt, paste0(expt, "_bam.Rda")))
}

load(file.path(data_dir, "corr_by_transcript.Rda"))

# aggregate footprint counts ----------------------------------------------

cts_by_transcript <- lapply(expts,
                            function(expt) {
                              tmp <- aggregate(cbind(count, correct_250) ~ transcript,
                                               data=get(paste0(expt, "_bam")),
                                               FUN=sum, na.rm=T, na.action=na.pass)
                              tmp$expt <- expt
                              tmp <- tmp[match(scer_lengths$transcript,
                                               tmp$transcript),]
                              return(tmp)
                            })
names(cts_by_transcript) <- expts
scer_lengths$fixed_raw <- cts_by_transcript$fixedLinker_fixedPrimer$count
scer_lengths$fixed_corrected <- cts_by_transcript$fixedLinker_fixedPrimer$correct_250
scer_lengths$random_raw <- cts_by_transcript$randomLinker_randomPrimer$count
scer_lengths$random_corrected <- cts_by_transcript$randomLinker_randomPrimer$correct_250

fixed_training <- readLines(file.path(data_dir, "fixedLinker_fixedPrimer",
                                      "training_set.txt"))
random_training <- readLines(file.path(data_dir, "randomLinker_randomPrimer",
                                       "training_set.txt"))
all_training <- unique(c(fixed_training, random_training))

scer_lengths_filtered <- subset(scer_lengths,
                                !(scer_lengths$transcript %in% all_training))

scer_lengths_filtered$delta_corr <- with(corr_by_transcript,
                                         corrected_corr - raw_corr)

scer_lengths_filtered$mean_raw <- rowMeans(scer_lengths_filtered[, c("fixed_raw", "random_raw")])
scer_lengths_filtered$mean_raw <- with(scer_lengths_filtered, mean_raw/num_codons)

# generate plot -----------------------------------------------------------

figure_4F_bottom <- ggplot(scer_lengths_filtered, aes(x=mean_raw, y=delta_corr)) +
  geom_point(size=0.05, alpha=0.05) +
  theme_classic(base_size=6) + geom_hline(yintercept=0, color="grey25", size=0.2) +
  geom_smooth(method="loess", formula=y~x, size=0.2, se=F) +
  scale_x_log10() + coord_cartesian(xlim=c(0.5, 1200)) +
  xlab("mean coverage") + ylab(expression(Delta*"(correlation)")) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
figure_4F_top <- ggplot(scer_lengths_filtered, aes(x=mean_raw)) +
  geom_density(fill="grey") + scale_x_log10() + coord_cartesian(xlim=c(0.5, 1200)) +
  xlab("") + ylab("density") + theme_classic(base_size=6)

figure_4F <- figure_4F_top / figure_4F_bottom +
  plot_layout(heights=c(0.2, 0.8))

ggsave(filename=file.path(figures_dir, "figure_4F.pdf"),
       plot=figure_4F_bottom, device="pdf", width=0.75, height=0.75)
