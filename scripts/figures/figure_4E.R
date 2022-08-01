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

for(expt in expts) {
  load(file.path(data_dir, expt, paste0(expt, "_bam.Rda")))
}

# aggregate footprint counts ----------------------------------------------

# codon counts
scer_lengths$num_codons <- with(scer_lengths, cds_length/3)

cts_by_codon <- data.frame(transcript=unlist(mapply(rep, scer_lengths$transcript,
                                                    scer_lengths$num_codons,
                                                    simplify=F)),
                           cod_idx=unlist(lapply(scer_lengths$num_codons, seq.int)),
                           row.names=NULL)

expts_codon_cts <- lapply(expts,
                          function(expt) {
                            tmp <- aggregate(cbind(count, correct_250) ~ transcript + cod_idx,
                                             data=get(paste0(expt, "_bam")),
                                             FUN=sum, na.rm=T, na.action=na.pass)
                            colnames(tmp) <- c("transcript", "cod_idx",
                                               paste0(expt, "_raw"),
                                               paste0(expt, "_corrected"))
                            tmp <- tmp[match_rows(cts_by_codon, tmp,
                                                  c("transcript", "cod_idx")),
                                       c(3:4)]
                            return(tmp)
                          })
cts_by_codon <- cbind(cts_by_codon, do.call(cbind, expts_codon_cts))
cts_by_codon[is.na(cts_by_codon)] <- 0
colnames(cts_by_codon) <- c("transcript", "cod_idx",
                            "fixed_raw", "fixed_corrected",
                            "random_raw", "random_corrected")

# transcript counts
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

cts_by_codon_filtered <- subset(cts_by_codon,
                                !(cts_by_codon$transcript %in% all_training))
cts_by_codon_filtered$transcript <- droplevels(cts_by_codon_filtered$transcript)
scer_lengths_filtered <- subset(scer_lengths,
                                !(scer_lengths$transcript %in% all_training))

# compute correlation by transcript ---------------------------------------

corr_by_transcript <- lapply(split(cts_by_codon_filtered,
                                   cts_by_codon_filtered$transcript),
                             function(x) {
                               raw_corr <- cor(x$fixed_raw, x$random_raw)
                               corrected_corr <- cor(x$fixed_corrected, x$random_corrected)
                               return(c(raw_corr, corrected_corr))
                             })
corr_by_transcript <- data.frame(do.call(rbind, corr_by_transcript))
colnames(corr_by_transcript) <- c("raw_corr", "corrected_corr")

scer_lengths_filtered$delta_corr <- with(corr_by_transcript,
                                         corrected_corr - raw_corr)
scer_lengths_filtered$mean_raw <- rowMeans(scer_lengths_filtered[, c("fixed_raw", "random_raw")])
scer_lengths_filtered$mean_raw <- with(scer_lengths_filtered, mean_raw/num_codons)

save(corr_by_transcript, file=file.path(data_dir, "corr_by_transcript.Rda"))

# generate plot -----------------------------------------------------------

figure_4E_left <- ggplot(corr_by_transcript, aes(x=raw_corr, y=corrected_corr)) +
  geom_point(size=0.05, alpha=0.05) +
  geom_abline(slope=1, intercept=0, col="blue", size=0.2) +
  theme_classic(base_size=6) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  xlab(expression(rho*"(raw counts)")) + ylab(expression(rho*"(corrected counts)"))

figure_4E_right <- ggplot(scer_lengths_filtered, aes(x=mean_raw, y=delta_corr)) +
  geom_point(size=0.05, alpha=0.05) +
  theme_classic(base_size=6) + geom_hline(yintercept=0, color="grey25", size=0.2) +
  geom_smooth(method="loess", formula=y~x, size=0.2, se=F) +
  scale_x_log10() + coord_cartesian(xlim=c(0.5, 1200)) +
  xlab("mean coverage") + ylab(expression(Delta*"(correlation)")) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

figure_4E <- figure_4E_left + figure_4E_right

ggsave(filename=file.path(figures_dir, "figure_4E.pdf"),
       plot=figure_4E, device="pdf", width=1.5, height=1)
