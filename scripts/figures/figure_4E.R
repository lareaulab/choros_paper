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

fixed_training <- readLines(file.path(data_dir, "fixedLinker_fixedPrimer",
                                      "training_set.txt"))
random_training <- readLines(file.path(data_dir, "randomLinker_randomPrimer",
                                       "training_set.txt"))
all_training <- unique(c(fixed_training, random_training))

cts_by_codon_filtered <- subset(cts_by_codon,
                                !(cts_by_codon$transcript %in% all_training))
cts_by_codon_filtered$transcript <- droplevels(cts_by_codon_filtered$transcript)

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

save(corr_by_transcript, file=file.path(data_dir, "corr_by_transcript.Rda"))

# generate plot -----------------------------------------------------------

figure_4E <- ggplot(corr_by_transcript, aes(x=raw_corr, y=corrected_corr)) +
  geom_point(size=0.1, alpha=0.2) + geom_abline(slope=1, intercept=0, col="blue") +
  theme_classic(base_size=6) +
  xlab("correlation of raw counts") + ylab("correlation of corrected counts")

ggsave(filename=file.path(figures_dir, "figure_4E.pdf"),
       plot=figure_4E, device="pdf", width=0.75, height=0.75)
