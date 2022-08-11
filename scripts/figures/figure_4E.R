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

training_set <- readLines(file.path(data_dir, "training_set.txt"))
test_set <- readLines(file.path(data_dir, "test_set.txt"))

# subset to transcripts with enough coverage
min_coverage <- 50
min_nonzero <- 100
bam_objs <- sapply(expts, function(expt) { paste0(expt, "_bam") })
transcript_coverage <- lapply(names(expts),
                              function(expt) {
                                tmp <- calculate_transcript_density(get(bam_objs[expt]),
                                                                    scer_lengths_fname)
                                tmp <- tmp[match(scer_lengths$transcript,
                                                 names(tmp))]
                                tmp[is.na(tmp)] <- 0
                                return(tmp)
                              })
transcript_coverage <- rowMeans(do.call(cbind, transcript_coverage))
transcript_num_nonzero <- lapply(names(expts),
                                 function(expt) {
                                   tmp <- count_nonzero_codons(get(bam_objs[expt]))
                                   tmp <- tmp[match(scer_lengths$transcript,
                                                    names(tmp))]
                                   tmp[is.na(tmp)] <- 0
                                   return(tmp)
                                 })
transcript_num_nonzero <- rowMeans(do.call(cbind, transcript_num_nonzero))
good_transcripts <- data.frame(transcript=as.character(scer_lengths$transcript),
                               coverage=transcript_coverage,
                               num_nonzero=transcript_num_nonzero)
good_transcripts <- subset(good_transcripts,
                           coverage > min_coverage & num_nonzero > min_nonzero)

# compute correlation by transcript ---------------------------------------

corr_by_transcript <- lapply(split(cts_by_codon,
                                   cts_by_codon$transcript),
                             function(x) {
                               raw_corr <- cor(x$fixed_raw, x$random_raw,
                                               method="spearman")
                               corrected_corr <- cor(x$fixed_corrected, x$random_corrected,
                                                     method="spearman")
                               return(c(raw_corr, corrected_corr))
                             })
corr_by_transcript <- data.frame(do.call(rbind, corr_by_transcript))
colnames(corr_by_transcript) <- c("raw_corr", "corrected_corr")
corr_by_transcript$transcript <- rownames(corr_by_transcript)
corr_by_transcript$good <- corr_by_transcript$transcript %in% as.character(good_transcripts$transcript)
corr_by_transcript$which_set <- sapply(corr_by_transcript$transcript,
                                       function(x) {
                                         ifelse(x %in% training_set, "training",
                                                ifelse(x %in% test_set,
                                                       "test", "other"))
                                       })
corr_by_transcript <- plyr::join(corr_by_transcript, scer_lengths,
                                 by="transcript")
corr_by_transcript$delta_corr <- with(corr_by_transcript, corrected_corr - raw_corr)
corr_by_transcript$mean_raw <- rowMeans(corr_by_transcript[, c("fixed_raw", "random_raw")])
corr_by_transcript$mean_raw <- with(corr_by_transcript, mean_raw/num_codons)

save(corr_by_transcript,
     file=file.path(data_dir, "corr_by_transcript.Rda"))

# generate plot -----------------------------------------------------------

figure_4E <- ggplot(subset(corr_by_transcript, good & mean_raw > 50),
                    aes(x=raw_corr, y=corrected_corr, col=which_set)) +
  geom_abline(slope=1, intercept=0, col="blue") +
  geom_point(size=0.5) + # facet_grid(~factor(which_set, levels=c("training", "test"))) +
  theme_classic(base_size=8) + labs(col="") +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
  xlab(expression(rho*"(raw counts)")) + ylab(expression(rho*"(corrected counts)")) +
  scale_color_manual(values=alpha(setNames(c("orange", "purple", "grey45"),
                                     c("training", "test", "other")), 0.5))

ggsave(filename=file.path(figures_dir, "figure_4E.pdf"),
       plot=figure_4E, device="pdf", width=2.5, height=2)
