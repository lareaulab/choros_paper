rm(list=ls())

library(here)
library(choros)

data_dir <- file.path(here(), "data")
figures_dir <- file.path(here(), "figures")

ref_dir <- file.path(here(), "reference_data")
transcript_lengths_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
transcript_lengths <- load_lengths(transcript_lengths_fname)
transcript_lengths$num_codons <- with(transcript_lengths, cds_length/3)
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_seq <- read_fasta_as_codons(transcript_fa_fname, transcript_lengths_fname)

# load data ---------------------------------------------------------------

expts <- c("tunney", "schuller", "weinberg")
names(expts) <- c("Tunney", "Schuller", "Weinberg")

load(file.path(data_dir, "tunney_2018", "tunney", "tunney_bam.Rda"))
load(file.path(data_dir, "schuller_2017", "schuller", "schuller_bam.Rda"))
load(file.path(data_dir, "weinberg_2016", "weinberg", "weinberg_bam.Rda"))

# aggregate data ----------------------------------------------------------

codon_cts <- lapply(expts,
                    function(x) {
                      within(aggregate(cbind(count, correct_250) ~ transcript + cod_idx,
                                       data=get(paste0(x, "_bam")), FUN=sum,
                                       na.rm=T, na.action=na.pass),
                             expt <- x)
                    })
codon_cts <- do.call(rbind, codon_cts)
codon_cts$expt <- factor(codon_cts$expt, levels=expts)
levels(codon_cts$expt) <- names(expts)

transcript_cts <- lapply(expts,
                         function(x) {
                           within(aggregate(cbind(count, correct_250) ~ transcript,
                                            data=get(paste0(x, "_bam")), FUN=sum,
                                            na.rm=T, na.action=na.pass),
                                  expt <- x)
                         })
transcript_cts <- do.call(rbind, transcript_cts)
transcript_cts$expt <- factor(transcript_cts$expt, levels=expts)
levels(transcript_cts) <- names(expts)

codon_cts$A <- sapply(seq(nrow(codon_cts)),
                      function(x) {
                        tmp_transcript <- as.character(codon_cts$transcript[x])
                        tmp_cod_idx <- codon_cts$cod_idx[x]
                        return(transcript_seq[[tmp_transcript]][as.character(tmp_cod_idx-1)])
                      })
Asite_cts <- lapply(expts,
                    function(x) {
                      # normalize to mean transcript count
                      tmp_cts <- subset(codon_cts, expt == x)
                      norm_cts <- lapply(levels(tmp_cts$transcript),
                                         function(y) {
                                           num_codons <- subset(transcript_lengths,
                                                                transcript==y)$num_codons
                                           if(num_codons < 40) { return(NULL) }
                                           tmp_data <- subset(tmp_cts,
                                                              transcript == y &
                                                                cod_idx %in% seq(21, num_codons-20))
                                           count_norm <- sum(tmp_data$count)/(num_codons-40)
                                           correct_250_norm <- sum(tmp_data$correct_250)/(num_codons-40)
                                           return(within(tmp_data, {
                                             count <- count / count_norm
                                             correct_250 <- correct_250 / correct_250_norm
                                           }))
                                         })
                      norm_cts <- do.call(rbind, norm_cts)
                      within(aggregate(cbind(count, correct_250) ~ A,
                                       data=norm_cts, FUN=mean, na.rm=T),
                             expt <- x)
                    })
Asite_cts <- do.call(rbind, Asite_cts)
Asite_cts$expt <- factor(Asite_cts$expt, levels=expts)
levels(Asite_cts$expt) <- names(expts)

# generate plots ----------------------------------------------------------

figure_S5A <- ggplot(codon_cts, aes(x=count, y=correct_250)) +
  geom_hex(bins=100) + geom_abline(slope=1, intercept=0, color="grey25") +
  facet_grid(expt~.) + theme_classic(base_size=6) +
  xlab("raw counts") + ylab("corrected counts") +
  scale_fill_gradient(low="grey", high="blue") +
  scale_x_log10(limits=c(0.5, 35250)) +
  scale_y_log10(limits=c(0.5, 62000))

ggsave(filename=file.path(figures_dir, "figure_S5A.pdf"),
       plot=figure_S5A, device="pdf", width=2, height=3.5, units="in")

figure_S5B <- ggplot(transcript_cts, aes(x=count, y=correct_250)) +
  geom_hex(bins=50) + geom_abline(slope=1, intercept=0, color="grey25") +
  facet_grid(expt~.) + theme_classic(base_size=6) +
  xlab("raw counts") + ylab("corrected counts") +
  scale_fill_gradient(low="grey", high="blue") +
  scale_x_log10(limits=c(0.5, 1.6e6)) + scale_y_log10(limits=c(0.5, 1.6e6))

ggsave(filename=file.path(figures_dir, "figure_S5B.pdf"),
       plot=figure_S5B, device="pdf", width=2, height=3.5, units="in")

figure_S5C <- ggplot(Asite_cts, aes(x=count, y=correct_250)) +
  geom_point(size=0.5, color="grey25") +
  geom_abline(slope=1, intercept=0, color="grey25") +
  facet_grid(expt~.) + theme_classic(base_size=6) +
  xlab("raw counts") + ylab("corrected counts") +
  scale_fill_gradient(low="grey", high="blue") +
  scale_x_log10(limits=c(0.5, 25)) + scale_y_log10(limits=c(0.5, 50))

ggsave(filename=file.path(figures_dir, "figure_S5C.pdf"),
       plot=figure_S5C, device="pdf", width=2, height=3.5, units="in")

save(codon_cts, transcript_cts, Asite_cts,
     file=file.path(here(), "data", "figure_S5_data.Rda"))
