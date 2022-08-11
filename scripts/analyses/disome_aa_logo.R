rm(list=ls())

library(here)
library(choros)
library(ggplot2)
library(patchwork)
library(ggseqlogo)

dataset <- "meydan_2020"
dataset_dir <- file.path(here(), "data", dataset)

ref_dir <- file.path(here(), "reference_data")

# transcript lengths
scer_lengths_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
scer_lengths <- load_lengths(scer_lengths_fname)
scer_lengths$num_codons <- with(scer_lengths, cds_length/3)

# transcript sequence
scer_seq_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
scer_seq_codons <- read_fasta_as_codons(scer_seq_fname, scer_lengths_fname)

# codon properties
scer_codon_properties_fname <- file.path(ref_dir, "yeast_codon_properties.txt")
scer_codon_properties <- read.table(scer_codon_properties_fname, header=T)

# disome alignment
disome_bam_fname <- file.path(dataset_dir, "disome", "disome_bam.Rda")
load(disome_bam_fname)

# monosome alignment
monosome_bam_fname <- file.path(dataset_dir, "monosome", "monosome_bam.Rda")
load(monosome_bam_fname)

aggregate_counts <- function(bam_data, which_column) {
  count_footprints(bam_data, cts_by_codon, which_column=which_column,
                   features=c("transcript", "cod_idx"), integer_counts=F)
}

# generate disome density profiles ----------------------------------------

cts_by_codon <- data.frame(transcript=unlist(mapply(rep, scer_lengths$transcript,
                                                    scer_lengths$num_codons,
                                                    simplify=F)),
                           cod_idx=unlist(lapply(scer_lengths$num_codons, seq.int)),
                           row.names=NULL)

cts_by_codon <- within(cts_by_codon, {
  # monosomes
  monosome_raw <- aggregate_counts(monosome_bam, "count")
  monosome_corrected <- aggregate_counts(monosome_bam, "correct_250")
  # disomes
  disome_leading_raw <- aggregate_counts(within(disome_bam,
                                                cod_idx <- cod_idx_leading),
                                         "count")
  disome_leading_corrected <- aggregate_counts(within(disome_bam,
                                                      cod_idx <- cod_idx_leading),
                                               "correct_250")
})

# remove transcripts with <10 monosome counts
no_monosome <- aggregate(monosome_raw ~ transcript, cts_by_codon, FUN=sum)
no_monosome <- subset(no_monosome, monosome_raw < 10)
no_monosome <- as.character(no_monosome$transcript)
cts_by_codon <- subset(cts_by_codon,
                       !(transcript %in% no_monosome))

# remove transcripts with <30 codon positions
short_transcripts <- subset(scer_lengths, cds_length/3 <= 30)
short_transcripts <- as.character(short_transcripts$transcript)
cts_by_codon <- subset(cts_by_codon,
                       !(transcript %in% short_transcripts))

# normalize disome counts by mean monosome density ------------------------

# cts_by_codon_normalized <- lapply(split(cts_by_codon, cts_by_codon$transcript),
#                                   function(tmp_cts) {
#                                     num_codons <- nrow(tmp_cts)
#                                     mean_monosome_raw <- mean(tmp_cts$monosome_raw[seq(11, num_codons-10)])
#                                     mean_monosome_corrected <- mean(tmp_cts$monosome_corrected[seq(11, num_codons-10)])
#                                     tmp_cts <- within(tmp_cts, {
#                                       disome_norm_raw <- disome_leading_raw / mean_monosome_raw
#                                       disome_norm_corrected <- disome_leading_corrected / mean_monosome_corrected
#                                     })
#                                     return(tmp_cts)
#                                   })
cts_by_codon_normalized <- lapply(split(cts_by_codon, cts_by_codon$transcript),
                                  function(tmp_cts) {
                                    tmp_cts <- within(tmp_cts, {
                                      disome_norm_raw <- disome_leading_raw
                                      disome_norm_corrected <- disome_leading_corrected
                                    })
                                  })

# generate observed and expected values -----------------------------------

cts_by_codon_normalized <- lapply(cts_by_codon_normalized,
                                  function(tmp_cts) {
                                    num_codons <- nrow(tmp_cts)
                                    mean_raw <- mean(tmp_cts$disome_norm_raw)
                                    mean_corrected <- mean(tmp_cts$disome_norm_corrected)
                                    median_raw <- mean(tmp_cts$disome_norm_raw)
                                    median_corrected <- mean(tmp_cts$disome_norm_corrected)
                                    tmp_cts <- within(tmp_cts, {
                                      mean_raw <- mean_raw
                                      median_raw <- median_raw
                                      mean_corrected <- mean_corrected
                                      median_corrected <- median_corrected
                                      binary_mean_raw <- as.numeric(disome_leading_raw > mean_raw)
                                      binary_mean_corrected <- as.numeric(disome_leading_corrected > mean_corrected)
                                      binary_median_raw <- as.numeric(disome_leading_raw > median_raw)
                                      binary_median_corrected <- as.numeric(disome_leading_corrected > median_corrected)
                                      mean_binary_mean_raw <- mean(binary_mean_raw)
                                      mean_binary_mean_corrected <- mean(binary_mean_corrected)
                                      mean_binary_median_raw <- mean(binary_median_raw)
                                      mean_binary_median_corrected <- mean(binary_median_corrected)
                                    })
                                    return(tmp_cts)
                                  })
cts_by_codon_normalized <- do.call(rbind, cts_by_codon_normalized)

# pull amino acid sequences -----------------------------------------------

# exclude positions whose neighborhood would overlap w/ start/stop codons
cts_by_codon_normalized$num_codons <- scer_lengths$num_codons[match(cts_by_codon_normalized$transcript,
                                                                    scer_lengths$transcript)]
which_codons <- with(cts_by_codon_normalized, ((cod_idx-20) > 1) & ((cod_idx + 8) < num_codons))
cts_by_codon_normalized <- cts_by_codon_normalized[which_codons,]

# pull codons for [-20, 8] around leading codon A site
codon_neighborhood <- t(sapply(seq(nrow(cts_by_codon_normalized)),
                               function(x) {
                                 tmp_transcript <- as.character(cts_by_codon_normalized$transcript[x])
                                 tmp_cod_idx <- cts_by_codon_normalized$cod_idx[x] - 1
                                 # start codon of scer_seq_codons is "0"
                                 tmp_start <- tmp_cod_idx - 20
                                 tmp_end <- tmp_cod_idx + 8
                                 tmp_codons <- scer_seq_codons[[tmp_transcript]][as.character(seq(tmp_start, tmp_end))]
                                 return(tmp_codons)
                               }))

# convert codons to amino acids
aa_neighborhood <- sapply(seq(ncol(codon_neighborhood)),
                          function(x) {
                            tmp_codon <- as.character(codon_neighborhood[, x])
                            tmp_aa <- as.character(scer_codon_properties$aa)[match(tmp_codon,
                                                                                   as.character(scer_codon_properties$codon))]
                            return(tmp_aa)
                          })

# remove disomes that overlap non-terminal stop codons
bad_disomes <- lapply(seq(ncol(codon_neighborhood)),
                      function(x) {
                        which(codon_neighborhood[,x] %in% c("TAG", "TAA", "TGA"))
                      })
bad_disomes <- unique(unlist(bad_disomes)) # n=182
cts_by_codon_normalized <- cts_by_codon_normalized[-bad_disomes,]
codon_neighborhood <- codon_neighborhood[-bad_disomes,]
aa_neighborhood <- aa_neighborhood[-bad_disomes,]

# compute enrichment scores -----------------------------------------------

calculate_enrichment <- function(dat, neighborhood,
                                 expected, observed) {
  # dat: data.frame; data to pull counts from
  # neighborhood: data.frame; neighborhood around leading A site to calculate enrichment scores
  # expected: character; column in dat with expected counts
  # observed: character; column in dat with observed counts
  enrichment <- data.frame(neighborhood,
                           dat[, match(c(expected, observed), colnames(dat))])
  # by position, compute enrichment per amino acid
  ## enrichment: log2(sum(observed)/sum(expected))
  enrichment <- sapply(grep("X", colnames(enrichment), value=T),
                       function(x) {
                         observed <- aggregate(formula(paste(observed, x, sep="~")),
                                               enrichment, FUN=sum)
                         observed[,2] <- observed[,2] / sum(observed[,2])
                         expected <- aggregate(formula(paste(expected, x, sep="~")),
                                               enrichment, FUN=sum)
                         # expected <- plyr::count(enrichment, vars=x)
                         expected[,2] <- expected[,2] / sum(expected[,2])
                         expected <- expected[match(observed[,1], expected[,1]),]
                         enrichment <- setNames(log2(observed[,2]/expected[,2]),
                                                observed[, 1])
                         return(enrichment)
                       })
  return(enrichment)
}

# plot_enrichment <- function(enrichment) {
#   ## enrichment: matrix; position weighted enrichments for raw counts
#   ### rownames correspond to letter for logo
#   logo_plot <- ggseqlogo(enrichment, method="custom") +
#     ylab("log2 enrichment") + theme(axis.text.x=element_blank()) +
#     guides(fill="none")
#   return(logo_plot)
# }

comparisons <- rbind(c("mean", "disome_norm"),
                     c("median", "disome_norm"),
                     c("mean_binary_mean", "binary_mean"),
                     c("mean_binary_median", "binary_median"))
comparisons <- as.data.frame(comparisons)
colnames(comparisons) <- c("expected", "observed")

raw_enrichments <- lapply(seq(nrow(comparisons)),
                          function(x) {
                            tmp_expected <- paste0(comparisons$expected[x], "_raw")
                            tmp_observed <- paste0(comparisons$observed[x], "_raw")
                            calculate_enrichment(cts_by_codon_normalized,
                                                 aa_neighborhood,
                                                 tmp_expected, tmp_observed)
                          })
corrected_enrichments <- lapply(seq(nrow(comparisons)),
                                function(x) {
                                  tmp_expected <- paste0(comparisons$expected[x], "_corrected")
                                  tmp_observed <- paste0(comparisons$observed[x], "_corrected")
                                  calculate_enrichment(cts_by_codon_normalized,
                                                       aa_neighborhood,
                                                       tmp_expected, tmp_observed)
                                })

save(cts_by_codon, cts_by_codon_normalized,
     codon_neighborhood, aa_neighborhood,
     raw_enrichments, corrected_enrichments,
     file=file.path(dataset_dir, "disome", "disome_aa_logo.Rda"))

(ggseqlogo(raw_enrichments[[3]], method="custom") + ggtitle("raw counts")) /
  (ggseqlogo(corrected_enrichments[[3]], method="custom") + ggtitle("corrected counts")) /
  (ggseqlogo(corrected_enrichments[[3]] - raw_enrichments[[3]], method="custom") + ggtitle("difference"))
