rm(list=ls())

library(here)
library(choros)

figures_dir <- file.path(here(), "figures")

# reference files
ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
transcript_lengths <- load_lengths(transcript_length_fname)
paralogs_fname <- file.path(ref_dir, "scer_100_80.paralogs.id.txt")

# data files
raw_data_dir <- file.path(here(), "data", "chou_2017", "raw_data")
sample_names <- grep("SRR", list.files(raw_data_dir), value=T)
sample_names <- sample_names[!(sample_names=="SRR_accessions")]
offsets_fname <- file.path(here(), "data", "chou_2017", "Asite_rules.txt")

# define functions --------------------------------------------------------

## taken from ~/archive/mean_variance/mean_variance_plots.Rmd

read_cts_by_codon <- function(fname) {
  ### read in cts_by_codon table
  # fname: character; file path to cts_by_codon table
  raw_dat <- readLines(fname)
  transcript_names <- sapply(raw_dat,
                             function(x) {
                               strsplit(x, split="\t")[[1]][1]
                             })
  codon_counts <- lapply(raw_dat,
                         function(x) {
                           as.numeric(strsplit(x, split="\t")[[1]][-1])
                         })
  names(codon_counts) <- transcript_names
  return(codon_counts)
}

sum_transcripts <- function(codon_cts) {
  ### return data frame of counts per transcript (row: transcript; col: sample)
  # codon_cts: list of per-codon footprint counts
  transcripts <- names(codon_cts[[1]])
  transcript_cts <- sapply(codon_cts,
                           function(x) {
                             sapply(transcripts,
                                    function(y) {
                                      sum(x[[match(y, names(x))]])
                                    })
                           })
  colnames(transcript_cts) <- names(codon_cts)
  return(transcript_cts)
}

reshape_cts_by_codon <- function(codon_counts_list) {
  ### convert list codon counts lists into data.frame
  # codon_counts_list: named list of cts_by_codon tables per experiment
  expt_names <- names(codon_counts_list)
  transcript_names <- names(codon_counts_list[[1]])
  transcript_lengths <- lengths(codon_counts_list[[1]])
  codon_counts <- data.frame(sapply(codon_counts_list, unlist))
  rownames(codon_counts) <- paste(unlist(mapply(rep.int, transcript_names, transcript_lengths)),
                                  unlist(lapply(transcript_lengths, seq)), sep="_")
  colnames(codon_counts) <- expt_names
  return(codon_counts)
}

compute_stats <- function(counts_table) {
  ### compute mean and variance of counts per codon
  # counts_table: data.frame of codon counts (row) and expt (col)
  stats <- data.frame(mean=rowMeans(counts_table),
                      var=apply(counts_table, 1, var))
  stats$disp <- with(stats, (var-mean)/mean^2)
  rownames(stats) <- rownames(counts_table)
  return(stats)
}

plot_mean_var <- function(stats, plot_title=NULL) {
  ### return mean-variance plots
  # stats: data.frame with 3 columns: mean, variance, dispersion
  # plot_title: character
  ggplot(stats, aes(mean, var)) + geom_point(alpha=0.1) +
    theme_classic(base_size=6) + geom_smooth(colour="red", method="gam", formula=y~s(x, bs="cs")) +
    geom_smooth(method=lm, se=T) + geom_abline(intercept=0, slope=1, colour="green") +
    ggtitle(plot_title) + xlab("mean") + ylab("variance")
}

make_plots_mean_var <- function(stats, plot_title=NULL, x_max=100) {
  ### return mean-variance plots
  # stats: data.frame with 3 columns: mean, variance, dispersion
  # plot_title: character
  # x_max: numeric; maximum for x-axis
  plot_all <- plot_mean_var(stats, plot_title)
  plot_subset <- plot_mean_var(subset(stats, mean<x_max),
                               paste("mean count <", x_max))
  gridExtra::grid.arrange(plot_all, plot_subset, ncol=2)
}

plot_mean_disp <- function(stats, plot_title=NULL) {
  ### return mean-dispersion plots
  # stats: data.frame with 3 columns: mean, variance, dispersion
  # plot_title: character
  ggplot(stats, aes(mean, disp)) + geom_point(alpha=0.1) +
    theme_classic(base_size=6) + geom_smooth(colour="red", method="gam", formula=y~s(x, bs="cs")) +
    geom_smooth(method=lm, se=T) +
    ggtitle(plot_title) + xlab("mean") + ylab("dispersion")
}

cts_fname <- "cts_by_codon.size.27.31.txt"

# generate plot -----------------------------------------------------------

expt_dir <- "~/archive/mean_variance/chou_2017"
reps <- paste0("ribo_", 1:14)[-10]
codon_counts <- lapply(reps,
                       function(x) {
                         read_cts_by_codon(file.path(expt_dir, x, "process", cts_fname))
                       })
names(codon_counts) <- reps
transcript_counts <- sum_transcripts(codon_counts)
codon_counts <- reshape_cts_by_codon(codon_counts)
chou_all_codons <- compute_stats(codon_counts)
chou_all_nCodons_zeroCounts <- sum(with(chou_all_codons, mean==0 & var==0))
print(paste0("codons with no footprint counts: ", chou_all_nCodons_zeroCounts, " (",
             round(chou_all_nCodons_zeroCounts/nrow(chou_all_codons)*100, 2), "%)"))
chou_all_transcripts <- compute_stats(transcript_counts)
chou_all_nTranscripts_zeroCounts <- sum(with(chou_all_transcripts, mean==0 & var==0))
print(paste0("transcripts with no footprint counts: ", chou_all_nTranscripts_zeroCounts, " (",
             round(chou_all_nTranscripts_zeroCounts/nrow(chou_all_transcripts)*100, 2), "%)"))

figure_S3 <- ggplot(chou_all_codons, aes(mean, var)) + theme_classic(base_size=6) +
  geom_hex(bins=100) + scale_fill_gradient(low="grey", high="blue") +
  geom_smooth(color="red", alpha=0.5, method="lm", formula=y~x+I(x^2)) +
  geom_smooth(color="green", alpha=0.5, method="lm", formula=y~x) +
  geom_smooth(color="purple", alpha=0.5, method="lm", formula=y~1) +
  xlab("mean") + ylab("variance")

ggsave(filename=file.path(figures_dir, "figure_S3.pdf"),
       plot=figure_S3, device="pdf", width=5, height=5, units="in")
