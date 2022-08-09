rm(list=ls())

library(here)
library(choros)

set.seed(10)

# reference files
ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
transcript_lengths <- load_lengths(transcript_length_fname)
paralogs_fname <- file.path(ref_dir, "scer_100_80.paralogs.id.txt")

# experiment-specific files
dataset <- "wu_2019"
expt <- "CHXTIG_WT"
dataset_dir <- file.path(here(), "data", dataset)
offsets_fname <- file.path(dataset_dir, "Asite_rules_shortmers.txt")
results_dir <- file.path(dataset_dir, expt)
if(!dir.exists(results_dir)) { dir.create(results_dir) }
bam_alignment_fname <- file.path(dataset_dir, "raw_data", "WT_CHXTIG",
                                 grep("transcript.bam",
                                      list.files(file.path(dataset_dir, "raw_data", "WT_CHXTIG")),
                                      value=T))

# choros parameters
min_prop <- 0.9
f5_length <- 3
f3_length <- 3
num_genes <- 250
min_coverage <- 5
min_nonzero <- 100
regression_model <- formula(count ~ transcript + A + P + E + d5*f5 + d3*f3 + gc)

# load footprint alignments -----------------------------------------------

bam_obj <- paste0(expt, "_bam")
bam_fname <- file.path(results_dir, paste0(bam_obj, ".Rda"))
if(!file.exists(bam_fname)) {
  assign(bam_obj,
         load_bam(bam_alignment_fname, transcript_fa_fname,
                  transcript_length_fname, offsets_fname,
                  f5_length=f5_length, f3_length=f3_length))
  save(list=bam_obj, file=bam_fname)
} else {
  load(bam_fname)
}

# compute size/frame subsets ----------------------------------------------

d5_d3_obj <- paste0(expt, "_d5_d3")
d5_d3_fname <- file.path(results_dir, paste0(expt, "_d5_d3.Rda"))
subsets_obj <- paste0(expt, "_subsets")
subsets_fname <- file.path(results_dir, paste0(expt, "_subsets.Rda"))
if(!file.exists(subsets_fname)) {
  assign(d5_d3_obj, count_d5_d3(get(bam_obj)))
  save(list=d5_d3_obj, file=d5_d3_fname)
  assign(subsets_obj, choose_subsets(get(d5_d3_obj)))
  save(list=subsets_obj, file=subsets_fname)
} else {
  load(subsets_fname)
}

# choose transcripts for training/test sets -------------------------------

training_set <- readLines(file.path(results_dir, "training_set.txt"))
test_set <- readLines(file.path(results_dir, "test_set.txt"))

# initialize training data for regression ---------------------------------

training_obj <- paste0(expt, "_training")
training_fname <- file.path(results_dir, paste0(training_obj, ".Rda"))
if(!file.exists(training_fname)) {
  assign(training_obj,
         init_data(transcript_fa_fname, transcript_length_fname,
                   d5_d3_subsets=get(subsets_obj),
                   f5_length=f5_length, f3_length=f3_length,
                   which_transcripts=training_set))
  assign(training_obj,
         within(get(training_obj),
                count <- count_footprints(get(bam_obj), get(training_obj), "count")))
  save(list=training_obj, file=training_fname)
} else {
  load(training_fname)
}

# compute regression ------------------------------------------------------

fit_obj <- paste0(expt, "_fit")
fit_fname <- file.path(results_dir, paste0(fit_obj, ".Rda"))
coef_obj <- paste0(expt, "_coef")
coef_fname <- file.path(results_dir, paste0(coef_obj, ".Rda"))
if(!file.exists(coef_fname)) {
  assign(fit_obj, glm.nb(regression_model, data=get(training_obj)))
  save(list=fit_obj, file=fit_fname)
  assign(coef_obj, parse_coefs(get(fit_obj)))
  save(list=coef_obj, file=coef_fname)
} else {
  load(coef_fname)
}

# correct counts ----------------------------------------------------------

correction_name <- paste0("correct_", num_genes)

assign(bam_obj,
       within(get(bam_obj),
              assign(correction_name,
                     correct_bias(get(bam_obj), fit_coefs=get(coef_obj)))))
save(list=bam_obj, file=bam_fname)

assign(training_obj,
       within(get(training_obj),
              assign(correction_name,
                     count_footprints(get(bam_obj), get(training_obj),
                                      correction_name))))
save(list=training_obj, file=training_fname)

# evaluate bias -----------------------------------------------------------

to_evaluate <- c("count", correction_name)
test_data <- subset(get(bam_obj), transcript %in% test_set)

# calculate codon position importance
codon_corr_obj <- paste0(expt, "_codon_corr")
codon_corr_fname <- file.path(results_dir, paste0(codon_corr_obj, ".Rda"))
if(!file.exists(codon_corr_fname)) {
  assign(codon_corr_obj,
         lapply(to_evaluate,
                function(x) {
                  evaluate_bias(test_data, which_column=x,
                                transcript_fa_fname, transcript_length_fname,
                                type="codon")
                }))
  save(list=codon_corr_obj, file=codon_corr_fname)
}

# calulate nt position importance
nt_corr_obj <- paste0(expt, "_nt_corr")
nt_corr_fname <- file.path(results_dir, paste0(nt_corr_obj, ".Rda"))
if(!file.exists(nt_corr_fname)) {
  assign(nt_corr_obj,
         lapply(to_evaluate,
                function(x) {
                  evaluate_bias(test_data, which_column=x,
                                transcript_fa_fname, transcript_length_fname,
                                type="nt")
                }))
  save(list=nt_corr_obj, file=nt_corr_fname)
}

q(save="no")
