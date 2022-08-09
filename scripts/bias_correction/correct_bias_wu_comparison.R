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
dataset_dir <- file.path(here(), "data", dataset)
comparison <- "CHXTIG_WT_v_3AT"
results_dir <- file.path(dataset_dir, comparison)
if(!dir.exists(results_dir)) { dir.create(results_dir) }
expts <- c("CHXTIG_WT", "CHXTIG_3AT")

# choros parameters
min_prop <- 0.9
f5_length <- 3
f3_length <- 3
num_genes <- 250
min_coverage <- 5
min_nonzero <- 100
regression_model <- formula(count ~ transcript + A*expt + P + E + d5*f5 + d3*f3 + gc)

# initialize training data for regression ---------------------------------

for(expt in expts) { 
  load(file.path(dataset_dir, expt, paste0(expt, "_training.Rda")))
}

training_obj <- paste0(comparison, "_training")
training_fname <- file.path(results_dir, paste0(training_obj, ".Rda"))
assign(training_obj,
       rbind(within(CHXTIG_WT_training, expt <- "WT"), 
             within(CHXTIG_3AT_training, expt <- "3AT")))
assign(training_obj,
       within(get(training_obj),
              expt <- factor(expt, levels=c("WT", "3AT"))))

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

q(save="no")
