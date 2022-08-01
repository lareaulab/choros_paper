rm(list=ls())

library(here)
library(choros)

ref_dir <- file.path(here(), "reference_data")
transcript_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
transcript_length_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
transcript_lengths <- load_lengths(transcript_length_fname)

f5_length <- 2
f3_length <- 3
num_genes <- 250

data_dir <- file.path(here(), "data", "meydan_2020")
offsets_5prime_fname <- file.path(data_dir, "Asite_rules_disome_5prime.txt")
offsets_3prime_fname <- file.path(data_dir, "Asite_rules_disome_3prime.txt")

expt <- "disome"
results_dir <- file.path(data_dir, expt)
if(!dir.exists(results_dir)) { dir.create(results_dir) }
raw_bam_fname <- file.path(data_dir, "raw_data", expt,
                           paste0(expt, "_trimUMI_footprints.transcript.bam"))

coef_fname <- file.path(data_dir, "monosome", "monosome_coef.Rda")
load(coef_fname)

# read in alignments ------------------------------------------------------

bam_obj <- paste0(expt, "_bam")
bam_fname <- file.path(results_dir, paste0(bam_obj, ".Rda"))
if(!file.exists(bam_fname)) {
  assign(bam_obj,
         load_bam(raw_bam_fname, transcript_fa_fname, transcript_length_fname,
                  f5_length=f5_length, f3_length=f3_length, read_type="disome",
                  offsets_5prime_fname=offsets_5prime_fname,
                  offsets_3prime_fname=offsets_3prime_fname))
  save(list=bam_obj, file=bam_fname)
} else {
  load(bam_fname)
}

# correct counts ----------------------------------------------------------

correction_name <- paste0("correct_", num_genes)
assign(bam_obj,
       within(get(bam_obj),
              assign(correction_name,
                     correct_bias(get(bam_obj),
                                  fit_coefs=monosome_coef))))
save(list=bam_obj, file=bam_fname)

q(save="no")
