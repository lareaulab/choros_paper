rm(list=ls())

library(here)
library(choros)

ref_dir <- file.path(here(), "reference_data")
scer_fa_fname <- file.path(ref_dir, "scer.transcripts.20cds20.fa")
scer_lengths_fname <- file.path(ref_dir, "scer.transcripts.20cds20.lengths.txt")
scer_lengths <- load_lengths(scer_lengths_fname)
scer_lengths$num_codons <- scer_lengths$cds_length / 3

dataset_dir <- file.path(here(), "data", "chou_2017")
offsets_fname <- file.path(dataset_dir, "Asite_rules.txt")
samples <- grep("SRR_", invert=T, value=T,
                grep("^SRR", list.files(file.path(dataset_dir, "raw_data")), value=T))

results_dir <- file.path(dataset_dir, "processed_data")
if(!dir.exists(results_dir)) { dir.create(results_dir) }


# generate diagnostic plot + A-site offset rules --------------------------

for(x in samples) {
  raw_bam_fname <- file.path(dataset_dir, "raw_data", x, paste0(x, ".transcript.bam"))
  diagnostic_plot_obj <- paste0(x, "_diagnostic_plot")
  diagnostic_plot_fname <- file.path(results_dir, paste0(diagnostic_plot_obj, ".Rda"))
  if(!file.exists(diagnostic_plot_fname)) {
    assign(diagnostic_plot_obj,
           plot_diagnostic(raw_bam_fname, scer_lengths_fname, plot_title=x))
    save(list=diagnostic_plot_obj, file=diagnostic_plot_fname)
  }
}

rm(list=grep("diagnostic_plot", ls(), value=T))

# aggregate counts --------------------------------------------------------

codon_cts <- data.frame(transcript=unlist(mapply(rep,
                                                 x=scer_lengths$transcript,
                                                 times=scer_lengths$num_codons)),
                        cod_idx=unlist(lapply(scer_lengths$num_codons, seq_len)))

for(x in samples) {
  print(paste(x, "loading bam file", sep=": "))
  raw_bam_fname <- file.path(dataset_dir, "raw_data", x, paste0(x, ".transcript.bam"))
  bam_obj <- paste0(x, "_bam")
  bam_fname <- file.path(results_dir, paste0(bam_obj, ".Rda"))
  if(!file.exists(bam_fname)) {
    assign(bam_obj,
           load_bam(raw_bam_fname, scer_fa_fname,
                    scer_lengths_fname, offsets_fname))
    save(list=bam_obj, file=bam_fname)
  } else {
    load(bam_fname)
  }
  print(paste(x, "counting footprints", sep=": "))
  codon_cts <- within(codon_cts,
                      assign(x,
                             count_footprints(get(bam_obj), codon_cts,
                                              features=c("transcript", "cod_idx"))))
  rm(list=bam_obj)
}

save(codon_cts, file=file.path(results_dir, "codon_cts.Rda"))