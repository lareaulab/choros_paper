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
dataset <- "lecanda_2016"
expts <- setNames(c("fixedLinker_fixedPrimer", "randomLinker_randomPrimer"),
                  c("fixed", "random"))
dataset_dir <- file.path(here(), "data", dataset)
offsets_fname <- file.path(dataset_dir, "Asite_rules.txt")
results_dirs <- sapply(expts, function(expt) file.path(dataset_dir, expt))
for(tmp_dir in results_dirs) { if(!dir.exists(tmp_dir)) { dir.create(tmp_dir) }}
bam_alignment_fnames <- sapply(expts,
                               function(expt) {
                                 file.path(dataset_dir, "raw_data", expt,
                                           grep("transcript.bam",
                                                list.files(file.path(dataset_dir, "raw_data", expt)),
                                                value=T))
                               })

# choros parameters
min_prop <- 0.9
f5_length <- 3
f3_length <- 3
num_genes <- 250
min_coverage <- 5
min_nonzero <- 100
regression_model <- formula(count ~ transcript + A + P + E + d5*f5 + d3*f3 + gc)

# generate diagnostic plot ------------------------------------------------

diagnostic_plot_objs <- sapply(expts, function(expt) { paste0(expt, "_diagnostic_plot") })
diagnostic_plot_fnames <- sapply(names(expts),
                                 function(expt) {
                                   file.path(results_dirs[expt],
                                             paste0(diagnostic_plot_objs[expt], ".Rda"))
                                 })
for(expt in names(expts)) {
  if(!file.exists(diagnostic_plot_fnames[expt])) {
    assign(diagnostic_plot_objs[expt],
           plot_diagnostic(bam_alignment_fnames[expt], transcript_length_fname))
    save(list=diagnostic_plot_objs[expt], file=diagnostic_plot_fnames[expt])
  }
}

# load footprint alignments -----------------------------------------------

bam_objs <- sapply(expts, function(expt) { paste0(expt, "_bam") })
bam_fnames <- sapply(names(expts),
                     function(expt) {
                       file.path(results_dirs[expt],
                                 paste0(bam_objs[expt], ".Rda"))
                     })
for(expt in names(expts)) {
  if(!file.exists(bam_fnames[expt])) {
    assign(bam_objs[expt],
           load_bam(bam_alignment_fnames[expt], transcript_fa_fname,
                    transcript_length_fname, offsets_fname,
                    f5_length=f5_length, f3_length=f3_length))
    save(list=bam_objs[expt], file=bam_fnames[expt])
  } else {
    load(bam_fnames[expt])
  }
}

# compute size/frame subsets ----------------------------------------------

d5_d3_objs <- sapply(expts, function(expt) { paste0(expt, "_d5_d3") })
d5_d3_fnames <- sapply(names(expts),
                       function(expt) {
                         file.path(results_dirs[expt],
                                   paste0(expts[expt], "_d5_d3.Rda"))
                       })
subsets_objs <- sapply(expts, function(expt) { paste0(expt, "_subsets") })
subsets_fnames <- sapply(names(expts),
                         function(expt) {
                           file.path(results_dirs[expt],
                                     paste0(expts[expt], "_subsets.Rda"))
                         })

for(expt in names(expts)) {
  if(!file.exists(subsets_fnames[expt])) {
    assign(d5_d3_objs[expt], count_d5_d3(get(bam_objs[expt])))
    save(list=d5_d3_objs[expt], file=d5_d3_fnames[expt])
    assign(subsets_objs[expt], choose_subsets(get(d5_d3_objs[expt])))
    save(list=subsets_objs[expt], file=subsets_fnames[expt])
  } else {
    load(subsets_fnames[expt])
  }
}

# choose transcripts for training/test sets -------------------------------

# per transcript: calculate mean codon coverage and number nonzero codon positions
transcript_coverage <- lapply(names(expts),
                              function(expt) {
                                tmp <- calculate_transcript_density(get(bam_objs[expt]),
                                                                    transcript_length_fname)
                                tmp <- tmp[match(transcript_lengths$transcript,
                                                 names(tmp))]
                                tmp[is.na(tmp)] <- 0
                                return(tmp)
                              })
transcript_coverage <- rowMeans(do.call(cbind, transcript_coverage))
transcript_num_nonzero <- lapply(names(expts),
                                 function(expt) {
                                   tmp <- count_nonzero_codons(get(bam_objs[expt]))
                                   tmp <- tmp[match(transcript_lengths$transcript,
                                                    names(tmp))]
                                   tmp[is.na(tmp)] <- 0
                                   return(tmp)
                                 })
transcript_num_nonzero <- rowMeans(do.call(cbind, transcript_num_nonzero))

# check if number good transcripts sufficient for training
good_transcripts <- data.frame(transcript=as.character(transcript_lengths$transcript),
                               coverage=transcript_coverage,
                               num_nonzero=transcript_num_nonzero)
good_transcripts <- subset(good_transcripts,
                           coverage > min_coverage & num_nonzero > min_nonzero)

# identify paralogs
paralogs <- readLines(paralogs_fname)
paralogs <- strsplit(paralogs, split="\t")
paralogs <- subset(paralogs, lengths(paralogs) > 1)
paralogs <- data.frame(transcript=unlist(paralogs),
                       group=unlist(mapply(rep, seq(length(paralogs)), times=lengths(paralogs))))
good_transcripts$paralog_group <- paralogs$group[match(good_transcripts$transcript,
                                                       paralogs$transcript)]

# take top transcript by read coverage per paralog group
paralog_groups <- split(good_transcripts, good_transcripts$paralog_group)
paralog_groups <- lapply(paralog_groups,
                         function(x) {
                           if(nrow(x) == 1) {
                             return(x)
                           } else {
                             x <- x[order(x$coverage, decreasing=T),]
                             return(x[1, ])
                           }
                         })
paralog_groups <- do.call(rbind, paralog_groups)
good_transcripts <- rbind(subset(good_transcripts,
                                 is.na(good_transcripts$paralog_group)),
                          paralog_groups)

# take top n=2*num_genes for training/test sets
good_transcripts <- good_transcripts[order(good_transcripts$coverage, decreasing=T),]
if(nrow(good_transcripts) > 2*num_genes) {
  good_transcripts <- good_transcripts[seq(2*num_genes),]
} else {
  num_genes <- round(nrow(good_transcripts)/2)
  print(paste("Insufficient genes for training; new number of genes for training:", num_genes))
}

# split transcripts evenly into training/test sets
training_set <- sample.int(nrow(good_transcripts), size=num_genes)
training_set <- as.character(good_transcripts$transcript)[training_set]
test_set <- subset(good_transcripts, !(as.character(transcript) %in% training_set))
test_set <- as.character(test_set$transcript)

writeLines(training_set, con=file.path(dataset_dir, "training_set.txt"))
writeLines(test_set, con=file.path(dataset_dir, "test_set.txt"))

q(save="no")