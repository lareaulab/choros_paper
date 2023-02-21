rm(list=ls())

if(!("simRiboSeq") %in% installed.packages()) {
  devtools::install_github("amandamok/simRiboSeq", build_vignettes=T)
}

library(here)
library(simRiboSeq)
library(Rsamtools)

expt_dir <- file.path(here(), "data", "simulated_data", "raw_data")
if(!dir.exists(expt_dir)) { dir.create(expt_dir, recursive=T) }

# transcriptome: sequences
scer_fasta_fname <- file.path(here(), "reference_data",
                              "scer.transcripts.20cds20.fa")
scer_fasta <- read_fasta(scer_fasta_fname,
                         utr5_length=20, utr3_length=20)

# transcriptome: CDS lengths
scer_lengths_fname <- file.path(here(), "reference_data",
                                "scer.transcripts.20cds20.lengths.txt")
scer_lengths <- read.table(scer_lengths_fname,
                           col.names=c("transcript", "utr5_length",
                                       "cds_length", "utr3_length"))
scer_lengths <- setNames(scer_lengths$cds_length, scer_lengths$transcript)
scer_lengths <- scer_lengths[match(names(scer_fasta), names(scer_lengths))]

# footprint abundances
weinberg_bam_fname <- file.path(here(), "data", "weinberg_2016",
                                "raw_data", "weinberg",
                                "weinberg_footprints.transcript.bam")
weinberg_bam <- data.frame(scanBam(BamFile(weinberg_bam_fname),
                                   param=ScanBamParam(tag=c("ZW", "MD"),
                                                      what="rname"))[[1]])
weinberg_abundances <- aggregate(tag.ZW ~ rname, weinberg_bam, FUN=sum)
weinberg_abundances <- setNames(weinberg_abundances$tag.ZW,
                                weinberg_abundances$rname)
weinberg_abundances <- weinberg_abundances[match(names(scer_fasta), names(weinberg_abundances))]

# set simulation parameters -----------------------------------------------

min_size <- 25
max_size <- 31
num_reads <- 50e6
chunk_size <- 1e6
num_chunks <- num_reads / chunk_size
num_cores <- 20

data("d5_bias", "d3_bias", "f5_bias", "f3_bias", "codon_TE", "rt_bias",
     package="simRiboSeq")
no_f5_bias <- setNames(rep(1, length(f5_bias)), names(f5_bias))
no_f3_bias <- setNames(rep(1, length(f3_bias)), names(f3_bias))

scer_rho <- simulate_rho(scer_fasta, scer_lengths, weinberg_abundances)
scer_pi <- simulate_pi(scer_fasta, utr5_length=6, utr3_length=7, codon_TE=codon_TE)

# simulation: no ligation bias --------------------------------------------

simulation <- "no_bias"
parts_dir <- file.path(expt_dir, simulation, "parts")
dir.create(parts_dir, recursive=T)

for(i in 1:num_chunks) {
  print(paste("Part", i, "of", num_chunks))
  chunk_name <- paste0(simulation, "_part", i)
  assign(chunk_name,
         value=simulate_footprints(scer_fasta, num_ribosomes=chunk_size,
                                   rhos=scer_rho, pis=scer_pi,
                                   delta_5=d5_bias, delta_3=d3_bias,
                                   min_size=min_size, max_size=max_size,
                                   lig_bias=no_f3_bias, circ_bias=no_f5_bias,
                                   rt_bias=rt_bias, mc.cores=num_cores))
  write_footprints_fastq(get(chunk_name),
                         file.path(parts_dir, paste0(chunk_name, ".fq")))
  rm(list=chunk_name)
}

# simulation: 3' bias -----------------------------------------------------

simulation <- "f3_bias"
parts_dir <- file.path(expt_dir, simulation, "parts")
dir.create(parts_dir, recursive=T)

for(i in 1:num_chunks) {
  print(paste("Part", i, "of", num_chunks))
  chunk_name <- paste0(simulation, "_part", i)
  assign(chunk_name,
         value=simulate_footprints(scer_fasta, num_ribosomes=chunk_size,
                                   rhos=scer_rho, pis=scer_pi,
                                   delta_5=d5_bias, delta_3=d3_bias,
                                   min_size=min_size, max_size=max_size,
                                   lig_bias=f3_bias, circ_bias=no_f5_bias,
                                   rt_bias=rt_bias, mc.cores=num_cores))
  write_footprints_fastq(get(chunk_name),
                         file.path(parts_dir, paste0(chunk_name, ".fq")))
  rm(list=chunk_name)
}


# simulation: 5' bias -----------------------------------------------------

simulation <- "f5_bias"
parts_dir <- file.path(expt_dir, simulation, "parts")
dir.create(parts_dir, recursive=T)

for(i in 1:num_chunks) {
  print(paste("Part", i, "of", num_chunks))
  chunk_name <- paste0(simulation, "_part", i)
  assign(chunk_name,
         value=simulate_footprints(scer_fasta, num_ribosomes=chunk_size,
                                   rhos=scer_rho, pis=scer_pi,
                                   delta_5=d5_bias, delta_3=d3_bias,
                                   min_size=min_size, max_size=max_size,
                                   lig_bias=no_f3_bias, circ_bias=f5_bias,
                                   rt_bias=rt_bias, mc.cores=num_cores))
  write_footprints_fastq(get(chunk_name),
                         file.path(parts_dir, paste0(chunk_name, ".fq")))
  rm(list=chunk_name)
}


# simulation: 3' and 5' bias ----------------------------------------------

simulation <- "f3_f5_bias"
parts_dir <- file.path(expt_dir, simulation, "parts")
dir.create(parts_dir, recursive=T)

for(i in 1:num_chunks) {
  print(paste("Part", i, "of", num_chunks))
  chunk_name <- paste0(simulation, "_part", i)
  assign(chunk_name,
         value=simulate_footprints(scer_fasta, num_ribosomes=chunk_size,
                                   rhos=scer_rho, pis=scer_pi,
                                   delta_5=d5_bias, delta_3=d3_bias,
                                   min_size=min_size, max_size=max_size,
                                   lig_bias=f3_bias, circ_bias=f5_bias,
                                   rt_bias=rt_bias, mc.cores=num_cores))
  write_footprints_fastq(get(chunk_name),
                         file.path(parts_dir, paste0(chunk_name, ".fq")))
  rm(list=chunk_name)
}

q(save="no")