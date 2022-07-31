#!/bin/bash

ref_dir=${HOME}/choros_paper/reference_data
script_dir=${HOME}/choros_paper/scripts
raw_data_dir=${HOME}/choros_paper/data/simulated_data/raw_data

cd ${raw_data_dir}

for expt in no_bias f3_bias f5_bias f3_f5_bias
do
  echo ${expt}
  cd ${expt}

  # 1. pool chunks
  echo "... pooling chunks"
  cat parts/*fq > ${expt}_raw.fq

  # 2. trim adapter
  echo "... 3' trimming adapter"
  cutadapt -a CTGTAGGCACCATCAAT -m 15 --trimmed-only -o ${expt}_trim3.fq \
  ${expt}_raw.fq > ${expt}_trim3.cutadapt

  # 3. fitler out reads aligning to rRNA
  echo "... removing rRNA reads"
  bowtie -p 10 -v 2 -S --un ${expt}_not_rRNA.fq ${ref_dir}/ScerRRNA \
  ${expt}_trim3.fq > ${expt}_rRNA.sam 2> ${expt}_rRNA.bowtiestats

  # 4. filter out reads aligning to non-coding genes
  echo "... removing ncRNA reads"
  bowtie -p 10 -v 2 -S --un ${expt}_not_rRNA_ncRNA.fq ${ref_dir}/rna_coding \
  ${expt}_not_rRNA.fq > ${expt}_ncRNA.sam 2> ${expt}_ncRNA.bowtiestats

  # 5. align to transcriptome
  echo "... aligning to transcriptome"
  bowtie -p 10 -v 2 -S --un ${expt}_unmapped.fq -a --norc \
  ${ref_dir}/scer.transcripts.20cds20 ${expt}_not_rRNA_ncRNA.fq \
  > ${expt}_footprints.sam 2> ${expt}_footprints.bowtiestats

  # 6. calculate multi-mapping weights
  echo "... calculating multimapping weights"
  rsem-calculate-expression --sam ${expt}_footprints.sam --seed-length 15 \
  ${ref_dir}/scer.transcripts.20cds20 ${expt} > ${expt}.rsem.stdout \
  2> ${expt}.rsem.stderr

  cd ..
  echo 
done
