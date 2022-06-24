#!/bin/bash

ref_dir=${HOME}/choros_paper/reference_data
raw_data_dir=${HOME}/choros_paper/data/chou_2017/raw_data

cd ${raw_data_dir}
# sample_names.txt copied from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100626
# SraRunTable.txt downloaded from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA392312&o=acc_s%3Aa

# pull SRR accessions for WT ribosome profiling samples
grep Ribo_BY4741_ sample_names.txt | cut -f1 > GSM_accessions
grep -f GSM_accessions SraRunTable.txt | cut -f1 -d',' > SRR_accessions

# download and pre-process files
for sample in $(cat SRR_accessions)
do
  echo ${sample}
  mkdir ${sample}
  cd ${sample}

  # 1. download fastq file
  echo "... downloading fastq"
  fasterq-dump ${sample}

  # 2. trim adapter
  echo "... trimming adapter"
  cutadapt -a TGGAATTCTCGGGTGCCAAGG --trimmed-only -m 20 -o ${sample}_trimmed.fq \
  ${sample}.fastq > ${sample}_cutadapt.out

  # 3. filter out reads aligning to rRNA
  echo "... removing rRNA reads"
  bowtie -p 10 -v 2 -S --un ${sample}_not_rRNA.fq ${ref_dir}/ScerRRNA \
  ${sample}_trimmed.fq > ${sample}_rRNA.sam 2> ${sample}_rRNA.bowtiestats

  # 4. filter out reads aligning to non-coding genes
  echo "... removing ncRNA reads"
  bowtie -p 10 -v 2 -S --un ${sample}_not_rRNA_ncRNA.fq ${ref_dir}/rna_coding \
  ${sample}_not_rRNA.fq > ${sample}_ncRNA.sam 2> ${sample}_ncRNA.bowtiestats

  # 5. align to transcriptome
  echo "... aligning to transcriptome"
  bowtie -p 10 -v 2 -S --un ${sample}_unmapped.fq -a --norc \
  ${ref_dir}/scer.transcripts.20cds20 ${sample}_not_rRNA_ncRNA.fq \
  > ${sample}_footprints.sam 2> ${sample}_footprints.bowtiestats

  # 6. calculate multi-mapping weights
  echo "... calculating multimapping weights"
  rsem-calculate-expression --sam ${sample}_footprints.sam \
  ${ref_dir}/scer.transcripts.20cds20 ${sample} > ${sample}.rsem.stdout \
  2> ${sample}.rsem.stderr

  echo 
done
