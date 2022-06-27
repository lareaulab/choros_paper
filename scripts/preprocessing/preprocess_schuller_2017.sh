#!/bin/bash

ref_dir=${HOME}/choros_paper/reference_data
raw_data_dir=${HOME}/choros_paper/data/schuller_2017/raw_data

cd ${raw_data_dir}
# SRR_Acc_List.txt downloaded from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA352976&o=acc_s%3Aa

# pull SRR accessions for WT ribosome profiling (not high salt wash)
grep "Wild-Type" SraRunTable.txt | grep -v "high salt" | cut -f1 -d',' > SRR_accessions

expt="schuller"
mkdir ${expt}

# 1. download fastq files
for sample in $(cat SRR_accessions)
do
  echo "... downloading" ${sample}
  fasterq-dump ${sample}
done

# 2. merge fastq files
echo "... merging lanes"
cd ${expt}
cat ../*fastq > ${expt}.fq

# 3. trim adapter
echo "... trimming adapter"
cutadapt -a CTGTAGGCACCATCAAT --trimmed-only -m 20 -o \
${expt}_trimmed.fq ${expt}.fq > ${expt}_cutadapt.out

# 4. filter out reads aligning to rRNA
echo "... removing rRNA reads"
bowtie -p 10 -v 2 -S --un ${expt}_not_rRNA.fq ${ref_dir}/ScerRRNA \
${expt}_trimmed.fq > ${expt}_rRNA.sam 2> ${expt}_rRNA.bowtiestats

# 5. filter out reads aligning to non-coding genes
echo "... removing ncRNA reads"
bowtie -p 10 -v 2 -S --un ${expt}_not_rRNA_ncRNA.fq ${ref_dir}/rna_coding \
${expt}_not_rRNA.fq > ${expt}_ncRNA.sam 2> ${expt}_ncRNA.bowtiestats

# 6. align to transcriptome
echo "... aligning to transcriptome"
bowtie -p 10 -v 2 -S --un ${expt}_unmapped.fq -a --norc \
${ref_dir}/scer.transcripts.20cds20 ${expt}_not_rRNA_ncRNA.fq \
> ${expt}_footprints.sam 2> ${expt}_footprints.bowtiestats

# 7. calculate multi-mapping weights
echo "... calculating multimapping weights"
rsem-calculate-expression --sam ${expt}_footprints.sam \
${ref_dir}/scer.transcripts.20cds20 ${expt} > ${expt}.rsem.stdout \
2> ${expt}.rsem.stderr
