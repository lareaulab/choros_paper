#!/bin/bash

ref_dir=${HOME}/choros_paper/reference_data
script_dir=${HOME}/choros_paper/scripts
raw_data_dir=${HOME}/choros_paper/data/wu_2019/raw_data

cd ${raw_data_dir}
# sample_names.txt copied from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115162
# SraRunTable.txt downloaded from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA473984&o=acc_s%3Aa

# pull SRR accessions for samples
## WT CHX
grep WT_CHX sample_names.txt | cut -f1 > WT_CHX_GSM_accessions
grep -f WT_CHX_GSM_accessions SraRunTable.txt | cut -f1 -d',' > WT_CHX_SRR_accessions
## WT CHXTIG
grep "Wild-Type CHXTIG" sample_names.txt | cut -f1 > WT_CHXTIG_GSM_accessions
grep -f WT_CHXTIG_GSM_accessions SraRunTable.txt | cut -f1 -d',' > WT_CHXTIG_SRR_accessions
## 3AT CHX
grep AT sample_names.txt | grep -v TIG | cut -f1 > 3AT_CHX_GSM_accessions
grep -f 3AT_CHX_GSM_accessions SraRunTable.txt | cut -f1 -d',' > 3AT_CHX_SRR_accessions
## 3AT TIG
grep AT sample_names.txt | grep TIG | cut -f1 > 3AT_CHXTIG_GSM_accessions
grep -f 3AT_CHXTIG_GSM_accessions SraRunTable.txt | cut -f1 -d',' > 3AT_CHXTIG_SRR_accessions

for expt in WT_CHX 3AT_CHX
do
  echo ${expt}
  mkdir ${expt}
  cd ${expt}

  # 1. download fastq files
  for SRR_id in $(cat ../${expt}_SRR_accessions)
  do
    echo "... downloading" ${SRR_id}
    fasterq-dump ${SRR_id}
  done

  # 2. pool replicates
  echo "... pooling replicates"
  cat SRR*fastq > ${expt}_raw.fq

  # 3. trim 3' adapter
  echo "... trimming 3' adapter"
  cutadapt -a CTGTAGGCACCATCAAT -m 15 --trimmed-only -o ${expt}_trim3.fq \
  ${expt}_raw.fq > ${expt}_trim3.cutadapt

  # 3. filter out reads aligning to rRNA
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

for expt in WT_CHXTIG 3AT_CHXTIG
do
  echo ${expt}
  mkdir ${expt}
  cd ${expt}

  # 1. download fastq files
  for SRR_id in $(cat ../${expt}_SRR_accessions)
  do
    echo "... downloading" ${SRR_id}
    fasterq-dump ${SRR_id}
  done

  # 2. pool replicates
  echo "... pooling replicates"
  cat SRR*fastq > ${expt}_raw.fq

  # 3. trim 3' adapter
  echo "... trimming 3' adapter"
  cutadapt -a CACTCGGGCACCAAGGA --trimmed-only -m 25 -o ${expt}_trim3.fq \
  ${expt}_raw.fq > ${expt}_trim3.cutadapt

  # 4. trim random nucleotides for deduplication
  echo "... trimming random linker"
  umi_tools extract --extract-method=regex \
  --bc-pattern="(?P<umi_1>.{4})(.+)(?P<umi_2>.{6})$" \
  -I ${expt}_trim3.fq -S ${expt}_trim3_trimUMI.fq -L ${expt}_trim3_trimUMI.umitools

  # 5. remove contaminant reads: rRNA
  echo "... removing rRNA"
  bowtie -v 2 -p 10 -S --un ${expt}_trim3_not_rrna.fq ${ref_dir}/ScerRRNA \
  ${expt}_trim3_trimUMI.fq > ${expt}_trim3_rrna.sam 2> \
  ${expt}_trim3_rrna.bowtiestats

  # 6. remove contaminant reads: ncRNA
  echo "... removing ncRNA"
  bowtie -v 2 -p 10 -S --un ${expt}_trim3_not_rrna_ncrna.fq \
  ${ref_dir}/rna_coding ${expt}_trim3_not_rrna.fq > \
  ${expt}_trim3_ncrna.sam 2> ${expt}_trim3_ncrna.bowtiestats

  # 7. deduplicate: report all best alignments (-a --best --strata)
  echo "... deduplicate: report all best alignments"
  bowtie -v 2 -p 10 -S --norc -a --best --strata --un \
  ${expt}_trim3_unaligned.fq ${ref_dir}/scer.transcripts.20cds20 \
  ${expt}_trim3_not_rrna_ncrna.fq > ${expt}_trim3_not_rrna_ncrna_best.sam 2> \
  ${expt}_trim3_not_rrna_ncrna_best.bowtiestats

  # 8. deduplicate: sort by read name
  echo "... deduplicate: sort by read name"
  samtools view -h -F 0x04 ${expt}_trim3_not_rrna_ncrna_best.sam | \
  samtools sort -n -O SAM -o ${expt}_trim3_not_rrna_ncrna_best_sorted.sam

  # 9. deduplicate: choose positionally-first alignment per reads
  echo "... deduplicate: choose positionally-first alignment"
  Rscript ${script_dir}/filter_sam.R \
  -i ${expt}_trim3_not_rrna_ncrna_best_sorted.sam \
  -o ${expt}_trim3_not_rrna_ncrna_best_sorted_filtered.sam

  # 10. deduplicate: sort and index
  echo "deduplicate: sort and index"
  samtools sort ${expt}_trim3_not_rrna_ncrna_best_sorted_filtered.sam \
  -o ${expt}_trim3_not_rrna_ncrna_unique.bam -O BAM
  samtools index ${expt}_trim3_not_rrna_ncrna_unique.bam

  # 11. deduplicate: deduplicate UMIs
  echo "deduplicate: deduplicate UMIs"
  umi_tools dedup --in-sam --out-sam --read-length --method=adjacency \
  -I ${expt}_trim3_not_rrna_ncrna_unique.bam \
  -S ${expt}_trim3_deduplicated.sam \
  -L ${expt}_trim3_deduplicated.umitools \
  -E ${expt}_trim3_deduplicated.error

  # 12. deduplicate: convert into fastq
  echo "deduplicate: convert into fastq"
  samtools fastq ${expt}_trim3_deduplicated.sam > ${expt}_trim3_deduplicated.fq

  # 13. align to transcriptome
  echo "align to transcriptome"
  bowtie -v 2 -p 10 -S --norc -a ${ref_dir}/scer.transcripts.20cds20 \
  ${expt}_trim3_deduplicated.fq > ${expt}_trim3_footprints.sam 2> \
  ${expt}_trim3_footprints.bowtiestats

  # 14. compute multimapping weights
  echo "compute multimapping weights"
  rsem-calculate-expression --sam ${expt}_trim3_footprints.sam --seed-length 15 \
  ${ref_dir}/scer.transcripts.20cds20 ${expt}_trim3_footprints > \
  ${expt}_trim3_footprints.rsem.stdout 2> ${expt}_trim3_footprints.rsem.stderr

  cd ..
done
