#!/bin/bash

ref_dir=${HOME}/choros_paper/reference_data
script_dir=${HOME}/choros_paper/scripts
raw_data_dir=${HOME}/choros_paper/data/lecanda_2016/raw_data

cd ${raw_data_dir}
# sample_names.txt copied from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84746
# SraRunTable.txt downloaded from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA330982&o=acc_s%3Aa

# pull SRR accessions for fixedLinker_fixedPrimer and randomLinker_randomPrimer samples
grep Nonrandom sample_names.txt | grep -v mRNA | cut -f1 > \
fixedLinker_fixedPrimer_GSM_accessions
grep -f fixedLinker_fixedPrimer_GSM_accessions SraRunTable.txt | cut -f1 -d',' > \
fixedLinker_fixedPrimer_SRR_accessions
grep 4+3N sample_names.txt | grep -v mRNA | cut -f1 > \
randomLinker_randomPrimer_GSM_accessions
grep -f randomLinker_randomPrimer_GSM_accessions SraRunTable.txt | cut -f1 -d',' > \
randomLinker_randomPrimer_SRR_accessions

### process fixedLinker_fixedPrimer sample

expt="fixedLinker_fixedPrimer"
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
cutadapt -a CTGTAGGCACCATCAAT -m 20 --trimmed-only -o \
${expt}_trim3.fq ${expt}_raw.fq > ${expt}_trim3.cutadapt

# 4. remove contaminant reads: rRNA
echo "... removing rRNA"
bowtie -v 2 -p 10 -S --un ${expt}_trim3_not_rrna.fq ${ref_dir}/ScerRRNA \
${expt}_trim3.fq > ${expt}_trim3_rrna.sam 2> ${expt}_trim3_rrna.bowtiestats

# 5. remove contaminant reads: ncRNA
echo "... removing ncRNA"
bowtie -v 2 -p 10 -S --un ${expt}_trim3_not_rrna_ncrna.fq \
${ref_dir}/rna_coding ${expt}_trim3_not_rrna.fq > \
${expt}_trim3_ncrna.sam 2> ${expt}_trim3_ncrna.bowtiestats

# 6. align to transcriptome 
echo "... aligning to transcriptome"
bowtie -v 2 -p 10 -S --norc -a ${ref_dir}/scer.transcripts.20cds20 \
${expt}_trim3_not_rrna_ncrna.fq > ${expt}_trim3_footprints.sam 2> \
${expt}_trim3_footprints.bowtiestats

# 7. compute multimapping weights
echo "... compute multimapping weights"
rsem-calculate-expression --sam ${expt}_trim3_footprints.sam \
${ref_dir}/scer.transcripts.20cds20 ${expt}_trim3_footprints > \
${expt}_trim3_footprints.rsem.stdout 2> ${expt}_trim3_footprints.rsem.stderr

cd ..

### process randomLinker_randomPrimer sample

expt="randomLinker_randomPrimer"
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
cutadapt -a CTGTAGGCACCATCAAT -m 20 --trimmed-only -o \
${expt}_trim3.fq ${expt}_raw.fq > ${expt}_trim3.cutadapt

# 4. trim random nucleotides for deduplication
echo "... trimming random linker"
umi_tools extract --extract-method=regex \
--bc-pattern="(?P<umi_1>.{3})(.+)(?P<umi_2>.{4})$" \
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
umi_tools dedup --output-stats=deduplicated --in-sam --out-sam --read-length \
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
rsem-calculate-expression --sam ${expt}_trim3_footprints.sam \
${ref_dir}/scer.transcripts.20cds20 ${expt}_trim3_footprints > \
${expt}_trim3_footprints.rsem.stdout 2> ${expt}_trim3_footprints.rsem.stderr

