#!/bin/bash

ref_dir=${HOME}/choros_paper/reference_data
raw_data_dir=${HOME}/choros_paper/data/weinberg_2016/raw_data

cd ${raw_data_dir}

expt="weinberg"
mkdir ${expt}

# 1. download fastq files
SRR_id="SRR1049521"
fasterq-dump ${SRR_id}

cd ${expt}

# 2. trim adapter
echo "... trimming adapter"
cutadapt -a TCGTATGCCGTCTTCTGCTTG --trimmed-only -m 20 -o \
${expt}_trim3.fq ../${SRR_id}.fastq > ${expt}_cutadapt.out

# 3. trim 5' random nucleotides for deduplication
echo "... trimming random linker"
umi_tools extract --extract-method=regex --bc-pattern="(?P<umi_1>.{8})(.+)$" \
-I ${expt}_trim3.fq -S ${expt}_trim3_trimUMI.fq -L ${expt}_trim3_trimUMI.umitools

# 3. filter out reads aligning to rRNA
echo "... removing rRNA reads"
bowtie -p 10 -v 2 -S --un ${expt}_not_rRNA.fq ${ref_dir}/ScerRRNA \
${expt}_trim3_trimUMI.fq > ${expt}_rRNA.sam 2> ${expt}_rRNA.bowtiestats

# 4. filter out reads aligning to non-coding genes
echo "... removing ncRNA reads"
bowtie -p 10 -v 2 -S --un ${expt}_not_rRNA_ncRNA.fq ${ref_dir}/rna_coding \
${expt}_not_rRNA.fq > ${expt}_ncRNA.sam 2> ${expt}_ncRNA.bowtiestats

# 5. deduplicate: report all best alignments (-a --best --strata)
echo "... deduplicate: report all best alignments"
bowtie -v 2 -p 10 -S --norc -a --best --strata --un \
${expt}_unaligned.fq ${ref_dir}/scer.transcripts.20cds20 \
../${expt}_not_rrna_ncrna.fq > ${expt}_not_rrna_ncrna_best.sam 2> \
${expt}_not_rrna_ncrna_best.bowtiestats

# 6. deduplicate: sort by read name
echo "... deduplicate: sort by read name"
samtools view -h -F 0x04 ${expt}_not_rrna_ncrna_best.sam | \
samtools sort -n -O SAM -o ${expt}_not_rrna_ncrna_best_sorted.sam

# 7. deduplicate: choose positionally-first alignment per reads
echo "... deduplicate: choose positionally-first alignment"
Rscript ${script_dir}/filter_sam.R \
-i ${expt}_not_rrna_ncrna_best_sorted.sam \
-o ${expt}_not_rrna_ncrna_best_sorted_filtered.sam

# 8. deduplicate: sort and index
echo "... deduplicate: sort and index"
samtools sort ${expt}_not_rrna_ncrna_best_sorted_filtered.sam \
-o ${expt}_not_rrna_ncrna_unique.bam -O BAM
samtools index ${expt}_not_rrna_ncrna_unique.bam

# 9. deduplicate: deduplicate UMIs
echo "... deduplicate: deduplicate UMIs"
umi_tools dedup --output-stats=deduplicated --in-sam --out-sam --read-length \
-I ${expt}_not_rrna_ncrna_unique.bam \
-S ${expt}_deduplicated.sam \
-L ${expt}_deduplicated.umitools \
-E ${expt}_deduplicated.error

# 10. deduplicate: convert into fastq
echo "... deduplicate: convert into fastq"
samtools fastq ${expt}_deduplicated.sam > ${expt}_deduplicated.fq

# 11. align to transcriptome
echo "... aligning to transcriptome"
bowtie -v 2 -p 10 -S --norc -a ${ref_dir}/scer.transcripts.20cds20 \
${expt}_deduplicated.fq > ${expt}_footprints.sam 2> \
${expt}_footprints.bowtiestats

# 12. compute multimapping weights
echo "... compute multimapping weights"
rsem-calculate-expression --sam ${expt}_footprints.sam --seed-length 15 \
${ref_dir}/scer.transcripts.20cds20 ${expt}_footprints > \
${expt}_footprints.rsem.stdout 2> ${expt}_footprints.rsem.stderr

