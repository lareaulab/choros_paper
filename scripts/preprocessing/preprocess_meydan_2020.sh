#!/bin/bash

ref_dir=${HOME}/choros_paper/reference_data
script_dir=${HOME}/choros_paper/scripts
raw_data_dir=${HOME}/choros_paper/data/meydan_2020/raw_data

cd ${raw_data_dir}

### download fastq files
# sample name in Table S3
# SRR accession from SraRunTable

mkdir monosome
cd monosome
fasterq-dump SRR10302098
fasterq-dump SRR10302100
cd ..

mkdir disome
cd disome
fasterq-dump SRR10302102
fasterq-dump SRR10302104
fasterq-dump SRR10302108
cd ..

mkdir rnaseq
cd rnaseq
fasterq-dump SRR10302118
fasterq-dump SRR10302120
cd ..

### process ribosome profiling data
for expt in monosome disome
do
  echo ${expt}
  cd ${expt}

  # 1. merge replicates
  cat SRR*fastq > ${expt}_raw.fq

  # 2. trim UMIs
  echo "... trimming UMI"
  umi_tools extract --extract-method=regex \
  --bc-pattern="(?P<umi_1>.{2})(.+)(?P<umi_2>.{5})$" \
  -I ${expt}_raw.fq -S ${expt}_trimUMI.fq -L ${expt}_trimUMI.umitools

  # 3. remove contaminant reads: rRNA
  echo "... removing rRNA reads"
  bowtie -v 2 -p 10 -S --un ${expt}_trimUMI_not_rrna.fq ${ref_dir}/ScerRRNA \
  ${expt}_trimUMI.fq > ${expt}_trimUMI_rrna.sam 2> \
  ${expt}_trimUMI_rrna.bowtiestats

  # 4. remove contaminant reads: ncRNA
  echo "... removing ncRNA reads"
  bowtie -v 2 -p 10 --un ${expt}_trimUMI_not_rrna_ncrna.fq ${ref_dir}/rna_coding \
  ${expt}_trimUMI_not_rrna.fq > ${expt}_trimUMI_ncrna.sam 2> \
  ${expt}_trimUMI_ncrna.bowtiestats

  # 5. deduplicate: report all best alignments (-a --best --strata)
  echo "... deduplication: report all best alignments"
  bowtie -v 2 -p 10 -S --norc -a --best --strata --un \
  ${expt}_trimUMI_unaligned.fq ${ref_dir}/scer.transcripts.20cds20 \
  ${expt}_trimUMI_not_rrna_ncrna.fq > ${expt}_trimUMI_not_rrna_ncrna_best.sam 2> \
  ${expt}_trimUMI_not_rrna_ncrna_best.bowtiestats

  # 6. deduplicate: sort by read name
  echo "... deduplication: sort by read name"
  samtools view -h -F 0x04 ${expt}_trimUMI_not_rrna_ncrna_best.sam | \
  samtools sort -n -O SAM -o ${expt}_trimUMI_not_rrna_ncrna_best_sorted.sam

  # 7. deduplicate: choose positionally-first alignment per reads
  echo "... deduplication: choose positionally-first alignments"
  Rscript ${HOME}/footprint-bias/scripts/filter_sam.R \
  -i ${expt}_trimUMI_not_rrna_ncrna_best_sorted.sam \
  -o ${expt}_trimUMI_not_rrna_ncrna_best_sorted_filtered.sam

  # 8. deduplicate: sort and index
  echo "... deduplicate: sort and index"
  samtools sort ${expt}_trimUMI_not_rrna_ncrna_best_sorted_filtered.sam \
  -o ${expt}_trimUMI_not_rrna_ncrna_unique.bam -O BAM

  samtools index ${expt}_trimUMI_not_rrna_ncrna_unique.bam

  # 9. deduplicate: deduplicate UMIs
  echo "... deduplicate: deduplicate UMIs"
  umi_tools dedup --output-stats=deduplicated --in-sam --out-sam --read-length \
  -I ${expt}_trimUMI_not_rrna_ncrna_unique.bam \
  -S ${expt}_trimUMI_deduplicated.sam \
  -L ${expt}_trimUMI_deduplicated.umitools \
  -E ${expt}_trimUMI_deduplicated.error

  # 10. deduplicate: convert into fastq
  echo "... deduplicate: convert deduplicated reads into fastq"
  samtools fastq ${expt}_trimUMI_deduplicated.sam > ${expt}_trimUMI_deduplicated.fq

  # 11. align to transcriptome
  echo "... aligning to transcriptome"
  bowtie -v 2 -p 10 -S --norc -a ${ref_dir}/scer.transcripts.20cds20 \
  ${expt}_trimUMI_deduplicated.fq > ${expt}_trimUMI_footprints.sam 2> \
  ${expt}_trimUMI_footprints.bowtiestats

  # 12. compute multimapping weights
  echo "... computing multimapping weights"
  rsem-calculate-expression --sam ${expt}_trimUMI_footprints.sam \
  ${ref_dir}/scer.transcripts.20cds20 ${expt}_trimUMI_footprints > \
  ${expt}_trimUMI_footprints.rsem.stdout 2> ${expt}_trimUMI_footprints.rsem.stderr

  cd ..
  
done

### rnaseq data

expt="rnaseq"
cd ${expt}
echo ${expt}

# 1. merge replicates
cat SRR*fastq > ${expt}.fq

# 2. remove contaminant reads: rRNA
bowtie -v 2 -p 10 -S --un ${expt}_not_rrna.fq ${ref_dir}/ScerRRNA \
${expt}.fq > ${expt}_rrna.sam 2> ${expt}_rrna.bowtiestats

# 3. remove contaminant reads: ncRNA
bowtie -v 2 -p 10 -S --un ${expt}_not_rrna_ncrna.fq ${ref_dir}/rna_coding \
${expt}_not_rrna.fq > ${expt}_ncrna.sam 2> ${expt}_ncrna.bowtiestats

# 4. align to transcriptome
bowtie -v 2 -p 10 -S -a ${ref_dir}/scer.transcripts.20cds20 \
${expt}_not_rrna_ncrna.fq > ${expt}_transcripts.sam 2> ${expt}_transcripts.bowtiestats

# 5. compute multimapping weights
rsem-calculate-expression --sam ${expt}_transcripts.sam \
${ref_dir}/scer.transcripts.20cds20 ${expt}_transcripts > \
${expt}_transcripts.rsem.stdout 2> ${expt}_transcripts.rsem.stderr


