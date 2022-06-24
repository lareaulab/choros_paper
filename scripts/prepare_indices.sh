#!/bin/bash

ref_dir=${HOME}/choros_paper/reference_data

cd ${ref_dir}

bowtie-build ScerRRNA.fa ScerRRNA

bowtie-build rna_coding.fa rna_coding

bowtie-build scer.transcripts.20cds20.fa scer.transcripts.20cds20
rsem-prepare-reference scer.transcripts.20cds20.fa scer.transcripts.20cds20
