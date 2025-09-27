#!/bin/bash
set -euxo pipefail

HUMAN_VERSION=$1
MOUSE_VERSION=$2
ASSEMBLIES=$3

########################################################################################################################
# CHM13 v2.0
########################################################################################################################
mkdir -p "${ASSEMBLIES}/CHM13v2"

# Fasta
wget -O "${ASSEMBLIES}/CHM13v2/CHM13v2.fa.gz" \
  "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"

## RefSeq
#wget -O "${ASSEMBLIES}/CHM13v2/refseq/CHM13v2.refseq.gff.gz" \
#  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz"

# GRCh38 GENCODE annotation and sequence for lift-over
wget -O "${ASSEMBLIES}/CHM13v2/gencode/liftoff/GRCh38.primary_assembly.genome.fa.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${HUMAN_VERSION}/GRCh38.primary_assembly.genome.fa.gz"

wget -O "${ASSEMBLIES}/CHM13v2/gencode/liftoff/gencode.v${HUMAN_VERSION}.primary_assembly.annotation.gff3.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${HUMAN_VERSION}/gencode.v${HUMAN_VERSION}.primary_assembly.annotation.gff3.gz"

########################################################################################################################
# GRCm39
########################################################################################################################
mkdir -p "${ASSEMBLIES}/GRCm39"

wget -O "${ASSEMBLIES}/GRCm39/GRCm39.primary_assembly.genome.fa.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${MOUSE_VERSION}/GRCm39.primary_assembly.genome.fa.gz"

wget -O "${ASSEMBLIES}/GRCm39/gencode/gencode.v${MOUSE_VERSION}.primary_assembly.annotation.gff3.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${MOUSE_VERSION}/gencode.v${MOUSE_VERSION}.primary_assembly.annotation.gff3.gz"

wget -O "${ASSEMBLIES}/GRCm39/gencode/gencode.v${MOUSE_VERSION}.primary_assembly.annotation.gtf.gz" \
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_${MOUSE_VERSION}/gencode.v${MOUSE_VERSION}.primary_assembly.annotation.gtf.gz"

## RefSeq
#wget -O "${ASSEMBLIES}/GRCm39/refseq/GRCm39.refseq.gff.gz" \
#  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gff.gz"

########################################################################################################################
# Postprocessing
########################################################################################################################

for file in "${ASSEMBLIES}"/*/*.fa.gz; do
  saveto="${file/.fa.gz/.fa.bgz}"
  gunzip "${file}" --stdout | bgzip --threads "$(nproc)" -l 6 -o "${saveto}" /dev/stdin
  rm "${file}"
  echo "Indexing" "${saveto}"
  samtools faidx "${saveto}"
done
