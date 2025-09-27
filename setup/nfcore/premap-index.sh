#!/bin/bash
set -euxo pipefail

FOLDER="$1"
PREMAP="$2"
BACKLINK="$(pwd)"

cd "$FOLDER"

# Required:
#  FOLDER:
#  - sequence.fa
#  - annotation.gff3
#  PREMAP:
#  - pre-mapping.bed.gz
#  - pre-mapping.fa.gz

cp "${PREMAP}/pre-mapping.fa.gz" .
gunzip pre-mapping.fa.gz

cp "${PREMAP}/pre-mapping.bed.bgz" .
gunzip pre-mapping.bed.bgz --stdout > pre-mapping.bed
rm pre-mapping.bed.bgz

# Generate STAR index
STAR --runMode genomeGenerate --runThreadN "$(nproc)" \
  --genomeDir STAR-premap-index --genomeFastaFiles sequence.fa pre-mapping.fa \
  --sjdbGTFfile annotation.gff3 --sjdbOverhang 149 --limitGenomeGenerateRAM 128849018880 # 120GB RAM

cd "$BACKLINK"
