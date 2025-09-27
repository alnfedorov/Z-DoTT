#!/bin/bash
set -euxo pipefail

FOLDER="$1"
BACKLINK="$(pwd)"

cd "$FOLDER"

# Required:
#  - sequence.fa
#  - annotation.gff3

samtools faidx sequence.fa

# Generate STAR index
STAR --runMode genomeGenerate --runThreadN "$(nproc)" \
  --genomeDir STAR-index --genomeFastaFiles sequence.fa \
  --sjdbGTFfile annotation.gff3 --sjdbOverhang 149

# Generate Salmon index
# (https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/)
mkdir Salmon-index
grep "^>" sequence.fa | cut -d " " -f 1 > Salmon-index/decoys.txt
sed -i.bak -e 's/>//g' Salmon-index/decoys.txt

gffread -F -w transcriptome.fa -g sequence.fa annotation.gff3

cat transcriptome.fa sequence.fa > Salmon-index/gentrome.fa
pigz Salmon-index/gentrome.fa && rm -f transcriptome.fa

salmon index -p "$(nproc)" \
  -t Salmon-index/gentrome.fa.gz \
  -d Salmon-index/decoys.txt \
  -i sindex
rm -rf Salmon-index && mv sindex Salmon-index

cd "$BACKLINK"
