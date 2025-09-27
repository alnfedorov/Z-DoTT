#!/bin/bash
set -euxo pipefail

ASSEMBLIES="$1"
PREMAPPING="$2"
INDEXES="$3"

# Drop old indexes
rm -rf "${INDEXES}/GRCm39" "${INDEXES}/CHM13v2"

# Cat the annotation data
mkdir -p "${INDEXES}/CHM13v2"
zcat "${ASSEMBLIES}/CHM13v2/CHM13v2.fa.bgz" > "${INDEXES}/CHM13v2/sequence.fa"
zcat "${ASSEMBLIES}/CHM13v2/gencode/CHM13v2.liftoff+gencode-47.gff3.gz" > "${INDEXES}/CHM13v2/annotation.gff3"

mkdir -p "${INDEXES}/GRCm39"
zcat "${ASSEMBLIES}/GRCm39/GRCm39.primary_assembly.genome.fa.bgz" > "${INDEXES}/GRCm39/sequence.fa"
zcat "${ASSEMBLIES}/GRCm39/gencode/gencode.vM36.primary_assembly.annotation.gff3.gz" > "${INDEXES}/GRCm39/annotation.gff3"

for assembly in CHM13v2 GRCm39; do
  # Add viruses
  for virus in EMCV SARSCov2;
  do
    cat "${ASSEMBLIES}/${virus}/sequence.fa" >> "${INDEXES}/${assembly}/sequence.fa"
    cat "${ASSEMBLIES}/${virus}/sequence.gff3" >> "${INDEXES}/${assembly}/annotation.gff3"
  done

  # Generate pre-mapping index
  bash premap-index.sh "${INDEXES}/${assembly}" "${PREMAPPING}/${assembly}"

  # Run the main indexing
  bash main-index.sh "${INDEXES}/${assembly}"

  # Link configs
  for fname in nextflow.config params.yaml;
  do
    ln -s "$(pwd)/${fname}" "${INDEXES}/${assembly}/${fname}"
  done

  # Rename gff3 to gff to satisfy nf-core/rnaseq
  mv "${INDEXES}/${assembly}/annotation.gff3" "${INDEXES}/${assembly}/annotation.gff"
done
