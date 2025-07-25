#!/bin/bash
set -euxo pipefail

VERSION=47
ASSEMBLIES=$1

# Create a list of chromosomes to map the annotation from
rm -f "${ASSEMBLIES}/CHM13v2/chroms.txt"
for chromosome in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y;
do
  echo chr${chromosome},chr${chromosome} >> "${ASSEMBLIES}/CHM13v2/chroms.txt"
done

# Unplaced contigs
rm -f "${ASSEMBLIES}/CHM13v2/contigs.txt"
for chromosome in KI270728.1 KI270727.1 KI270442.1 KI270729.1 GL000225.1 KI270743.1 GL000008.2 GL000009.2 KI270747.1 \
                  KI270722.1 GL000194.1 KI270742.1 GL000205.2 GL000195.1 KI270736.1 KI270733.1 GL000224.1 GL000219.1 \
                  KI270719.1 GL000216.2 KI270712.1 KI270706.1 KI270725.1 KI270744.1 KI270734.1 GL000213.1 GL000220.1 \
                  KI270715.1 GL000218.1 KI270749.1 KI270741.1 GL000221.1 KI270716.1 KI270731.1 KI270751.1 KI270750.1 \
                  KI270519.1 GL000214.1 KI270708.1 KI270730.1 KI270438.1 KI270737.1 KI270721.1 KI270738.1 KI270748.1 \
                  KI270435.1 GL000208.1 KI270538.1 KI270756.1 KI270739.1 KI270757.1 KI270709.1 KI270746.1 KI270753.1 \
                  KI270589.1 KI270726.1 KI270735.1 KI270711.1 KI270745.1 KI270714.1 KI270732.1 KI270713.1 KI270754.1 \
                  KI270710.1 KI270717.1 KI270724.1 KI270720.1 KI270723.1 KI270718.1 KI270317.1 KI270740.1 KI270755.1 \
                  KI270707.1 KI270579.1 KI270752.1 KI270512.1 KI270322.1 GL000226.1 KI270311.1 KI270366.1 KI270511.1 \
                  KI270448.1 KI270521.1 KI270581.1 KI270582.1 KI270515.1 KI270588.1 KI270591.1 KI270522.1 KI270507.1 \
                  KI270590.1 KI270584.1 KI270320.1 KI270382.1 KI270468.1 KI270467.1 KI270362.1 KI270517.1 KI270593.1 \
                  KI270528.1 KI270587.1 KI270364.1 KI270371.1 KI270333.1 KI270374.1 KI270411.1 KI270414.1 KI270510.1 \
                  KI270390.1 KI270375.1 KI270420.1 KI270509.1 KI270315.1 KI270302.1 KI270518.1 KI270530.1 KI270304.1 \
                  KI270418.1 KI270424.1 KI270417.1 KI270508.1 KI270303.1 KI270381.1 KI270529.1 KI270425.1 KI270396.1 \
                  KI270363.1 KI270386.1 KI270465.1 KI270383.1 KI270384.1 KI270330.1 KI270372.1 KI270548.1 KI270580.1 \
                  KI270387.1 KI270391.1 KI270305.1 KI270373.1 KI270422.1 KI270316.1 KI270340.1 KI270338.1 KI270583.1 \
                  KI270334.1 KI270429.1 KI270393.1 KI270516.1 KI270389.1 KI270466.1 KI270388.1 KI270544.1 KI270310.1 \
                  KI270412.1 KI270395.1 KI270376.1 KI270337.1 KI270335.1 KI270378.1 KI270379.1 KI270329.1 KI270419.1 \
                  KI270336.1 KI270312.1 KI270539.1 KI270385.1 KI270423.1 KI270392.1 KI270394.1;
do
  echo ${chromosome} >> "${ASSEMBLIES}/CHM13v2/contigs.txt"
done

# Unzip fasta files, otherwise Liftoff will not work
gunzip "${ASSEMBLIES}/CHM13v2/CHM13v2.fa.gz" &
gunzip "${ASSEMBLIES}/CHM13v2/gencode/liftoff/GRCh38.primary_assembly.genome.fa.gz" &
wait

liftoff \
  -g "${ASSEMBLIES}/CHM13v2/gencode/liftoff/gencode.v${VERSION}.primary_assembly.annotation.gff3.gz" \
  -o "${ASSEMBLIES}/CHM13v2/gencode/CHM13v2.gencode-${VERSION}.gff3" \
  -u "${ASSEMBLIES}/CHM13v2/gencode/CHM13v2.gencode-${VERSION}.unmapped-gff3.txt" \
  -chroms "${ASSEMBLIES}/CHM13v2/chroms.txt" \
  -unplaced "${ASSEMBLIES}/CHM13v2/contigs.txt" \
  -exclude_partial -a 0.95 -s 0.95 \
  -p "$(nproc)" \
  -copies -sc 0.95 -polish -cds \
  "${ASSEMBLIES}/CHM13v2/CHM13v2.fa" \
  "${ASSEMBLIES}/CHM13v2/gencode/liftoff/GRCh38.primary_assembly.genome.fa"

# Cleanup
rm -rf intermediate_files \
  "${ASSEMBLIES}/CHM13v2/gencode/liftoff/GRCh38.primary_assembly.genome.fa.fai" \
  "${ASSEMBLIES}/CHM13v2/gencode/liftoff/gencode.v${VERSION}.primary_assembly.annotation.gff3.gz_db" \
  "${ASSEMBLIES}/CHM13v2/chroms.txt" \
  "${ASSEMBLIES}/CHM13v2/contigs.txt" \
  "${ASSEMBLIES}/CHM13v2/CHM13v2.fa.mmi" \
  "${ASSEMBLIES}/CHM13v2/CHM13v2.fa.fai" \
  "${ASSEMBLIES}/CHM13v2/gencode/CHM13v2.gencode-${VERSION}.gff3"

mv "${ASSEMBLIES}/CHM13v2/gencode/CHM13v2.gencode-${VERSION}.gff3_polished" \
   "${ASSEMBLIES}/CHM13v2/gencode/CHM13v2.liftoff+gencode-${VERSION}.gff3"

gzip -f "${ASSEMBLIES}/CHM13v2/gencode/CHM13v2.liftoff+gencode-${VERSION}.gff3" &

# Compress fasta files back
bgzip -@ "$(nproc)" "${ASSEMBLIES}/CHM13v2/CHM13v2.fa" && \
  samtools faidx "${ASSEMBLIES}/CHM13v2/CHM13v2.fa.gz" &
bgzip -@ "$(nproc)" "${ASSEMBLIES}/CHM13v2/gencode/liftoff/GRCh38.primary_assembly.genome.fa" && \
  samtools faidx "${ASSEMBLIES}/CHM13v2/gencode/liftoff/GRCh38.primary_assembly.genome.fa.gz" &
wait

