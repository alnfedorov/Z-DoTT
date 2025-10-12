import gzip
import logging
from collections import defaultdict

from pybedtools import BedTool, Interval
from pysam import FastxFile
from tqdm import tqdm

import resources
from lib import bed, fasta, logs
from lib.assembly import Assembly, RefSeqAnnotome, HasGencodeAnnotome, HasGencodeLiftoffAnnotome, HasRefSeqAnnotome, \
    HasRepeatMasker
from workflows.assemblies import zdmm, zdhs

logs.setup()
resources.RESULTS.mkdir(parents=True, exist_ok=True)

intervals = defaultdict(list[Interval])
collected = defaultdict(int)

for assembly in zdmm(), zdhs():  # type: Assembly
    # Fetch GENCODE genes
    if isinstance(assembly, HasGencodeAnnotome):
        gencode = assembly.gencode()
    elif isinstance(assembly, HasGencodeLiftoffAnnotome):
        gencode = assembly.gencode_liftoff()
    else:
        raise TypeError(assembly)

    for transcript in gencode.rnas.values():
        if transcript.attrs.type in {'rRNA_pseudogene', 'rRNA'}:
            collected[assembly.name, "GENCODE"] += 1
            seqid, strand = transcript.loc.seqid, transcript.loc.strand.symbol()
            intervals[assembly].append(Interval(seqid, transcript.loc.start, transcript.loc.end, strand=strand))

    # Refseq genes
    if isinstance(assembly, HasRefSeqAnnotome):
        refseq: RefSeqAnnotome = assembly.refseq()
        for transcript in refseq.rnas.values():
            if transcript.attrs.biotype == 'rRNA':
                collected[assembly.name, "RefSeq"] += 1
                seqid, strand = transcript.loc.seqid, transcript.loc.strand.symbol()
                intervals[assembly].append(Interval(seqid, transcript.loc.start, transcript.loc.end, strand=strand))

    # Repmasker annotations
    if isinstance(assembly, HasRepeatMasker):
        repcls = assembly.repeat_masker_classes()
        for i in assembly.repeat_masker():
            classification = repcls.classify(i.name)
            if classification is None:
                raise ValueError(f"Unknown repeat: {i.name} ({i})")

            _, _, cls = classification
            if cls in {"rRNA"}:
                collected[assembly.name, "RepeatMasker"] += 1
                intervals[assembly].append(i)

    # Curated sequences
    path = RESOURCES / f"{assembly.name}.bed"
    if path.exists():
        logging.info(f"Loading curated rRNA annotations for {assembly.name} from {path}")
        for i in BedTool(RESOURCES / f"{assembly.name}.bed"):
            collected[assembly.name, "Curated"] += 1
            intervals[assembly].append(Interval(i.chrom, i.start, i.end, strand=i.strand))

print("Collected sequences from the annotation:")
for k in collected:
    print(f"{k[0]}: {collected[k]} ({k[1]})")

RNA, inside_assembly = [], {}
for assembly, ints in intervals.items():
    merged = bed.merge_stranded([BedTool(ints)]).sort()
    for i in merged:
        seq = fasta.sequence(assembly._fasta, i.chrom, i.start, i.end, strand=i.strand)
        name = f"{assembly.name}_{i.chrom}_{i.start}_{i.end}_{i.strand}"
        assert name not in inside_assembly
        inside_assembly[name] = i
        RNA.append((name, seq))

# NCBI 45s / RFAM 5s & 5.8s / silva 28s & 18s
for x in "NCBI", "rfam", "silva":
    file = RESOURCES / f"rRNA.{x}.fa.gz"
    assert file.exists()

    with FastxFile(file.as_posix(), persist=True) as f:
        for seq in f:
            name, seq = seq.name, seq.sequence
            name = name.replace(' ', '_').replace('-', '_')
            RNA.append((name, seq))

# Deduplicate
RNA = sorted(RNA, key=lambda x: len(x[1]), reverse=True)
before, todrop = len(RNA), set()

for i in tqdm(range(len(RNA))):
    if i in todrop:
        continue

    reference = RNA[i][1]
    for j in range(i + 1, len(RNA)):
        seq = RNA[j][1]
        if j not in todrop and seq in reference:
            todrop.add(j)
optimized = [x for ind, x in enumerate(RNA) if ind not in todrop]

after = len(optimized)
print(f'Collapsed: {before - after} sequences ({(before - after) / before:.2%}%)')

for assembly in GRCm39, CHM13v2:
    # Partition RNA sequences into those that are inside the assembly and those that are outside
    # Sequences that are inside the assembly are saved as BED intervals (even if they were 'optimized' away).
    # Outside optimized sequences are saved as BED intervals too and also written to a FASTA file.

    bed, outside = [], []
    included = set()

    # First pass: collect all sequences that are inside the assembly
    for name, seq in RNA:
        if name.startswith(assembly.name):
            interval = inside_assembly[name]
            # Intentionally drop the strand information
            bed.append(Interval(interval.chrom, interval.start, interval.end))
            included.add(name)

    # Second pass: collect all sequences that are outside the assembly
    for name, seq in optimized:
        if name in included:
            continue
        assert not name.startswith(assembly.name)
        outside.append((name, seq))
        bed.append(Interval(name, 0, len(seq)))

    saveto = RESULTS / assembly.name
    saveto.mkdir(exist_ok=True)

    # Write RNA intervals to BED file
    BedTool(bed).sort().merge().sort().saveas(saveto / "pre-mapping.bed.gz")

    # Write pre-mapping sequences to FASTA file
    with gzip.open(saveto / "pre-mapping.fa.gz", 'wt') as stream:
        for name, seq in outside:
            stream.write(f">{name}\n{seq}\n")
