import pickle
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
from BCBio import GFF
from Bio import SeqIO
from biobit.core.loc import PerStrand
from matplotlib import colors
from pybedtools import BedTool, Interval

from assemblies import HSV1
from . import visual


def primers(saveto: Path):
    primers = {
        "UL25/26": {
            "fwd": "ATCGCCTGCGTGCAGCTGGAGC", "rev": "CGAGGGCAGGTAGTTGGTGATC"
        },
        "UL25/26 AS": {
            "fwd": "CGTAGGTGACGATAGTGCC", "rev": "CTGTTGTACCTGATCACCAAC"
        },
        "UL54 AS": {
            "fwd": "CCAGGCCGAGGTCAATTAG", "rev": "ACCAGAGGCCATATCCGAC",
        },
        "UL54": {
            "fwd": 'ACCAGAGGCCATATCCGACAC', "rev": 'GTCCAGGCCGAGGTCAATTAG'
        },
        "RS1": {
            "fwd": "TGATCACGCGGCTGCTGTACAC", "rev": "GGTGATGAAGGAGCTGCTGTTG",
        },
        "RS1 AS": {
            "fwd": "ATGAAGGAGCTGCTGTTGC", "rev": "AGCAGTACGCCCTGATCAC"
        },
        "UL48": {
            "fwd": "CGAAGCGCTCTCTCGTTTCTTC", "rev": "TACAGGGCCGAGCAGAAGTTG",
        },
        "UL48 AS": {
            "fwd": "GTACAGGGCCGAGCAGAAGTTG", "rev": "CGAAGCGCTCTCTCGTTTCTT"
        },
        "UL39": {
            "fwd": "GGTCGTGGTGGGAAATGTTC", "rev": "CAGGTAGCAGCTGGAGGTGTAG"
        },
        "UL39 AS": {
            "fwd": "GTGTACATCTGGAAGACGGAC", "rev": "TTACGACTGTCTGATCCACAG"
        },
        "UL19/UL20": {
            "fwd": "TCCTTAGCACGATCGAGGTG", "rev": "GACAGGGTGTTGCAATACGAC"
        },
        "UL19/UL20 AS": {
            "fwd": "GCTGGTGGACCTCGAACTGTAC", "rev": "CTTTCTGGAGCTCGGGTTGTC"
        },
        "UL29": {
            "fwd": "CCGGCTCCTGTCGCGCGAGGA", "rev": "GATCGCAGTACCGCAGAATC"
        },
        "UL29 AS": {
            "fwd": "TGATCGCAGTACCGCAGAATC", "rev": "TCCTCGCGCGACAGGAGCCGG"
        },
    }

    with open(HSV1.fasta) as stream:
        contig = list(SeqIO.parse(stream, 'fasta'))
    assert len(contig) == 1
    contig = contig[0]
    assert set(contig) == {"A", "C", "G", "T"}
    forward = str(contig.seq)
    reverse = str(contig.reverse_complement().seq)

    bed = []
    for gene, seq in primers.items():
        for name, seq in seq.items():
            found = False
            for it in re.finditer(seq, forward):
                found = True
                bed.append(Interval(HSV1.contig, it.start(), it.end(), name=f"{gene}-{name}", strand="+"))
            for it in re.finditer(seq, reverse):
                found = True
                start, end = HSV1.size - it.end(), HSV1.size - it.start()
                bed.append(Interval(HSV1.contig, start, end, name=f"{gene}-{name}", strand="-"))
            if not found:
                print(f"Failed to match {gene}-{name}: {seq}")

    print()
    primers = []
    for region in bed:
        primers.append(region)

    saveto.parent.mkdir(parents=True, exist_ok=True)
    BedTool(primers).sort().saveas(saveto)


def miRNA(saveto: Path):
    # http://crdd.osdd.net/servers/virmirna/browse.php?by=HSV1&TYPE=Virus
    miRNA = {
        "hsv1-miR-h8*": "gcccccggucccuguauaua",
        "hsv1-miR-h8": "uauauagggucaggggguuc",
        "hsv1-miR-h7*": "uuuggaucccgaccccucuuc",
        "hsv1-miR-h7": "aaaggggucugcaaccaaagg",
        "hsv1-miR-h6-5p": "gguggaaggcaggggggugua",
        "hsv1-miR-h6-3p": "cacuucccguccuuccauccc",
        "hsv1-miR-h6": "cacuucccguccuuccauccc",
        "hsv1-miR-h5-5p": "ggggggguucgggcaucucuac",
        "hsv1-miR-h5-3p": "gucagagauccaaacccuccgg",
        "hsv1-miR-h4-5p": "gguagaguuugacaggcaagca",
        "hsv1-miR-h4-3p": "cuugccugucuaacucgcuagu",
        "hsv1-miR-h4": "gguagaguuugacaggcaagca",
        "hsv1-miR-h3*": "cuccugaccgcggguuccgagu",
        "hsv1-miR-h3": "cugggacugugcgguugggac",
        "hsv1-miR-h27": "cagaccccuuucuccccccucuu",
        "hsv1-miR-h2-3p": "ccugagccagggacgagugcgacu",
        "hsv1-miR-h2*": "ucgcacgcgcccggcacagacu",
        "hsv1-miR-h2": "ccugagccagggacgagugcgacu",
        "hsv1-miR-h18": "cccgcccgccggacgccgggacc",
        "hsv1-miR-h17": "uggcgcuggggcgcgaggcgg",
        "hsv1-miR-h16": "ccaggaggcugggaucgaaggc",
        "hsv1-miR-h15": "ggccccgggccgggccgccacg",
        "hsv1-miR-h14-5p": "agucgcacucgucccuggcucagg",
        "hsv1-miR-h14-3p": "ucugugccgggcgcgugcgac",
        "hsv1-miR-h13": "uuagggcgaagugcgagcacugg",
        "hsv1-miR-h12": "uugggacgaagugcgaacgcuu",
        "hsv1-miR-h11": "uuaggacaaagugcgaacgc",
        "hsv1-miR-h1*": "uacaccccccugccuuccacccu",
        "hsv1-miR-h1": "gauggaaggacgggaagugga",
        "hsv1-miR-92944": "uggcucggugagcgacgguc",
        "hsv1-miR-69853": "acgacugggccgaggucccg",
        "hsv1-miR-59899": "caccggcgugucggugcug",
        "hsv1-miR-56054b": "ugccucggucuaccggugc",
        "hsv1-miR-55593": "guucaccaagcugcugcuga",
        "hsv1-miR-47858b": "gcgucaacagcgugccgcaga",
        "hsv1-miR-470": "cuugccugucaaacucuaccacc",
        "hsv1-miR-425": "uagcaguuagacaggcaag",
        "hsv1-miR-40262b": "guuuccggagcuggccuac",
        "hsv1-miR-38539": "acgccgcgcaccgcguguug",
        "hsv1-miR-34208b": "aagugcuccagggcgaagau",
        "hsv1-miR-2932": "uuugguugcagaccccuuucuc",
        "hsv1-miR-132792": "agucugaggucgaauccgagac",
        "hsv1-miR-132144": "aacgacgggagcggcugcgg",
        "hsv1-miR-114783": "cgaggcgcuggccuccgccgac",
        "hsv1-miR-107969": "ucacgggaccaaagcgcuucg",
        "hsv1-miR-107631": "gccaucgccuuccuuccaaa",
        "hsv1-miR-106947b": "gggcgggccuguuguuugucuu",
        "hsv1-miR-105617b": "ucagcgcgccaacgaguugg",
    }
    miRNA = {k: v.upper().replace("U", "T") for k, v in miRNA.items()}

    with open(HSV1.fasta, 'r') as stream:
        hsv1 = list(SeqIO.parse(stream, 'fasta'))[0].seq.upper()
    hsv1 = {"+": str(hsv1), "-": str(hsv1.complement())}

    hits = {k: {"+": [], "-": []} for k in miRNA}
    for st in range(0, HSV1.size):
        for mir, seq in miRNA.items():
            en = st + len(seq)
            if hsv1["+"][st: en] == seq:
                hits[mir]["+"].append((st, en))

            if hsv1["-"][st: en][::-1] == seq:
                hits[mir]["-"].append((st, en))

    bed = []
    for mir, strands in hits.items():
        for strand, coords in strands.items():
            for st, en in coords:
                bed.append(Interval(
                    HSV1.contig, st, en, name=mir, strand=strand
                ))

    saveto.parent.mkdir(parents=True, exist_ok=True)
    BedTool(bed).sort().saveas(saveto)


def repeats(saveto: Path):
    # Parse GTF
    contigs = list(GFF.parse(HSV1.ncbi))
    assert len(contigs) == 1
    contig = contigs[0]

    bed = []
    for repeat in contig.features:
        if repeat.type in {"sequence_feature", "inverted_repeat"}:
            note = repeat.qualifiers['Note']
            assert len(note) == 1

            name = note[0]

            loc = repeat.location
            start, end = loc.start, loc.end
            # strand = "+" if loc.strand == 1 else "-"

            bed.append(Interval(contig.id, start, end, name, strand='.', otherfields=[str(start), str(end), "0,0,0"]))

    saveto.parent.mkdir(parents=True, exist_ok=True)
    BedTool(bed).sort().saveas(saveto, trackline='track itemRgb=On')


def orfs(saveto: PerStrand[Path]):
    file = HSV1.ROOT / "doi.org-10.1038-s41467-020-15992-5/supplements-ORFs.csv"
    ORFS = pd.read_csv(file)[["Type", "Location", "Strand"]]

    # Create BED file
    per_strand = {"+": [], "-": []}
    for _, row in ORFS.iterrows():
        exons = row['Location'].split(':')[1]
        exons = [x.split('-') for x in exons.split('|')]
        exons = sorted((int(s), int(e)) for s, e in exons)
        assert all(s <= e for s, e in zip(exons[:-1], exons[1:]))

        blockCount = str(len(exons))
        blockSizes = ",".join(map(str, [e - s for s, e in exons]))
        blockStarts = ",".join(map(str, [s - exons[0][0] for s, e in exons]))

        start, end = exons[0][0], exons[-1][1]
        thickStart, thickEnd = str(start), str(end)
        per_strand[row['Strand']].append(Interval(
            HSV1.contig, start, end, row['Type'], strand=row['Strand'],
            otherfields=[thickStart, thickEnd, "0,0,0", blockCount, blockSizes, blockStarts]
        ))

    for data, path in (per_strand['+'], saveto.forward), (per_strand['-'], saveto.reverse):
        path.parent.mkdir(parents=True, exist_ok=True)
        BedTool(data).sort().saveas(path)


def genes(saveto: PerStrand[Path]):
    PALETTE = {k: colors.to_rgb(v) for k, v in HSV1.utils.palette.items()}
    PALETTE = {k: f"{int(v[0] * 255)},{int(v[1] * 255)},{int(v[2] * 255)}" for k, v in PALETTE.items()}

    # "Immediate Early": "0,162,168",
    # "Early": "186,59,252",
    # "Late": "250,166,35",
    # "Latency": "186,64,64",
    # "Unknown": "180,180,180"

    # Parse GTF
    contigs = list(GFF.parse(HSV1.gff3))
    assert len(contigs) == 1
    contig = contigs[0]

    fwd, rev = [], []
    for gene in contig.features:
        if gene.type != 'gene':
            print(f"Skipped: {gene.type}")
        for RNA in gene.sub_features:
            assert RNA.type in {'mRNA', "ncRNA", 'transcript'}

            loc = RNA.location
            start, end = loc.start, loc.end
            assert start >= gene.location.start and end <= gene.location.end

            exons = sorted((x.location.start, x.location.end) for x in RNA.sub_features if x.type == 'exon')
            assert len(exons) >= 1, RNA

            assert exons[0][0] == start and exons[-1][1] == end, gene
            blockCount, blockSizes, blockStarts = len(exons), [], []
            for exstart, exend in exons:
                blockSizes.append(exend - exstart)
                blockStarts.append(exstart - start)

            # Almost all genes lack annotated 5`UTRs
            blockCount = str(blockCount)
            blockSizes = ",".join(map(str, blockSizes))
            blockStarts = ",".join(map(str, blockStarts))

            if loc.strand == 1:
                strand, writeto = "+", fwd
            else:
                assert loc.strand == -1
                strand, writeto = "-", rev

            name = RNA.qualifiers['Name'][0] if "Name" in RNA.qualifiers else RNA.id
            group = HSV1.utils.classify(name)

            color = PALETTE[group]
            thickStart, thickEnd = str(start), str(end)
            writeto.append(Interval(
                contig.id, start, end, name, strand=strand,
                otherfields=[thickStart, thickEnd, color, blockCount, blockSizes, blockStarts]
            ))

    # Save canonical transcripts
    for data, saveto in (fwd, saveto.forward), (rev, saveto.reverse):
        saveto.parent.mkdir(parents=True, exist_ok=True)
        BedTool(data).sort().saveas(saveto, trackline='track itemRgb=On')


def circos(genes: PerStrand[Path], saveto: Path):
    # Prepare circos annotation - collapse 'almost' identical non-canonical transcripts
    SKIP = {
        "UL5.5-#1", "UL6-*2", "UL6-*1", "UL40.5-#1", "US3-#1", "UL6-#2", "UL8.5-*1",
        "UL30.5-#1", "UL37.6-#1", "UL37.5-#1", "UL37-#1", "UL5-*1", "UL41-#1",
        "TRL2.5-*1", "IRL2.5-*1", "UL7-*1", "UL40.5-*1", "UL51-#1",
        "UL5.6", "UL4.5-*1", "UL4.5"
    }
    SKIP |= HSV1.nagnag_splicing

    fwd, rev = BedTool(genes.forward), BedTool(genes.reverse)
    per_strand = {}
    for genes, strand in [fwd, "+"], [rev, "-"]:
        # Skip NAGNAG
        genes = [x for x in genes if x.name not in SKIP]
        result = genes

        # Save height arrangement in the score column
        result = visual.assort_genes_by_height(result, limits=(0.8, 0), offset=0, step=0.17)
        # Renormalize to 0, 1
        maxscore, minscore = max(float(x.score) for x in result), min(float(x.score) for x in result)
        for x in result:
            x.score = str((float(x.score) - minscore) / (maxscore - minscore))
        per_strand[strand] = result
    fwd, rev = per_strand.pop("+"), per_strand.pop("-")

    for i in rev:
        i.score = str(-float(i.score))

    result = defaultdict(list)
    for genes, strand in (fwd, "+"), (rev, "-"):
        # save introns-exons-5`ends
        introns, exons, fivepr, threepr = [], [], [], []
        for gene in genes:
            _, _, sizes, starts = gene.fields[8:]
            color = HSV1.utils.palette[HSV1.utils.classify(gene.name)]

            # exons
            excoords = []
            for length, start in zip(sizes.split(","), starts.split(',')):
                length, start = int(length), int(start) + gene.start
                end = start + length
                excoords.append((start, end))
                exons.append((start, end, gene.name, float(gene.score), color))
            # introns
            for prevex, nextex in zip(excoords[:-1], excoords[1:]):
                start, end = prevex[1], nextex[0]
                introns.append((start, end, gene.name, float(gene.score), color))

            # 5`end
            if gene.strand == "+":
                start, end = excoords[0]
            else:
                start, end = excoords[-1]
            fivepr.append((start, end, gene.name, float(gene.score), '#000000'))

            # 3`end
            if gene.strand == "+":
                start, end = excoords[-1]
            else:
                start, end = excoords[0]
            threepr.append((start, end, gene.name, float(gene.score), '#000000'))

        for x, name in (introns, "introns"), (exons, "exons"), (fivepr, "5`end"), (threepr, "3`end"):
            result[name, strand] = x

    saveto.parent.mkdir(parents=True, exist_ok=True)
    with open(saveto, 'wb') as stream:
        pickle.dump(dict(result), stream)
