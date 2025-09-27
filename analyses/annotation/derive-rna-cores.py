import pickle

from biobit.core.loc import Strand, Interval
from biobit.toolkit.annotome.transcriptome import Location
from joblib import Parallel, delayed
from pybedtools import BedTool
from pybedtools.cbedtools import defaultdict

import ld
import utils
from stories.annotation import RNACore, IntronCore


def job(assembly: str):
    results = defaultdict(list)

    print(f"Processing {assembly}...")
    assembly = utils.assembly.get(name=assembly)

    # Boundaries include only boundaries that are "relevant" for the current RNA isolation protocol.
    # They might/will not exactly recapitulate RNAs that are processed in this script.
    # = we shouldn't assume that a given TSS/TES/Splice site is recorded in the index
    assembly_boundaries = ld.load.transcription_boundaries(assembly.name)
    annotation = ld.load.gencode(assembly.name)

    # Aberrant templates are derived on a per-gene basis
    for gene in annotation.genes.values():
        # Ignore MT - it's easier to present/analyse separately if needed
        if gene.loc.seqid == 'chrM':
            continue

        boundaries = assembly_boundaries[gene.loc.seqid, gene.loc.strand]
        os_boundaries = assembly_boundaries[gene.loc.seqid, gene.loc.strand.flipped()]

        # Iterate over all transcripts and find the canonical one
        canonical = []
        for tid in gene.transcripts:
            rna = annotation.rnas[tid]
            if ld.filters.is_transcription_boundary(rna) and 'Ensembl canonical' in rna.attrs.tags:
                canonical.append(rna)
        if len(canonical) == 0:
            continue

        assert len(canonical) == 1, canonical
        refrna = canonical[0]

        # Shift gene start/end to either the closest annotated TES/TSS or up to the reference RNA
        start = min(refrna.loc.start, boundaries.closest(gene.loc.start, 'right')[0])
        end = max(refrna.loc.end, boundaries.closest(gene.loc.end, 'left')[0])
        assert start <= refrna.loc.start < refrna.loc.end <= end, (start, end, refrna.loc.start, refrna.loc.end)

        loc = Location(gene.loc.seqid, gene.loc.strand, start, end)

        # Deduce upstream and downstream directions
        match loc.strand:
            case Strand.Forward:
                upstream, downstream = 'left', 'right'
            case Strand.Reverse:
                upstream, downstream = 'right', 'left'
            case _:
                raise ValueError(loc.strand)

        # Derive reference exons:
        # - First exon: from the closest boundary to the end of the first exon.
        #               The closest boundary doesn't have to be the transcript's TSS.
        # - Last exon: from the start of the last exon to the closest boundary.
        #              The closest boundary doesn't have to be the transcript's TES.
        refexons = refrna.exons
        if len(refexons) > 1:
            match loc.strand:
                case Strand.Forward:
                    first, _ = boundaries.window(refexons[0].end, upstream, maxsize=refexons[0].len())
                    last, _ = boundaries.window(refexons[-1].start, downstream, maxsize=refexons[-1].len())
                case Strand.Reverse:
                    last, _ = boundaries.window(refexons[0].end, downstream, maxsize=refexons[0].len())
                    first, _ = boundaries.window(refexons[-1].start, upstream, maxsize=refexons[-1].len())
                case _:
                    raise ValueError(loc.strand)

            refexons = Interval.merge([first, *refrna.exons[1:-1], last])

        # Upstream ("Read-in") regions
        pos = loc.start if loc.strand == Strand.Forward else loc.end
        read_in, key = boundaries.window(pos, upstream)
        while key == 'tss':
            newpos = read_in.start if loc.strand == Strand.Forward else read_in.end
            read_in, key = boundaries.window(newpos, upstream)

        # Divergent transcription region
        pos = read_in.end if loc.strand == Strand.Forward else read_in.start
        divergent, _ = os_boundaries.window(pos, upstream)

        # Downstream ("Read-through") regions
        pos = loc.end if loc.strand == Strand.Forward else loc.start
        read_through, key = boundaries.window(pos, downstream)
        while key == 'tes':
            newpos = read_through.end if loc.strand == Strand.Forward else read_through.start
            read_through, key = boundaries.window(newpos, downstream)

        # Introns
        introns = None
        if len(refexons) > 1:
            introns = []
            for i in range(1, len(refexons)):
                start, end = refexons[i - 1].end, refexons[i].start
                assert start < end, (start, end)
                length = end - start

                # Split into chunks = from start to the closest boundary, from the end to the closest boundary
                first, _ = boundaries.window(start, 'right', maxsize=length)
                last, _ = boundaries.window(end, 'left', maxsize=length)
                blocks = Interval.merge([first, last])

                if loc.strand == Strand.Forward:
                    donor, acceptor = refexons[i - 1], refexons[i]
                else:
                    donor, acceptor = refexons[i], refexons[i - 1]
                introns.append(IntronCore(donor, acceptor, tuple(blocks)))

        results[loc.seqid, loc.strand].append(RNACore(
            gene.ind, loc.seqid, loc.strand, refexons, introns, read_through, read_in, divergent
        ))

    # BED for visualization
    records = []
    for (seqid, strand), templates in results.items():
        for template in templates:  # type: RNACore
            records.append(utils.bed.blocks.make(
                seqid, template.exons, orientation=strand.symbol(), name=f"Reference [{template.gid}]"
            ))
            if template.read_through:
                records.append(utils.bed.blocks.make(
                    seqid, [template.read_through], orientation=strand.symbol(), name=f"Read-through [{template.gid}]"
                ))
            if template.read_in:
                records.append(utils.bed.blocks.make(
                    seqid, [template.read_in], orientation=strand.symbol(), name=f"Read-in [{template.gid}]"
                ))
            if template.divergent:
                records.append(utils.bed.blocks.make(
                    seqid, [template.divergent], orientation=strand.flipped().symbol(),
                    name=f"Divergent [{template.gid}]"
                ))
            if template.introns:
                for ind, intron in enumerate(template.introns):
                    blocks = [intron.donor, *intron.intron, intron.acceptor]
                    records.append(utils.bed.blocks.make(
                        seqid, blocks, orientation=strand.symbol(), name=f"Intron {ind} [{template.gid}]"
                    ))

    # Save as pkl
    saveto = ld.paths.rna_cores.pkl[assembly.name]
    saveto.parent.mkdir(parents=True, exist_ok=True)
    with open(saveto, 'wb') as stream:
        pickle.dump(dict(results), stream, protocol=pickle.HIGHEST_PROTOCOL)

    # Save as BED
    utils.bed.tbindex(BedTool(records).sort(), saveto.with_suffix('.bed.gz'))


Parallel(n_jobs=-1)(delayed(job)(assembly) for assembly in ['GRCm39', 'CHM13v2'])
