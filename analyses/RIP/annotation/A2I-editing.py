import gzip
import pickle
from collections import defaultdict

from biobit.core.loc import Orientation
from intervaltree import intervaltree
from joblib import Parallel, delayed
from tqdm import tqdm

from stories import A2I
from stories.RIP import annotation


# Load and index editing sites
def index_editing_sites(assembly: str) -> dict[tuple[str, Orientation], intervaltree.IntervalTree]:
    bed = A2I.tracks.all_passed / f"{assembly}.bed.gz"

    sites = defaultdict(intervaltree.IntervalTree)
    with gzip.open(bed, 'rt') as stream:
        for line in stream:
            seqid, start, end, _, _, strand = line.strip().split('\t')
            sites[seqid, Orientation(strand)].addi(int(start), int(end))
    return dict(sites)


EDITING_SITES = {assembly: index_editing_sites(assembly) for assembly in ["CHM13v2", "GRCm39"]}


def job(config: annotation.Config):
    assembly = {x.assembly for x in config.comparisons}
    assert len(assembly) == 1, assembly

    all_esites = EDITING_SITES[assembly.pop()]

    results = {}
    for partition in tqdm(config.elements):
        if (partition.contig, partition.orientation) not in all_esites:
            results[partition.ind] = 0
            continue

        esites = all_esites[partition.contig, partition.orientation]
        overlap = set()

        # Calculate scores for individual peaks
        for peak in partition.peaks:
            overlap |= esites.overlap(peak.start, peak.end)

        for invrep in partition.invrep:
            for seqrng in invrep.seqranges():
                overlap |= esites.overlap(seqrng.start, seqrng.end)

        results[partition.ind] = len(overlap)

    config.a2i.parent.mkdir(parents=True, exist_ok=True)
    with open(config.a2i, 'wb') as stream:
        pickle.dump(results, stream)


Parallel(n_jobs=-1, verbose=100, backend='sequential')(delayed(job)(config) for config in annotation.Config.load())
