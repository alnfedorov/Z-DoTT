import pickle
from collections import defaultdict

from biobit.core.loc import Orientation, Interval
from joblib import Parallel, delayed

import ld
import utils.assembly
from stories.RIP import pcalling, annotation, clustering

POOL = Parallel(n_jobs=-1, verbose=100)


def job(cmp: pcalling.Config, seqid: str, strand: Orientation, segments: list[Interval]):
    seqsize = utils.assembly.seqsizes(cmp.organism)[seqid]
    tracks = clustering.invrep_scoring.ExperimentTracks(cmp) \
        .open(seqid, seqsize, strand.to_strand())

    scores = {}
    for segment in segments:
        scores[segment] = float(
            (tracks.signal(segment.start, segment.end) - tracks.control(segment.start, segment.end)).max()
        )
    return cmp.project, cmp.ind, seqid, strand, scores


def index(config: annotation.Config):
    # Parse all segments of interest
    allsegments = defaultdict(set)
    for e in config.elements:
        allsegments[e.contig, e.orientation].update(e.peaks)
        for repeat in e.invrep:
            for segment in repeat.seqranges():
                allsegments[e.contig, e.orientation].add(segment)
    allsegments = dict(allsegments)

    print(f"Processing {config.ind}: {len(config.comparisons) * len(allsegments)} tasks")
    processed = POOL(
        delayed(job)(cmp, seqid, strand, segments)
        for cmp in config.comparisons
        for (seqid, strand), segments in allsegments.items()
    )

    values = defaultdict(dict)
    for project, cmp, contig, strand, score in processed:
        assert (project, cmp) not in values[contig, strand]
        values[contig, strand][project, cmp] = score
    return values


signal = {config.ind: index(config) for config in annotation.Config.load()}

saveto = ld.cache.signal.pkl
saveto.parent.mkdir(parents=True, exist_ok=True)
with open(saveto, 'wb') as stream:
    pickle.dump(signal, stream)
