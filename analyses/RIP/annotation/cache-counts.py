import pickle
from collections import defaultdict

from biobit import io
from biobit.core.loc import Interval
from biobit.core.ngs import Layout
from biobit.toolkit import countit, nfcore

import ld
import utils.assembly
from stories.RIP import annotation


def reader(path: str, layout: Layout) -> io.bam.Reader:
    assert isinstance(layout, Layout.Paired), f"Unsupported layout: {layout} ({path})"
    return io.bam.Reader(path, inflags=3, exflags=2572, minmapq=0)


def job(config: annotation.Config):
    elements = config.elements

    builder = countit.rigid.Engine.builder().set_threads(-1)

    # Add all elements to the engine
    covered = defaultdict(list)
    _elements = []
    for e in elements:
        segments = e.peaks.copy()
        segments.extend(block for r in e.invrep for block in r.seqranges())
        segments = Interval.merge(segments)
        covered[e.contig].extend(segments)
        _elements.append((e.ind, [(e.contig, e.orientation, segments)]))
    builder.add_elements(_elements)

    # Add all sources to the engine
    added = set()
    sources = []
    for cmp in config.comparisons:
        for exp in cmp.signal + cmp.control:
            if (cmp.project, exp.ind) in added:
                continue

            source, layout = nfcore.rnaseq.extract.bam(exp, factory=reader)
            sources.append((
                (cmp.project, exp.ind),
                source,
                layout
            ))
            added.add((cmp.project, exp.ind))
            print(f"Added {cmp.project} {exp.ind}")
    print(f"Added {len(added)} sources")

    # # Add all processing regions to the engine
    # partitions = [
    #     (seqid, interval)
    #     for seqid, segments in covered.items()
    #     for interval in Interval.merge_within(segments, distance=16_000)
    # ]
    seqsizes = utils.assembly.seqsizes(config.comparisons[0].organism)
    partitions = [(seqid, (0, seqsizes[seqid])) for seqid in covered]

    builder.add_partitions(partitions)

    engine = builder.build()

    # Run and parse the results
    print(f"Running {config.ind}")
    result = engine.run(sources, countit.rigid.resolution.AnyOverlap())
    del engine
    print(f"Finished {config.ind}")

    counts, _ = countit.utils.result_to_pandas(result)
    return counts.set_index('source').T


counts = {config.ind: job(config) for config in annotation.Config.load()}

saveto = ld.cache.counts.pkl
saveto.parent.mkdir(parents=True, exist_ok=True)
with open(saveto, 'wb') as stream:
    pickle.dump(counts, stream)
