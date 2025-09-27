import pandas as pd
import pybedtools
from biobit.core.loc import Interval

import ld
import utils


def job(assembly: str):
    gencode = utils.assembly.get(name=assembly).gencode.load()
    gname2ind = {}
    for gene in gencode.genes.values():
        gname2ind[gene.attrs.name] = gene.ind

    # Curated peak groups
    df = ld.cache.regions.load(assembly)
    mask = df['Ensembl ID'].isna() & (~df['Gene name'].isna())
    df.loc[mask, 'Ensembl ID'] = df.loc[mask, 'Gene name'].apply(lambda x: gname2ind.get(x, None))
    df['Ensembl ID'] = df['Ensembl ID'].apply(lambda x: x if pd.isna(x) else x.split('.')[0])

    reformatted = []
    for intervals in df['intervals']:
        intervals = intervals.replace(",", "").replace("chr", "")

        intervals = [coord.split('-') for coord in intervals.split(';')]
        intervals = [Interval(int(start), int(end)) for start, end in intervals]
        intervals = Interval.merge(intervals)
        intervals = [f"{it.start}-{it.end}" for it in intervals]
        reformatted.append(";".join(intervals))

    df['intervals'] = reformatted
    df = df.sort_values(['seqid', 'intervals'])
    ld.cache.regions.save(assembly, df)

    # Sanity checks
    index = {ind.split(".")[0]: gene for ind, gene in gencode.genes.items()}
    for _, row in df.iterrows():
        gid, gname, strand, location = row['Ensembl ID'], row['Gene name'], row['strand'], row['location']
        if pd.isna(gid):
            continue
        gene = index[gid]

        if gene.attrs.name != gname:
            print(f"Incorrect Ensembl ID & Gene name: {gid} -> {gene.attrs.name} vs {gname}")

        if location == 'Antisense RNA' or location == 'Divergent transcript':
            if gene.loc.strand == strand:
                print(f"Incorrect strand for {gname}")
        elif gene.loc.strand != strand:
            print(f"Incorrect strand for {gname}")

        try:
            # Check locations
            starts, ends = [], []
            for tid in gene.transcripts:
                rna = gencode.rnas[tid]
                if "Ensembl canonical" in rna.attrs.tags:
                    starts.append(rna.loc.start)
                    ends.append(rna.loc.end)
            assert len(starts) == 1, f"N={len(starts)} canonical transcripts for {gname}"
            gstart, gend = starts[0], ends[0]
        except Exception as e:
            print(e)
            continue

        start, *_, end = row['intervals'].split('-')
        start, end = int(start), int(end)
        if end - start > 500_000:
            print(f"Too large interval for {gname}: {end - start}")

        distance = min(abs(start - gstart), abs(start - gend), abs(end - gstart), abs(end - gend))
        if Interval(start, end).intersects(Interval(gstart, gend)):
            distance = 0
        if distance > 100_000:
            print(f"Distance between {gname} and {location} is suspiciously large: {distance}")

        match (location, strand):
            case ('Antisense RNA', '+'):
                is_error = start <= gstart
            case ('Antisense RNA', '-'):
                is_error = start >= gend
            case ('Elongated 3` end', '+'):
                is_error = start <= gend
            case ('Elongated 3` end', '-'):
                is_error = end >= gstart
            case _:
                is_error = False

        if is_error:
            print(f"Incorrect localization for {gname} ({location})")

    duplocations = df.groupby(['Ensembl ID', 'Gene name', 'location']).filter(lambda x: len(x) >= 2)
    if len(duplocations) > 0:
        print("Duplications:")
        print(duplocations)

    # Save as BED
    intervals = []
    for _, row in df.iterrows():
        for coord in row['intervals'].split(';'):
            start, end = coord.split('-')
            start, end = int(start), int(end)
            intervals.append(pybedtools.Interval(
                row['seqid'], start, end, name=f"{row['location']}({row['Gene name']})", strand=row['strand']
            ))
    pybedtools.BedTool(intervals).sort().saveas(ld.cache.regions.root / f"{assembly}.bed")


for name in ("GRCm39", "CHM13v2"):
    job(name)
