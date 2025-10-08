def from_ucsc(ucsc: str) -> str | None:
    return {
        "chr1": "chr1",
        "chr2": "chr2",
        "chr3": "chr3",
        "chr4": "chr4",
        "chr5": "chr5",
        "chr6": "chr6",
        "chr7": "chr7",
        "chr8": "chr8",
        "chr9": "chr9",
        "chr10": "chr10",
        "chr11": "chr11",
        "chr12": "chr12",
        "chr13": "chr13",
        "chr14": "chr14",
        "chr15": "chr15",
        "chr16": "chr16",
        "chr17": "chr17",
        "chr18": "chr18",
        "chr19": "chr19",
        "chr20": "chr20",
        "chr21": "chr21",
        "chr22": "chr22",
        "chrM": "chrM",
        "chrX": "chrX",
        "chrY": "chrY",
    }.get(ucsc, None)


def from_ensembl(ensembl: str) -> str:
    return {
        '1': 'chr1',
        '2': 'chr2',
        '3': 'chr3',
        '4': 'chr4',
        '5': 'chr5',
        '6': 'chr6',
        '7': 'chr7',
        '8': 'chr8',
        '9': 'chr9',
        '10': 'chr10',
        '11': 'chr11',
        '12': 'chr12',
        '13': 'chr13',
        '14': 'chr14',
        '15': 'chr15',
        '16': 'chr16',
        '17': 'chr17',
        '18': 'chr18',
        '19': 'chr19',
        '20': 'chr20',
        '21': 'chr21',
        '22': 'chr22',
        'MT': 'chrM',
        'X': 'chrX',
        'Y': 'chrY'
    }[ensembl]


def from_refseq(ncbi: str) -> str:
    return {
        "NC_060925.1": "chr1",
        "NC_060926.1": "chr2",
        "NC_060927.1": "chr3",
        "NC_060928.1": "chr4",
        "NC_060929.1": "chr5",
        "NC_060930.1": "chr6",
        "NC_060931.1": "chr7",
        "NC_060932.1": "chr8",
        "NC_060933.1": "chr9",
        "NC_060934.1": "chr10",
        "NC_060935.1": "chr11",
        "NC_060936.1": "chr12",
        "NC_060937.1": "chr13",
        "NC_060938.1": "chr14",
        "NC_060939.1": "chr15",
        "NC_060940.1": "chr16",
        "NC_060941.1": "chr17",
        "NC_060942.1": "chr18",
        "NC_060943.1": "chr19",
        "NC_060944.1": "chr20",
        "NC_060945.1": "chr21",
        "NC_060946.1": "chr22",
        "NC_060947.1": "chrX",
        "NC_060948.1": "chrY",
    }[ncbi]


def all() -> tuple[str, ...]:
    return (
        "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13",
        "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM",
    )


def sizes() -> dict[str, int]:
    return {
        "chr1": 248387328,
        "chr2": 242696752,
        "chr3": 201105948,
        "chr4": 193574945,
        "chr5": 182045439,
        "chr6": 172126628,
        "chr7": 160567428,
        "chr8": 146259331,
        "chr9": 150617247,
        "chr10": 134758134,
        "chr11": 135127769,
        "chr12": 133324548,
        "chr13": 113566686,
        "chr14": 101161492,
        "chr15": 99753195,
        "chr16": 96330374,
        "chr17": 84276897,
        "chr18": 80542538,
        "chr19": 61707364,
        "chr20": 66210255,
        "chr21": 45090682,
        "chr22": 51324926,
        "chrX": 154259566,
        "chrY": 62460029,
        "chrM": 16569,
    }
