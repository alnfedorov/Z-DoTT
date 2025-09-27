import pickle
from dataclasses import dataclass, field
from pathlib import Path

from stories.RIP import pcalling

PKL = Path(__file__).with_suffix(".pkl")


@dataclass(frozen=True, slots=True)
class PeaksConfig:
    min_length: int
    coverage_gaps_tolerance: int

    # Min number of experiments in which a given peak must be supported to be deemed "replicated"
    strong_min_replication: int
    # Min number of experiments in which a given peak must be supported to be deemed "relaxed"
    # Relaxed peaks with editing sites or reverse-complementarity to an adjacent peak are considered replicated
    relaxed_min_editing_sites: int
    relaxed_min_replication: int

    # Paths in BED format
    universe: Path = field(init=False)  # All detected peaks
    prefiltered: Path = field(init=False)  # Peaks replicated in at least N experiments
    filtered: Path = field(init=False)  # Peaks passing the final filter
    covered: Path = field(init=False)  # All genomic regions covered in at least 1 experiment
    not_covered: Path = field(init=False)  # All genomic regions not covered in all experiments
    # Curated
    curated_include: Path = field(init=False)  # include regions
    curated_exclude: Path = field(init=False)  # exclude regions


@dataclass(frozen=True, slots=True)
class dsRNAConfig:
    # Min alignment score for each inverted repeat
    min_score: int = 36
    # Min length of complementary segments
    min_length: int = 36
    # Min total "gap-less" overlap with peaks
    min_roi_overlap: int = 12
    # Max distance for chaining peaks into a single dsRNA prediction region (ROI)
    max_distance: int = 15_000
    # Sequence offset for automatically derived ROIs
    offset: int = 15_000
    # Minimum number of experiments in which a given dsRNA must be supported
    min_replication: int = 2

    # Paths
    repeto: Path = field(init=False)  # Peak groups for repeto
    cache: Path = field(init=False)  # Repeto cache
    predicted: Path = field(init=False)  # All predicted dsRNA
    filtered: Path = field(init=False)  # All dsRNAs passing the filtering
    insulators: Path = field(init=False)  # Curated insulators


@dataclass(frozen=True, slots=True)
class ClusteringConfig:
    # Max distance between bounding boxes of peaks / filtered dsRNAs within a group
    max_distance: int = 5_000

    curated: Path = field(init=False)
    results: Path = field(init=False)


@dataclass(frozen=True, slots=True)
class Config:
    ind: str
    comparisons: tuple[pcalling.Config, ...]

    root: Path
    peaks: PeaksConfig
    dsRNA: dsRNAConfig
    clusters: ClusteringConfig

    @property
    def organism(self) -> set[str]:
        return set.union(*(set(e.organism) for e in self.comparisons))

    @property
    def host(self) -> str:
        hosts = {x.host for x in self.comparisons}
        assert len(hosts) == 1, hosts
        return hosts.pop()

    @property
    def assembly(self) -> str:
        assemblies = {x.assembly for x in self.comparisons}
        assert len(assemblies) == 1, assemblies
        return assemblies.pop()

    def __post_init__(self):
        # Peaks
        object.__setattr__(self.peaks, "universe", self.root / "peaks" / "universe.bed.gz")
        object.__setattr__(self.peaks, "prefiltered", self.root / "peaks" / "prefiltered.bed.gz")
        object.__setattr__(self.peaks, "filtered", self.root / "peaks" / "filtered.bed.gz")
        object.__setattr__(self.peaks, "covered", self.root / "peaks" / "covered.bed.gz")
        object.__setattr__(self.peaks, "not_covered", self.root / "peaks" / "not_covered.bed.gz")
        object.__setattr__(self.peaks, "curated_include", self.root / "curated" / "include.bed")
        object.__setattr__(self.peaks, "curated_exclude", self.root / "curated" / "exclude.bed")
        # dsRNA
        object.__setattr__(self.dsRNA, "repeto", self.root / "dsRNA" / "repeto.bed.gz")
        object.__setattr__(self.dsRNA, "cache", self.root / "dsRNA" / "cache")
        object.__setattr__(self.dsRNA, "predicted", self.root / "dsRNA" / "predicted.pkl")
        object.__setattr__(self.dsRNA, "filtered", self.root / "dsRNA" / "filtered.pkl")
        object.__setattr__(self.dsRNA, "insulators", self.root / "curated" / "insulators.bed")
        # Clusters
        object.__setattr__(self.clusters, "curated", self.root / "curated" / "clusters.bed")
        object.__setattr__(self.clusters, "results", self.root / "clusters.pkl")

    @staticmethod
    def load() -> list["Config"]:
        with open(PKL, 'rb') as stream:
            return pickle.load(stream)
