import logging
from collections import Counter
from itertools import chain
from pathlib import Path
from typing import Any

from .resources import GRCm39, HSV1, IAV
from .resources.assembly import Assembly, HasGencodeAnnotation, HasRefSeqAnnotation


class ChimericMouseAssembly(Assembly, HasGencodeAnnotation, HasRefSeqAnnotation):
    """
    Chimeric mouse assembly based on GRCm39 with additional HSV-1 and IAV sequences.
    """

    def __init__(self, root: Path):
        super().__init__()
        self.root = root
        try:
            ...
        except FileNotFoundError as err:
            logging.error(f"Failed to load ChimericMouseAssembly: {err}\nDid you forget to setup the assembly?")
            raise err

    @property
    def name(self) -> str:
        return "chimeric-mouse-assembly"

    @property
    def organisms(self) -> frozenset[str]:
        return frozenset(["Mus musculus", "HSV-1", "IAV"])

    def seqsizes(self) -> dict[str, int]:
        cnts = Counter(chain(HSV1.SEGMENTS.keys(), IAV.SEGMENTS.keys(), GRCm39.seqid.sizes().keys()))
        repeats = [k for k, v in cnts.items() if v > 1]
        if repeats:
            raise ValueError(f"Found repeated sequence IDs: {repeats}")
        seqsizes = {**HSV1.SEGMENTS, **IAV.SEGMENTS, **GRCm39.seqid.sizes()}
        return seqsizes

    def annotations(self) -> frozenset[str]:
        raise NotImplementedError()

    def files(self, annotation: str) -> dict[str, Path]:
        raise NotImplementedError()

    def load(self, annotation: str) -> Any:
        raise NotImplementedError()

    @staticmethod
    def setup(root: Path):
        # What do I need to do here?
        # 1. FASTA:
        # - Download GRCm39 primary assembly
        # - Concat HSV-1, IAV, GRCm39 into one FASTA
        # - bgzip it and index with samtools faidx
        # 2. Annotations:
        # - Download GENCODE VM28 annotation
        # - Download RefSeq annotation
        # 3. Link from resources:
        # - REDIPortal annotation
        # - Repeatmasker annotation
        raise NotImplementedError()
