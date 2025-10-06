from collections.abc import Iterable

from assemblies import GRCm39, CHM13v2, HSV1, IAV


def get(*, name: str | None = None, organism: str | None = None):
    if name:
        return {
            GRCm39.name: GRCm39, CHM13v2.name: CHM13v2, HSV1.name: HSV1, IAV.name: IAV
        }[name]
    else:
        assert organism is not None
        return {
            GRCm39.organism: GRCm39, CHM13v2.organism: CHM13v2, HSV1.organism: HSV1, IAV.organism: IAV
        }[organism]


def seqsizes(organisms: Iterable[str]) -> dict[str, int]:
    result = {}
    for org in organisms:
        assembly = get(organism=org)
        if assembly.name in {"GRCh38", "GRCm39", "CHM13v2"}:
            seqids = assembly.seqid.sizes()
        else:
            seqids = assembly.segments

        for name, length in seqids.items():
            assert name not in result
            result[name] = length

    return result
