import logging
import multiprocessing
import shutil
import subprocess

import resources
from lib import setup
from workflows.assemblies import GRCh38, CHM13v2

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    reference, target = GRCh38(), CHM13v2()

    # 1. Create a target directory for liftoff files
    cache = resources.RESULTS / "liftoff"
    if cache.exists():
        shutil.rmtree(cache)
    cache.mkdir(parents=True, exist_ok=True)

    # 2. Create a list of chromosomes to map the annotation from
    # Chromosomes are matched based on identity of names
    refseqs = set(reference.seqsizes().keys())
    tgtseqs = set(target.seqsizes().keys())
    common = refseqs.intersection(tgtseqs)
    logging.info(f"Common chromosomes: {sorted(common)}")

    chroms = cache / "chroms.txt"
    with chroms.open("w") as f:
        for seq in sorted(common):
            f.write(f"{seq},{seq}\n")
    logging.info(f"Wrote chromosome mapping to {chroms}")

    # 3. Create a list of unplaced contigs
    # Contigs are those sequences in the reference that are not in the target
    exclusive = refseqs - tgtseqs
    logging.info(f"Exclusive reference contigs: {exclusive}")

    contigs = cache / "contigs.txt"
    with contigs.open("w") as f:
        for seq in sorted(exclusive):
            f.write(f"{seq}\n")
    logging.info(f"Wrote contig list to {contigs}")

    # 4. Unzip fasta file for each assembly
    reffa, tgtfa = cache / "reference.fa", cache / "target.fa"
    for fasta, unzipto in [
        (reference.files('FASTA'), reffa), (target.files('FASTA'), tgtfa)
    ]:
        assert len(fasta) == 1, f"Expected exactly one FASTA file for {fasta}"
        fasta = fasta[0]
        assert fasta.name.endswith(".gz") or fasta.name.endswith(".bgz"), f"Expected gzipped FASTA file got {fasta}"
        assert fasta.exists(), f"Expected FASTA file {fasta} does not exist"

        with unzipto.open("w") as f:
            subprocess.check_call(["gunzip", "-k", fasta, "--stdout"], stdout=f)
        logging.info(f"Unzipped {fasta} to {unzipto}")

    # 5. Run liftoff
    refgff = reference.files("GENCODE")
    assert len(refgff) == 1, f"Expected exactly one GENCODE file for {refgff}"
    refgff = refgff[0]

    savemapped = cache / "mapped.gff3"
    saveunmapped = cache / "unmapped.txt"

    cmd = [
        "liftoff",
        "-g", refgff, "-o", savemapped, "-u", saveunmapped,
        "-chroms", chroms, "-unplaced", contigs,
        "-exclude_partial", "-a", "0.95", "-s", "0.95",
        "-copies", "-sc", "0.95", "-polish", "-cds",
        "-p", multiprocessing.cpu_count(),
        tgtfa, reffa
    ]
    cmd = " ".join(map(str, cmd))
    logging.info(f"Running liftoff:\n{cmd}")
    setup.run_in_pixi("liftoff", cwd=cache, cmd=cmd)  # Run inside the pixi environment

    # 6. Save results
    saveto = target.files("GENCODE-Liftoff")
    assert len(saveto) == 1, f"Expected exactly one GENCODE-Liftoff file for {saveto}"
    saveto = saveto[0]

    # Move the mapped file to the final location
    polished = cache / "mapped.gff3_polished"
    assert polished.exists(), f"Expected polished GFF3 file {polished} does not exist"
    if saveto.name.endswith(".gz"):
        with open(saveto, "wb") as f:
            subprocess.check_call(["gzip", "--stdout", polished], stdout=f)
    else:
        polished.rename(saveto)
    logging.info(f"Saved lifted annotation to {saveto}")

    saveunmapped.rename(resources.RESULTS / "liftoff_unmapped.txt")
    logging.info(f"Wrote unmapped annotation to {saveunmapped}")

    # 7. Cleanup
    shutil.rmtree(cache)
