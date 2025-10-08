import resources
from lib import setup
from workflows.assemblies import HSV1, IAV, CHM13v2

if __name__ == "__main__":
    iav, hsv1, chm13v2 = IAV(), HSV1(), CHM13v2()

    # Concatenate, compress, and index a chimeric FASTA reference
    fasta = [x.files('FASTA') for x in (iav, hsv1, chm13v2)]
    assert all(len(x) == 1 for x in fasta), f"Expected exactly one FASTA file per reference, got {fasta}"

    uncompressed = resources.FASTA.with_suffix(".tmp")
    setup.autocat([x[0] for x in fasta], uncompressed)
    setup.bgzip(uncompressed, saveto=resources.FASTA)
    setup.faidx(resources.FASTA)
    uncompressed.unlink()

    # Concatenate and compress a GFF3 annotation
    gff3 = [iav.files("GFF3"), hsv1.files("GFF3"), chm13v2.files("GENCODE-Liftoff")]
    assert all(len(x) == 1 for x in gff3), f"Expected exactly one GFF3 file per reference, got {gff3}"

    uncompressed = resources.GFF3.with_suffix(".tmp")
    setup.autocat([x[0] for x in gff3], uncompressed)
    setup.gzip(uncompressed, saveto=resources.GFF3)
    uncompressed.unlink()
