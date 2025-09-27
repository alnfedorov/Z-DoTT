import ld

for func, args in [
    (ld.features.repeats, [ld.annotation.repeats]),
    (ld.features.miRNA, [ld.annotation.miRNA]),
    (ld.features.orfs, [ld.annotation.orfs]),
    (ld.features.genes, [ld.annotation.genes]),
    (ld.features.primers, [ld.annotation.primers]),
]:
    func(*args)

ld.features.circos(ld.annotation.genes, ld.annotation.circos)
