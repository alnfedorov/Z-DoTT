import pickle

import RNA
from biobit.core.loc import Interval
from joblib import Parallel, delayed

import ld
from assemblies import HSV1
from utils import fasta


def job(seqstart, seqend, strand, saveto):
    seq = fasta.sequence(HSV1.fasta, HSV1.contig, seqstart, seqend, strand=strand)
    seq = seq.replace("T", "U")
    assert len(seq) >= 1

    fold, _ = RNA.fold(seq)

    # Derive stems
    stack, stems = [], []
    for pos, symbol in enumerate(fold):
        if symbol == '(':
            stack.append(pos)
        elif symbol == ')' and stems and stems[-1][1][1] == pos:
            prevleft, prevright = stems[-1]
            newleft = stack.pop()

            # Continue previous loop
            if newleft == prevleft[0] - 1:
                stems[-1] = (
                    (newleft, prevleft[1]), (prevright[0], pos + 1)
                )
            # Bulge -> start new loop
            else:
                stems.append((
                    (newleft, newleft + 1), (pos, pos + 1)
                ))
        elif symbol == ')':
            left = stack.pop()
            left = (left, left + 1)
            right = (pos, pos + 1)
            stems.append((left, right))

    # Sanity checks
    assert fold.count("(") + fold.count(")") == sum(left[1] - left[0] for left, _ in stems) * 2
    for left, right in stems:
        assert left[1] - left[0] == right[1] - right[0], (left, right)
        assert set(fold[left[0]: left[1]]) == {'('}
        assert set(fold[right[0]: right[1]]) == {')'}

    print(
        f'{seqstart}-{seqend}({strand}) -> {saveto}\n'
        f'{seq}\n'
        f'{fold}\n'
        f'stems={len(stems)}'
    )

    stemit = []
    for left, right in stems:
        stemit.append((Interval(*left), Interval(*right)))

    with open(saveto, 'wb') as stream:
        pickle.dump((seqstart, seqend, strand, seq, fold, stemit), stream)


ld.MFE.root.mkdir(parents=True, exist_ok=True)
Parallel(n_jobs=-1, verbose=100)(
    delayed(job)(start, end, strand, ld.MFE.root / f"{title}.raw-fold.pkl")
    for start, end, strand, title in [
        (108260, 115846, '-', 'UL54-5'),
        (127170, 131430, '-', 'RS1'),
        (35004, 41614, '-', 'UL20'),
        (35004, 40767, '-', 'UL19'),
        (101286, 108176, '+', 'UL46-5'),
        (103415, 108176, '+', 'UL47-5'),
    ]
)
