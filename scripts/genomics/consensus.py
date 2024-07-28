from collections import defaultdict
from math import log

class MotifMatrix:
    """
    Represents a motif matrix as a list of sequences
    with equal length.
    """
    def __init__(self, *seqs):
        self.seqs = self._validate_seqs(*seqs)
        self.nrows = len(self.seqs)
        self.ncols = len(self.seqs[0])
    
    def __str__(self) -> str:
        s = ''
        for i, seq in enumerate(self.seqs, start=1):
            s += f'{i}\t'
            for base in seq:
                s += f'{base}  '
            s += '\n'
        return s

    def _validate_seqs(self, *seqs) -> list[str]:
        """
        Check if all input sequences have the same length.
        """
        processed_seqs = []
        l = len(seqs[0])
        for seq in seqs:
            if len(seq) != l:
                raise ValueError("All sequences must have the same length")
            processed_seqs.append(seq.upper())
        return processed_seqs

    def _transpose(self) -> dict[int, str]:
        """
        Transforms the columns of the matrix into rows.
        """
        seqs_t = defaultdict(list)
        for seq in self.seqs:
            for i in range(self.nrows):
                seqs_t[i].append(seq[i])
        return seqs_t

    def counts(self) -> dict[str, int]:
        """
        Return the nucleotide counts of each column in the matrix.
        """
        M = {n: [] for n in 'ACGT'}
        for row in self._transpose().values():
            for base in 'ACGT':
                M[base].append(row.count(base))
        return dict(M)

    def profile(self) -> dict[str, float]:
        """
        Return the nucleotide proportion of each column in the matrix.
        """
        M = {n: [] for n in 'ACGT'}
        for row in self._transpose().values():
            for base in 'ACGT':
                M[base].append(row.count(base)/self.nrows)
        return dict(M)

    def consensus(self) -> str:
        consensus_str = ''
        pmatrix = self.profile()
        for i in range(len(pmatrix['A'])):
            max_prob = 0
            base_identity = None
            for base in 'ACGT':
                if pmatrix[base][i] > max_prob:
                    max_prob = pmatrix[base][i]
                    base_identity = base
            consensus_str += base_identity
        return consensus_str

    @property
    def entropy(self) -> float:
        score = 0
        for row in self.profile().values():
            for n in row:
                for n == 0: continue
                score += -(n * log(n, 2))
        return score


def test_motif_matrix() -> None:
    motifs = [
    "TCGGGGGTTTTT",
    "CCGGTGACTTAC",
    "ACGGGGATTTTC",
    "TTGGGGACTTTT",
    "AAGGGGACTTCC",
    "TTGGGGACTTCC",
    "TCGGGGATTCAT",
    "TCGGGGATTCCT",
    "TAGGGGAACTAC",
    "TCGGGTATAACC"
    ]
    M = MotifMatrix(*motifs)
    print(M.consensus())

test_motif_matrix()