from collections import defaultdict
from math import log
import numpy as np


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
        n = len(seqs[0])
        for seq in seqs:
            if len(seq) != n:
                raise ValueError("All sequences must have the same length")
            processed_seqs.append(seq.upper())
        return processed_seqs

    # def score(self, pattern: str) -> int:
    #     """
    #     Compute the hamming distance between `pattern` and
    #     each string in a list of sequnces.
    #     """
    #     k = len(pattern)
    #     dist = 0
    #     for seq in self.seqs:
    #         hd_min = float('inf')
    #         for i in range(len(seq)-k+1):
    #             kmer = seq[i:i+k]
    #             hd_curr = self.hamming_distance(kmer, pattern)
    #             if hd_min > hd_curr:
    #                 hd_min = hd_curr
    #         dist += hd_min
    #     return dist


class CountMatrix:
    """
    Represents a count matrix derived from a list
    of sequences with equal length.
    """
    def __init__(self, motifs: MotifMatrix, pseudocounts: bool=True) -> None:
        self.motifs = motifs
        self.pseudocounts = pseudocounts
        self.C = self._generate_counts()

    def __str__(self) -> str:
        fmt = ''
        for base, counts in self.C.items():
            c_str = '\t'.join(map(str, counts))
            fmt += f'[{base}]\t{c_str}\n'
        return fmt

    def _generate_counts(self) -> dict[str, int]:
        """
        Return the nucleotide counts of each column in the matrix.
        If pseudocounts is set to True, each position in the matrix
        is incremented by 1.
        """
        if self.pseudocounts:
            C = {base: np.ones(self.motifs.ncols, dtype=np.int16) for base in 'ACGT'}
        else:
            C = {base: np.zeros(self.motifs.ncols, dtype=np.int16) for base in 'ACGT'}

        for i in range(self.motifs.nrows):
            for j in range(self.motifs.ncols):
                base = self.motifs.seqs[i][j]
                C[base][j] += 1
        return C


class ProfileMatrix(CountMatrix):
    """
    Represents a profile matrix dervied from a count matrix.
    """
    def __init__(self, motifs: MotifMatrix, **kwargs) -> None:
        super().__init__(motifs, pseudocounts=kwargs.get('pseudocounts', True))
        self.P = self._generate_profile()

    def _generate_profile(self) -> dict[str, float]:
        """
        Retrun the relative proportion of each nucleotide at each
        index of a motif matrix.
        """
        # Vectorized function to round to 4 decimals.
        round_4 = lambda x: round(x, 4)
        round_4_vec = np.vectorize(round_4)

        P = self.C.copy()
        denom = self.motifs.nrows
        if self.pseudocounts:
            denom += 4
        for base, counts in self.C.items():
            P[base] = round_4_vec(counts / denom)
        return P

    def __str__(self) -> str:
        fmt = ''
        for base, counts in self.P.items():
            p_str = '\t'.join(map(str, counts))
            fmt += f'[{base}]\t{p_str}\n'
        return fmt

    @property
    def consensus(self) -> str:
        consensus_str = ''
        m = len(self.P['A'])
        for i in range(m):
            max_prob = 0
            base_identity = None
            for base in 'ACGT':
                if self.P[base][i] > max_prob:
                    max_prob = self.P[base][i]
                    base_identity = base
            consensus_str += base_identity
        return consensus_str

    @property
    def entropy(self) -> float:
        score = 0
        for row in self.P.values():
            for n in row:
                if n == 0: continue
                score += -(n * log(n, 2))
        return round(score, 4)

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
    cm = CountMatrix(M)
    pm = ProfileMatrix(M)
    print(pm)
    print(pm.consensus)
    print(pm.entropy)

test_motif_matrix()