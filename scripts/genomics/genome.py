from pathlib import Path
import matplotlib.pyplot as plt
from collections import defaultdict
from .sequtils import *


class Genome:
    def __init__(self, seq: str, header: str=None, source: str=None) -> None:
        self.seq = seq.upper()
        self.header = header
        self.source = source
    
    @property
    def composition(self) -> dict[str, float]:
        """
        Compute the relative abundance of each nucleotide
        in the genome sequence.

        Return: a dict with bases a keys and proportions as values.
        """
        l = len(self.seq)
        return {base: (self.seq.count(base) / l * 100) for base in 'ACGT'}

    def find_clumps(self, k: int, L: int, t: int) -> list[str]:
        """
        Return a list of k-mers that occur at least t-times within subsets with length L.

        k: (int) k-mer length.
        L: (int) subset length.
        t: (int) minimumum number of occurrences.
        """
        patterns = []
        n = len(self.seq)
        for i in range(n-L+1):
            window = self.seq[i:i+L]
            freq_map = frequency_table(window, k)
            for kmer in freq_map:
                if freq_map[kmer] >= t:
                    patterns.append(kmer)
        patterns = list(set(patterns))
        return patterns

    def skew_array(self) -> list[int]:
        """
        Compute the running difference between cytosine cout and guanine count
        in a given genome string.
        """
        counts = [0]
        for base in self.seq:
            val = counts[-1]
            if base == "G":
                val += 1
            elif base == "C":
                val -= 1
            counts.append(val)
        # The counts array should be longer than the genome length by one element.
        assert len(counts) == len(self.seq) + 1
        return counts

    def skew_diagram(self) -> None:
        """
        Plot values from a skew array with x-values corresponding to genome position
        and y-values corresponding to skew value.
        """
        skew_arr = self.skew_array()
        positions = [*range(len(skew_arr))]

        # Positions of orisite and tersite, respectively.
        min_count, max_count = min(skew_arr), max(skew_arr)
        ori = (positions[skew_arr.index(min_count)], min_count)
        ter =  (positions[skew_arr.index(max_count)], max_count)

        # Plot skew diagram and label ori/ter sites.
        plt.plot(positions, skew_arr)
        plt.annotate(
            text="ori",
            xy=ori,
            xytext=(ori[0], ori[1]*1.1),
            ha="center"
        )
        plt.annotate(
            text="ter",
            xy=ter,
            xytext=(ter[0], ter[1]*1.02),
            ha="center"
        )
        plt.title(f"Skew diagram of {self.source}")
        plt.ylabel("skew")
        plt.xlabel("position")
        plt.xticks([0, 1e6, 2e6, 3e6, 4e6, 5e6], map(int, [0, 1e6, 2e6, 3e6, 4e6, 5e6]))
        plt.show()

    def frequent_words_mismatches_reverse_complements(self, k: int, d: int) -> list[str]:
        kmer_count = defaultdict(int)
        for i in range(len(self.seq)-k+1):
            kmer_count[self.seq[i:i+k]] += 1
            kmer_count[reverse_complement(self.seq[i:i+k])] += 1

        mismatch_count = defaultdict(int)
        for kmer, count in kmer_count.items():
            for nb in neighbors(kmer, d):
                mismatch_count[nb] += count

        max_count = max(mismatch_count.values())
        return sorted([kmer for kmer, count in mismatch_count.items() if count == max_count])

def load_genome(filepath: Path, **kwargs) -> Genome:
    try:
        with open(filepath, "r") as file:
            header = None
            genome_str = ""
            for line in file.readlines():
                line = line.rstrip()
                if line.startswith(">"):
                    header = line
                else:
                    genome_str += line.upper()
            return Genome(genome_str, header=header, **kwargs)
    except FileNotFoundError as e:
        print(f'An error has occurred: {e}')