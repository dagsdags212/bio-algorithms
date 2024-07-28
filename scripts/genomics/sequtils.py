def frequency_table(text: str, k: int) -> dict[str, int]:
    """
    Compute for all unique k-mers and count their occurrences. 
    Returns a dictionary with k-mers as keys and counts as values.
    """
    freq_map = {}
    for i in range(len(text)-k+1):
        kmer = text[i:i+k]
        if kmer in freq_map:
            freq_map[kmer] += 1
        else:
            freq_map[kmer] = 1
    return freq_map

def reverse_complement(pattern: str) -> str:
    """
    Returns the reverse complement of a given DNA string.
    """
    comp_map = {"A": "T", "T": "A", "G": "C", "C": "G"}
    comp = ""
    for n in pattern:
        comp = comp_map[n] + comp
    return comp

def hamming_distance(p: str, q: str) -> int:
    """
    Compute for the distance between two strings.
    """
    assert len(p) == len(q)
    dist = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            dist += 1
    return dist

def neighbors(pattern: str, d: int):

    suffix = lambda pattern: pattern[1:]
    first_symbol = lambda pattern: pattern[0]

    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    
    neighborhood = set()
    suffix_neighbors = neighbors(suffix(pattern), d)

    for nb in suffix_neighbors:
        if hamming_distance(suffix(pattern), nb) < d:
            for base in 'ACGT':
                neighborhood.add(base + nb)
        else:
            neighborhood.add(first_symbol(pattern) + nb)
    return neighborhood
