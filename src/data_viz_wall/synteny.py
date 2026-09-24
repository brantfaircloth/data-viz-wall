"""Chromosome evolution: how an ancestral karyotype is taken apart.

A genome is a set of chromosomes, each an ordered, signed list of conserved
blocks.  Inversions reverse a stretch and flip its orientation, reciprocal
translocations swap the tails of two chromosomes, and fusions and fissions
change the chromosome count itself.  None of these destroy sequence -- they
only rearrange it -- so every block keeps the identity it had in the ancestor
and the descendant karyotype can be painted in ancestral colours.

What decays is contiguity.  Watch the number of maximal syntenic runs climb
and their N50 fall, and watch whether breakpoints fall in fresh places or
reuse the same fragile sites over and over.
"""

from itertools import pairwise

import numpy as np

# A block is a signed integer: |s| - 1 is its identity, the sign its strand.
# Two blocks stay syntenic while the next equals the previous plus one, which
# holds for an inverted run (-5, -4, -3) exactly as it does for a forward one.


def ancestral_karyotype(n_blocks, n_chrom, emit):
    """Lay the ancestral blocks out in order across the ancestral chromosomes."""
    cuts = np.linspace(0, n_blocks, n_chrom + 1).astype(int)
    genome = [list(range(a + 1, b + 1)) for a, b in pairwise(cuts)]
    emit("ancestor", n_blocks=n_blocks, n_chrom=n_chrom, cuts=cuts.tolist())
    return genome


def inversion(genome, rng, breaks, emit):
    """Reverse a segment in place and flip the strand of everything inside it."""
    c = weighted_chromosome(genome, rng)
    i, j = sorted(rng.choice(len(genome[c]) + 1, size=2, replace=False))
    if i == j:
        return genome

    record(genome[c], i, breaks)
    record(genome[c], j, breaks)
    genome[c][i:j] = [-b for b in reversed(genome[c][i:j])]
    emit("inversion", chrom=c, start=int(i), end=int(j), span=int(j - i))
    return genome


def translocation(genome, rng, breaks, emit):
    """Reciprocal exchange: two chromosomes trade everything past a cut."""
    if len(genome) < 2:
        return genome
    a, b = rng.choice(len(genome), size=2, replace=False)
    i = int(rng.integers(1, len(genome[a])))
    j = int(rng.integers(1, len(genome[b])))

    record(genome[a], i, breaks)
    record(genome[b], j, breaks)
    genome[a][i:], genome[b][j:] = genome[b][j:], genome[a][i:]
    emit("translocation", chroms=[int(a), int(b)], cuts=[i, j])
    return genome


def fusion(genome, rng, emit, floor=6):
    """Two chromosomes become one; the karyotype contracts."""
    if len(genome) <= floor:
        return genome
    a, b = rng.choice(len(genome), size=2, replace=False)
    merged = genome[a] + genome[b]
    genome = [c for i, c in enumerate(genome) if i not in (a, b)] + [merged]
    emit("fusion", chroms=[int(a), int(b)], length=len(merged), n_chrom=len(genome))
    return genome


def fission(genome, rng, breaks, emit, ceiling=28):
    """One chromosome becomes two; the karyotype expands."""
    if len(genome) >= ceiling:
        return genome
    c = weighted_chromosome(genome, rng)
    if len(genome[c]) < 4:
        return genome
    i = int(rng.integers(2, len(genome[c]) - 1))

    record(genome[c], i, breaks)
    genome = [x for k, x in enumerate(genome) if k != c] + [genome[c][:i], genome[c][i:]]
    emit("fission", chrom=int(c), at=i, n_chrom=len(genome))
    return genome


def weighted_chromosome(genome, rng):
    """Breakpoints land per unit length, so long chromosomes are hit more often."""
    lengths = np.array([len(c) for c in genome], dtype=float)
    return int(rng.choice(len(genome), p=lengths / lengths.sum()))


def record(chromosome, index, breaks):
    """Remember which ancestral adjacency was severed, to detect reuse."""
    if 0 < index < len(chromosome):
        left, right = abs(chromosome[index - 1]), abs(chromosome[index])
        breaks[(min(left, right), max(left, right))] += 1


def syntenic_runs(genome):
    """Maximal stretches still colinear with the ancestor."""
    runs = []
    for chromosome in genome:
        length = 1
        for previous, current in pairwise(chromosome):
            if current == previous + 1:
                length += 1
            else:
                runs.append(length)
                length = 1
        runs.append(length)
    return sorted(runs, reverse=True)


def summarize(genome, breaks, events, emit):
    """Contiguity decay and breakpoint reuse, the two things comparative maps report."""
    runs = syntenic_runs(genome)
    total = sum(runs)

    cumulative, n50 = 0, runs[-1]
    for run in runs:
        cumulative += run
        if cumulative >= total / 2:
            n50 = run
            break

    reuse = sum(1 for n in breaks.values() if n > 1) / max(len(breaks), 1)
    emit(
        "contiguity",
        events=events,
        n_chrom=len(genome),
        runs=len(runs),
        n50=n50,
        largest=runs[0],
        reuse=reuse,
    )
    return runs
