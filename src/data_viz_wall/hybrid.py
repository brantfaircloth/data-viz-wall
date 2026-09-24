"""A hybrid zone held in place by selection against hybrids.

Two populations meet along a line of demes.  Every individual carries L
unlinked ancestry markers on two haplotypes, so we can ask of anyone in the
zone both what fraction of their genome came from each parent (hybrid index)
and how much of it is still one copy of each (interclass heterozygosity).

Migration mixes ancestry outward; selection removes the mixed individuals.
Where the two balance, the cline stops spreading and settles at a width that
depends only on dispersal and the strength of selection -- Barton and Hewitt's
tension zone, which sits wherever it likes and tells you nothing about the
environment, only about the incompatibility.
"""

import numpy as np

# anc[deme, individual, locus, haplotype] == 1 for ancestry from the right parent.


def secondary_contact(demes, size, loci, emit):
    """Two pure populations, previously separate, brought into contact."""
    anc = np.zeros((demes, size, loci, 2), dtype=np.uint8)
    anc[demes // 2 :] = 1
    emit("contact", demes=demes, size=size, loci=loci)
    return anc


def migrate(anc, m, rng, emit):
    """Stepping stone: neighbouring demes exchange a fraction m of their residents."""
    demes, size = anc.shape[0], anc.shape[1]
    k = rng.binomial(size, m)
    for d in range(demes - 1):
        here = rng.choice(size, size=k, replace=False)
        there = rng.choice(size, size=k, replace=False)
        swap = anc[d, here].copy()
        anc[d, here] = anc[d + 1, there]
        anc[d + 1, there] = swap
    emit("migrate", swapped=int(k), sigma=float(np.sqrt(2 * m)))
    return anc


def heterozygosity(anc):
    """Fraction of loci carrying one copy of each ancestry -- the cost of being a hybrid."""
    return (anc[..., 0] != anc[..., 1]).mean(axis=2)


def reproduce(anc, s, rng, emit):
    """Weighted by fitness, then one haplotype drawn from each parent at every locus."""
    demes, size, loci, _ = anc.shape
    fitness = 1.0 - s * heterozygosity(anc)

    # Gumbel-max: exact weighted sampling of parents, vectorised across demes.
    keys = np.log(fitness)[:, None, :] + rng.gumbel(size=(demes, 2 * size, size))
    parents = keys.argmax(axis=2).reshape(demes, size, 2)
    emit("select", mean_fitness=float(fitness.mean()), load=float(1 - fitness.min()))

    rows = np.arange(demes)[:, None]
    child = np.empty_like(anc)
    for side in (0, 1):
        gametes = anc[rows, parents[..., side]]
        draw = rng.integers(0, 2, size=(demes, size, loci, 1))
        child[..., side] = np.take_along_axis(gametes, draw, axis=3)[..., 0]
    emit("recombine", loci=loci)
    return child


def anchor(anc, emit):
    """The ends are continuous source populations, still pure."""
    anc[0] = 0
    anc[-1] = 1
    emit("anchor", left=0, right=1)
    return anc


def cline(anc, emit):
    """Ancestry frequency deme by deme, and the width of the step between them."""
    p = anc.mean(axis=(1, 2, 3))
    per_locus = anc.mean(axis=(1, 3))

    slope = np.abs(np.diff(p)).max()
    width = 1.0 / slope if slope > 0 else float(len(p))
    centre = float(np.interp(0.5, p, np.arange(len(p)))) if p[0] < p[-1] else len(p) / 2

    emit("cline", p=p.tolist(), width=float(width), centre=centre, loci=per_locus.T.tolist())
    return p, width


def sample_hybrids(anc, span, emit):
    """Hybrid index against interclass heterozygosity for everyone near the centre."""
    middle = anc.shape[0] // 2
    core = anc[max(middle - span, 0) : middle + span + 1]

    index = core.mean(axis=(2, 3)).ravel()
    inter = heterozygosity(core).ravel()

    emit(
        "hybrids",
        index=index.round(3).tolist(),
        inter=inter.round(3).tolist(),
        f1_like=int(((index > 0.4) & (index < 0.6) & (inter > 0.7)).sum()),
    )
    return index, inter


def expected_width(m, s, loci):
    """Single-locus tension zone width, sigma * sqrt(6 / s_locus), with s spread
    evenly over the L markers.  It is an upper bound here: statistical coupling
    between loci makes the realised cline narrower than independent loci would be."""
    return np.sqrt(2 * m) * np.sqrt(6.0 * loci / s)
