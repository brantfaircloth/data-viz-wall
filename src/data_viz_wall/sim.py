"""Kingman's coalescent: a genealogy assembled backward in time.

Starting from n sampled chromosomes in the present, we walk into the past.
While k lineages remain, any of the C(k,2) pairs may find a common ancestor
in the previous generation with probability 1/(2Ne), so the waiting time to
the next coalescence is exponential with rate k(k-1)/(4Ne) generations.
Mutations are then sprinkled on the finished tree under the infinite-sites
model at rate theta/2 per unit branch length.
"""

import numpy as np

# ---------------------------------------------------------------- genealogy


def simulate_tree(n, ne, rng, emit):
    """Coalesce n samples back to their most recent common ancestor."""
    time = np.zeros(2 * n - 1)  # node -> age in generations
    children = {}  # internal node -> (left, right)
    lineages = list(range(n))  # what is still uncoalesced
    t = 0.0

    while len(lineages) > 1:
        k = len(lineages)

        # Exponential waiting time: every pair is a candidate ancestor.
        rate = k * (k - 1) / (4 * ne)
        t += rng.exponential(1.0 / rate)
        emit("wait", t=t, k=k, rate=rate)

        # Pick the lucky pair uniformly; all pairs are exchangeable.
        i, j = rng.choice(k, size=2, replace=False)
        left, right = lineages[i], lineages[j]

        parent = n + (n - len(lineages))  # next free internal node id
        time[parent] = t
        children[parent] = (left, right)
        lineages = [x for x in lineages if x not in (left, right)] + [parent]
        emit("coalesce", t=t, k=k - 1, parent=parent, children=[left, right])

    emit("mrca", t=t, root=lineages[0], tmrca=t)
    return time, children, lineages[0]


# ---------------------------------------------------------------- mutations


def drop_mutations(time, children, root, n, ne, theta, rng, emit):
    """Poisson mutations along each branch; each one defines a segregating site."""
    descendants = leaf_counts(children, root, n)
    total_length = 0.0
    spectrum = np.zeros(n, dtype=int)

    for node, (left, right) in sorted(children.items()):
        for child in (left, right):
            length = time[node] - time[child]
            total_length += length

            # Infinite sites: every mutation lands on a fresh position.
            count = rng.poisson(theta * length / (4 * ne))
            if count:
                spectrum[descendants[child]] += count
                emit(
                    "mutate",
                    branch=child,
                    parent=node,
                    count=int(count),
                    carriers=descendants[child],
                    offsets=sorted(rng.random(count).tolist()),
                )

    emit("spectrum", counts=spectrum[1:n].tolist(), total_length=total_length)
    return spectrum, total_length


def leaf_counts(children, root, n):
    """How many samples sit below each node (derived allele count of a branch)."""
    counts = {}

    def walk(node):
        if node < n:
            counts[node] = 1
        else:
            left, right = children[node]
            counts[node] = walk(left) + walk(right)
        return counts[node]

    walk(root)
    return counts


# ---------------------------------------------------------------- statistics


def summarize(spectrum, n, emit):
    """Watterson's theta, nucleotide diversity, and Tajima's D."""
    counts = spectrum[1:n]
    s = int(counts.sum())
    i = np.arange(1, n)

    a1 = float((1.0 / i).sum())
    a2 = float((1.0 / i**2).sum())
    theta_w = s / a1 if a1 else 0.0
    pi = float((counts * i * (n - i)).sum()) / (n * (n - 1) / 2)

    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n**2 + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / a1**2
    var = c1 / a1 * s + c2 / (a1**2 + a2) * s * (s - 1)
    d = (pi - theta_w) / np.sqrt(var) if var > 0 else 0.0

    emit("stats", s=s, theta_w=theta_w, pi=pi, tajimas_d=float(d))
    return theta_w, pi, d
