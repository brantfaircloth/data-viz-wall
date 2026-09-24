from collections import Counter
from itertools import pairwise

import numpy as np
import pytest

from data_viz_wall import hybrid, patterns, sim, synteny
from data_viz_wall.server import leaf_order, run_replicate


def collect():
    events = []
    return events, lambda kind, **d: events.append(dict(d, kind=kind))


def test_tree_shape():
    _events, emit = collect()
    rng = np.random.default_rng(0)
    time, children, root = sim.simulate_tree(12, 1000, rng, emit)
    assert len(children) == 11
    assert root == 22
    assert max(sim.leaf_counts(children, root, 12).values()) == 12
    assert all(time[p] > time[c] for p, kids in children.items() for c in kids)


def test_tmrca_expectation():
    """E[TMRCA] = 4Ne(1 - 1/n) generations."""
    rng = np.random.default_rng(1)
    n, ne = 10, 1000
    tmrca = [sim.simulate_tree(n, ne, rng, lambda *a, **k: None)[0][-1] for _ in range(400)]
    assert np.mean(tmrca) == pytest.approx(4 * ne * (1 - 1 / n), rel=0.15)


def test_total_branch_length_expectation():
    """E[L] = 4Ne * sum(1/i) for i in 1..n-1."""
    rng = np.random.default_rng(2)
    n, ne = 8, 1000
    lengths = []
    for _ in range(300):
        events, emit = collect()
        time, children, root = sim.simulate_tree(n, ne, rng, emit)
        sim.drop_mutations(time, children, root, n, ne, 10.0, rng, emit)
        lengths.append(events[-1]["total_length"])
    expected = 4 * ne * sum(1 / i for i in range(1, n))
    assert np.mean(lengths) == pytest.approx(expected, rel=0.15)


def test_watterson_recovers_theta():
    rng = np.random.default_rng(3)
    n, ne, theta = 20, 1000, 25.0
    est = []
    for _ in range(200):
        _events, emit = collect()
        time, children, root = sim.simulate_tree(n, ne, rng, emit)
        spectrum, _ = sim.drop_mutations(time, children, root, n, ne, theta, rng, emit)
        est.append(sim.summarize(spectrum, n, emit)[0])
    assert np.mean(est) == pytest.approx(theta, rel=0.12)


def test_tajimas_d_is_centered():
    rng = np.random.default_rng(4)
    ds = []
    for _ in range(200):
        _events, emit = collect()
        time, children, root = sim.simulate_tree(15, 1000, rng, emit)
        spectrum, _ = sim.drop_mutations(time, children, root, 15, 1000, 20.0, rng, emit)
        ds.append(sim.summarize(spectrum, 15, emit)[2])
    assert abs(np.mean(ds)) < 0.3


def test_leaf_order_is_a_permutation():
    rng = np.random.default_rng(5)
    _, children, root = sim.simulate_tree(16, 500, rng, lambda *a, **k: None)
    assert sorted(leaf_order(children, root, 16)) == list(range(16))


def test_events_carry_source_lines():
    events, order = run_replicate(10, 1000, 15.0, np.random.default_rng(6))
    kinds = {e["kind"] for e in events}
    assert {"wait", "coalesce", "mrca", "spectrum", "stats"} <= kinds
    assert all(e["line"] > 0 for e in events)
    assert len(order) == 10


def test_laplacian_on_a_plane_is_zero():
    x = np.tile(np.arange(16.0), (16, 1))
    inner = patterns.laplacian(x)[1:-1, 1:-1]
    assert np.allclose(inner, 0.0)


def test_laplacian_matches_analytic_sinusoid():
    n = 64
    x = np.arange(n) * 2 * np.pi / n
    a = np.sin(x)[None, :] * np.ones((n, 1))
    q = 2 * np.pi / n
    assert np.allclose(patterns.laplacian(a), -(q**2) * a, atol=1e-4)


def test_react_reaches_a_patterned_steady_state():
    rng = np.random.default_rng(7)
    _events, emit = collect()
    u, v = patterns.seed(96, rng, emit)
    for _ in range(300):
        u, v = patterns.react(u, v, 0.035, 0.065, 20, emit)
    assert np.isfinite(v).all()
    assert 0.0 <= v.min() and v.max() < 0.5
    assert v.std() > 0.05  # structure, not a uniform field


def test_drift_is_monotone_between_endpoints():
    _events, emit = collect()
    a, b = patterns.REGIMES[0], patterns.REGIMES[1]
    fs = [patterns.drift(a, b, p / 20, emit)[0] for p in range(21)]
    assert fs[0] == pytest.approx(a[1]) and fs[-1] == pytest.approx(b[1])
    assert all(x <= y + 1e-12 for x, y in pairwise(fs))


def test_measure_recovers_a_known_wavelength():
    _events, emit = collect()
    n, cycles = 128, 8
    x = np.arange(n) * 2 * np.pi * cycles / n
    field = np.sin(x)[None, :] * np.ones((n, 1))
    assert patterns.measure(field, emit) == pytest.approx(n / cycles, rel=0.02)


def apply_many(steps, seed=11, blocks=200, chroms=6):
    rng = np.random.default_rng(seed)
    _events, emit = collect()
    genome = synteny.ancestral_karyotype(blocks, chroms, emit)
    breaks = Counter()
    for _ in range(steps):
        draw = rng.random()
        if draw < 0.7:
            genome = synteny.inversion(genome, rng, breaks, emit)
        elif draw < 0.85:
            genome = synteny.translocation(genome, rng, breaks, emit)
        elif draw < 0.92:
            genome = synteny.fusion(genome, rng, emit)
        else:
            genome = synteny.fission(genome, rng, breaks, emit)
    return genome, breaks


def test_rearrangement_conserves_every_block():
    """Rearrangements move sequence; they never create or destroy it."""
    genome, _ = apply_many(500)
    assert sorted(abs(b) for c in genome for b in c) == list(range(1, 201))


def test_karyotype_stays_within_bounds():
    genome, _ = apply_many(800, seed=12)
    assert 6 <= len(genome) <= 28
    assert all(len(c) > 0 for c in genome)


def test_inversion_flips_the_segment():
    _events, emit = collect()
    genome = [[1, 2, 3, 4, 5, 6]]
    rng = np.random.default_rng(3)
    for _ in range(20):
        genome = synteny.inversion(genome, rng, Counter(), emit)
    inner = [abs(b) for b in genome[0]]
    assert sorted(inner) == [1, 2, 3, 4, 5, 6]
    assert any(b < 0 for b in genome[0])


def test_syntenic_runs_start_whole_and_only_fragment():
    _events, emit = collect()
    intact = synteny.ancestral_karyotype(200, 6, emit)
    assert len(synteny.syntenic_runs(intact)) == 6
    genome, _ = apply_many(400)
    assert len(synteny.syntenic_runs(genome)) > 6


def test_inverted_run_counts_as_one_block():
    assert synteny.syntenic_runs([[-3, -2, -1]]) == [3]
    assert synteny.syntenic_runs([[1, 2, -9, -8]]) == [2, 2]


def test_contiguity_summary_is_consistent():
    genome, breaks = apply_many(300)
    _events, emit = collect()
    runs = synteny.summarize(genome, breaks, 300, emit)
    report = _events[-1]
    assert sum(runs) == 200
    assert report["largest"] == runs[0]
    assert report["n50"] <= report["largest"]
    assert 0.0 <= report["reuse"] <= 1.0


def run_zone(generations, demes=40, size=30, loci=10, m=0.15, s=0.35, seed=21):
    rng = np.random.default_rng(seed)
    _events, emit = collect()
    anc = hybrid.secondary_contact(demes, size, loci, emit)
    for _ in range(generations):
        anc = hybrid.migrate(anc, m, rng, emit)
        anc = hybrid.reproduce(anc, s, rng, emit)
        anc = hybrid.anchor(anc, emit)
    return anc


def test_ancestry_stays_binary_and_sources_stay_pure():
    anc = run_zone(40)
    assert set(np.unique(anc)) <= {0, 1}
    assert anc[0].sum() == 0
    assert anc[-1].mean() == 1.0


def test_cline_is_monotone_and_spans_both_parents():
    _events, emit = collect()
    p, width = hybrid.cline(run_zone(60), emit)
    assert p[0] == 0.0 and p[-1] == 1.0
    assert 0 < width < len(p)
    assert p[len(p) // 4] < p[3 * len(p) // 4]


def test_selection_narrows_the_cline():
    """Stronger selection against hybrids means a tighter zone."""
    _events, emit = collect()
    weak = hybrid.cline(run_zone(120, s=0.05), emit)[1]
    strong = hybrid.cline(run_zone(120, s=0.60), emit)[1]
    assert strong < weak


def test_migration_widens_the_cline():
    _events, emit = collect()
    low = hybrid.cline(run_zone(120, m=0.05), emit)[1]
    high = hybrid.cline(run_zone(120, m=0.35), emit)[1]
    assert high > low


def test_width_stays_under_the_uncoupled_expectation():
    """Coupling between markers holds the cline tighter than independent loci."""
    _events, emit = collect()
    width = hybrid.cline(run_zone(200), emit)[1]
    assert width < hybrid.expected_width(0.15, 0.35, 10)


def test_hybrid_index_and_heterozygosity_are_proportions():
    _events, emit = collect()
    index, inter = hybrid.sample_hybrids(run_zone(80), 3, emit)
    assert index.min() >= 0.0 and index.max() <= 1.0
    assert inter.min() >= 0.0 and inter.max() <= 1.0
    assert 0.2 < index.mean() < 0.8  # the centre really is intermediate


def test_pure_individuals_carry_no_interclass_heterozygosity():
    anc = hybrid.secondary_contact(10, 5, 8, lambda *a, **k: None)
    assert hybrid.heterozygosity(anc).max() == 0.0
