"""Gray-Scott reaction-diffusion: how two morphogens paint a vertebrate.

Substrate u is consumed by the autocatalytic activator v in the reaction
u + 2v -> 3v.  u is fed in at rate F, v is drained at rate F + k, and both
diffuse -- but v spreads half as fast as u.  Short-range activation with
long-range inhibition is exactly Turing's instability, and the small set of
patterns it settles into is the set vertebrate skin actually uses: spots,
stripes, labyrinths, reticulation.  Walking (F, k) slowly across parameter
space walks the pattern from one animal to the next.
"""

import numpy as np

DU, DV, DT = 0.16, 0.08, 1.0  # diffusion of u, diffusion of v, timestep

# (name, F, k, low colour, high colour)
REGIMES = [
    ("reticulated python", 0.0290, 0.0570, "#1a1206", "#d8b168"),
    ("leopard rosettes", 0.0350, 0.0650, "#241a08", "#f0c04a"),
    ("giraffe reticulation", 0.0390, 0.0580, "#2a1608", "#c9762f"),
    ("pufferfish vermiculation", 0.0540, 0.0630, "#07161c", "#57d6c8"),
    ("zebra flank", 0.0140, 0.0450, "#0c0c0f", "#ece7dd"),
    ("guineafowl spotting", 0.0367, 0.0649, "#0a0f18", "#9fb6d8"),
]


def laplacian(a):
    """Five-point stencil on a periodic grid -- the skin is a torus."""
    return np.roll(a, 1, 0) + np.roll(a, -1, 0) + np.roll(a, 1, 1) + np.roll(a, -1, 1) - 4.0 * a


def react(u, v, f, k, steps, emit):
    """Advance the pair of PDEs by `steps` explicit Euler steps."""
    for step in range(steps):
        uvv = u * v * v  # the autocatalytic term, u + 2v -> 3v

        u += DT * (DU * laplacian(u) - uvv + f * (1.0 - u))
        v += DT * (DV * laplacian(v) + uvv - (f + k) * v)

        if step == steps - 1:
            emit("react", f=f, k=k, activator=float(v.mean()))
    return u, v


def seed(size, rng, emit):
    """Start from a saturated substrate disturbed by a few drops of activator."""
    u = np.ones((size, size))
    v = np.zeros((size, size))
    for _ in range(24):
        r, c = rng.integers(0, size - 8, size=2)
        u[r : r + 8, c : c + 8] = 0.50
        v[r : r + 8, c : c + 8] = 0.25
    v += 0.02 * rng.random((size, size))
    emit("seed", size=size)
    return u, v


def drift(start, target, phase, emit):
    """Ease (F, k) from one animal's corner of parameter space toward the next."""
    s = phase * phase * (3.0 - 2.0 * phase)
    f = start[1] + (target[1] - start[1]) * s
    k = start[2] + (target[2] - start[2]) * s
    emit("drift", toward=target[0], phase=phase, f=f, k=k)
    return f, k


def measure(v, emit):
    """Characteristic wavelength: the peak of the radially averaged power spectrum."""
    power = np.abs(np.fft.fftshift(np.fft.fft2(v - v.mean()))) ** 2
    n = v.shape[0]
    y, x = np.indices(power.shape) - n // 2
    radius = np.hypot(y, x).astype(int)

    profile = np.bincount(radius.ravel(), power.ravel()) / np.bincount(radius.ravel())
    q = int(np.argmax(profile[1 : n // 2])) + 1  # cycles across the grid
    emit("wavelength", cells=float(n / q), power=profile[1 : n // 2].tolist())
    return n / q
