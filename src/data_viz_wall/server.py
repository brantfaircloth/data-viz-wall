"""Simulate, then replay the simulations' own execution to the wall's browsers."""

import asyncio
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import websockets
from websockets.datastructures import Headers
from websockets.http11 import Response

from . import hybrid, patterns, sim, synteny

HERE = Path(__file__).parent
STATIC = HERE / "static"

MIME = {".html": "text/html", ".css": "text/css", ".js": "text/javascript", ".py": "text/plain"}
ROUTES = {
    "/": "index.html",
    "/coalescent": "viz.html",
    "/coalescent-code": "code.html",
    "/pattern": "pattern.html",
    "/pattern-code": "code.html",
    "/synteny": "synteny.html",
    "/synteny-code": "code.html",
    "/hybrid": "hybrid.html",
    "/hybrid-code": "code.html",
}


def recorder():
    """An emit() that stamps every event with the source line that raised it."""
    events = []

    def emit(kind, **data):
        data["kind"] = kind
        data["line"] = sys._getframe(1).f_lineno
        events.append(data)

    return events, emit


class Feed:
    """A set of clients watching one simulation."""

    def __init__(self):
        self.clients = set()
        self.history = []

    async def register(self, ws):
        self.clients.add(ws)
        for message in self.history:
            await ws.send(message)
        try:
            await ws.wait_closed()
        finally:
            self.clients.discard(ws)

    def send(self, payload, keep=True):
        message = payload if isinstance(payload, bytes) else json.dumps(payload)
        if keep:
            self.history.append(message)
        websockets.broadcast(self.clients, message)


# ----------------------------------------------------------------- coalescent


def run_replicate(n, ne, theta, rng):
    events, emit = recorder()
    time, children, root = sim.simulate_tree(n, ne, rng, emit)
    spectrum, _ = sim.drop_mutations(time, children, root, n, ne, theta, rng, emit)
    sim.summarize(spectrum, n, emit)
    return events, leaf_order(children, root, n)


def leaf_order(children, root, n):
    """Left-to-right sample order that draws the tree without crossing branches."""
    order, stack = [], [root]
    while stack:
        node = stack.pop()
        if node < n:
            order.append(int(node))
        else:
            stack.extend(reversed(children[node]))
    return order


class CoalescentWall(Feed):
    def __init__(self, n, ne, theta, speed, hold, seed):
        super().__init__()
        self.n, self.ne, self.theta = n, ne, theta
        self.speed, self.hold = speed, hold
        self.rng = np.random.default_rng(seed)

    async def loop(self):
        while True:
            events, order = run_replicate(self.n, self.ne, self.theta, self.rng)
            self.history = []
            self.send(
                {"kind": "reset", "n": self.n, "ne": self.ne, "theta": self.theta, "order": order}
            )
            for event in events:
                self.send(event)
                await asyncio.sleep(self.pause(event))
            await asyncio.sleep(self.hold)

    def pause(self, event):
        """Slow down as lineages become scarce; the last coalescence is the drama."""
        if event["kind"] == "coalesce":
            return self.speed * (1.0 + 2.0 / max(event["k"], 1))
        if event["kind"] in ("wait", "mutate"):
            return self.speed * 0.3
        if event["kind"] == "mrca":
            return self.speed * 4.0
        return self.speed


# -------------------------------------------------------------------- Turing


class PatternWall(Feed):
    """Gray-Scott on a torus, drifting between named coat-pattern regimes."""

    def __init__(self, size, steps, fps, dwell, morph, seed):
        super().__init__()
        self.size, self.steps = size, steps
        self.interval = 1.0 / fps
        self.dwell, self.morph = dwell, morph
        self.rng = np.random.default_rng(seed)

    async def loop(self):
        events, emit = recorder()
        u, v = patterns.seed(self.size, self.rng, emit)
        current = patterns.REGIMES[0]
        target = self.pick(current)
        f, k = current[1], current[2]
        elapsed = 0.0
        ticks = 0

        while True:
            events.clear()
            elapsed += self.interval
            ticks += 1
            if elapsed < self.dwell:
                emit("hold", regime=current[0], f=f, k=k)
            elif elapsed < self.dwell + self.morph:
                phase = (elapsed - self.dwell) / self.morph
                f, k = patterns.drift(current, target, phase, emit)
            else:
                current, target, elapsed = target, self.pick(target), 0.0
                f, k = current[1], current[2]

            u, v = patterns.react(u, v, f, k, self.steps, emit)
            if ticks % 40 == 0:
                patterns.measure(v, emit)

            self.announce(current, target, elapsed, events)
            self.send(self.frame(v), keep=False)
            await asyncio.sleep(self.interval)

    def pick(self, current):
        others = [r for r in patterns.REGIMES if r[0] != current[0]]
        return others[int(self.rng.integers(len(others)))]

    def announce(self, current, target, elapsed, events):
        """Only the latest state is kept, so a late client starts up to date."""
        state = {
            "kind": "regime",
            "current": current[0],
            "target": target[0],
            "low": current[3],
            "high": current[4],
            "target_low": target[3],
            "target_high": target[4],
            "size": self.size,
            "elapsed": elapsed,
            "dwell": self.dwell,
            "morph": self.morph,
        }
        self.history = [json.dumps(state)]
        self.send(state, keep=False)
        for event in events:
            self.send(event, keep=False)

    def frame(self, v):
        return np.clip(v * (255.0 / 0.4), 0, 255).astype(np.uint8).tobytes()


# ------------------------------------------------------------------- synteny


class SyntenyWall(Feed):
    """An ancestral karyotype taken apart by inversions, translocations, fusions."""

    MOVES = (("inversion", 0.70), ("translocation", 0.15), ("fusion", 0.07), ("fission", 0.08))

    def __init__(self, blocks, chroms, events, pace, hold, seed):
        super().__init__()
        self.blocks, self.chroms = blocks, chroms
        self.events, self.pace, self.hold = events, pace, hold
        self.rng = np.random.default_rng(seed)

    async def loop(self):
        while True:
            events, emit = recorder()
            breaks = Counter()
            genome = synteny.ancestral_karyotype(self.blocks, self.chroms, emit)
            self.history = []
            self.send(
                {
                    "kind": "reset",
                    "blocks": self.blocks,
                    "chroms": self.chroms,
                    "total_events": self.events,
                }
            )
            self.flush(events, genome)

            for step in range(1, self.events + 1):
                events.clear()
                genome = self.move(genome, breaks, emit)
                synteny.summarize(genome, breaks, step, emit)
                self.flush(events, genome)
                await asyncio.sleep(self.pace)
            await asyncio.sleep(self.hold)

    def move(self, genome, breaks, emit):
        draw = self.rng.random()
        for name, weight in self.MOVES:
            if draw < weight:
                if name == "inversion":
                    return synteny.inversion(genome, self.rng, breaks, emit)
                if name == "translocation":
                    return synteny.translocation(genome, self.rng, breaks, emit)
                if name == "fusion":
                    return synteny.fusion(genome, self.rng, emit)
                return synteny.fission(genome, self.rng, breaks, emit)
            draw -= weight
        return genome

    def flush(self, events, genome):
        for event in events:
            self.send(event, keep=False)
        state = {"kind": "karyotype", "chroms": genome}
        self.history = [json.dumps(state)]
        self.send(state, keep=False)


# -------------------------------------------------------------------- hybrid


class HybridWall(Feed):
    """Two populations in secondary contact, and the cline that settles between them."""

    def __init__(self, demes, size, loci, m, s, pace, generations, hold, seed):
        super().__init__()
        self.demes, self.size, self.loci = demes, size, loci
        self.m, self.s = m, s
        self.pace, self.generations, self.hold = pace, generations, hold
        self.rng = np.random.default_rng(seed)

    async def loop(self):
        while True:
            events, emit = recorder()
            anc = hybrid.secondary_contact(self.demes, self.size, self.loci, emit)
            self.history = []
            self.send(
                {
                    "kind": "reset",
                    "demes": self.demes,
                    "size": self.size,
                    "loci": self.loci,
                    "m": self.m,
                    "s": self.s,
                    "expected": hybrid.expected_width(self.m, self.s, self.loci),
                    "total": self.generations,
                }
            )
            for generation in range(1, self.generations + 1):
                events.clear()
                anc = hybrid.migrate(anc, self.m, self.rng, emit)
                anc = hybrid.reproduce(anc, self.s, self.rng, emit)
                anc = hybrid.anchor(anc, emit)
                hybrid.cline(anc, emit)
                hybrid.sample_hybrids(anc, 3, emit)

                self.send({"kind": "generation", "n": generation}, keep=False)
                for event in events:
                    self.send(event, keep=False)
                self.history = [json.dumps(e) for e in events if e["kind"] in ("cline", "hybrids")]
                await asyncio.sleep(self.pace)
            await asyncio.sleep(self.hold)


# -------------------------------------------------------------------- serving


def static_files(theme):
    """Serve the pages. In light mode wall.css carries its palette overrides."""

    async def handler(connection, request):
        path = request.path.split("?")[0]
        if path.startswith("/ws"):
            return None
        name = ROUTES.get(path, path.lstrip("/"))
        source = HERE / name if name.endswith(".py") else STATIC / name
        if not source.is_file():
            return connection.respond(404, "not found\n")
        body = source.read_bytes()
        if name == "wall.css" and theme == "light":
            body += (STATIC / "light.css").read_bytes()

        headers = Headers(
            {
                "Content-Type": MIME.get(source.suffix, "application/octet-stream"),
                "Content-Length": str(len(body)),
                "Cache-Control": "no-store",
            }
        )
        return Response(200, "OK", headers, body)

    return handler


def query(path):
    return dict(p.split("=", 1) for p in path.partition("?")[2].split("&") if "=" in p)


async def serve(feeds, host, port, theme="dark"):
    async def handler(ws):
        feed = feeds.get(query(ws.request.path).get("feed", "coalescent"))
        if feed:
            await feed.register(ws)

    async with websockets.serve(handler, host, port, process_request=static_files(theme)):
        await asyncio.gather(*(f.loop() for f in feeds.values()))
