# Wall displays

Four ambient simulations, each built for two large monitors: the render page on
one, the simulator's own Python source on the other with the currently executing
line highlighted. Nothing is pre-recorded and nothing loops a video. Every event
that reaches a browser carries the line number of the `emit()` call that raised
it (`sys._getframe(1).f_lineno`), so the highlight on the source monitor is the
real call site.

```
uv run data-viz-wall
```

| display | source monitor | render monitor |
| --- | --- | --- |
| Kingman coalescent | `/coalescent-code` | `/coalescent` |
| Gray-Scott coat patterns | `/pattern-code` | `/pattern` |
| Karyotype evolution | `/synteny-code` | `/synteny` |
| Hybrid zone | `/hybrid-code` | `/hybrid` |

`/` is a landing page linking to all four. Open a pair fullscreen, one per
monitor, at `http://127.0.0.1:8770`. All four simulations always run; the pages
reconnect on their own, so the server can be restarted without touching the
browsers.

## Running it

Requires [uv](https://docs.astral.sh/uv/) and Python 3.11 or newer. Nothing else
to install and no data files to fetch. `uv run` builds the environment on first
use.

```
git clone <this repo> && cd data-viz-wall
uv run data-viz-wall
```

The server prints all eight routes and then stays in the foreground; Ctrl-C
stops it. One process serves every simulation, so you launch it once regardless
of which display you intend to put on the wall.

Then, on the machine driving the monitors:

1. Open `http://127.0.0.1:8770` in a browser and pick a simulation from the
   landing page. That is the render page. Drag it to the secondary monitor and
   fullscreen it (F11, or ⌃⌘F on macOS).
2. Open the matching `-code` route in a second window, drag it to the primary
   monitor, and fullscreen that too.
3. Leave both alone. The pages reconnect on their own if the server restarts, so
   you can change options and relaunch without touching the browsers.

To show a different simulation, navigate both windows to another pair of routes;
nothing needs restarting. Putting two different simulations on the two monitors
works too, say the coalescent render beside the Turing render, since every feed
runs whether or not anyone is subscribed.

### Server options

| option | default | effect |
| --- | --- | --- |
| `--host` | `127.0.0.1` | Bind address. Set to `0.0.0.0` to drive the wall from another machine on the network. |
| `--port` | `8770` | Port for both the pages and the WebSocket. |
| `--seed` | none | Seeds every simulation, making a run reproducible. Leave unset for a different wall each launch. |
| `--hold` | `8.0` | Seconds to rest on a finished result before the next run. Shared by the coalescent, synteny, and hybrid displays. |
| `--theme` | `dark` | `dark` or `light`. Repaints every page, including everything drawn on canvas. |

```
uv run data-viz-wall --host 0.0.0.0 --port 9000     # drive the wall remotely
uv run data-viz-wall --seed 42                      # reproducible run
uv run data-viz-wall --theme light                  # for a bright room
```

### Light and dark

Large panels in a bright room often handle a white ground better than a black
one, so every page ships both palettes. `--theme light` swaps them for the whole
wall, landing pages and source monitors included.

Both themes are one stylesheet. Every colour, including the ones painted onto
canvas, is a CSS custom property that the drawing code reads back at runtime, so
`light.css` overrides the token values and nothing else changes. The one thing
the flag does not touch is the Gray-Scott palettes: those are pigment colours
belonging to each named regime, and they come from `patterns.py` rather than the
theme.

Changing the theme means restarting the server, since the stylesheet is chosen
at serve time. The open browser windows pick it up on reconnect after a reload.

Per-simulation options are listed with each display below. They are all one flat
namespace on a single command, so options for displays you are not watching are
simply ignored.

---

## 1. Kingman coalescent

**For a visitor.** Every dot along the bottom is a chromosome sampled from a
living population. The tree above them is the family history that connects
them, drawn backward into the past: follow any two dots upward and where their
lines meet is the ancestor they share. The surprising part is the shape. Most
pairings happen quickly, low in the tree, and then the last two lineages sit
there for an enormous stretch before finally meeting. That is not artistic
license. With only two lineages left there is just one pair that can find a
common ancestor, so that final wait takes about as long as everything below it
combined. The orange ticks are mutations landing on branches. One on a long
branch high in the tree ends up in many of the samples; one near the bottom
sits in a single individual. That relationship between where a mutation falls
and how many individuals carry it is the whole basis of reading population
history out of DNA.

**Technically.** n sampled lineages coalesce under Kingman's coalescent: while
k remain, the waiting time to the next event is exponential with rate
k(k−1)/4Ne generations and the coalescing pair is uniform over all C(k,2)
pairs. Mutations are then Poisson along each branch at rate θ/4Ne per
generation under an infinite-sites model, so every mutation defines a new
segregating site whose derived allele count is the number of tips below its
branch. The panel reports S, nucleotide diversity π, Watterson's θ_W = S/a₁,
and Tajima's D with the standard variance normalisation. Leaf order comes from
an in-order traversal sent ahead of the replay so branches never cross without
revealing the topology early, and the time axis is drawn as t^0.62 so the dense
recent coalescences do not pile up on the baseline while the long final wait
still dominates the frame.

### Options

| option | default | effect |
| --- | --- | --- |
| `-n`, `--samples` | 24 | Sampled chromosomes, the tips of the tree. Above ~60 the tips crowd on a 4K monitor. |
| `--ne` | 10,000 | Effective population size. Scales the time axis; TMRCA runs near 4Ne(1 − 1/n) generations. |
| `--theta` | 20.0 | 4Ne·μ for the locus. Sets how many mutations land, so how dense the orange ticks and how tall the spectrum. |
| `--speed` | 0.65 | Base seconds per event. Coalescences stretch further as lineages grow scarce, so the last few are slowest. |
| `--hold` | 8.0 | Seconds on the finished tree before the next replicate. |

```
uv run data-viz-wall -n 40 --theta 40 --speed 0.9    # bigger, denser, slower
uv run data-viz-wall -n 12 --speed 0.3 --hold 3      # quick replicates for a talk
```

Kingman 1982, *Stochastic Processes and their Applications* 13:235.
Watterson 1975, *Theoretical Population Biology* 7:256.
Tajima 1989, *Genetics* 123:585.

---

## 2. Gray-Scott coat patterns

**For a visitor.** Nobody drew those spots. Imagine two substances in
developing skin: one switches pigment on and makes more of itself, the other
switches it off and spreads faster. Wherever the first gets a head start it
reinforces itself into a spot, while the second races outward and stops another
spot forming right next door. Run that rule everywhere at once and evenly
spaced spots appear on their own, or stripes, or a maze, depending on nothing
but the relative rates. Alan Turing worked this out in 1952, before anyone
could look for the molecules involved. The animal is not following a blueprint
with a spot marked at each location; it is running a rule, and the pattern
falls out. The two numbers on the right are those rates, and when they drift
the pattern reorganises into a different animal with no redrawing step in
between.

**Technically.** The Gray-Scott system on a 320² periodic grid, integrated with
explicit Euler steps and a five-point Laplacian: u is consumed by the
autocatalytic reaction u + 2v → 3v, fed at rate F, while v is removed at rate
F + k, with D_u = 0.16 and D_v = 0.08. Short-range activation with faster,
longer-range inhibition is Turing's instability, and (F, k) selects which
branch of the pattern zoo the system settles into. The scheduler holds one
named regime, then eases (F, k) toward the next with a smoothstep over 25
seconds, so the pattern reorganises continuously rather than cutting. The PDE
runs in numpy, 10 steps per frame in about 7 ms, and frames go to the browser
as raw uint8 over the same WebSocket, with the colour LUT and upscaling done
client-side so palettes can cross-fade during a morph. The panel tracks the walk
through (F, k) space and the characteristic wavelength taken from the peak of
the radially averaged 2-D power spectrum.

One honest caveat: the regime names are interpretive. The (F, k) values were
chosen because they look like a python or a guineafowl, not fitted to measured
pigment data. The general claim that reaction-diffusion explains vertebrate
skin patterning is well supported; the specific labels on screen are not a
result.

### Options

| option | default | effect |
| --- | --- | --- |
| `--size` | 320 | Grid edge in cells. Higher means finer pattern and more CPU; cost grows as the square. |
| `--steps` | 10 | Euler steps per rendered frame. Raise it to make the pattern evolve faster without raising the frame rate. |
| `--fps` | 24 | Frames per second sent to the browser. |
| `--dwell` | 45.0 | Seconds held on each named regime. |
| `--morph` | 25.0 | Seconds easing (F, k) to the next regime. Short values snap; long values are hypnotic. |

Roughly `size² × steps × fps` sets the CPU cost. The defaults run about 7 ms per
frame on one core; `--size 512 --steps 16` is still comfortable on a modern
desktop and looks noticeably crisper on a large monitor.

```
uv run data-viz-wall --size 512 --steps 16          # crisper, for a 4K panel
uv run data-viz-wall --dwell 15 --morph 10          # cycle the animals quickly
```

Turing 1952, *Philosophical Transactions of the Royal Society B* 237:37.
Kondo & Asai 1995, *Nature* 376:765.
Kondo & Miura 2010, *Science* 329:1616.
Pearson 1993, *Science* 261:189.

---

## 3. Karyotype evolution

**For a visitor.** The top bar is the genome of an ancestor, the bottom the same
genome tens of millions of years later. Nothing was added and nothing was lost;
it was only shuffled. Colours mark where each piece sat in the ancestor, so
every fragment can still be traced home, which is exactly how a human genome
gets aligned to a chicken's. What degrades is neighbourliness. Early on the
colours sit in clean blocks because long stretches remain in ancestral order,
and as the simulation runs those blocks fragment. The scatter plot on the right
starts as a clean diagonal staircase and dissolves; that diagonal is what "these
two genomes are in the same order" looks like, and watching it break is the
story. How far apart two species' maps have been shuffled is a rough clock on
how long they have been separated.

**Technically.** An ancestral karyotype of 600 conserved blocks on 10
chromosomes is dismantled by inversions (70%), reciprocal translocations (15%),
fusions (7%), and fissions (8%), one event per 1.1 s. A chromosome is an ordered
list of signed integers, where the sign is orientation, so two blocks remain
syntenic exactly when the next equals the previous plus one. That holds for an
inverted run (−5, −4, −3) exactly as for a forward one, and makes maximal
syntenic runs fall straight out of the representation with no special case.
Breakpoints are placed per unit chromosome length, not per chromosome. Severed
ancestral adjacencies are tallied so breakpoint reuse can be reported. The panel
gives run count, run N50, largest run, and reuse; the render page paints both
karyotypes in ancestral colours with ribbons between them and tweens block
positions over 700 ms when an event lands.

Breakpoints here fall at random, so reuse stays near zero. That is the random
breakage null that the fragile breakage model argues against with observed reuse
well above it. Biasing placement toward a fixed set of fragile sites is a few
lines in `weighted_chromosome` and `record` if you want the wall to contrast the
two models.

### Options

| option | default | effect |
| --- | --- | --- |
| `--blocks` | 600 | Conserved blocks in the ancestral genome. More blocks means finer ribbons and a denser Oxford plot. |
| `--karyotype` | 10 | Ancestral chromosome count. Fusions and fissions move it between a floor of 6 and a ceiling of 28. |
| `--rearrangements` | 240 | Events per run, after which the ancestor is rebuilt. Raise it to watch synteny decay further. |
| `--pace` | 1.1 | Seconds per rearrangement. |
| `--hold` | 8.0 | Seconds on the finished karyotype before restarting. |

```
uv run data-viz-wall --rearrangements 600 --pace 0.5   # run contiguity into the ground
uv run data-viz-wall --blocks 1200 --karyotype 24      # a finer, more fragmented genome
```

Nadeau & Taylor 1984, *PNAS* 81:814.
Pevzner & Tesler 2003, *PNAS* 100:7672.
Murphy et al. 2005, *Science* 309:613.

---

## 4. Hybrid zone

**For a visitor.** Two populations meet along a line of habitat and interbreed
where they touch. Migration keeps carrying ancestry from each side into the
other, but hybrids leave fewer offspring, so mixed ancestry keeps getting
removed. Those two forces balance, and the transition between the populations
stops spreading and settles at a fixed width: the white curve, with the orange
bracket marking it. The triangle on the right is one dot per individual living
in the middle of the zone. Left corner is pure one parent, right corner pure the
other, and the apex is a first-generation hybrid with one full chromosome set
from each. Dots fill in below the apex as later generations recombine the two
ancestries into finer and finer mixtures. The width of this zone says nothing
about the habitat; it is set by how far individuals move and how badly hybrids
do.

**Technically.** Individual-based on 80 demes of 50 diploids, each carrying 20
unlinked ancestry markers on two haplotypes, so hybrid index and interclass
heterozygosity are properties of simulated individuals rather than
reconstructions from allele frequencies. Each generation: neighbouring demes
exchange a binomial fraction of residents; fitness is 1 − s × interclass
heterozygosity; parents are drawn by exact weighted sampling using the
Gumbel-max trick, vectorised across all demes at once; gametes take one
haplotype per locus at free recombination; and the two flanking demes are reset
to pure as continuous source populations, without which the zone wanders under
drift. Cline width is 1/max|Δp| across the transect.

The realised width sits below the single-locus expectation σ√(6L/s), roughly 8
demes against 10 at the defaults, and that is the right direction: spreading
selection over L markers makes each one weak, but statistical coupling among
them holds the cline tighter than independent loci would be. The panel labels
that value "uncoupled expectation" rather than a prediction, and a test pins the
inequality.

### Options

| option | default | effect |
| --- | --- | --- |
| `--demes` | 80 | Demes across the transect. Widens the habitat the cline sits in. |
| `--deme-size` | 50 | Individuals per deme. Small values let drift roughen the cline. |
| `--loci` | 20 | Unlinked ancestry markers per individual. More markers means a smoother hybrid index and a better-filled triangle. |
| `--migration` | 0.15 | Fraction exchanged with each neighbour per generation. Dispersal σ = √(2m) demes. **Widens the cline.** |
| `--selection` | 0.35 | Fitness cost of a fully heterozygous hybrid. **Narrows the cline.** |
| `--generations` | 600 | Generations per run before secondary contact restarts. |
| `--tick` | 0.14 | Seconds per generation. About 11 ms of that is computation at the defaults. |
| `--hold` | 8.0 | Seconds on the settled cline before restarting. |

`--migration` and `--selection` are the pair worth demonstrating live: relaunch
with one changed and the zone visibly settles at a different width, while its
position stays wherever it happens to be.

```
uv run data-viz-wall --selection 0.8                   # a hard, narrow zone
uv run data-viz-wall --selection 0.05 --migration 0.3  # a broad, leaky one
uv run data-viz-wall --deme-size 15                    # drift roughens the cline
```

Barton & Hewitt 1985, *Annual Review of Ecology and Systematics* 16:113.
Barton 1983, *Evolution* 37:454.
Barton & Gale 1993, in *Hybrid Zones and the Evolutionary Process*, Oxford.

---

## Architecture

`server.py` holds one `Feed` per simulation. Each runs its own asyncio loop,
broadcasting JSON events to whichever clients subscribed with `?feed=<name>`,
plus raw binary frames in the Gray-Scott case. `code.html` is shared by all four
source monitors: it infers its feed from its own path, fetches the matching
module as plain text, and highlights lines as events arrive, rate-limited so a
fast emitter cannot strobe. Each feed keeps a small history so a client that
connects late paints something immediately instead of waiting for the next
event.

## Tests

`uv run pytest`. The tests check the estimators, not just the plumbing:
E[TMRCA] = 4Ne(1 − 1/n) and E[L] = 4Ne·Σ1/i over replicates, θ_W recovering the
input θ, Tajima's D centred under neutrality, the Laplacian against an analytic
sinusoid, the wavelength measurement against a known one, block conservation
across thousands of rearrangements, inverted runs counting as single syntenic
blocks, and the cline responding the right way to both migration and selection.

## Caveats worth repeating to visitors

None of these display real data. They are live simulations, which is why there
is a source monitor beside each one. In particular, the coalescent tree is the
generic expectation under a neutral model with no selection, no population
structure, and constant size, not any species' actual genealogy.
