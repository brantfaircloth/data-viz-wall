import asyncio

import click

from .server import CoalescentWall, HybridWall, PatternWall, SyntenyWall, serve


@click.command()
@click.option("-n", "--samples", default=24, show_default=True, help="Sampled chromosomes.")
@click.option("--ne", default=10_000, show_default=True, help="Effective population size.")
@click.option("--theta", default=20.0, show_default=True, help="4*Ne*mu for the locus.")
@click.option("--speed", default=0.65, show_default=True, help="Base seconds per event.")
@click.option("--hold", default=8.0, show_default=True, help="Seconds to admire the result.")
@click.option("--size", default=320, show_default=True, help="Turing grid edge in cells.")
@click.option("--steps", default=10, show_default=True, help="Euler steps per rendered frame.")
@click.option("--fps", default=24, show_default=True, help="Turing frames per second.")
@click.option("--dwell", default=45.0, show_default=True, help="Seconds held on each regime.")
@click.option("--morph", default=25.0, show_default=True, help="Seconds morphing between them.")
@click.option("--blocks", default=600, show_default=True, help="Conserved synteny blocks.")
@click.option("--karyotype", default=10, show_default=True, help="Ancestral chromosomes.")
@click.option("--rearrangements", default=240, show_default=True, help="Events per run.")
@click.option("--pace", default=2.0, show_default=True, help="Seconds per rearrangement.")
@click.option("--demes", default=80, show_default=True, help="Demes across the transect.")
@click.option("--deme-size", default=50, show_default=True, help="Individuals per deme.")
@click.option("--loci", default=20, show_default=True, help="Unlinked ancestry markers.")
@click.option("--migration", default=0.15, show_default=True, help="Exchanged per neighbour.")
@click.option("--selection", default=0.35, show_default=True, help="Cost of a full hybrid.")
@click.option("--generations", default=600, show_default=True, help="Generations per run.")
@click.option("--tick", default=0.14, show_default=True, help="Seconds per generation.")
@click.option("--seed", default=None, type=int, help="Random seed.")
@click.option(
    "--theme",
    type=click.Choice(["dark", "light"]),
    default="dark",
    show_default=True,
    help="Palette for every page.",
)
@click.option("--host", default="127.0.0.1", show_default=True)
@click.option("--port", default=8770, show_default=True)
def main(
    samples,
    ne,
    theta,
    speed,
    hold,
    size,
    steps,
    fps,
    dwell,
    morph,
    blocks,
    karyotype,
    rearrangements,
    pace,
    demes,
    deme_size,
    loci,
    migration,
    selection,
    generations,
    tick,
    seed,
    theme,
    host,
    port,
):
    """Serve the wall: a source monitor and a rendering monitor for each display."""
    feeds = {
        "coalescent": CoalescentWall(samples, ne, theta, speed, hold, seed),
        "pattern": PatternWall(size, steps, fps, dwell, morph, seed),
        "synteny": SyntenyWall(blocks, karyotype, rearrangements, pace, hold, seed),
        "hybrid": HybridWall(
            demes, deme_size, loci, migration, selection, tick, generations, hold, seed
        ),
    }
    for label, route in (
        ("coalescent", "coalescent-code"),
        ("", "coalescent"),
        ("turing", "pattern-code"),
        ("", "pattern"),
        ("synteny", "synteny-code"),
        ("", "synteny"),
        ("hybrid zone", "hybrid-code"),
        ("", "hybrid"),
    ):
        click.echo(f"{label:>11}  http://{host}:{port}/{route}")
    asyncio.run(serve(feeds, host, port, theme))
