# Path B — Canvas 2D + JS scene layer (fallback)

Taken only if the spike ([07](07-spike-and-decision-log.md)) fails on bundle
size, adapter performance, or SVG export quality.

Prerequisites: [04-plan-common-trunk.md](04-plan-common-trunk.md) complete —
identical prerequisites to Path A, which is why the trunk is worth doing before
deciding.

## Architecture

```
Python: pysam / VirtualBAM / classification / annotation      [trunk, unchanged]
   ↓ payload (trunk T1/T4)
JS core:  coordinate projection (linear + exon-slice)
          row packing / layout
          LOD policy
          interval index + hit-testing
   ↓
Renderers:  Canvas2D (interactive)
            SVG serializer (figures)
            WebGL (opt-in, later, same interface)
   ↓
Hosts:  anywidget (live)  |  static self-contained HTML (share)
```

## Why Canvas 2D rather than WebGL

- Feature counts per viewport are ~10^4, not 10^6. Canvas 2D handles 50–100k
  shapes/frame with `fillStyle` batching and `Path2D` reuse.
- Alignment views are text-heavy — read labels, mismatch letters, base-level
  zoom. `fillText` is trivial; WebGL needs SDF/bitmap font atlases.
- The binding constraint is data transfer, not fill rate.

Keep WebGL as an opt-in layer behind the same renderer interface for specific
cases (deep coverage heatmaps, large dotplots).

## B1 — JS core

The real work, and the part Path A gets for free.

- Coordinate projection, including the nonlinear exon-slice transform with
  per-region width scaling. Path B has an advantage here: the transform is
  yours, so `normalize_interval_width` is straightforward rather than a
  subclass of someone else's block layout.
- Port row packing from `bamtrack.py::layout_interval` (L68–151) to JS. It must
  run per-zoom, client-side — this is the change that makes interactivity
  possible.
- LOD policy: reads → binned coverage above a `bases_per_pixel` threshold.
- Interval index (interval tree or sorted+binary-search) for culling and
  hit-testing.

## B2 — Canvas renderer

Mirrors `svg.py`'s primitive set: `rect`, `line`, `text`,
`text_with_background`, `arrow`, `block_arrow`, clipped groups. Batch by
`fillStyle`. Handle devicePixelRatio explicitly.

Hit-testing via the B1 interval index; fall back to an offscreen picking buffer
if that proves awkward.

## B3 — SVG serializer

Second renderer against the same scene interface, emitting SVG for figures.

**Risk: the two renderers drifting apart.** Mitigate with snapshot tests that
render identical scenes through both and compare. This is a permanent tax that
Path A does not pay, since JBrowse shares Canvas2D functions between screen and
SVG export.

Much of `svg.py` can be ported directly — it is already shaped as a primitive
emitter.

## B4 — Off-main-thread rendering

`OffscreenCanvas` + Web Workers for tile rendering once the single-threaded
version works. Do not do this first.

## B5 — Hosts

- **Static self-contained HTML** — smaller fixed bundle than Path A (no React,
  no MST, no JBrowse core), which is a genuine advantage for the
  many-small-reports case.
- **anywidget** for notebook exploration with a live kernel.
- **Headless figure generation** — reuse B3 server-side by running the scene
  layer in node, or keep the Python SVG renderer from trunk T3.

## B6 — UI chrome

Everything JBrowse would have provided must be built: track selector, location
box, zoom controls, tooltips, feature detail panels, drag-pan, keyboard
navigation. **This is easy to underestimate and is the largest hidden cost of
Path B.**

## Honest comparison

| | Path A | Path B |
|---|---|---|
| Time to first interactive view | weeks | months |
| UI chrome | free | build it all |
| SVG export | free, no drift | build it, permanent drift risk |
| Multi-region view | free | build it |
| Worker architecture | free | build it |
| Exon width normalisation | subclass someone else's layout | straightforward, it's your transform |
| Bundle size | multi-MB | small |
| Dependency risk | JBrowse major versions | none |
| Team skills needed | React + MST | plain TS |
| Bespoke visual control | plugin API limits | total |

Path B wins on bundle size, dependency independence, and control. It loses
badly on time-to-value. It is the right choice only if Path A's constraints
prove genuinely blocking, not merely annoying.
