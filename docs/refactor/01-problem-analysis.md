# Problem analysis

The performance ceiling is not "SVG is slow". It is three independent problems
that happen to compound. Naming them separately matters, because they have
different fixes and different urgencies.

## Problem 1 — No data/view separation

Geometry is baked directly into presentation strings. `svg.py` functions yield
markup; there is no intermediate representation of "what is where".

Row packing (`bamtrack.py::layout_interval`) calls `scale.topixels()`, so the
layout itself is a function of the current zoom level, and it is computed in
Python.

**Consequence:** pan, zoom, re-sort, and re-filter each require a full Python
round-trip. Interactivity is structurally impossible, not merely slow. No
amount of rendering-substrate work fixes this.

**This is the root problem.** It must be solved under every candidate
architecture, including adopting JBrowse 2.

## Problem 2 — Eager materialisation

`templates.py::plot_sorted_support_as_tabs` renders every (feature ×
classification) SVG string before anything is displayed. All tabs exist in the
DOM simultaneously; `display:none` does not make them free, because the
elements are still live in the tree.

**Consequence:** cost scales with the *total* number of tabs rather than the
number being looked at. This is what makes single-cell cluster splits and tag
splits unusable.

**Note:** this problem largely disappears for free under Path A, since JBrowse
renders only visible blocks of visible tracks.

## Problem 3 — Presentation-encoded payload

A read serialised as SVG markup costs ~500–2000 bytes. The same read as
columnar binary — start, end, row, flags, varint CIGAR, tag indices — costs
~40–80 bytes. That is 20–50× before compression.

**Consequence:** self-contained shareable files hit a size wall far earlier
than they need to. 1.7 MB for one subsampled gene is the symptom.

**This is also the fix for self-containment.** Once the payload is compact, a
single HTML file can carry hundreds of thousands of reads within the practical
browser envelope.

## Ordering

1. Fix Problem 1 first. It unblocks the others and is path-independent.
2. Problem 3 falls out of Problem 1 almost for free — once you have a data
   representation, serialising it compactly is the obvious next step.
3. Problem 2 is fixed by whichever renderer you adopt, or by lazy tab rendering
   as an interim patch on the current SVG path.

## Interim mitigation (optional, ~1–2 weeks)

If the rewrite will take a while and the current pain is acute, these can be
applied to the existing SVG path without architectural change:

- Lazy per-tab rendering (render on first tab activation, cache).
- Level-of-detail: collapse reads to a coverage track when
  `bases_per_pixel` exceeds a threshold.
- Viewport culling of off-screen elements.

This is a bridge, not a destination — the ceiling stays around 10–20k live
nodes. Do not let it absorb effort that belongs in the common trunk.
