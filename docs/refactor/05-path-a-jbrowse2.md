# Path A — JBrowse 2 plugin (preferred)

Prerequisites: [04-plan-common-trunk.md](04-plan-common-trunk.md) complete,
spike passed ([07](07-spike-and-decision-log.md)).

Capability findings and confidence levels:
[03-jbrowse2-assessment.md](03-jbrowse2-assessment.md).

## Shape of the deliverable

```
itv-jbrowse-plugin/          (TypeScript)
    adapters/
        ITVReadsAdapter          — serves a readset from the inlined payload
        ITVCoverageAdapter       — serves precomputed per-split coverage
        ITVAnnotationAdapter     — BED/GTF features from the payload
    displays/
        ITVIsoformPairedDisplay  — isoform model + its aligned reads
    views/
        ITVExonSliceView         — LGV subclass with per-region width scaling
        ITVOverviewView          — optional, nice-to-have
    index.ts

itv/  (Python, in this repo)
    report.py                — payload → self-contained HTML
    figure.py                — payload → SVG via headless render
```

## A1 — Plugin skeleton and in-memory adapter

Start from `GMOD/jbrowse-plugin-template`. First milestone: a custom adapter
extending `BaseFeatureDataAdapter`, `getFeatures()` returning an RxJS
observable over reads decoded from an inlined payload, rendered in
`@jbrowse/react-app2`.

This single milestone validates the riskiest integration point and is spike
question 4.

## A2 — Classification splits

**Recommended idiom: N tracks in a categorised track selector, not tabs.**

Generate one track per classification split from Python, each backed by its own
`ITVReadsAdapter` instance over the corresponding `VirtualBAM`-derived readset,
grouped with track `category` so they collapse into a tree.

Paired with each read track, an `ITVCoverageAdapter` track for that same split
— this covers the "coverage track per split" requirement directly, since the
payload already carries precomputed per-readset coverage.

Alternatives considered:
- `jexl` `filter` + `color` callbacks over a single track. Simpler config, but
  it pushes filtering to render time and does not give per-split coverage as
  naturally. Good for quick ad-hoc colouring, not for the primary split axis.
- **Multiple JBrowse instances inside existing ITV tabs.** Viable fallback but
  expensive — each instance carries its own worker pool, session state and
  render loop. Prefer one instance with many tracks. Only reach for this if
  track counts become unmanageable in the selector UI.
- Own tab chrome above `react-app2`, flipping track visibility rather than
  remounting. Keeps one instance and one worker pool while preserving the tab
  metaphor users know. **This is the compromise to reach for if track
  categories test badly with users.**

Worth validating with an actual user (spike question 5): categories may be
*better* than tabs here, because comparing cell types side by side beats
flipping between them.

## A3 — Exon slices

Baseline: native multi-region LGV, driven from Python by emitting the region
list into `displayedRegions`. No custom code.

If per-region width normalisation is required (ITV's
`normalize_interval_width`), implement `ITVExonSliceView` as an LGV subclass
applying a per-region scaling factor in the coordinate transform. Assessment is
that this is tractable — it is a scaling factor per block — but it touches
block-layout machinery and is the most uncertain piece of Path A. See spike
question 2.

## A4 — Isoform-paired annotation

No native JBrowse concept. Implement as a custom **display type** that reads
`linked_annotation_id` from the payload and draws the transcript model together
with its aligned reads in one track.

Defer until A1–A3 are working; it is the most bespoke piece and benefits from
familiarity with the display API.

## A5 — Static export

The piece JBrowse provides least help with. Two products:

### A5a — Self-contained single report

Vite build producing one HTML file: app bundle + ITV plugin + payload injected
as a global, gzip+base64, decoded in-page via `DecompressionStream`, served to
the adapter from memory.

No network, no sidecar, no `file://` CORS problem. Size = fixed bundle overhead
+ compressed payload.

### A5b — Multi-report bundle with shared runtime

For "one report per gene" over a gene list, inlining the bundle N times is
wasteful. Instead:

```
report_bundle/
    index.html            — landing page: searchable gene list + stats
    assets/itv.<hash>.js  — shared runtime, loaded once, browser-cached
    data/GENE.itv.bin     — per-gene payload
```

Navigation: sidebar gene list, client-side routing, payload fetched on
selection.

**Constraint to design around:** `fetch()` of `data/*.bin` is blocked under
`file://` in Chrome. Options, in order of preference:
1. Ship with a one-line launcher (`python -m http.server`) and document it.
2. Inline all payloads into `index.html` when the total is small enough.
3. Single-file-per-gene (A5a) for the handful of genes that matter, plus a
   lightweight static index page linking to them — this works under `file://`
   because plain navigation is not subject to CORS.

Option 3 is probably the best default for sharing with collaborators: it
degrades gracefully, and the per-file bundle overhead only hurts if the gene
list is long.

## A6 — Batch figure generation

`@jbrowse/img` / `jb2export` renders SVG headlessly via node + React SSR, with
no browser required. Python wraps it as a subprocess call, giving back the
scripted figure generation ITV has today.

Open question: whether `jb2export` can load a custom plugin. If not, write a
small node script using the same SSR approach against the ITV plugin. Either
way the mechanism exists.

## A7 — Overview view

Confirmed nice-to-have. The intended form is a **dense genome- or
chromosome-scale track** showing where there is anything worth exploring —
closer to a coverage track than to a table.

Cheapest implementation: emit the T5 stats bundle as a quantitative track
(reads/base density, max depth) served by a custom adapter, displayed in a
zoomed-out LGV. Clicking or navigating drills into the detailed view. No custom
view type needed — it is just another track.

Only build `ITVOverviewView` as a custom view type if the quantitative-track
form proves insufficient.

## Sequencing

| Step | Depends on | Notes |
|---|---|---|
| A1 adapter | trunk T4 | validates the hinge |
| A2 splits | A1 | the main daily-use feature |
| A3 exon slices | A1 | native first, subclass only if needed |
| A5a single-file | A1 | unblocks sharing early |
| A6 batch figures | A1 | restores existing capability |
| A4 isoform display | A1–A3 | most bespoke |
| A5b multi-report | A5a | |
| A7 overview | trunk T5 | nice-to-have, last |

## Ongoing costs to accept

- React + MobX-State-Tree competence required.
- Plugin maintenance across JBrowse major versions (v1→v4 so far).
- Multi-MB fixed bundle overhead per self-contained file.
- One version of a plugin per page.
