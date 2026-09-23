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

**Chosen idiom: one track per split, in one ordered track list, with a button
row driving `showTrack`/`hideTrack` to give tab-like behaviour.**

An LGV holds a single ordered list of mixed-type tracks. So for each
classification split emit a *group* of tracks:

```
[split: celltype=Tcell]
    ITVAnnotationAdapter  ->  isoforms / BED entries relevant to this split
    ITVCoverageAdapter    ->  coverage for this split
    ITVReadsAdapter       ->  reads for this split
[split: celltype=Bcell]
    ...
```

Order within a group is controlled, so each split's annotation sits directly
above or below its own reads. A button row above the view calls
`showTrack`/`hideTrack` over group members — one JBrowse instance, one worker
pool, and the perceived result is ITV's current tabs.

**This reuses ITV's existing work directly.** The current code already splits
both reads *and* BED annotation entries per tab; that Python logic feeds the
payload unchanged and simply emits more tracks instead of more SVG documents.

Alternatives considered:
- **Categorised track selector** instead of a button row. Same track
  structure, different chrome; free either way. Worth showing users alongside
  the buttons (spike Q5) — comparing splits side by side may beat flipping
  between them.
- `jexl` `filter` + `color` callbacks over a single track. Simpler config, but
  it pushes filtering to render time and does not give per-split coverage as
  naturally. Good for quick ad-hoc colouring, not for the primary split axis.
- **Multiple JBrowse instances inside existing ITV tabs.** Rejected as the
  default — each instance carries its own worker pool, session state and
  render loop, reimporting the heaviness the rewrite exists to remove. Keep as
  a last resort if track counts become unmanageable.

## A2b — Stacked "sediment" coverage layers

ITV supports coverage split into cumulative stacked layers by BAM tag, by
binned alignment start/end position, by read-end peak, or by strand
(`add_tagged_coverage`, `add_binned_coverage`, `add_peak_coverage`,
`add_stranded_coverage`, all funnelling into `_add_multi_coverage`).

JBrowse multi-quantitative tracks have no stacked/cumulative mode — only
multi-row and overlapping, across five plot types.

**Planned approach — no custom rendering needed.** `_add_multi_coverage`
already emits each layer *pre-summed* (`cumulative_coverage[ix]`). Carry those
pre-summed series into the payload as ordinary subtracks and render as an
**overlapping filled XY plot** with opaque fills, tallest drawn first.
Visually identical to the current sediment bands, and it inherits Export SVG
for free.

Requires control over subtrack draw order and fill opacity. Subtrack colours
and row ordering are configurable (`facet` `domain`, per-track colour, manual
reorder); *z-order in overlapping mode* is the specific unknown — spike Q6.

Fallback if z-order is not controllable: a custom display type with its own
renderer. Precedent exists in the third-party `cancerit/proportionalmultibw`
plugin, which does stacked proportional bigwigs. A well-supported plugin
point, so the fallback is real work but not risky work.

## A3 — Exon slices

Baseline: native multi-region LGV, driven from Python by emitting the region
list into `displayedRegions`. No custom code.

If per-region width normalisation is required (ITV's
`normalize_interval_width`), implement `ITVExonSliceView` as an LGV subclass
applying a per-region scaling factor in the coordinate transform. Tractable —
it is a scaling factor per block — though it touches block-layout machinery.
See spike Q2.

**ITV's existing slice mode is experimental and scales elements badly** (it
manipulates an already-rendered document via `<svg viewBox>`/`<use href>`
rather than transforming coordinates). Doing this properly is an improvement,
not a reproduction, so Q2 is scope-affecting rather than blocking.

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
`file://` — in **Firefox and Chrome alike**. Firefox has treated each
`file://` document as a unique origin since version 68
(`privacy.file_unique_origin`, CVE-2019-11730). Targeting Firefox first does
not avoid this. Options, in order of preference:
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
| A2b sediment coverage | A2 | try pre-summed overlapping XY before custom renderer |
| A4 isoform display | A1–A3 | most bespoke |
| A5b multi-report | A5a | |
| A7 overview | trunk T5 | nice-to-have, last |

## Browser support

**Firefox is the primary target**; Chrome is secondary but must work.

- Put Firefox in the screenshot/CI baseline from A1 onwards rather than
  retrofitting it later.
- JBrowse layers WebGPU/WebGL2 shaders above a Canvas2D baseline, so Firefox
  degrades to the baseline rather than failing. That baseline is also the path
  that feeds Export SVG, so the export-critical code is the most portable code.
- The `file://` sidecar-fetch restriction applies to Firefox too — see A5b.

## Ongoing costs to accept

- React + MobX-State-Tree competence required.
- Plugin maintenance across JBrowse major versions (v1→v4 so far).
- Multi-MB fixed bundle overhead per self-contained file.
- One version of a plugin per page.
