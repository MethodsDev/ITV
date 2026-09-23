# Current ITV architecture (as of 2026-09)

Baseline description of the package being replaced.

**All line references below are relative to the `typing_and_doc` branch**,
which `refactor/interactive-rewrite` is based on. They do not resolve against
`main` — the two diverge by ~2,100 lines, and `_plot_exons_slices` does not
exist on `main` at all. Re-verify these references after any merge of
`typing_and_doc` into `main`.

Measurements taken from the `examples/` outputs committed on that branch.

## Pipeline

```
pysam / VirtualBAM
      ↓
Python track objects  (BAMTrack, BEDTrack, GraphTrack, ...)
      ↓
layout()  — row packing, height calculation   [depends on zoom]
      ↓
render()  — generators yielding SVG markup strings
      ↓
one big SVG blob
      ↓
embedded in HTML / ipywidgets
```

## Key modules

| File | Role |
|---|---|
| `src/integrative_transcriptomics_viewer/svg.py` | `SVG` backend: generator functions yielding SVG markup strings. `Renderer` does coordinate offsetting and clip-group nesting. |
| `genomeview.py` | `Document` → `ViewRow` → `GenomeView` → `Track` hierarchy. `Scale` holds the single linear genomic→pixel projection. |
| `bamtrack.py` | Read tracks. Row packing in `layout_interval()` (L68–151). `VirtualBAM` (L510) is an in-memory read set with `fetch`/`pileup`/`count`. |
| `convenience.py` | 2517 lines. `Configuration` class — the whole user-facing API: annotation indexing, feature lookup, tab organisation, plotting entry points. |
| `convenience.py::_plot_exons_slices` (L995–1088) | Nonlinear transcriptome projection, implemented with `<svg viewBox>` + `<use href>` tricks. |
| `templates.py` | Jinja2 tab assembly. `plot_sorted_support_as_tabs` (L31–98) materialises a full SVG string per (feature × classification) **before** display. |
| `export.py` | Save to SVG/PNG/PDF via resvg/cairosvg. Defines `SvgSplitter`, which `genomeview.py:118` invokes with `max_height=10000` to chop oversized output. |

## Measurements

Shipped examples, subsampled BAMs, single gene:

| File | Size | DOM nodes | `<svg>` elements |
|---|---|---|---|
| `examples/AIF1_classification_as_tabs.html` | 1.7 MB | ~10,000 | 24 |
| `examples/ACTB_exons_plot.html` | 1.5 MB | ~7,600 | 1 |

Node counts are `<rect>` + `<path>` + `<line>` occurrences.

Extrapolation: a per-cluster tab split at ~20 clusters and realistic depth
lands at 10^5–10^6 DOM nodes, which is past the point where browser style and
layout recalculation stays interactive.

## Capabilities that must not be lost

These are the reasons people choose ITV over an off-the-shelf browser. Any
replacement must preserve them.

1. **Publication-quality vector export.** SVG/PDF output suitable for
   Illustrator/Inkscape. This is the single strongest current feature.
2. **Transcriptome-coordinate / exon-slice projection.** Intron-collapsed views
   with optional per-exon width normalisation (`normalize_interval_width`).
3. **Classification-driven splitting.** Reads grouped by cell type, cell
   barcode, SQANTI-like read class, or arbitrary BAM tag — each group getting
   its own read track *and* its own coverage track. Note that ITV already
   splits **both reads and BED annotation entries** per group, which makes the
   logic directly reusable by a track-emitting backend.
4. **Stacked "sediment" coverage.** Coverage split into cumulative stacked
   layers by BAM tag (`add_tagged_coverage`), binned alignment start/end
   position (`add_binned_coverage`), nearest read-end peak
   (`add_peak_coverage`) or strand (`add_stranded_coverage`). All funnel into
   `_add_multi_coverage`, which accumulates into a running total and draws
   layers back-to-front so they read as sediment bands. Importantly the layers
   are emitted **already summed**, which matters for reimplementation.
5. **Isoform-paired annotation.** Pairing a specific transcript/isoform
   reference model with the reads that align to it.
6. **Scriptable batch figure generation** from Python, for a list of genes.
7. **`VirtualBAM`** — in-memory, region-restricted, filtered read sets. This is
   an asset for the rewrite, not a liability: it is already the compact
   pre-split payload that a new renderer would want.

## Where the cost actually is

- Layout is a function of zoom (`scale.topixels()` inside `layout_interval`),
  and it runs in Python. So any zoom change requires a Python round-trip.
- Every tab is rendered before any tab is shown.
- A read costs ~500–2000 bytes as SVG markup versus ~40–80 bytes in a columnar
  binary encoding (start, end, row, flags, varint CIGAR, tag indices) — a
  20–50× payload difference before compression.
