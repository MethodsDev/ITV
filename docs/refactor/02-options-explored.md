# Options explored

Three independent axes. Choices on one do not force choices on another.

- [Axis 1: rendering substrate](#axis-1--rendering-substrate)
- [Axis 2: data delivery and self-containment](#axis-2--data-delivery-and-self-containment)
- [Axis 3: host / runtime](#axis-3--host--runtime)

---

## Axis 1 — Rendering substrate

### A. Keep SVG, add virtualization

Viewport culling, lazy tab rendering, LOD, DOM recycling.

**Pros** — smallest diff; `svg.py` and all track classes survive; vector export
stays free; CSS/DOM tooltips and hover for free; no new language.

**Cons** — ceiling ~10–20k live nodes; 100–300 ms/frame on style recalc; smooth
60 fps panning unreachable.

**Verdict:** good 2-week stopgap, dead end as an architecture. See
[01-problem-analysis.md](01-problem-analysis.md#interim-mitigation-optional-12-weeks).

### B. Canvas 2D + JS scene layer

Python emits a compact scene payload; JS owns layout, culling, drawing,
hit-testing.

**Pros** — 50–100k shapes/frame comfortable with `fillStyle` batching and
`Path2D` reuse, well past what is legible in a viewport; `fillText` makes
base-level zoom, read labels and mismatch letters trivial; `OffscreenCanvas` +
workers for off-main-thread tiles; hit-testing via JS interval tree or offscreen
picking buffer; debuggable; a second serializer from the same scene graph emits
SVG for figures.

**Cons** — you hand-write the drawing layer (mirrors `svg.py` closely);
no DOM accessibility; manual HiDPI handling.

**Verdict:** the fallback path. See [06-path-b-canvas.md](06-path-b-canvas.md).

### C. WebGL (PixiJS / regl / deck.gl)

**Pros** — 10^6+ instanced quads at 60 fps; free smooth zoom via transform
matrices.

**Cons** — text requires SDF/bitmap font atlases, painful in a text-heavy
alignment view; GPU/driver variance; harder headless CI screenshots; much more
code.

**Verdict:** not a starting point. The binding constraint is data transfer, not
fill rate. Reserve for specific layers (deep coverage heatmaps, large dotplots)
behind the same scene interface.

### D. Adopt an existing browser

#### D1. igv.js
**Pros** — canvas-based, BAM via HTTP range requests, battle-tested, familiar
interaction model.
**Cons** — genome-centric track model. The transcriptome-coordinate projection,
exon slicing, classification tab splits and `VirtualBAM` read sets are exactly
what it does not express. Custom track plugins constrained. Self-contained
single-file export is not its model.
**Verdict:** rejected — the differentiators do not survive.

#### D2. JBrowse 2 — **selected as Path A**
Plugin architecture, canvas + GPU rendering, embeddable React components.
Full assessment in [03-jbrowse2-assessment.md](03-jbrowse2-assessment.md).
**Verdict:** try first.

#### D3. Gosling.js / HiGlass
**Pros** — declarative grammar over PixiJS/WebGL; built for multiscale
genomics; linked views; nonlinear and circular layouts are first-class, which
is unusually close to the exon-slice need; Python API exists.
**Cons** — the declarative grammar fights bespoke visuals; per-read CIGAR detail
with custom classification colouring is not its sweet spot; scalable paths want
pre-tiled formats (beddb/multivec), reintroducing a build step and hurting
self-containment.
**Verdict:** second choice if JBrowse fails, ahead of Path B, but the
per-read-detail mismatch is significant.

### E. Precomputed tile pyramid (the "Google Maps" answer)

Python pre-renders PNG or vector tiles at N zoom levels; browser blits.

**Pros** — rendering cost moves offline; guaranteed smooth pan; trivial client.

**Cons** — tile count explodes as zoom × tabs × tracks. Critically, it buys
*smoothness but not interactivity* — hovering a read, filtering by tag,
re-sorting or toggling a track all require regeneration. Google Maps gets away
with this because the map is read-only; ITV users interrogate reads.

**Verdict:** rejected as a primary architecture. Possible hybrid: raster tiles
for a zoomed-out overview band, live vector once the window is small enough.
Overlaps with the "nice to have" overview view — see
[05-path-a-jbrowse2.md](05-path-a-jbrowse2.md).

---

## Axis 2 — Data delivery and self-containment

### 2.1 Inline compact payload in the HTML
Columnar binary or Arrow IPC → gzip → base64, decompressed in-page with
`DecompressionStream`.

**Pros** — genuinely single-file; no CORS/`file://` problems; no server.
**Cons** — practical envelope ~20–30 MB compressed before browsers struggle.
At 40–80 bytes/read that is still hundreds of thousands of reads.
**Verdict:** primary export tier.

### 2.2 HTML + sidecar data directory (or `.zip` read in-page)
**Pros** — no size ceiling; shared runtime across many reports.
**Cons** — not single-file; **`fetch()` of sidecar files is blocked under
`file://` in both Firefox and Chrome** (Firefox since v68,
`privacy.file_unique_origin`, CVE-2019-11730), so recipients need a local web
server or hosting. This is the main practical gotcha, and it is not avoidable
by choosing a browser. Plain link navigation between local HTML files is not
affected.
**Verdict:** secondary tier for large sessions and for multi-report bundles.

### 2.3 In-browser BAM readers over HTTP range requests
`@gmod/bam`, `@gmod/tabix`, `@gmod/bbi`, or htslib/noodles via WASM.
**Pros** — unlimited pan/zoom, no precomputation.
**Cons** — requires a range-serving host; the opposite of self-contained.
**Verdict:** out of scope given the confirmed per-gene use case. Note
`VirtualBAM` already produces region-restricted read sets, which is easier to
serialise than shipping htslib-WASM.

### 2.4 Live kernel mode (anywidget / ipywidgets comm)
**Pros** — unlimited input size, no precomputation, full Python available.
**Cons** — dies outside a running kernel, so it cannot be the sharing path.
**Verdict:** exploration mode only; complements but does not replace 2.1.

**Chosen combination:** 2.1 for single-gene reports, 2.2 for multi-report
bundles and oversized sessions, 2.4 later for notebook exploration.

---

## Axis 3 — Host / runtime

| Option | Pros | Cons |
|---|---|---|
| Static self-contained HTML | shareable by email, archivable, no deps | size ceiling; no live Python |
| Static HTML + sidecar | no size ceiling; shared runtime | `file://` CORS; needs a server |
| Jupyter widget (anywidget) | live Python, unlimited data | needs running kernel; no sharing |
| Headless CLI render | scriptable batch figures | no interactivity (by design) |

All four are wanted eventually. The common trunk makes all four reachable from
one payload format.

---

## Cross-cutting conclusion

Axis 1 gets the most attention and matters least. Axis 2 determines whether the
"shareable" goal is met. And the prerequisite for everything is the data/view
separation described in [01-problem-analysis.md](01-problem-analysis.md), which
is why [04-plan-common-trunk.md](04-plan-common-trunk.md) comes before any path
choice.
