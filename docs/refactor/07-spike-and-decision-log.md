# Spike plan and decision log

## Go / no-go spike

**Timebox: ~1 week.** Run after the common trunk
([04](04-plan-common-trunk.md)) is at least through T4, so there is a real
payload to test against.

Purpose: decide between [Path A](05-path-a-jbrowse2.md) and
[Path B](06-path-b-canvas.md) on evidence rather than preference.

### Q1 — Does `exportSvg` work from `@jbrowse/react-app2`, and is the output clean?

**Why it matters:** publication vector export is a hard requirement and the
main objection raised against JBrowse.

**Test:** embed `react-app2`, load an alignments track, trigger `exportSvg`,
open the result in Inkscape or Illustrator.

**Pass:** SVG opens, text is real text (not paths), elements are grouped
sensibly, colours and theme survive, file is editable.

**Fail:** rasterised content embedded, unusable grouping, or the action is not
reachable from the embedded component.

---

### Q2 — Can multi-region LGV do per-region width normalisation?

**Why it matters:** ITV's `normalize_interval_width` gives each exon equal
display width. If unreachable, either the feature is lost or a custom view
subclass is needed.

**Test:** open a multi-region view; attempt to give regions independent
bases-per-pixel. Read the LGV block-layout source to judge how contained a
subclass would be.

**Pass:** achievable natively, *or* the subclass looks like a contained change
to the coordinate transform.

**Fail:** width scaling is entangled through block layout, rendering and
export such that a subclass would fork significant machinery.

**Prior assessment:** likely a pass — it is a scaling factor per block — but
this is the most uncertain piece of Path A.

---

### Q3 — What is the real gzipped bundle for a minimal app + one custom plugin?

**Why it matters:** determines whether single-file self-contained reports are
practical, especially for many small per-gene reports.

**Test:** Vite production build of `react-app2` + a trivial plugin. Measure
gzipped JS.

**Pass:** small enough that bundle overhead is acceptable next to payload size
for the intended report sizes.

**Fail:** so large that per-gene single-file reports are dominated by runtime
overhead — in which case fall back to the shared-runtime bundle
([A5b](05-path-a-jbrowse2.md#a5b--multi-report-bundle-with-shared-runtime))
or reconsider Path B.

**Known:** npm unpacked sizes are ~30 MB (`react-app2`) and ~26 MB (LGV2) at
v4.3.0, but those include sourcemaps and multiple build targets. The
production number is what matters and is not yet measured.

---

### Q4 — Can a custom adapter serve in-memory reads at realistic depths?

**Why it matters:** this is the hinge for the whole self-contained design.

**Test:** implement `ITVReadsAdapter` over a decoded payload for a realistic
gene at production depth with 10–20 classification splits. Measure time to
first render, pan/zoom responsiveness, memory.

**Pass:** interactive pan/zoom, no main-thread stalls beyond a few hundred ms
on load.

**Fail:** the RxJS observable path or feature-object overhead makes per-read
serving too slow, and the columnar payload has to be re-boxed per feature at
prohibitive cost.

---

### Q5 — Does the track-category UX actually beat tabs for comparing cell types?

**Why it matters:** determines whether ITV's tab metaphor needs recreating at
all. **Ask a user, not yourself.**

**Test:** put the A2 prototype in front of someone who uses ITV for cell-type
comparison. Watch what they do.

**Pass:** categories are acceptable or preferred — comparing splits side by
side beats flipping between them.

**Fail:** users need the tab metaphor, in which case build own chrome above
`react-app2`
([A2 alternative](05-path-a-jbrowse2.md#a2--classification-splits)) rather than
multiple JBrowse instances.

---

### Decision rule

- **Q1, Q3, Q4 pass** → take Path A. Q2 and Q5 failures are absorbed by the
  documented fallbacks; they change scope, not direction.
- **Q3 or Q4 fail hard** → Path B, with the common trunk already banked.
- **Q1 fails** → Path B. Vector export is non-negotiable and is the main reason
  people use ITV.

---

## Decision log

Append-only. Record what was decided, when, and why.

### 2026-09-23 — Scope confirmed: per-feature, not whole-transcriptome
Primary use case is viewing pre-defined genes/transcripts with their isoforms;
most inter-genic space is empty. Arbitrary whole-transcriptome pan/zoom is not
required.

**Consequence:** in-browser BAM readers over HTTP range requests
([option 2.3](02-options-explored.md#23-in-browser-bam-readers-over-http-range-requests))
are out of scope. Inlined payloads are sufficient, which keeps
self-containment achievable.

### 2026-09-23 — Static vector export is a hard requirement
Figures for papers and posters. Any architecture that cannot produce clean
editable SVG is disqualified.

**Consequence:** drives spike Q1; made JBrowse viable once Export SVG and the
shared Canvas2D/SVG path were confirmed.

### 2026-09-23 — JBrowse 2 selected as the first path to try
Interactivity and usability match the target. Export SVG, multi-region LGV,
headless `jb2export`, and programmatic `createViewState`/`navToLocString`
control all confirmed to exist.

**Consequence:** [Path A](05-path-a-jbrowse2.md) is primary,
[Path B](06-path-b-canvas.md) is the fallback.

### 2026-09-23 — Common trunk before path choice
The Python compute/presentation separation is required under both paths.

**Consequence:** [04-plan-common-trunk.md](04-plan-common-trunk.md) runs first;
the path decision is deferred until after it, making the choice reversible and
evidence-based.

### 2026-09-23 — Overview view downgraded to nice-to-have
Intended form clarified as a dense genome/chromosome-scale track, not a
spreadsheet.

**Consequence:** implementable as a quantitative track over the T5 stats
bundle; no custom view type needed initially. Scheduled last.

### 2026-09-23 — Multiple JBrowse instances rejected as the default for tabs
Considered as a way to reuse ITV's existing JS tabs. Each instance carries its
own worker pool, session state and render loop.

**Consequence:** prefer one instance with many categorised tracks; own tab
chrome above `react-app2` is the compromise if users need tabs. Multiple
instances stay as a last resort.

<!-- Next entry goes here -->
