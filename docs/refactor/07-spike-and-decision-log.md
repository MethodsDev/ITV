# Spike plan and decision log

## Go / no-go spike

Purpose: decide between [Path A](05-path-a-jbrowse2.md) and
[Path B](06-path-b-canvas.md) on evidence rather than preference.

**The spike runs in two stages.** Q1 and Q3 depend on nothing and answer the
existential question, so they run *before* the common trunk. The rest need a
real payload and run after.

### Stage 1 — before the trunk (~1 day, no ITV code)

| Q | Question | Why it can run now |
|---|---|---|
| [Q1](#q1--does-exportsvg-work-from-jbrowsereact-app2-and-is-the-output-clean) | Is Export SVG Illustrator-clean? | stock `react-app2`, any BAM |
| [Q3](#q3--what-is-the-real-gzipped-bundle-for-a-minimal-app--one-custom-plugin) | Real gzipped bundle size? | a Vite build and a trivial plugin |

**Q1 is existential.** Clean vector export is the single hard requirement and
the main objection raised against JBrowse. If it fails, Path A is dead and
months of trunk work would otherwise have been done under a false assumption.
Answer it before committing to anything.

Q3 is cheap to answer at the same time and determines the packaging strategy
([A5a vs A5b](05-path-a-jbrowse2.md#a5--static-export)), which is useful input
to the trunk's serialisation design (T4).

### Stage 2 — after the trunk (~1 week)

Requires [04-plan-common-trunk.md](04-plan-common-trunk.md) through at least
T4, so there is a real payload to test against.

| Q | Question |
|---|---|
| [Q2](#q2--can-multi-region-lgv-do-per-region-width-normalisation) | Per-region width normalisation |
| [Q4](#q4--can-a-custom-adapter-serve-in-memory-reads-at-realistic-depths) | Custom adapter performance on a real payload |
| [Q5](#q5--does-the-track-category-ux-actually-beat-tabs-for-comparing-cell-types) | Track groups vs tabs, with a real user |
| [Q6](#q6--can-stacked-sediment-coverage-be-reproduced-with-built-in-renderers) | Stacked sediment coverage |
| [Q7](#q7--does-everything-behave-in-firefox) | Firefox parity |

---

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

**Prior assessment:** likely a pass — it is a scaling factor per block.

**The bar is lower than it first appeared.** ITV's current exon-slice mode is
experimental: it works, but scales elements badly, because it applies
`<svg viewBox>` / `<use href>` tricks over an already-rendered document
instead of doing a real coordinate transform. So this question is *not* "can
we reproduce current output" — a genuine per-region transform under either
path should be an improvement on what exists today. Treat Q2 as
scope-affecting rather than blocking.

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

---

### Q6 — Can stacked "sediment" coverage be reproduced with built-in renderers?

**Why it matters:** ITV splits coverage into cumulative stacked layers by tag,
binned start/end position, read-end peak, or strand. JBrowse
multi-quantitative tracks have **no stacked/cumulative mode** — only
multi-row and overlapping.

**Test:** emit pre-summed layers (which `_add_multi_coverage` already produces
as `cumulative_coverage[ix]`) as subtracks of a multi-quantitative track.
Render as overlapping filled XY plot with opaque fills. Try to force draw
order tallest-first.

**Pass:** z-order and fill opacity are controllable, and the result is visually
equivalent to current ITV sediment bands — including after Export SVG.

**Fail:** subtrack z-order in overlapping mode is not controllable, so layers
occlude each other wrongly.

**On failure:** custom display type with its own renderer. Well-supported
plugin point, precedent in `cancerit/proportionalmultibw`. Adds scope, does
not change path.

---

### Q7 — Does everything behave in Firefox?

**Why it matters:** Firefox is the primary target browser.

**Test:** run Q1 (Export SVG), the multi-region view, and the Q4 adapter in
Firefox as well as Chrome.

**Pass:** no behavioural difference in export output, rendering, or adapter
performance.

**Fail:** Firefox-specific rendering or export defects.

**Note:** JBrowse's GPU shaders sit above a Canvas2D baseline, so absent
WebGPU the baseline is used — and the baseline is also what drives SVG export.
No Firefox-specific blocker was found during research, so this is a
confirmation check rather than an expected risk.

---

### Decision rule

**After stage 1:**
- **Q1 fails** → stop. Path B. Clean vector export is non-negotiable and is
  the main reason people choose ITV. Do the trunk with Path B in mind.
- **Q1 passes, Q3 large** → continue to the trunk on Path A, but plan for the
  shared-runtime bundle
  ([A5b](05-path-a-jbrowse2.md#a5b--multi-report-bundle-with-shared-runtime))
  rather than single-file-per-gene, and feed that back into T4.
- **Both pass** → continue to the trunk on Path A.

**After stage 2:**
- **Q4 and Q7 pass** → commit to Path A. Q2, Q5 and Q6 failures are absorbed
  by documented fallbacks; they change scope, not direction.
- **Q4 fails hard** → Path B, with the common trunk already banked.
- **Q7 fails** → Path B. Firefox is the primary target.

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

### 2026-09-23 — Tabs implemented as track groups + show/hide buttons
An LGV holds one ordered list of mixed-type tracks (alignments, quantitative,
annotation), freely interleaved and reorderable, with programmatic
`showTrack`/`hideTrack`. So each classification split becomes a *group* of
tracks — annotation + coverage + reads — and a button row toggles group
visibility.

**Consequence:** ITV's current per-tab splitting of both reads and BED
annotation entries is reused unchanged; it emits more tracks instead of more
SVG documents. Supersedes the earlier "categorised track selector, not tabs"
recommendation — the track structure is the same either way, and the button
row is closer to the UX users already have. One instance, one worker pool.
See [A2](05-path-a-jbrowse2.md#a2--classification-splits).

### 2026-09-23 — Firefox is the primary target browser
Chrome secondary but must work.

**Consequence:** Firefox joins the screenshot/CI baseline from the first
milestone; added spike Q7. No Firefox-specific blocker found — JBrowse's GPU
shaders sit above a Canvas2D baseline, and that baseline is also what drives
SVG export.

### 2026-09-23 — Correction: the `file://` fetch restriction is not Chrome-specific
Earlier notes framed sidecar-file `fetch()` blocking as a Chrome behaviour.
It applies to **Firefox as well**, which has treated each `file://` document
as a unique origin since version 68 (`privacy.file_unique_origin`,
CVE-2019-11730).

**Consequence:** choosing Firefox as the primary target does *not* relax the
multi-report packaging constraint. The index-page-plus-single-file-reports
approach remains the offline-safe option, because plain link navigation
between local HTML files is unaffected.

### 2026-09-23 — Sediment coverage: no native JBrowse mode, workaround identified
JBrowse multi-quantitative tracks support five plot types under multi-row or
overlapping layouts, but **no stacked/cumulative area mode**.

**Consequence:** plan to exploit the fact that
`bamtrack.py::_add_multi_coverage` already emits *pre-summed* layers — render
them as overlapping filled XY subtracks, opaque, tallest first. Added spike Q6
to confirm z-order is controllable. Custom display type is the documented
fallback. See [A2b](05-path-a-jbrowse2.md#a2b--stacked-sediment-coverage-layers).

### 2026-09-23 — Spike split into two stages; Q1 and Q3 move ahead of the trunk
Q1 (Export SVG quality) and Q3 (real bundle size) depend on no ITV code —
a stock `react-app2` build answers both in about a day.

Q1 is existential: clean vector export is the one hard requirement, and if it
fails, Path A is dead. Running it after the trunk would risk months of work
done under a false assumption. Q3 is free to run alongside and its answer
feeds the trunk's serialisation design (T4), since packaging strategy
determines the payload size budget.

**Consequence:** stage 1 (Q1, Q3) runs before
[04-plan-common-trunk.md](04-plan-common-trunk.md); stage 2 (Q2, Q4, Q5, Q6,
Q7) runs after T4. Decision rule split accordingly.

### 2026-09-23 — Exon-slice mode is experimental; current output is not the quality bar
The `plot_exons` slice mode works but produces badly scaled elements, hurting
quality and readability. Cause: it applies `<svg viewBox>` / `<use href>`
tricks over an already-rendered document rather than performing a real
coordinate transform.

**Consequence:** keep the capability through the rewrite, where a proper
renderer should fix it, but do not treat current output as the target to
match. Spike Q2 downgraded from blocking to scope-affecting — any genuine
per-region transform is expected to beat what exists.

<!-- Next entry goes here -->
