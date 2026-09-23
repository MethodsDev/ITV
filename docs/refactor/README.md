# ITV interactive rewrite — documentation index

Working notes for re-architecting ITV from a static SVG generator into an
interactive, efficient, shareable viewer.

**Status:** exploration complete, plan drafted, no code written yet.
**Branch:** `refactor/interactive-rewrite`
**Last updated:** 2026-09-23

## How to use these documents

Each file stands alone. Load only what you need for the task at hand.

| File | Read it when |
|---|---|
| [00-current-architecture.md](00-current-architecture.md) | You need to know how ITV works today and where the cost is. |
| [01-problem-analysis.md](01-problem-analysis.md) | You want the root-cause diagnosis (three separate problems, not "SVG is slow"). |
| [02-options-explored.md](02-options-explored.md) | You are revisiting a substrate / data-delivery / hosting decision. Full pros-cons matrix. |
| [03-jbrowse2-assessment.md](03-jbrowse2-assessment.md) | You are working on Path A. Capability findings, verified vs unverified. |
| [04-plan-common-trunk.md](04-plan-common-trunk.md) | **Start here for implementation.** Path-independent work, do this first. |
| [05-path-a-jbrowse2.md](05-path-a-jbrowse2.md) | Common trunk is done and the spike passed. |
| [06-path-b-canvas.md](06-path-b-canvas.md) | Common trunk is done and the spike failed. |
| [07-spike-and-decision-log.md](07-spike-and-decision-log.md) | You are running or recording the go/no-go spike. |

## One-paragraph summary

ITV today bakes geometry into SVG markup strings in Python, materialises every
tab eagerly, and ships presentation-encoded payloads. A single subsampled gene
produces a 1.7 MB HTML with ~10k DOM nodes. The fix is to separate the Python
compute layer from presentation and emit a compact data payload — work that is
required under every candidate architecture and should therefore be done first
(the "common trunk"). After that, two paths: adopt JBrowse 2 and write a plugin
(Path A, try first), or build a Canvas 2D scene layer in-house (Path B,
fallback). A one-week spike with five falsifiable questions decides between
them.

## Decision state

- [x] Scope confirmed: primary use case is **per pre-defined gene/transcript**, not arbitrary whole-transcriptome browsing. Whole-genome overview is a nice-to-have.
- [x] Static vector export for papers/posters is a **hard requirement**.
- [x] **Firefox is the primary target browser**; Chrome secondary but must work.
- [x] Path A (JBrowse 2) is the preferred first attempt.
- [x] Tabs to be implemented as track groups + `showTrack`/`hideTrack` buttons, one JBrowse instance.
- [ ] **Spike stage 1** (Q1 Export SVG quality, Q3 bundle size) — no ITV code needed, ~1 day. **Do this first.**
- [ ] Common trunk ([04](04-plan-common-trunk.md)) through T4.
- [ ] Spike stage 2 (Q2, Q4, Q5, Q6, Q7).
- [ ] Path chosen.

## Execution order

```
spike stage 1  →  common trunk T0-T5  →  spike stage 2  →  Path A or B
   (Q1, Q3)                               (Q2,Q4,Q5,Q6,Q7)
```

Q1 is existential — clean vector export is the one hard requirement. It needs
no ITV code, so answer it before investing in the trunk.

## Known gaps requiring work under Path A

| Gap | Planned approach | Spike Q |
|---|---|---|
| Per-exon width normalisation | LGV subclass scaling each region | Q2 |
| Self-contained single-file HTML | Vite bundle + inlined payload + in-memory adapter | Q3 |
| Stacked "sediment" coverage | Pre-summed layers as overlapping filled XY | Q6 |
| Isoform-paired annotation | Custom display type | — |
