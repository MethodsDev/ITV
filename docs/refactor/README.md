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
- [x] Path A (JBrowse 2) is the preferred first attempt.
- [ ] Spike run — see [07-spike-and-decision-log.md](07-spike-and-decision-log.md).
- [ ] Path chosen.
