# ITV interactive rewrite — documentation index

Working notes for re-architecting ITV from a static SVG generator into an
interactive, efficient, shareable viewer.

**Status:** exploration complete, plan drafted, no code written yet.
**Branch:** `refactor/interactive-rewrite`, based on **`typing_and_doc`** (not
`main` — code line references only resolve against `typing_and_doc`).
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

## Execution units

Work in one unit per session. Each unit lists the only docs that need loading —
resist loading all nine, the files are split precisely so you don't have to.

| # | Unit | Load | Model | Mode |
|---|---|---|---|---|
| 1 | Spike Q3 — bundle size | `03`, `07` | Sonnet | normal |
| 2 | Spike Q1 — Export SVG quality | `03`, `07` | Sonnet to set up; **verdict is human** | normal |
| 3 | T0 — baseline fixtures | `00`, `04` | Sonnet | normal |
| 4 | **T1 — payload schema** | `00`, `04`; skim `03`, `06` | **Opus** | **plan mode** |
| 5 | **T2 — extract compute layer** | `04`, `CLAUDE.md` | **Opus throughout** | plan → execute |
| 6 | T3 — SVG renderer on payload | `04`, T1 schema doc | Sonnet | normal |
| 7 | T4 — serialisation | `04` | Sonnet | normal |
| 8 | T5 — overview stats bundle | `04` | Sonnet | normal |
| 9 | Spike Q4/Q6 | `03`, `05`, `07` | Opus | normal |

**Units 4 and 5 are where to spend effort.** T1 is load-bearing for both paths
and expensive to get wrong. T2 is hazardous because of the import-time kwargs
metaprogramming — see
[04 § Hazard](04-plan-common-trunk.md#hazard-the-import-time-kwargs-metaprogramming).
Everything else is recoverable, and `/fast` is worth enabling for units 3, 6
and 7.

Two things that are **not** delegable:
- **Q1's verdict.** Whether an exported SVG is clean enough for a figure is a
  human judgment on a rendered artefact opened in Inkscape/Illustrator.
- **Q5.** Ask an actual ITV user, don't self-assess the UX.

### Starting a unit

Each session: state the unit, point at the docs, state the constraint. E.g.

> Working on the ITV interactive rewrite, **unit 4 (T1 — payload schema)**.
> Read `docs/refactor/04-plan-common-trunk.md` and
> `docs/refactor/00-current-architecture.md` first; skim `03` and `06` for
> what each renderer will need from the payload. Design the schema so it can
> feed both a JBrowse custom adapter and a Canvas scene layer. Use plan mode
> and argue the tradeoffs before writing anything.

> Working on the ITV interactive rewrite, **unit 5 (T2)**. Read
> `docs/refactor/04-plan-common-trunk.md`, especially the hazard section on
> the import-time kwargs metaprogramming, plus `CLAUDE.md`. Move in small
> steps and re-import the package after each one.

### Finishing a unit

Append an entry to
[07-spike-and-decision-log.md](07-spike-and-decision-log.md) for any decision
made or reversed, and tick the box in [Decision state](#decision-state). The
log is append-only — it is what lets a session months from now reconstruct
*why*, not just *what*.

## Known gaps requiring work under Path A

| Gap | Planned approach | Spike Q |
|---|---|---|
| Per-exon width normalisation | LGV subclass scaling each region | Q2 |
| Self-contained single-file HTML | Vite bundle + inlined payload + in-memory adapter | Q3 |
| Stacked "sediment" coverage | Pre-summed layers as overlapping filled XY | Q6 |
| Isoform-paired annotation | Custom display type | — |
