# Plan — common trunk (path-independent)

**Do this first.** Every task here is required under Path A *and* Path B. None
of it is wasted regardless of which path is chosen, which makes it the
risk-reducing move: the path decision can be deferred until the spike results
are in.

Goal: separate ITV's Python compute layer from SVG string generation, and emit
a compact, versioned data payload.

## T0 — Freeze a behavioural baseline

Before touching anything, capture what current ITV produces so the rewrite can
be checked against it.

- Pick ~10 representative cases: single gene, multi-isoform gene, exon-slice
  view, classification tab split, coverage-only, read-level zoom.
- Render each with current ITV to SVG; commit to `tests/baseline/`.
- Record the Python call that produced each one.

These become regression fixtures. Without them the rewrite has no ground truth.

## T1 — Define the payload schema

The single most load-bearing artefact. **Version it from commit one.**

Sketch (to be refined):

```
ITVPayload
  schema_version: int
  assembly:     { name, refNames[], (optional inlined reference seq slice) }
  regions[]:    { refName, start, end, label, display_weight? }   # exon slices
  features[]:   per annotation/BED/GTF entry
  readsets[]:   one per classification split
      key:        "celltype=Tcell" | "sqanti=FSM" | ...
      label, color
      reads:      columnar arrays
          start[], end[], flags[], cigar (varint-packed), mapq[]
          tag_ids[] -> tag dictionary
      coverage:   binned array (precomputed per readset)
      coverage_layers[]?      # stacked "sediment" mode, see note below
          key, color, values[]    # ALREADY cumulative-summed
      annotation_ids[]        # BED/GTF entries belonging to this split
      linked_annotation_id?   # isoform-paired annotation, see note below
  stats:        per-feature summary for the overview view
```

Design constraints:
- Columnar, not row-of-objects — 40–80 bytes/read target.
- Must round-trip to *both* a JBrowse custom adapter and a Canvas scene layer.
- Coverage is precomputed per readset, not derived client-side, because ITV
  already computes it and the splits are known ahead of time.
- `coverage_layers` carries ITV's stacked "sediment" modes
  (`add_tagged_coverage`, `add_binned_coverage`, `add_peak_coverage`,
  `add_stranded_coverage`). **Keep the values cumulative-summed**, exactly as
  `_add_multi_coverage` already produces them — Path A's planned renderer
  depends on that, since JBrowse has no stacked mode and the workaround is to
  hand it pre-summed series. Record the layer order explicitly; draw order
  matters.
- `annotation_ids` per readset preserves the fact that ITV splits BED/GTF
  entries alongside reads, so a backend can place each split's annotation
  next to its own reads.
- `linked_annotation_id` carries ITV's isoform-paired-annotation relationship
  (capability 5 in [00-current-architecture.md](00-current-architecture.md));
  it must be in the payload even if neither renderer consumes it on day one.

Deliverable: `docs/refactor/payload-schema.md` + a Python dataclass module +
JSON Schema or equivalent for the JS side.

## T2 — Extract the compute layer

Pull the data-producing logic out of `convenience.py` (2517 lines) into modules
that return payload objects instead of `Document`s.

Keep, essentially unchanged:
- `VirtualBAM` and the read filtering/classification machinery
- `annotation_matching.py`, `cellbarcode.py`, `classification.py`
- GTF/BED indexing and feature lookup in `Configuration`
- exon merging and interval computation from `_plot_exons_helper_get_info`

Change:
- `plot_*` methods split into `build_payload_*` (returns data) and a thin
  rendering wrapper (returns a `Document`, for backwards compatibility).

**Backwards compatibility matters here.** ITV is in use. The existing `plot_*`
API should keep working throughout the trunk work, implemented on top of the
new payload path once T3 lands.

## T3 — Reimplement the current SVG renderer on top of the payload

Prove the schema is sufficient by rewiring the *existing* SVG output to consume
the payload rather than the track objects.

This is the key validation step. If the current SVG output can be reproduced
from the payload, the payload is complete. If it cannot, the schema is wrong
and it is far cheaper to learn that now than after a JS renderer is built on it.

Success criterion: T0 baseline fixtures reproduce (modulo cosmetic diffs).

## T4 — Serialisation

- Binary/columnar encoder (Arrow IPC is the obvious candidate — good Python and
  JS support, columnar by nature, avoids inventing a format).
- Gzip + base64 wrapper for inlining, with `DecompressionStream` on the JS side.
- Measure: bytes/read achieved, and total size for a realistic gene with 10–20
  classification splits.

## T5 — Overview stats bundle

Small precomputed per-feature table: max depth, span, n_reads, reads/base
density, per-classification breakdown.

Confirmed as a **nice-to-have**, so keep it cheap. The intended presentation is
a genome- or chromosome-scale density track (like a coverage track), not a
spreadsheet — see [05-path-a-jbrowse2.md](05-path-a-jbrowse2.md) for how that
lands in JBrowse.

Emitting it is trunk work because it is pure Python aggregation; rendering it
is path-specific.

## Exit criteria

The trunk is done when:

1. `build_payload_*` functions exist for every current `plot_*` entry point.
2. The existing SVG renderer runs entirely off payloads.
3. T0 baseline fixtures reproduce.
4. A realistic multi-split gene serialises to a measured size, and that number
   is written down in [07-spike-and-decision-log.md](07-spike-and-decision-log.md).
5. The schema is versioned and documented.

At that point, run the spike ([07](07-spike-and-decision-log.md)) and choose a
path.
