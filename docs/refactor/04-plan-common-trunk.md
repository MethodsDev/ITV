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

**Source the cases from `examples/demo-ITV.ipynb`** (published as the
[online vignette](https://kinnex-documentation-external.readthedocs.io/en/latest/jupyter-notebooks/demo-ITV.html)).
It already exercises nearly every rendering capability, so there is no need to
invent cases. Inventory of what it covers:

| Notebook call | Capability exercised |
|---|---|
| `plot_interval(Interval(...), with_reads=False)` | explicit interval, coverage only |
| `plot_feature("APOE", padding_perc=-0.1)` | feature lookup, negative padding |
| `plot_feature("ACTB"/"LEPR"/"IL32")` | plain feature view, `view_width`, `coverage_height` |
| `plot_exons("ACTB"/"CAT"/...)` | **exon-slice projection** (experimental mode) |
| `plot_exons(..., tighter_track=True)` | compact track layout |
| `plot_exons(..., coverage_bin_size=100)` | `add_binned_coverage` → sediment layers |
| `plot_exons("PCMT1", coverage_peak_min_distance=20)` | `add_peak_coverage` → sediment layers |
| `plot_exons("STUB1", priming_orientation="3p")` | strand/priming-aware binning |
| `plot_exons("AIF1", coverage_tag="XI")` | `add_tagged_coverage` → sediment layers |
| `plot_exons("JUN"/"ACTB", with_reads=True)` | read-level rendering |
| `plot_splice_junctions("ACTB", as_widget=True)` | junction view, widget output |
| `plot_by_features_as_tab([...])` | feature tabs |
| `plot_by_classification_over_features(...)` | `BAMtagClassification` + `TaggedBAMAnnotationMatching`, and a custom classification subclass |
| `plot_by_classification_as_tabs("AIF1", ...)` | **classification tab split** — the heaviest case |
| `save(view, "ACTB_exons_plot.html")` | HTML/PNG export path |

Note this covers **all four sediment-coverage feeders** (binned, peak, tagged,
stranded), which matters for spike Q6.

Mechanics:
- Render each with current ITV to SVG; commit to `tests/baseline/`.
- Record the exact Python call alongside each output.
- There is **no test suite in this repo** — T0 means building the harness from
  scratch, not adding cases to an existing one.

**Comparison is approximate, with a human gate.** The new renderer is not
expected to be pixel- or node-identical, and for the exon-slice mode it should
be visibly *better* (see capability 2 in
[00-current-architecture.md](00-current-architecture.md)). So the fixtures are
a smoke test — do the right features appear, in the right places, at the right
scale, with the right groupings — not a diff gate. Flag divergences for human
review rather than failing on them.

**Data dependency:** the notebook needs `examples/data/` (~8 GB: GRCh38 refs,
GENCODE v39 GTF/BED, subsampled Kinnex bulk and single-cell BAMs). It is
untracked and sourced from GCS per `examples/data/data_files_needed.txt`.
Fixtures must not commit the data; record the manifest and expect the harness
to skip cleanly when it is absent.

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

### Hazard: the import-time kwargs metaprogramming

**Read this before touching `Configuration`.** `convenience.py` lines
2142–2517 (`# === ITV OPTIONS SPEC: begin/end ===`) contain ~375 lines of
import-time metaprogramming that T2 will disturb directly:

- `BuildViewRowOptions` is a frozen dataclass that is the **authoritative,
  hand-maintained** spec of every kwarg `_build_view_row` accepts. It is *not*
  derived from the signature — only defaults/annotations are pulled from the
  live function.
- `_itv__compute_option_exclusions()` statically analyses, via `ast`, which of
  the 11 `plot_*` methods actually consume or forward each kwarg, so a method
  is not advertised as accepting options it silently drops.
- `_itv__augment_docs_from_spec()` appends "Other Parameters" to each
  docstring at import time.
- `_itv__install_signatures_from_spec()` rewrites `func.__signature__` on
  every method with a `**kwargs` catch-all, so `help()` and Jupyter
  tab-completion show real parameters.
- `convenience.pyi` is **generated** — never hand-edit. Regenerate with
  `python tools/generate_convenience_stub.py`.

Why this matters for T2 specifically: splitting each `plot_*` into
`build_payload_*` plus a thin wrapper changes *exactly what the `ast` pass
analyses* — which function forwards which kwargs, and through how many hops.
Get it wrong and the advertised signatures silently drift from reality, or
import fails in a way that points at the metaprogramming rather than at the
refactor that caused it.

Practical rules for T2:
1. Move code in small steps, re-importing the package after each one. Failures
   surface at **import time**, not at call time.
2. After any change to which method forwards which kwargs, update
   `BuildViewRowOptions` by hand and regenerate the stub.
3. Decide early whether `build_payload_*` functions participate in the spec
   system at all. Keeping them outside it — plain explicit signatures, with
   the metaprogramming confined to the back-compat `plot_*` wrappers — is
   simpler and is the recommended default.
4. Treat "does `import integrative_transcriptomics_viewer` still work, and
   does `help(Configuration.plot_exons)` still show real parameters?" as a
   per-step check, not an end-of-task one.

This is why T2 is not a plan-then-delegate task — see the execution guidance
in [README.md](README.md#execution-order).

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

At that point, run **stage 2** of the spike
([07](07-spike-and-decision-log.md)) and choose a path.

## Before starting

Run **stage 1** of the spike first (Q1 Export SVG quality, Q3 bundle size).
Neither needs ITV code, both take about a day, and Q1 is existential — if
JBrowse's vector export is not good enough for figures, Path A is dead and the
trunk should be designed with Path B in mind instead. Q3's answer also feeds
T4, since the packaging strategy sets the payload size budget.
