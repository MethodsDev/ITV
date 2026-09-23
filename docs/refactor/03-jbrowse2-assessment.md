# JBrowse 2 capability assessment

Research date: 2026-09-23. JBrowse version at time of writing: **4.3.0**.

Every claim below is tagged:
- **[verified]** — confirmed against official docs or the npm registry.
- **[likely]** — strongly implied by docs but not directly confirmed.
- **[open]** — must be answered by the spike. See
  [07-spike-and-decision-log.md](07-spike-and-decision-log.md).

## Capability matrix against ITV requirements

| ITV requirement | JBrowse 2 answer | Confidence |
|---|---|---|
| Publication vector export | Built-in **Export SVG** on every genome view | verified |
| Custom renderer also exports to SVG | Shared Canvas2D functions drive screen *and* SVG export | verified |
| Headless / scripted figure generation | `@jbrowse/img` → `jb2export` CLI, node-only, no browser | verified |
| Programmatic view control | `createViewState()`, `navToLocString()`, `showTrack()` | verified |
| Exon slices side by side | Native multi-region LGV, space-delimited locations | verified |
| Per-exon width normalisation | Not standard; needs a custom view/display subclass | open |
| Reads split by cell type / read class | `jexl` `filter`+`color` callbacks, or one adapter per split | verified |
| Coverage track per split | Quantitative tracks; custom adapter can serve coverage | likely |
| Isoform-paired annotation with its reads | No native concept; custom display type | open |
| In-memory / inlined data source | Custom adapter extending `BaseFeatureDataAdapter` | verified |
| Self-contained single-file HTML | **No built-in support.** Must be built. | verified |
| Tabbed UI | No native concept. Track categories, or own chrome. | verified |

## Findings in detail

### Rendering architecture [verified]

From the [pluggable elements docs](https://jbrowse.org/jb2/docs/developer_guides/pluggable_elements/):

> Worker thread: Fetches features via RPC, returns compact typed arrays.
> Main thread (Canvas-backed): Draws via WebGPU/WebGL2/Canvas2D for features,
> variants, wiggle, synteny, MAF, Hi-C, and GWAS.
> **All canvas displays output identical pixels on-screen and in SVG export
> through shared Canvas2D functions.** Shaders provide optional GPU
> acceleration layered above the baseline implementation.

This is significant: it is exactly the dual-serializer design that Path B would
have had to build from scratch. A custom renderer written once gets vector
export with no second code path and no drift risk.

Note also the worker-thread/typed-array split — JBrowse has already solved the
"don't block the main thread" problem that Path B would inherit.

### Pluggable element types [verified]

Ten types, registered in dependency order: adapters, text search adapters,
display types, track types, connection types, **view types**, widget types,
**RPC method types**, internet account types, add-track workflow types.

Relevant to ITV:
- **Adapter** — serve `VirtualBAM` contents, one instance per classification.
- **Display type** — control how a track draws within a view; this is where
  isoform-paired-annotation rendering would live.
- **View type** — a custom subclass for per-exon width normalisation, and/or
  the overview view.
- **RPC method** — offload heavy computation to a worker.

Constraint: *"Only one version of a given plugin can be loaded on a page, even
if multiple products use it."*

### Custom adapters [verified]

Extend `BaseFeatureDataAdapter`, implement `getFeatures(region)` returning an
**RxJS observable** of features, plus optional `getRefNames` and
`freeResources`. Region carries `refName`, `start` (0-based half-open), `end`,
`assemblyName`, `originalRefName`.

An in-memory adapter reading from an inlined global is therefore small — on the
order of 50 lines. This is the hinge for self-contained export.

### Programmatic control [verified]

`createViewState({ assembly, tracks, defaultSession, plugins, location })`,
then `session.view.navToLocString(...)` and `showTrack(id)` imperatively.
`showTrack` with an unknown id returns `undefined` and reports via `session`.
`exportSvg` is an action on the LGV state model.

Caveat for the single-view components: *"The props are initial values, like an
input's `defaultValue`."* Changing assemblies or plugins requires remounting
with a new React `key`. The `react-app2` component does not have this
limitation.

### Headless figure generation [verified]

`@jbrowse/img` provides `jb2export`:

```
jb2export --fasta yourfile.fa --bam yourfile.bam \
          --loc chr1:1,000,000-1,001,000 --out file.svg
```

Notably it *does not require a browser, not even a headless one* — it runs a
node script and uses React SSR to produce the SVG. That makes batch figure
generation from Python a plain subprocess call, no Playwright/Puppeteer.

**[open]** Whether `jb2export` can load a custom plugin. If not, batch export
needs a small custom node script using the same SSR approach.

### Multi-region view [verified]

The LGV displays multiple discontinuous regions side by side in one view,
specified space-delimited in the location box (`chr1:1..100 chr2:1..100`) or
via `displayedRegions`.

This replaces `convenience.py::_plot_exons_slices` with a first-class feature.
Because it is *one view with one set of tracks*, data sharing across slices is
automatic — no view synchronisation to build.

**[open]** ITV's `normalize_interval_width` gives each exon equal display width
regardless of bp span. JBrowse multi-region very likely uses uniform
bases-per-pixel across regions. Assessment: this is a scaling factor per region
in the coordinate transform, so a custom view subclass should be tractable, but
it touches the block-layout machinery and is the single most uncertain piece of
Path A.

### Embedded components [verified]

| Feature | LGV2 / CGV2 | react-app2 | Full app |
|---|---|---|---|
| View types | one only | all | all |
| Multiple views | no | yes | yes |
| Plugin support | limited | yes | yes |
| Session management | none built-in | your app's state | full save/load/autosave |
| Feature details | dialog | drawer | oriented drawer |

Packages: `@jbrowse/react-linear-genome-view2`,
`@jbrowse/react-circular-genome-view2`, `@jbrowse/react-app2`.

**`react-app2` is the right choice for ITV** — plugins and multiple views are
both required.

### Bundle size [verified, but needs refinement]

npm registry, version 4.3.0, unpacked package size:

- `@jbrowse/react-app2` — 30,141,301 bytes, 45 dependencies
- `@jbrowse/react-linear-genome-view2` — 26,463,019 bytes, 28 dependencies

These are *unpacked package* sizes including sourcemaps and multiple build
targets. The real gzipped browser bundle is substantially smaller but still
multi-MB. **[open]** Measure the actual production bundle in the spike.

### Python integration [verified]

`jbrowse-jupyter` exists (Bioinformatics 2023) but is **Dash-based** and
supports only LGV and CGV. Not a good fit for ITV. Better options: an
anywidget/ESM wrapper around `react-app2`, or simply generating HTML files and
`display()`-ing them.

## Known gaps and risks

1. **Self-contained export is not provided.** This is the piece JBrowse helps
   with least and the one ITV most needs. Must be built: a Vite bundle of app +
   plugin, with data injected as a global and a custom adapter reading from
   memory.
2. **No tab concept.** Track categories in the selector are the idiomatic
   translation; a custom chrome above `react-app2` is the alternative.
3. **React + MobX-State-Tree is mandatory.** Real learning curve, real bus-factor
   change for the project.
4. **Plugin API churn.** Major versions v1→v4 so far. A plugin is an ongoing
   maintenance commitment.
5. **App-chrome polish.** A menu-bar bug was observed on
   `https://jbrowse.org/storybook/app/` — after adding a track via the top menu
   bar, the menu bar disappears. This is *app chrome*, not core view code, and
   embedded use supplies its own shell, so the exposure is limited. Worth
   knowing as a signal about the app layer's test coverage.

## Sources

- [Pluggable elements](https://jbrowse.org/jb2/docs/developer_guides/pluggable_elements/)
- [Creating custom adapters](https://jbrowse.org/jb2/docs/developer_guides/creating_adapter/)
- [Creating a custom view type](https://jbrowse.org/jb2/docs/developer_guides/creating_view/)
- [Embedded components](https://jbrowse.org/jb2/docs/embedded_components/)
- [Embedding tutorial](https://jbrowse.org/jb2/docs/tutorials/embed_linear_genome_view/)
- [jexl callbacks](https://jbrowse.org/jb2/docs/config_guides/jexl/)
- [Customizing feature colors](https://jbrowse.org/jb2/docs/config_guides/customizing_feature_colors/)
- [Basic usage / multi-region](https://jbrowse.org/jb2/docs/user_guides/basic_usage/)
- [Export SVG, v1.2.0 release notes](https://jbrowse.org/jb2/blog/2021/05/03/v1.2.0-release/)
- [@jbrowse/img](https://www.npmjs.com/package/@jbrowse/img)
- [JBrowse Jupyter paper](https://academic.oup.com/bioinformatics/article/39/1/btad032/6989625)
- [JBrowse 2 paper, Genome Biology 2023](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02914-z)
