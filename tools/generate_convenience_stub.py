#!/usr/bin/env python3
"""Generate the convenience.pyi stub from the runtime signature spec.

Usage:
    python tools/generate_convenience_stub.py [--output PATH]

The script imports integrative_transcriptomics_viewer.convenience, so the
package dependencies must be available in the current environment.
"""

from __future__ import annotations

import argparse
import inspect
import types as pytypes
import sys
from collections import OrderedDict
from pathlib import Path
from typing import Any, Iterable

import collections.abc as cabc
import importlib.util
import types as py_mod_types
import typing
from types import SimpleNamespace


def _install_stub_modules() -> None:
    """Install minimal stub modules for optional dependencies if missing."""

    if "cairosvg" not in sys.modules:
        import types as _types

        cairosvg_stub = _types.ModuleType("cairosvg")

        def _svg2png(*args: Any, **kwargs: Any) -> bytes:
            return b""

        cairosvg_stub.svg2png = _svg2png  # type: ignore[attr-defined]
        sys.modules["cairosvg"] = cairosvg_stub

    if "ipywidgets" not in sys.modules:
        import types as _types

        class _Widget:  # minimal placeholder
            def __init__(self, *args: Any, **kwargs: Any) -> None:
                pass

        ipywidgets_stub = _types.ModuleType("ipywidgets")
        ipywidgets_stub.VBox = _Widget
        ipywidgets_stub.HTML = _Widget
        ipywidgets_stub.Dropdown = _Widget
        ipywidgets_stub.Stack = _Widget
        ipywidgets_stub.Tab = _Widget
        ipywidgets_stub.jslink = lambda *args, **kwargs: None

        widgets_mod = _types.ModuleType("ipywidgets.widgets")
        widgets_mod.widget_selectioncontainer = SimpleNamespace(Tab=_Widget)
        widgets_mod.widget_box = SimpleNamespace(VBox=_Widget)

        ipywidgets_stub.widgets = widgets_mod  # type: ignore[attr-defined]

        sys.modules["ipywidgets"] = ipywidgets_stub
        sys.modules["ipywidgets.widgets"] = widgets_mod
        sys.modules["ipywidgets.widgets.widget_selectioncontainer"] = SimpleNamespace(Tab=_Widget)
        sys.modules["ipywidgets.widgets.widget_box"] = SimpleNamespace(VBox=_Widget)
        sys.modules["ipywidgets.embed"] = SimpleNamespace(
            embed_minimal_html=lambda *args, **kwargs: None,
            dependency_state=lambda *args, **kwargs: None,
        )

    if "intervaltree" not in sys.modules:
        import types as _types

        class _Interval:
            def __init__(self, begin: int = 0, end: int = 0, data: Any = None):
                self.begin = begin
                self.end = end
                self.data = data
                self.chrom = ""
                self.strand = True

            def __iter__(self):
                return iter((self.begin, self.end, self.data))

        class _IntervalTree(list):
            pass

        intervaltree_stub = _types.ModuleType("intervaltree")
        intervaltree_stub.Interval = _Interval  # type: ignore[attr-defined]
        intervaltree_stub.IntervalTree = _IntervalTree  # type: ignore[attr-defined]
        sys.modules["intervaltree"] = intervaltree_stub

    if "pysam" not in sys.modules:
        import types as _types

        pysam_stub = _types.ModuleType("pysam")
        sys.modules["pysam"] = pysam_stub


def _install_placeholder_submodules(package_name: str) -> None:
    """Provide placeholder submodules so relative imports succeed without loading real code."""

    class _PlaceholderModule(py_mod_types.ModuleType):
        def __getattr__(self, name: str) -> Any:
            placeholder = type(name, (), {})
            setattr(self, name, placeholder)
            return placeholder

    for name in [
        "utilities",
        "axis",
        "genomeview",
        "genomesource",
        "track",
        "bamtrack",
        "bedtrack",
        "cellbarcode",
        "templates",
        "bam_read_operations",
        "svg",
        "export",
        "quickconsensus",
        "graphtrack",
        "intervaltrack",
    ]:
        full_name = f"{package_name}.{name}"
        if full_name in sys.modules:
            continue
        module = _PlaceholderModule(full_name)
        sys.modules[full_name] = module
        parent = sys.modules.get(package_name)
        if parent is not None:
            setattr(parent, name, module)


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REPO_ROOT / "src" / "integrative_transcriptomics_viewer" / "convenience.pyi"


def _ensure_on_path() -> None:
    src = REPO_ROOT / "src"
    if str(src) not in sys.path:
        sys.path.insert(0, str(src))


def _format_annotation(
    annotation: Any,
    typing_imports: set[str],
) -> str:
    """Convert a Python annotation object into a stub-friendly string."""
    if annotation is inspect._empty:
        typing_imports.add("Any")
        return "Any"
    if annotation is Any:
        typing_imports.add("Any")
        return "Any"
    if annotation is None or annotation is type(None):
        return "None"
    if isinstance(annotation, str):
        return annotation

    origin = typing.get_origin(annotation)
    args = typing.get_args(annotation)

    union_types = {typing.Union}
    union_type = getattr(pytypes, "UnionType", None)
    if union_type is not None:
        union_types.add(union_type)

    if origin in union_types:
        parts = [_format_annotation(arg, typing_imports) for arg in args]
        if len(parts) == 2 and "None" in parts:
            other = next(p for p in parts if p != "None")
            typing_imports.add("Optional")
            return f"Optional[{other}]"
        typing_imports.add("Union")
        return f"Union[{', '.join(parts)}]"

    if origin in (typing.Literal, getattr(typing, "Literal", None)):
        typing_imports.add("Literal")
        literal_values = ", ".join(repr(arg) for arg in args)
        return f"Literal[{literal_values}]"

    if origin in (typing.Callable, cabc.Callable):
        typing_imports.add("Callable")
        if args == (Ellipsis,):
            return "Callable[..., Any]"
        if args:
            *param_types, ret = args
            if param_types == [Ellipsis]:
                params = "..."
            else:
                params = ", ".join(_format_annotation(arg, typing_imports) for arg in param_types)
            ret_str = _format_annotation(ret, typing_imports)
            return f"Callable[[{params}], {ret_str}]"
        typing_imports.add("Any")
        return "Callable[..., Any]"

    mapping_origins = {typing.Mapping, cabc.Mapping, dict}
    if origin in mapping_origins:
        typing_imports.update({"Mapping", "Any"})
        if args:
            key_str = _format_annotation(args[0], typing_imports)
            val_str = _format_annotation(args[1], typing_imports)
            return f"Mapping[{key_str}, {val_str}]"
        return "Mapping[str, Any]"

    sequence_origins = {typing.Sequence, cabc.Sequence, list, tuple, set, frozenset}
    if origin in sequence_origins:
        typing_imports.add("Sequence")
        if args:
            inner = _format_annotation(args[0], typing_imports)
            return f"Sequence[{inner}]"
        typing_imports.add("Any")
        return "Sequence[Any]"

    iterable_origins = {typing.Iterable, cabc.Iterable}
    if origin in iterable_origins:
        typing_imports.add("Iterable")
        inner = _format_annotation(args[0], typing_imports) if args else "Any"
        typing_imports.add("Any")
        return f"Iterable[{inner}]"

    if isinstance(annotation, type):
        if annotation.__module__ == "builtins":
            return annotation.__name__
        typing_imports.add("Any")
        return "Any"

    typing_imports.add("Any")
    return "Any"


def _render_signature(
    name: str,
    sig: inspect.Signature,
    typing_imports: set[str],
    option_types: dict[str, Any],
) -> list[str]:
    """Return the lines representing a single method stub."""
    params: Iterable[inspect.Parameter] = sig.parameters.values()
    params = list(params)

    if not params:
        header = f"    def {name}(self) -> Any: ..."
        typing_imports.add("Any")
        return [header]

    lines: list[str] = [f"    def {name}("]
    indent = " " * 12
    kw_only_inserted = False

    for param in params:
        if param.kind == inspect.Parameter.VAR_KEYWORD:
            piece = f"**{param.name}"
        elif param.kind == inspect.Parameter.VAR_POSITIONAL:
            piece = f"*{param.name}"
        else:
            if param.kind == inspect.Parameter.KEYWORD_ONLY and not kw_only_inserted:
                lines.append(f"{indent}*,")
                kw_only_inserted = True
            piece = param.name

        if param.kind == inspect.Parameter.VAR_KEYWORD:
            ann_obj = Any
        elif param.name in option_types:
            ann_obj = option_types[param.name]
        elif param.name == "bams_dict":
            ann_obj = typing.Mapping[str, Any]
        else:
            ann_obj = param.annotation

        if ann_obj is inspect._empty:
            if param.name != "self":
                typing_imports.add("Any")
            if param.name == "self":
                ann_str = None
            else:
                ann_str = "Any"
        else:
            ann_str = _format_annotation(ann_obj, typing_imports)

        if ann_str:
            piece = f"{piece}: {ann_str}"

        if param.kind != inspect.Parameter.VAR_POSITIONAL and param.kind != inspect.Parameter.VAR_KEYWORD:
            if param.default is not inspect._empty:
                piece = f"{piece} = ..."

        lines.append(f"{indent}{piece},")

    return_annotation = _format_annotation(sig.return_annotation, typing_imports)
    lines.append(f"    ) -> {return_annotation}: ...")
    return lines


def generate_stub(output: Path) -> None:
    _ensure_on_path()
    _install_stub_modules()

    package_name = "integrative_transcriptomics_viewer"
    if package_name not in sys.modules:
        pkg = py_mod_types.ModuleType(package_name)
        pkg.__path__ = [str(REPO_ROOT / "src" / package_name)]
        sys.modules[package_name] = pkg
        _install_placeholder_submodules(package_name)

    convenience_path = REPO_ROOT / "src" / package_name / "convenience.py"
    spec = importlib.util.spec_from_file_location(
        f"{package_name}.convenience",
        convenience_path,
        submodule_search_locations=[str(convenience_path.parent)],
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to create module spec for convenience.py")
    convenience = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = convenience
    spec.loader.exec_module(convenience)

    if hasattr(convenience, "_itv__install_signatures_from_spec"):
        convenience._itv__install_signatures_from_spec()

    config_cls = convenience.Configuration
    option_types: dict[str, Any] = {}
    try:
        import dataclasses

        row_sig = inspect.signature(config_cls._build_view_row)  # type: ignore[attr-defined]
        row_params = row_sig.parameters
        for field in dataclasses.fields(convenience.BuildViewRowOptions):  # type: ignore[attr-defined]
            param = row_params.get(field.name)
            annotation = None
            if param and param.annotation is not inspect._empty:
                annotation = param.annotation
            elif field.type is not None:
                annotation = field.type
            option_types[field.name] = annotation if annotation is not None else Any
    except Exception:
        option_types = {}
    methods: OrderedDict[str, inspect.Signature] = OrderedDict()
    for name, func in sorted(config_cls.__dict__.items()):
        if name.startswith("_"):
            continue
        if not inspect.isfunction(func):
            continue
        sig = inspect.signature(func)
        if not any(param.kind == inspect.Parameter.VAR_KEYWORD for param in sig.parameters.values()):
            continue
        methods[name] = sig

    typing_imports: set[str] = set()
    method_blocks: list[str] = []

    for name, sig in methods.items():
        method_lines = _render_signature(name, sig, typing_imports, option_types)
        method_blocks.append("\n".join(method_lines))

    typing_imports.discard("Any")  # ensure consistent ordering later; add back explicitly
    typing_imports_sorted = ["Any"] + sorted(typing_imports)

    header = [
        "# Auto-generated by tools/generate_convenience_stub.py; do not edit by hand.",
        "from __future__ import annotations",
        f"from typing import {', '.join(typing_imports_sorted)}",
        "",
        "",
        "class Configuration:",
    ]

    body = "\n\n".join(method_blocks) if method_blocks else "    pass"
    content = "\n".join(header) + "\n" + body + "\n"

    output.write_text(content, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Target .pyi file path.")
    args = parser.parse_args()
    generate_stub(args.output)


if __name__ == "__main__":
    main()
