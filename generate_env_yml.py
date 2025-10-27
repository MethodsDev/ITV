#!/usr/bin/env python3
from pathlib import Path
import re
from textwrap import dedent


REPO_ROOT = Path(__file__).resolve().parent
PACKAGE_INIT = REPO_ROOT / "src" / "integrative_transcriptomics_viewer" / "__init__.py"
ENVIRONMENT_FILE = REPO_ROOT / "environment.yml"
PACKAGE_REQUIREMENTS = REPO_ROOT / "requirements.txt"
DOCS_REQUIREMENTS = REPO_ROOT / "docs" / "requirements.txt"

PYTHON_DEPENDENCIES = [
    "cairosvg",
    "intervaltree",
    "ipywidgets>=8",
    "numpy",
    "pandas",
    "pyBigWig",
    "pysam",
]

CONDA_ONLY_DEPENDENCIES = [
    "git",
    "jupyterlab",
    "cython",
    "resvg",
    "py-open-fonts",
    "ipykernel",
    "samtools",
]

DOCS_DEPENDENCIES = [
    "sphinx>=5.0",
    "pydata-sphinx-theme",
    "ipywidgets>=8",
]


def parse_version(init_contents: str) -> str:
    """Extract the package version from the __init__ module."""
    version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
    match = re.search(version_re, init_contents, re.MULTILINE)
    if match is None:
        raise RuntimeError("Unable to determine package version from __init__.py")
    return match.group(1)


def write_lines(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def unique(sequence: list[str]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for item in sequence:
        if item not in seen:
            ordered.append(item)
            seen.add(item)
    return ordered


def main() -> None:
    version = parse_version(PACKAGE_INIT.read_text(encoding="utf-8"))

    environment_dependencies = unique(
        ["python>=3.9", "pip"] + CONDA_ONLY_DEPENDENCIES + PYTHON_DEPENDENCIES
    )

    env_contents = dedent(
        f"""\
        name: itv-{version}
        channels:
          - conda-forge
          - bioconda
          - r
        dependencies:
        """
    ).splitlines()

    env_contents.extend(f"  - {dep}" for dep in environment_dependencies)
    write_lines(ENVIRONMENT_FILE, env_contents)

    write_lines(PACKAGE_REQUIREMENTS, PYTHON_DEPENDENCIES)
    write_lines(DOCS_REQUIREMENTS, DOCS_DEPENDENCIES)


if __name__ == "__main__":
    main()
