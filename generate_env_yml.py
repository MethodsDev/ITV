import os

# Read the version number from the VERSION file
def get_version(string):
    """ Parse the version number variable __version__ from a script. """
    import re
    version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
    version_str = re.search(version_re, string, re.M).group(1)
    return version_str

version = get_version(open('src/itv/__init__.py').read())

env_yml_content = f"""
name: itv-{version}
channels:
  - conda-forge
  - bioconda
  - r
dependencies:
  - python>=3.9
  - pip
  - git
  - jupyterlab
  - pysam
  - cython
  - python
  - numpy
  - pyBigWig
  - biopython
  - pandas
  - resvg
  - cairosvg
  - ipywidgets
  - intervaltree
  - py-open-fonts
"""

with open("environment.yml", "w") as env_yml_file:
    env_yml_file.write(env_yml_content)
