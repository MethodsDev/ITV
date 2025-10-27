ITV Documentation
=================

.. toctree::
   :maxdepth: 1
   :titlesonly:

   usage
   convenience_configuration
   advanced_usage


What is ITV?
===================

Integrative Transcriptomics Viewer (ITV) visualizes genomic data straight from python. It is based on GenomeView. Features include:

* Easily extensible
* Integrates with `jupyter notebook <http://jupyter.readthedocs.io/en/latest/index.html>`_ / `jupyterlab <https://github.com/jupyterlab/jupyterlab>`_
* High-quality vector output to standard SVG format
* Includes built-in tracks to visualize:

    * BAMs (short and long reads)

       * Both single-ended and paired-ended views available
       * Includes a cython-optimized quick consensus module to visualize error-prone long-read data
       * Group BAM reads by tag or other features using python callbacks

    * Graphical data such as coverage tracks, wiggle files, etc

The output is suitable for static visualization in screen or print formats. ITV is not currently designed to produce interactive visualizations, although the python interface, through jupyter, provides an easy interface to quickly create new visualizations.


Installation
============

ITV requires python 3.3 or greater. The following shell command should typically suffice for installing the latest release:

.. code-block:: bash

    pip install integrative_transcriptomics_viewer

Or to install the bleeding edge from github:

.. code-block:: bash

    pip install -U git+https://github.com/MethodsDev/integrative_transcriptomics_viewer.git

To display `bigWig <https://genome.ucsc.edu/goldenpath/help/bigWig.html>`_ graphical tracks, the `pyBigWig <https://github.com/deeptools/pyBigWig>`_ python package must also be installed, eg ``pip install pyBigWig``.

For a hands-on walkthrough that mirrors the demo notebook shipped with the repository, head to :doc:`usage`.
