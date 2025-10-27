ITV API
============================================
.. rubric:: Reference specific configuration

Objects of type :class:`~integrative_transcriptomics_viewer.Configuration` store information related to a reference to allow easily plotting a number of standard views/reports for the BAMs provided. The :doc:`Configuration documentation <convenience_configuration/configuration>` covers constructor arguments and links to one-page references for each plotting helper.

.. toctree::
   :maxdepth: 1
   :titlesonly:
   :caption: Reference helpers

   convenience_configuration/configuration
   convenience_configuration/plotting
   annotation_matching
   cell_barcode
   classification


.. rubric:: Plotting
   A collection of methods to easily generate views data aligned to a given :class:`~integrative_transcriptomics_viewer.Configuration`.

.. rubric:: Cell barcode

See :doc:`cell_barcode` for the abstract :class:`~integrative_transcriptomics_viewer.cellbarcode.CellBarcode` API and implementations covering common barcode encodings.


.. rubric:: Read classification

Read-level classification helpers are documented on :doc:`classification`, which covers the abstract :class:`~integrative_transcriptomics_viewer.classification.Classification` interface and its built-in implementations.
