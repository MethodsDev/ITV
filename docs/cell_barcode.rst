Cell barcode helpers
====================

Helpers for extracting a cell barcode from a sequencing read build on the abstract :class:`~integrative_transcriptomics_viewer.cellbarcode.CellBarcode`, which defines :meth:`~integrative_transcriptomics_viewer.cellbarcode.CellBarcode.get_barcode` to return a barcode (or ``None``) for a given :class:`pysam.AlignedSegment`. Implementations differ by where the barcode is stored:

- :class:`~integrative_transcriptomics_viewer.cellbarcode.HaasStyleCellBarcode` parses read names that prefix the barcode and separate metadata with ``^``; the barcode is reverse complemented so it matches the original molecule orientation.
- :class:`~integrative_transcriptomics_viewer.cellbarcode.ONTCellBarcode` reads the ``BC`` BAM tag produced by Oxford Nanopore barcoding workflows.
- :class:`~integrative_transcriptomics_viewer.cellbarcode.StandardCellBarcode` extracts the ``CB`` BAM tag used by 10x Genomics and related pipelines.

.. toctree::
   :maxdepth: 1
   :caption: API reference

   cell_barcode/cellbarcode
   cell_barcode/haas_style_cell_barcode
   cell_barcode/ont_cell_barcode
   cell_barcode/standard_cell_barcode
