MDL helpers for generating views and reports
============================================
.. contents:: :local:


===
API
===

Reference specific configuration
================================

Objects of type :class:`~integrative_transcriptomics_viewer.Configuration` store information related to a reference to allow easily plotting a number of standard views/reports for the BAMs provided.

.. autoclass:: integrative_transcriptomics_viewer.Configuration
   :members:


Cell Barcode
============

Implementations for obtaining cell barcode information for the most common encodings.

.. autoclass:: integrative_transcriptomics_viewer.HaasStyleCellBarcode
   :members:

.. autoclass:: integrative_transcriptomics_viewer.ONTCellBarcode
   :members:

.. autoclass:: integrative_transcriptomics_viewer.StandardCellBarcode
   :members:


Read Classification
===================

Implementation for obtaining the classification of a read for the most common encodings.

.. autoclass:: integrative_transcriptomics_viewer.IsoQuantClassification
   :members:

