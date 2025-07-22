MDL helpers for generating views and reports
============================================
.. contents:: :local:


===
API
===

Reference specific configuration
================================

Objects of type :class:`~genomeview.Configuration` store information related to a reference to allow easily plotting a number of standard views/reports for the BAMs provided.

.. autoclass:: genomeview.Configuration
   :members:


Cell Barcode
============

Implementations for obtaining cell barcode information for the most common encodings.

.. autoclass:: genomeview.HaasStyleCellBarcode
   :members:

.. autoclass:: genomeview.ONTCellBarcode
   :members:

.. autoclass:: genomeview.StandardCellBarcode
   :members:


Read Classification
===================

Implementation for obtaining the classification of a read for the most common encodings.

.. autoclass:: genomeview.IsoQuantClassification
   :members:

