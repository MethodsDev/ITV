Read classification helpers
===========================

Read classification strategies inherit from :class:`~integrative_transcriptomics_viewer.classification.Classification`, whose :meth:`~integrative_transcriptomics_viewer.classification.Classification.get_classification` method returns zero or more labels for a read given a gene context. Available implementations include:

- :meth:`~integrative_transcriptomics_viewer.classification.Classification.get_classification`
  accepts a :class:`pysam.AlignedSegment` and a gene identifier and returns ``None`` when
  the read cannot be assigned, a single label, or multiple labels to flag ambiguous cases.
- :class:`~integrative_transcriptomics_viewer.classification.IsoQuantClassification`, which loads IsoQuant ``read_assignments.tsv`` files and reports the tool's assignments (tracking ambiguous reads when requested).
- :class:`~integrative_transcriptomics_viewer.classification.BAMtagClassification`, which returns values stored in an arbitrary BAM tag as strings.

.. toctree::
   :maxdepth: 1
   :caption: API reference

   classification/classification
   classification/isoquant_classification
   classification/bamtag_classification
