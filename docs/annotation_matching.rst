Annotation Matching helpers
===========================

Classification-to-annotation alignment is handled by classes that extend
:class:`~integrative_transcriptomics_viewer.annotation_matching.AnnotationMatching`.
They decide whether a classification label refers to a particular BED entry.

- :class:`~integrative_transcriptomics_viewer.annotation_matching.IsoquantSubstringAnnotationMatching`
  treats labels as matching when one string (or its prefix before ``|``) is a substring of the other—mirroring
  the IsoQuant naming convention.
- :class:`~integrative_transcriptomics_viewer.annotation_matching.TaggedBAMAnnotationMatching`
  compares the semicolon-separated labels emitted in BAM tags against the pipe-separated transcript IDs found
  in annotations, and reports a match when the two sets intersect.

.. toctree::
   :maxdepth: 1
   :caption: API reference

   annotation_matching/annotation_matching
   annotation_matching/isoquant_substring_annotation_matching
   annotation_matching/tagged_bam_annotation_matching
