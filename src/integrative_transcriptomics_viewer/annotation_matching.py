from abc import ABC, abstractmethod


class AnnotationMatching(ABC):
    """
    Abstract Class that returns whether two annotations are considered to be matching.
    Should be inherited from and have the .match(self, query_annotation, target_annotation) method implemented where both args are strings.
    """

    @abstractmethod
    def match(self, query_annotation: str, target_annotation: str) -> bool:
        """Return True when query and target annotations should be considered matching."""
        raise NotImplementedError


class IsoquantSubstringAnnotationMatching(AnnotationMatching):
    """
    Implementation of AnnotationMatching class that matches annotations if either is fully contained within the other, or if the first field when splitting by "|" is fully contained in the other.
    """

    def match(self, query_annotation: str, target_annotation: str) -> bool:
        if query_annotation in target_annotation or \
           target_annotation in query_annotation or \
           query_annotation.split("|")[0] in target_annotation or \
           target_annotation.split("|")[0] in query_annotation:
            return True
        else:
            return False


class TaggedBAMAnnotationMatching(AnnotationMatching):
    """

    """

    def match(self, query_annotation: str, target_annotation: str) -> bool:  # target = known, query = bam tag
        queries = set(query_annotation.split(";"))
        targets = set(target_annotation.split("|"))
        intersect = queries & targets
        if len(intersect) == 0:
            return False
        else:
            return True
