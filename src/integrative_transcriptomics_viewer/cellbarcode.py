from integrative_transcriptomics_viewer.bam_read_operations import get_read_tag

from abc import ABC, abstractmethod
from typing import Optional

import pysam
from Bio.Seq import Seq


class CellBarcode(ABC):
    """
    Abstract Class to return a cell barcode given a read.
    Should be inherited from and have the .get_barcode(self, read) method implemented where read is a pysam.AlignmentSegment object
    """

    @abstractmethod
    def get_barcode(self, read: pysam.AlignedSegment) -> Optional[str]:
        raise NotImplementedError


class HaasStyleCellBarcode(CellBarcode):
    """
    Implementation of CellBarcode class that assumes the read name is formated with the cell barcode at the start and uses "^" as a separator after.
    """

    def get_barcode(self, read: pysam.AlignedSegment) -> Optional[str]:
        query_name = read.query_name
        if query_name is None:
            return None
        barcode = query_name.split("^", 1)[0]
        return str(Seq(barcode).reverse_complement())


class ONTCellBarcode(CellBarcode):
    """
    Implementation of CellBarcode class to use the information stored in the "BC" BAM tag.
    """

    def get_barcode(self, read: pysam.AlignedSegment) -> Optional[str]:
        barcode = get_read_tag(read, "BC")
        return str(barcode) if barcode is not None else None


# 10X, PipSeq
class StandardCellBarcode(CellBarcode):
    """
    Implementation of CellBarcode class to use the information stored in the "CB" BAM tag.
    """

    def get_barcode(self, read: pysam.AlignedSegment) -> Optional[str]:
        barcode = get_read_tag(read, "CB")
        return str(barcode) if barcode is not None else None
