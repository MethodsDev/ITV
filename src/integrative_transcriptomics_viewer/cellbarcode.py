from integrative_transcriptomics_viewer.bam_read_operations import get_read_tag

from abc import ABC, abstractmethod
from typing import Optional

import pysam
from integrative_transcriptomics_viewer.genomesource import reverse_comp


class CellBarcode(ABC):
    """
    Abstract helper that extracts a cell barcode from a read.
    Subclasses implement :meth:`get_barcode` for their storage format.
    """

    @abstractmethod
    def get_barcode(self, read: pysam.AlignedSegment) -> Optional[str]:
        """
        Return the cell barcode carried by ``read``.

        Parameters
        ----------
        read : pysam.AlignedSegment
            Read from which to recover the barcode.

        Returns
        -------
        str or None
            The barcode string, or ``None`` when not present.
        """
        raise NotImplementedError


class HaasStyleCellBarcode(CellBarcode):
    """
    Implementation of CellBarcode class that assumes the read name is formated with the cell barcode at the start and uses "^" as a separator after.
    """

    def get_barcode(self, read: pysam.AlignedSegment) -> Optional[str]:
        """Extract the barcode from the read name, reverse complemented to original orientation."""
        query_name = read.query_name
        if query_name is None:
            return None
        barcode = query_name.split("^", 1)[0]
        return reverse_comp(barcode)


class ONTCellBarcode(CellBarcode):
    """
    Implementation of CellBarcode class to use the information stored in the "BC" BAM tag.
    """

    def get_barcode(self, read: pysam.AlignedSegment) -> Optional[str]:
        """Return the value stored in the ``BC`` tag if present."""
        barcode = get_read_tag(read, "BC")
        return str(barcode) if barcode is not None else None


# 10X, PipSeq
class StandardCellBarcode(CellBarcode):
    """
    Implementation of CellBarcode class to use the information stored in the "CB" BAM tag.
    """

    def get_barcode(self, read: pysam.AlignedSegment) -> Optional[str]:
        """Return the value stored in the ``CB`` tag if present."""
        barcode = get_read_tag(read, "CB")
        return str(barcode) if barcode is not None else None
