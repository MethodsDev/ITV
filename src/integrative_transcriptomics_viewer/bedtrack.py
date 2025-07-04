import logging
import math
import numpy

from integrative_transcriptomics_viewer.intervaltrack import Interval, IntervalTrack
from integrative_transcriptomics_viewer.utilities import match_chrom_format

DEFAULT_BED_FIELD_DEFS = {
    "chrom":0,
    "start":1,
    "end":2,
    "name":3,
    "strand":5,
    "coding_start":6,
    "coding_end":7,
    "color":8,
    "exon_lengths":10,
    "exon_starts":11,
}


class Transcript:
    def __init__(self, chrom, start, end, strand=True, name=None, coding_start=None,
                 coding_end=None, exons=None, color=None, **kwargs):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

        self.name = name

        self.coding_start = None
        if coding_start is not None:
            self.coding_start = int(coding_start)
        
        self.coding_end = None
        if coding_end is not None:
            self.coding_end = int(coding_end)

        if exons is None:
            exons = [(self.start, self.end)]
        self.exons = exons # list of (start, end) absolute coordinates

        self.color = color

        self.kwargs = kwargs

def tx_from_bedfields(bedfields, field_defs=None):
    """
    create a Transcript instance from bed fields (ie the result of
    bed_line.strip().split())
    """

    if field_defs is None:
        field_defs = DEFAULT_BED_FIELD_DEFS

    values = {}
    for field_name, field in field_defs.items():
        cur_value = None
        if len(bedfields) > field:
            cur_value = bedfields[field]
        values[field_name] = cur_value
    
    assert values["chrom"] is not None
    assert values["start"] is not None
    assert values["end"] is not None

    # if values["strand"] is None:
    #     values["strand"] = True
    # else:
    #     if values["strand"] == "+":
    #         values["strand"] = True
    #     else:
    #         values["strand"] = False

    if values["strand"] == "-":
        values["strand"] = False
    else:  # "+" or None
        values["strand"] = True
    
    values["exons"] = None
    if values["exon_starts"]:
        exon_lengths = [int(x) for x in values.pop("exon_lengths").split(",") if x]
        exon_starts = [int(x) for x in values.pop("exon_starts").split(",") if x]

        values["exons"] = [(int(values["start"])+exon_start, int(values["start"])+exon_start+exon_length)
                 for exon_start, exon_length in zip(exon_starts, exon_lengths)]

    return Transcript(**values)

def bed_fetch(path, chrom, start, end, field_defs=None):
    if type(path) is VirtualBEDTrack:
        yield from path.fetch(chrom, start, end, field_defs=field_defs)
        return

    try:
        yield from fetch_from_tabix(path, chrom, start, end, field_defs=field_defs)
        return
    except:
        pass

    try:
        yield from fetch_from_bigbed(path, chrom, start, end, field_defs=field_defs)
        return
    except ImportError:
        logging.warn("Unable to import pyBigWig")
    except:
        # raise
        pass

    try:
        yield from fetch_from_plainbed(path, chrom, start, end, field_defs=field_defs)
        return
    except:
        raise
        pass

    raise NotImplementedError("Not sure how to handle this file: {}".format(path))


def fetch_from_tabix(path, chrom, start, end, field_defs=None):
    import pysam

    bed = pysam.TabixFile(path)

    chrom = match_chrom_format(chrom, bed.contigs)
    for locus in bed.fetch(chrom, start, end):
        locus = locus.split()
        yield tx_from_bedfields(locus, field_defs=field_defs)

def fetch_from_bigbed(path, chrom, start, end, field_defs=None):
    import pyBigWig

    bed = pyBigWig.open(path)
    assert bed.isBigBed(), "Oops, for some reason I was expecting a bed file: {}".format(path)

    chrom = match_chrom_format(chrom, bed.chroms().keys())
    for cur_start, cur_end, bed_line in bed.entries(chrom, start, end):
        bed_line = bed_line.split()
        yield tx_from_bedfields([chrom, cur_start, cur_end] + bed_line, field_defs=field_defs)

def fetch_from_plainbed(path, chrom, start, end, field_defs=None):
    found_chrom = False
    for line in open(path):
        fields = line.strip().split()
        if fields[0] != chrom: continue
        found_chrom = True

        cur_start, cur_end = fields[1:3]
        if int(cur_end) < start or int(cur_start) > end: continue
        yield tx_from_bedfields(fields, field_defs=field_defs)

    if not found_chrom:
        warning = "Didn't find chromosome {}; make sure it's formatted correctly (eg 'chr1' vs '1')".format(chrom)
        logging.warn(warning)


class BEDTrack(IntervalTrack):
    def __init__(self, bed_path, name=None):
        """
        Args:
            name (str): name of the track
            bed_path (str): path of the bed file to display

        """
        super().__init__([], name=name)
        
        self.bed_path = bed_path
        self.intervals = self

        self.draw_locus_labels = True
        self.include_locus_fn = None

        self.row_height = 12
        self.thick_width = self.row_height
        self.thin_width = 5
        self.thinnest_width = 1

        self.min_exon_width = 1

        self.field_defs = None

    def fetch(self):
        """
        iterator over reads from the bed file
        """
        chrom = self.scale.chrom
        start, end = self.scale.start, self.scale.end
        
        for locus in bed_fetch(self.bed_path, chrom, start, end, field_defs=self.field_defs):
            if not self.include_locus_fn or self.include_locus_fn(locus):
                yield locus

    def __iter__(self):
        c = 0
        for i, tx in enumerate(self.fetch()):
            c += 1
            id_ = tx.name if tx.name else "feature"
            id_ = id_ + "_" + str(i)

            interval = Interval(id_, tx.chrom, tx.start, tx.end, tx.strand)

            interval.tx = tx
            if self.draw_locus_labels and tx.name:
                interval.label = tx.name
            yield interval


    def draw_interval(self, renderer, interval):
        if interval.id not in self.intervals_to_rows:
            return

        interval_pixel_width = self.scale.relpixels(interval.tx.end-interval.tx.start)
        if interval_pixel_width < 12:
            # could probably improve on this
            yield from super().draw_interval(renderer, interval)
            return

        row = self.intervals_to_rows[interval.id]
        top = row*(self.row_height+self.margin_y)
        top_thin = top + self.row_height/2 - self.thin_width/2
        midline = top + self.row_height/2 - self.thinnest_width/2
        
        color = self.color_fn(interval)
        temp_label = interval.label
        if interval.label is None:
            temp_label = interval.id
        
        tx = interval.tx

        # Draw the thin lines between "exons", along with arrows pointing in transcript direction
        for i in range(len(tx.exons)-1):
            cur_exon = tx.exons[i]
            next_exon = tx.exons[i+1]

            if cur_exon[1] < self.scale.start:
                if next_exon[0] < self.scale.start:
                    continue
                cur_start = self.scale.topixels(self.scale.start - 1)
            else:
                cur_start = self.scale.topixels(cur_exon[1])

            if next_exon[0] > self.scale.end:
                if cur_exon[1] > self.scale.end:
                    continue
                cur_end = self.scale.topixels(self.scale.end + 1)
            else:
                cur_end = self.scale.topixels(next_exon[0])

            direction = "right" if interval.strand else "left"
            n_arrows = int(round((cur_end-cur_start) / (self.row_height*0.75)))

            # print("  ", i, cur_exon, next_exon)

            arrows = (numpy.arange(1, n_arrows+1) / (n_arrows+1))# * 0.8 + 0.1

            # print(tx.exons[i])
            yield from renderer.line_with_arrows(cur_start, midline, cur_end, midline,
                direction=direction, color=color, arrows=arrows, filled=False,
                arrow_scale=self.thinnest_width*0.4, arrowKwdArgs={"stroke-width":self.thinnest_width*0.75})

        # print(3)

        # Draw the "exons", both thin (non-coding/UTR) and thick (coding)
        # print(interval, tx.exons)
        for which in ["thin", "thick"]:
            for cur_start, cur_end in tx.exons:

                # out of viewing bounds
                if cur_end < self.scale.start or cur_start > self.scale.end:
                    continue

                if which == "thick":
                    cur_y = top
                    cur_width = self.thick_width

                    if tx.coding_start is None: continue
                    if cur_end < tx.coding_start: continue
                    if cur_start > tx.coding_end: continue

                    cur_start = max(tx.coding_start, cur_start)
                    cur_end = min(tx.coding_end, cur_end)
                else:
                    if (tx.coding_start and 
                        (cur_start > tx.coding_start) and (cur_end < tx.coding_end)): continue

                    cur_y = top_thin
                    cur_width = self.thin_width

                cur_start = self.scale.topixels(cur_start)
                cur_end = self.scale.topixels(cur_end)
                width = cur_end - cur_start

                if width < self.min_exon_width:
                    cur_start -= self.min_exon_width / 2
                    width = self.min_exon_width

                yield from renderer.rect(cur_start, cur_y, width, cur_width, fill=color, 
                                         **{"stroke":"none", "id":temp_label})


        if interval.label is not None:
            end = self.scale.topixels(interval.end)

            yield from renderer.text(end+self.label_distance, top+self.row_height-2, interval.label, anchor="start")


class VirtualBEDTrack(BEDTrack):

    def __init__(self, transcripts=None, name=None):
        super().__init__(None, name=name)
        if transcripts is None:
            transcripts = []
        self.transcripts = transcripts

    def index(self, bed_path, chrom, start, end, field_defs=None):
        if self.transcripts is None:
            self.transcripts = []

        for transcript in bed_fetch(bed_path, chrom, start, end, field_defs):
            self.transcripts.append(transcript)


    def fetch(self, chrom=None, start=None, end=None, field_defs=None):
        yield from self.transcripts
        return

    def copy(self):
        return VirtualBEDTrack(transcripts=self.transcripts.copy(), name=self.name)

