import collections
import gzip
import math
import os
import pysam
import inspect
import ipywidgets as widgets
import re
import time
from functools import partial

from intervaltree import Interval, IntervalTree

import integrative_transcriptomics_viewer as itv
from integrative_transcriptomics_viewer import utilities
from integrative_transcriptomics_viewer.bam_read_operations import *
from integrative_transcriptomics_viewer.utilities import my_hook_compressed


def visualize_data(file_paths, chrom, start, end, reference_path=None, 
                   width=900, axis_on_top=False):
    """
    Creates a GenomeView document to display the data in the specified
    files (eg bam, bed, etc).

    Args:
        file_paths: this specifies the file paths to be rendered. It must be 
            either a list/tuple of the paths, or a dictionary mapping 
            {track_name:path}. (If you are using a python version prior to 3.6, 
            use collections.ordereddict to ensure the order remains the same.)
            Currently supports files ending in .bam, .cram, .bed, .bed.gz, 
            .bigbed, or .bigwig (or .bw). Most of these file types require a
            separate index file to be present (eg a .bam.bai or a .bed.gz.tbi 
            file must exist).
        chrom: chromosome (or contig) to be rendered
        start: start coordinate of region to be rendered
        end: end coordinate of region to be rendered
        reference_path: path to fasta file specifying reference genomic 
            sequence. This is required in order to display mismatches
            in bam tracks.
        width: the pixel width of the document
        axis_on_top: specifies whether the axis should be added at the bottom
            (default) or at the top
    """
    if reference_path is not None:
        source = itv.FastaGenomeSource(reference_path)
    else:
        source = None

    doc = itv.Document(width)
    
    view = itv.GenomeView(chrom, start, end, "+", source)
    doc.add_view(view)

    def add_axis():
        axis_track = itv.Axis("axis")
        view.add_track(axis_track)

    if axis_on_top:
        add_axis()

    if isinstance(file_paths, collections.abc.Mapping):
        names = file_paths.keys()
        file_paths = [file_paths[name] for name in names]
    else:
        names = [None] * len(file_paths)
        file_paths = file_paths
        
    for name, path in zip(names, file_paths):
        if path.lower().endswith(".bam") or path.lower().endswith(".cram"):
            if utilities.is_paired_end(path):
                cur_track = itv.PairedEndBAMTrack(path, name=name)
            else:
                cur_track = itv.SingleEndBAMTrack(path, name=name)
                if utilities.is_long_frag_dataset(path):
                    cur_track.min_indel_size = 5

        elif path.lower().endswith(".bed") or path.lower().endswith(".bed.gz") or path.lower().endswith(".bigbed") or path.lower().endswith(".bb"):
            cur_track = itv.BEDTrack(path, name=name)

        elif path.lower().endswith(".bigwig") or path.lower().endswith(".bw"):
            cur_track = itv.BigWigTrack(path, name=name)

        else:
            suffix =  os.path.basename(path)
            raise ValueError("Unknown file suffix: {}".format(suffix))

        view.add_track(cur_track)

    if not axis_on_top:
        add_axis()

    return doc


### newly added wrappers

class TighterSingleEndBAMTrack(itv.SingleEndBAMTrack):
    def __init__(self, *args, **kwdargs):
        super().__init__(*args, **kwdargs)
        self.row_height = 3
        self.margin_y = 2

def color_from_bed(interval):
    if interval.tx.color and len(interval.tx.color.split(",")) == 3:
        # hex_colors = interval.tx.color.split(",")
        # return "#{0:02x}{1:02x}{2:02x}".format(int(hex_colors[0]), int(hex_colors[1]), int(hex_colors[2]))
        return "rgb(" + interval.tx.color + ")"
    else:
        return itv.color_by_strand(interval)


def make_bed_track(bed, name=None): # , chrom=None, start=None, end=None):
    if isinstance(bed, itv.VirtualBEDTrack):
        bed_track = bed.copy()
    else:
        bed_track = itv.BEDTrack(bed, name)
    return bed_track


# adding chrom and strand accessors for Interval from data slot based on the usage made in the code below
@property
def interval_chrom(self):
    return self.data[:-1]

Interval.chrom = interval_chrom

@property
def interval_strand(self):
    return self.data[-1:]

Interval.strand = interval_strand


def interval_data_reduce(current_data, new_data):
    if current_data == new_data:
        return current_data
    else:
        return None



class Configuration:
    """
    Creates a helper object that is reference specific and can be provided annotations
    to display (BED files) and internally index (GTF files) to allow to easily generate 
    standard views and reports by only providing feature IDs/names and BAMs.
    Currently methods can return HTML, SVGs or ipywidgets but this will be standardized.

    Parameters
    ----------        
    genome_fasta : str
        Path to the FASTA sequence of the reference, preferably in an indexed format.
    gtf_annotation : str
        GTF file to index features from so they can later be queried by ID or name.
    bed_annotation : dict or list, optional
        Dictionary or list of BED annotation files, BGzipped and indexed. Can be an empty list.
    bed_color_fn : typing.Callable[[:py:class:`integrative_transcriptomics_viewer.intervaltrack.Interval`], str], optional
        Method to determine how the BED annotations will be colored (default uses BED color code).
    """

    def __init__(self, genome_fasta, bed_annotation, gtf_annotation = None, bed_color_fn=color_from_bed):
        self._fasta = genome_fasta
        self.genome_fasta = itv.genomesource.FastaGenomeSource(genome_fasta)
        self.bed_annotation = bed_annotation
        self.bed_color_fn = bed_color_fn

        self.gene_name_to_gene_id = {}
        self.gene_id_to_gene_name = {}
        self.gene_to_transcripts = {}
        self.transcript_to_gene = {}
        self.gene_to_exons = {}
        self.transcript_to_exons = {}
        self.id_to_coordinates = {}

        if gtf_annotation:
            self.index_gtf(gtf_annotation)


    def shallow_copy(self):
        return Configuration(self._fasta, None, bed_color_fn = self.bed_color_fn)


    def index_gtf(self, gtf_annotation):
        """
        Index features from a GTF to keep internally for easy querying. Stores Gene IDs/names, Transcript IDs/names, Exon IDs

        Parameters
        ----------
        gtf_annotation (str):
            Path of the gtf file to index.
        """

        gene_id_regex = re.compile(r'gene_id "([a-zA-Z0-9\.\-_]+)";')
        gene_name_regex = re.compile(r'gene_name "([a-zA-Z0-9\.\-_]+)";')
        transcript_id_regex = re.compile(r'transcript_id "([a-zA-Z0-9\._\^\-]+=?)";')
        exon_id_regex = re.compile(r'exon_id "([a-zA-Z0-9\.\-_]+)";')
        # transcript_name_regex = re.compile('transcript_name "([a-zA-Z0-9\.]+)";')
        # protein_id_regex = re.compile('protein_id "([a-zA-Z0-9\.]+)";')

        with my_hook_compressed(gtf_annotation, "r") as gtf_file:
            for entry in pysam.tabix_iterator(gtf_file, pysam.asGTF()):
                if entry.feature == "gene":
                    res = gene_id_regex.search(entry.attributes)
                    if res:
                        gene_id = res.group(1)
                    else:
                        print("missing gene_id in a gene entry, skipping entry")


                    res = gene_name_regex.search(entry.attributes)
                    if res:
                        gene_name = res.group(1)
                        self.gene_name_to_gene_id[gene_name] = gene_id
                        self.gene_id_to_gene_name[gene_id] = gene_name

                    self.id_to_coordinates[gene_id] = Interval(entry.start, entry.end, entry.contig + entry.strand)
                    self.gene_to_transcripts[gene_id] = []
                    self.gene_to_exons[gene_id] = IntervalTree()

                elif entry.feature == "transcript":
                    gene_id = None
                    res = gene_id_regex.search(entry.attributes)
                    if res:
                        gene_id = res.group(1)

                        if gene_id not in self.id_to_coordinates:
                            self.id_to_coordinates[gene_id] = Interval(entry.start, entry.end, entry.contig + entry.strand)
                            self.gene_to_transcripts[gene_id] = []
                            self.gene_to_exons[gene_id] = IntervalTree()

                    res = transcript_id_regex.search(entry.attributes)
                    if res:
                        transcript_id = res.group(1)
                    else:
                        print("missing transcript_id in a transcript entry, skipping entry:")
                        print(entry)

                    self.id_to_coordinates[transcript_id] = Interval(entry.start, entry.end, entry.contig + entry.strand)
                    self.transcript_to_gene[transcript_id] = gene_id
                    self.transcript_to_exons[transcript_id] = IntervalTree()
                    if gene_id:
                        self.gene_to_transcripts[gene_id].append(transcript_id)

                elif entry.feature == "exon":
                    # gene_id = None
                    res = gene_id_regex.search(entry.attributes)
                    if res:
                        gene_id = res.group(1)
                        self.gene_to_exons[gene_id].add(Interval(entry.start, entry.end, entry.contig + entry.strand))

                    transcript_id = None
                    res = transcript_id_regex.search(entry.attributes)
                    if res:
                        transcript_id = res.group(1)
                        self.transcript_to_exons[transcript_id].add(Interval(entry.start, entry.end, entry.contig + entry.strand))

                    exon_id = None
                    res = exon_id_regex.search(entry.attributes)
                    if res:
                        exon_id = res.group(1)
                        if exon_id not in self.id_to_coordinates:
                            self.id_to_coordinates[exon_id] = Interval(entry.start, entry.end, entry.contig + entry.strand)


    def update_bed(self, bed_annotation):
        """
        Replace current bed annotations with those provided.

        Parameters
        ----------
        bed_annotation : dict or list
            Dictionary or list of BED annotation files, BGzipped and indexed. Can be an empty list.
        """ 
        self.bed_annotation = bed_annotation


    def get_bed_entries(self, interval):
        """
        Parameters
        ----------
        interval : :class:`~intervaltree.Interval`
            Interval to fetch annotations within.

        Returns
        -------
        All BED entries found within the provided interval across all BEDs. : dict
        """
        all_known_annotations = {}
        if self.bed_annotation:
            if type(self.bed_annotation) is list:
                for bed_path in self.bed_annotation:
                    for bed_entry in itv.bedtrack.bed_fetch(bed_path, interval.chrom, interval.begin, interval.end):
                        all_known_annotations[bed_entry.name] = itv.VirtualBEDTrack(transcripts = [bed_entry])
            elif type(self.bed_annotation) is dict:
                for bed_name, bed_path in self.bed_annotation.items():
                    for bed_entry in itv.bedtrack.bed_fetch(bed_path, interval.chrom, interval.begin, interval.end):
                        all_known_annotations[bed_entry.name] = itv.VirtualBEDTrack(transcripts = [bed_entry])
            else:
                for bed_entry in itv.bedtrack.bed_fetch(self.bed_annotation, interval.chrom, interval.begin, interval.end):
                    all_known_annotations[bed_entry.name] = itv.VirtualBEDTrack(transcripts = [bed_entry])
        return all_known_annotations


    def add_bed_track_to_view(self, view, bed_track, vertical_layout=True, strand_specific=False, use_names=True):
        bed_track.color_fn = self.bed_color_fn
        bed_track.vertical_layout = vertical_layout
        bed_track.strand_specific = strand_specific
        view.add_track(bed_track)


    def add_bed_tracks_to_view(self, view, vertical_layout=True, strand_specific=False, use_names=True):
        """
        Transparently adds BED tracks as needed to a view for all BEDs in this configuration,
        regardless of whether they reference a file or are in-memory, and of self.bed_annotation format.

        Parameters
        ----------
        view : :py:class:`integrative_transcriptomics_viewer.GenomeView`
            The view to add BED tracks to.
        vertical_layout : bool, optional
            Control if a single entry should be drawn per row within the tracks. (Default: True)
        use_names : bool or dict, optional
            Control if BED names are displayed or not (requires self.bed_annotations to be a dict). If a dict is provided, it should have the same keys as self.bed_annotation
        """

        if self.bed_annotation:
            if type(self.bed_annotation) is list:
                for bed_path in self.bed_annotation:
                    bed_track = make_bed_track(bed_path)
                    self.add_bed_track_to_view(view, bed_track, vertical_layout, strand_specific, use_names)
            elif type(self.bed_annotation) is dict:
                for bed_name, bed_path in self.bed_annotation.items():
                    if use_names:
                        if type(use_names) is dict:
                            if use_names[bed_name] is not None:
                                view.add_track(itv.track.TrackLabel(use_names[bed_name]))
                        else:
                            view.add_track(itv.track.TrackLabel(bed_name))
                    # else:
                    #     view.add_track(itv.track.TrackLabel(""))
                    virtual_bed = make_bed_track(bed_path, name="")
                    self.add_bed_track_to_view(view, virtual_bed, vertical_layout, strand_specific, use_names)
            else:
                virtual_bed = make_bed_track(self.bed_annotation)
                self.add_bed_track_to_view(view, virtual_bed, vertical_layout, strand_specific, use_names)


    def build_view_row(self, start, end, chrom, strand, bams_dict,
                       padding_perc = 0.1,
                       add_track_label = "auto",
                       add_reads_label = "auto",
                       add_coverage_label = "auto",
                       with_reads = True,
                       with_axis = True,
                       with_coverage = True,
                       with_bed = True,
                       with_bed_label = False,
                       coverage_bin_size = 0,
                       coverage_height = 100,
                       coverage_tag = None,
                       coverage_tag_fn = None,
                       coverage_by_strand = False,
                       priming_orientation = "3p",
                       strand_specific_bed = False,
                       vertical_layout_reads = False,
                       max_read_depth = 100,
                       include_secondary = False,
                       quick_consensus = False,
                       row = None, 
                       view_width = None,
                       view_margin_y = None,
                       fill_coverage = True,
                       coverage_track_max_y = None,
                       tighter_track = False,
                       **kwargs):

        """
        Core method to generate a customizable view in a new row or append it to an existing row.

        Parameters
        ----------
        start : int
            Start position of the view
        end : int
            End position of the view
        chrom : str
            Chromosome/contig of the reference
        strand : str
            Strand specificity if any (currently unused)
        bams_dict : dict
            Dictionnary of (virtual) BAMs to display (empty dict supported)

        Other Parameters
        ----------------
        padding_perc : float, optional
            Multiply the (end-start) window size by this multiplier to "pad" the view with for context. (Default: 0.1)
        add_track_label : bool, optional
            Label to add at the top of the view, default is the coordinates of the window displayed. Can be set to None. (Default: "auto")
        add_reads_label : bool, optional
            Label of each BAM track, default is the dictionnary key. Can be set to None. [future: add support for dict] (Default: "auto")
        add_coverage_label : bool, optional
            Label of BAM coverage track, default is the dictionnary key. Can be set to None. [future: add support for dict] (Default: "auto")
        with_reads : bool, optional
            Control if reads are drawn. (Default: True)
        with_axis : bool, optional
            Control if an axis is drawn. (Default: True)
        with_coverage : bool, optional
            Control if a coverage track for each BAM is drawn. (Default: True)
        with_bed : bool, optional
            Control if bed annotations are drawn. (Default: True)
        with_bed_label : bool, optional
            Control if a label is displayed for each BED source. (Default: False)
        coverage_bin_size : int, optional
            Control the bin size for binned coverage (Default: 0, meaning no binning)
        priming_orientation : str, optional
            Define what side reads are supposed to be primed to based on chemistry for binned coverage direction. One of "5p" or "3p". (Default: "5p")
        vertical_layout_reads : bool, optional
            Control if reads are displayed in a vertical manner, meaning a single read per line, or not. (Default: False)
        include_secondary : bool, optional
            Control if secondary alignments are drawn. (Default: False)
        quick_consensus : bool, optional
            Control if quick consensus (hide SNPs that have a frequency lower than 0.2) mode is used. Slows plotting. (Default: False)
        row : :class:`~itv.ViewRow`, optional
            Row to which to append the view generated by this method. If None is provided, will create a new one. (Default: None)
        view_width : int, optional
            Width of the view in "pixels". (Default: None, which uses the global default of using that size of all views' width is defined, or divide the overall width evenly otherwise)
        view_margin_y : int, optional
            Size of the margin in "pixels" at the top and bottom of the view. (Default: None, which uses the global default of 5)
        fill_coverage : bool, optional
            Boolean to control if the coverage track should be a filled solid or just an outline. (Default: True)
        coverage_track_max_y : int, optional
            Specify a value that should be the max on the y-axis of the coverage track (useful to compare BAMs or regions with different coverage depth by forcing a common axis). (Default: None)
        tighter_track : bool, optional
            Boolean to control if the reads should be displayed in an alternate mode where they are thinner on the y-axis. (default: False)

        Returns
        -------
        The ViewRow to which the view was added : :py:class:`integrative_transcriptomics_viewer.ViewRow`
        """


        padding = math.ceil((end - start) * padding_perc)

        if row is None:
            row = itv.ViewRow("row")
        gene_view = itv.GenomeView(chrom, max(0, start - padding), end + padding, strand, self.genome_fasta)
        # gene_view = itv.GenomeView(chrom, start - padding, end + padding, "+", self.source)
        if add_track_label:
            if add_track_label == "auto":
                gene_view.add_track(itv.track.TrackLabel(chrom + " : " + str(start - padding) + " - " + str(end + padding)))
            else:
                gene_view.add_track(itv.track.TrackLabel(add_track_label))

        if with_bed:
            self.add_bed_tracks_to_view(gene_view,  strand_specific = strand_specific_bed, use_names = with_bed_label)

        if with_axis:
            gene_view.add_track(itv.Axis())
        for key, value in bams_dict.items():
            if isinstance(value, itv.VirtualBAM):
                opener_kwargs = {'opener_fn': lambda x: x}
            else:
                opener_kwargs = {}


            if with_coverage:
                coverage_label = ""
                if add_coverage_label:
                    if add_coverage_label == "auto":
                        coverage_label = key
                coverage_track = itv.BAMCoverageTrack(value, name=coverage_label, **opener_kwargs)
                coverage_track.bin_size = coverage_bin_size
                coverage_track.priming_orientation = priming_orientation
                coverage_track.height = coverage_height
                coverage_track.tag = coverage_tag
                coverage_track.tag_fn = coverage_tag_fn
                coverage_track.stranded_coverage = coverage_by_strand

                if fill_coverage:
                    coverage_track.fill_coverage = True

                if coverage_track_max_y:
                    if type(coverage_track_max_y) is dict:
                        coverage_track.max_y = coverage_track_max_y[key]
                    else:
                        coverage_track.max_y = coverage_track_max_y
                gene_view.add_track(coverage_track)
            if with_reads:
                reads_label = ""
                if add_reads_label:
                    if add_reads_label == "auto":
                        reads_label = key

                if tighter_track:
                    bam_track = TighterSingleEndBAMTrack(value, name=reads_label, **opener_kwargs)
                else:
                    bam_track = itv.SingleEndBAMTrack(value, name=reads_label, **opener_kwargs)
                if include_secondary:
                    coverage_track.include_secondary = True
                    bam_track.include_secondary = True
                bam_track.max_depth = max_read_depth
                bam_track.quick_consensus = quick_consensus
                bam_track.vertical_layout = vertical_layout_reads
                gene_view.add_track(bam_track)

        if view_width:
            gene_view.pixel_width = view_width
        if view_margin_y:
            gene_view.margin_y = view_margin_y

        row.add_view(gene_view)
        return row


    def add_single_view_row_to_plot(self, doc,
                                    interval = None, 
                                    data = None, **kwargs):
        """
        Helper that calls self.build_view_row given a window to make a view for. Requires to provide either "interval" or "data".

        Parameters
        ----------
        doc : :class:`~itv.Document`
            Document object to which to add the view row.
        
        interval : :class:`~intervaltree.Interval`
            an Interval(start, end, chrom + strand)
        data : :class:`~collections.abc.Sequence`
            a Sequence (list, tuple, etc) in the format (any, Interval(start, end), chromosome, strand)
        **kwargs
            anything that can be passed to :meth:`build_view_row()`
        """

        if data is not None:
            chrom = data[2]
            strand = data[3] if len(data) > 2 else True 
            start = data[1].begin
            end = data[1].end
        elif interval is not None:
            chrom = interval.chrom
            strand = interval.strand
            start = interval.begin
            end  = interval.end
        else:
            raise ValueError("Neither an Interval or data structure has been provided.")

        row = self.build_view_row(start, end, chrom, strand, **kwargs)
        doc.elements.append(row)
        return doc


    def add_multi_view_row_to_plot(self, doc,
                                   interval_list = None, 
                                   data_list = None, **kwargs):
        """
        Helper that calls self.build_view_row given a list of windows to make views for. Requires to provide either "interval_list" or "data_list".
        
        Parameters
        ----------
        doc : :class:`~itv.Document`
            Document object to which to add the view row
        interval_list : list
            a list of :class:`~intervaltree.Interval` (start, end, chrom + strand)
        data_list : list
            a list of lists in the format (any, :class:`~intervaltree.Interval` (start, end), chromosome, strand)
        **kwargs
            anything that can be passed to :meth:`build_view_row()`
        """

        row = itv.ViewRow("row")

        if interval_list is not None:
            for interval in interval_list:
                chrom = interval.chrom
                strand = interval.strand
                start = interval.begin
                end  = interval.end
                
                self.build_view_row(start, end, chrom, strand,
                                         row = row,  **kwargs)

        elif data_list is not None:
            for data in data_list:
                chrom = data[2]
                strand = data[3] if len(data) > 2 else True 
                start = data[1].begin
                end = data[1].end
                
                self.build_view_row(start, end, chrom, strand,
                                         row = row,  **kwargs)

        doc.elements.append(row)
        return doc


    def plot_interval(self, view_width = 1600, **kwargs):
        """
        Returns a Document with a view for the interval provided. Requires to provide either "interval" or "data".

        Parameters
        ----------
        interval : :class:`~intervaltree.Interval`
            An :class:`~intervaltree.Interval` (start, end, chrom + strand)
        data : list
            A list in the format (any, :class:`~intervaltree.Interval` (start, end), chromosome, strand)
        view_width : int, optional
            width of the view in "pixels". (Default: 1600)
        **kwargs
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        A new Document that contains a view of the interval : :py:class:`integrative_transcriptomics_viewer.Document`
        """

        doc = itv.Document(view_width)
        return self.add_single_view_row_to_plot(doc, **kwargs)



    def plot_intervals(self,
                       view_width = 1600,
                       interval_list = None, 
                       data_list = None, 
                       N_per_row = 1, 
                       **kwargs):
        """
        Returns a Document with views for all the intervals provided.

        Parameters
        ----------
        interval_list : list
            a list of :class:`~intervaltree.Interval` (start, end, chrom + strand)
        data_list : list
            a list of lists in the format (any, :class:`~intervaltree.Interval` (start, end), chromosome, strand)
        view_width : int, optional
            width of the view in "pixels". (Default: 1600)
        N_per_row : int, optional
            How many views to have per row. (Default: 1)
        **kwargs
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        A new Document that contains the views of the intervals : :py:class:`integrative_transcriptomics_viewer.Document`
        """

        doc = itv.Document(view_width)

        if interval_list is not None:
            for i in range(0, len(interval_list), N_per_row):
                doc = self.add_multi_view_row_to_plot(doc,
                                                      interval_list = interval_list[i:i+N_per_row], 
                                                      data_list = None,
                                                      **kwargs)

        elif data_list is not None:
            for i in range(0, len(data_list), N_per_row):  
                doc = self.add_multi_view_row_to_plot(doc,
                                                      interval_list = None, 
                                                      data_list = data_list[i:i+N_per_row],
                                                      **kwargs)
        return doc

    def get_feature_info(self, feature):
        """
        Parameters
        ----------
        feature :str
            The name or id of the feature to look for in the internal index.

        Returns
        -------
        A (feature_id, feature_type) tuple with the id of the feature and its type (gene/transcript/exon).(None, None) if the feature was not found. : tuple(str, str)
        """

        feature_id = None
        feature_type = None
        if feature in self.gene_name_to_gene_id:
            feature_id = self.gene_name_to_gene_id[feature]
            feature_type = "gene"
        elif feature in self.gene_id_to_gene_name:
            feature_id = feature
            feature_type = "gene"
        elif feature in self.transcript_to_exons:
            feature_id = feature
            feature_type = "transcript"
        elif feature in self.id_to_coordinates:
            feature_id = feature
            feature_type = "exon"
        
        return(feature_id, feature_type)


    def get_interval_from_feature(self, feature):
        """
        Return the Interval with coordinates of a given feature.

        Parameters
        ----------
        feature : str or tuple(str, str)
            the feature name, feature id, or (feature_id, feature_type) tuple, to look for in the internal index.

        Returns
        -------
        The Interval with the coordinates of the feature. : :class:`~intervaltree.Interval`
        """

        if isinstance(feature, Interval):
            return feature
        elif isinstance(feature, tuple):
            (feature_id, feature_type) = feature
        else:
            (feature_id, feature_type) = self.get_feature_info(feature)
        return self.id_to_coordinates[feature_id]


    def get_gene_name(self, feature):
        """
        Return the name of the gene a given feature is part of.

        Parameters
        ----------
        feature : str
            the feature name, feature id, or (feature_id, feature_type) tuple, to look for in the internal index.

        Returns
        -------
        The gene name of the feature. : str
        """

        if isinstance(feature, tuple):
            (feature_id, feature_type) = feature
        else:
            (feature_id, feature_type) = self.get_feature_info(feature)

        if feature_type == "transcript":
            if feature_id not in self.transcript_to_gene:
                print("transcript not associated with any gene")
                return feature_id
            feature_id = self.transcript_to_gene[feature_id]
            feature_type = "gene"

        if feature_type == "gene":
            if feature_id in self.gene_name_to_gene_id:
                return feature_id
            elif feature_id in self.gene_id_to_gene_name:
                return self.gene_id_to_gene_name[feature_id]
            else:
                print("no gene name associated with this gene id")
                return feature_id

        print("unknown feature type")
        return feature_id


    def plot_feature(self, feature, **kwargs):
        """
        Parameters
        ----------
        feature : str
            the feature id or name of the feature of interest.
        **kwargs:
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        A Document that contains a view of the feature. : :py:class:`integrative_transcriptomics_viewer.Document`
        """

        if isinstance(feature, tuple):
            (feature_id, feature_type) = feature
        else:
            (feature_id, feature_type) = self.get_feature_info(feature)
        return self.plot_interval(interval = self.id_to_coordinates[feature_id], **kwargs)

    # plot_feature for a list in tabs
    def plot_features(self,
                      features,
                      output_format="svg",
                      **kwargs):
        """
        Parameters
        ----------
        features : list(str)
            the list of features (id or name, can be mixed) of interest.
        output_format : str, optional
            format of the views, "svg" or "png". (Default: "svg")
        **kwargs:
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        An ipywidget where each tab is a view for one of the features provided. : :py:class:`ipywidgets.widgets.widget_selectioncontainer.Tab`

        """

        features_tab = widgets.Tab()
        tab_contents = []

        for feature in features:
            tab_contents.append(self.plot_feature(feature, **kwargs).get_widget(output_format))

        features_tab.children = tab_contents
        features_tab.titles = features

        return features_tab


    def plot_read(self, read_id, bams_dict, interval="auto", output_format="svg", silence_error=False, **kwargs):
        """
        Parameters
        ----------
        read_id : str
            the id of the read of interest.
        bams_dict : dict
            the dict of (Virtual) BAMs in which the read should be found
        interval : Interval or str, optional
            "auto" to automatically select the windows to plot based on where the read aligns. Can otherwise be set to an Interval of choice, but might result in empty views when the read does not align in that window. (Default: "auto")
        output_format : str, optional
            format of the views, "svg" or "png". (Default: "svg")
        silence_error : boolean, optional
            Boolean to control if read not being found in one of the BAMs is reported (useful when providing a dict of BAMs but the read is only present in one of them). (Default: False)
        **kwargs:
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        An ipywidget of all the windows in which a given read aligns. : :py:class:`ipywidgets.widgets.widget_selectioncontainer.Tab`
        """

        (regions, virtual_bams) = find_read_in_bam(read_id, bams_dict, silence_error)

        if regions == -1:
            return None


        all_widgets = []
        for region, virtual_bam in zip(regions, virtual_bams):
            if isinstance(interval, str) and interval == "auto":
                tmp = self.plot_interval(bams_dict={"read": virtual_bam}, interval=region, **kwargs).get_widget(output_format)
            else:
                tmp = self.plot_interval(bams_dict={"read": virtual_bam}, interval=interval, **kwargs).get_widget(output_format)
            if tmp is not None:
                all_widgets.append(tmp)

        return widgets.VBox(all_widgets)


    def plot_reads(self, read_ids, bams_dict, interval, output_format="svg", **kwargs):
        """
        Parameters
        ----------
        read_ids : list(str)
            ids of the reads of interest.
        bams_dict : dict
            (virtual) BAMs with the reads to display
        interval : :class:`~intervaltree.Interval` or str
            "auto" to automatically select the windows to plot based on where the read aligns. Can otherwise be set to an Interval of choice, but might result in empty views when the read does not align in that window. (Default: "auto")
        output_format : str, optional
            format of the views, "svg" or "png". (Default: "svg")
        silence_error : boolean, optional
            control if read not being found in one of the BAMs is reported (useful when providing a dict of BAMs but the read is only present in one of them). (Default: False)
        **kwargs:
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        An ipywidget of all the windows in which the given reads align. : :py:class:`ipywidgets.widgets.widget_selectioncontainer.Tab`
        """

        if not isinstance(interval, Interval):
            print("Error, Interval() required but not provided")
            return -1

        first = False
        all_widgets = []

        all_widgets.append(self.plot_interval(interval=interval, bams_dict={}, **kwargs).get_widget(output_format))

        for bam_name, bam_file in bams_dict.items():
            for read_id in read_ids:
                all_widgets.extend(self.plot_read(read_id,
                                                  {bam_name: bam_file},
                                                  interval,
                                                  output_format,
                                                  silence_error=True,
                                                  with_coverage = False,
                                                  with_axis = False,
                                                  with_bed = False,
                                                  add_track_label = False,
                                                  add_reads_label = False,
                                                  **kwargs).children)
        return widgets.VBox(all_widgets)


    def plot_exons_helper(self,
                          feature,
                          merge_exons = True,
                          **kwargs):

        (feature_id, feature_type) = self.get_feature_info(feature)
        exons_list = None

        if feature_type == "exon":
            return self.plot_feature(feature_id, **kwargs)
        elif feature_type == "transcript":
            exons_list = sorted(self.transcript_to_exons[feature_id])
        elif feature_type == "gene":
            if merge_exons:
                tmp_exons = self.gene_to_exons[feature_id].copy()
                tmp_exons.merge_overlaps(data_reducer = interval_data_reduce)
                exons_list = sorted(tmp_exons)
            else:
                exons_list = sorted(self.gene_to_exons[feature_id])
        else:
            raise TypeError("Error, could not find exons for feature\nfeature_id = ", feature_id, " ; feature_type = ", feature_type)

        return exons_list


    def plot_exons_slices(self, 
                          feature,
                          merge_exons = True, 
                          normalize_interval_width = False,
                          padding_perc = 0.1,
                          **kwargs):

        exons_list = self.plot_exons_helper(feature, merge_exons, **kwargs)
        if isinstance(exons_list, itv.Document):
            return exons_list
        # else res is an exons_list

        total_interval_size = 0
        left_bound = math.inf
        right_bound = -math.inf

        smallest_interval_size = math.inf
        for interval in exons_list:
            total_interval_size += interval.end - interval.begin
            smallest_interval_size = min(smallest_interval_size, interval.end - interval.begin)
            left_bound = min(left_bound, interval.begin)
            right_bound = max(right_bound, interval.end)

        padding = None
        if not normalize_interval_width:
            padding = math.ceil(smallest_interval_size * padding_perc)
            padding_perc = padding / (right_bound - left_bound + 1)  # recalculate so that we minimize the full svg size
        else:
            pass # HANDLE, maybe later in the code

        kwargs["vertical_layout"] = True  # necessary to make sure reads/BED entries cannot be on the same line and appear like they are the same because when they are not

        doc = self.plot_feature(feature = feature, padding_perc = padding_perc, **kwargs)
        doc.hide = True  # so that viewbox height and width are set to 0
        base_svg = doc._repr_svg__()
        doc_actual_height = doc.height
      

        per_base_size = None
        reserved_width = (len(exons_list) - 1) * doc.elements[0].space_between + doc.margin_x * 2
        if not normalize_interval_width:
            total_interval_size = total_interval_size + (padding * len(exons_list))
            per_base_size = (doc.width - reserved_width)/total_interval_size


        padding_doc = math.ceil((doc.elements[0].views[0].scale.end - doc.elements[0].views[0].scale.start) * padding_perc)
        per_base_size_doc = doc.width / (doc.elements[0].views[0].scale.end - doc.elements[0].views[0].scale.start + 2*padding_doc)

        left_bound = doc.elements[0].views[0].scale.start - padding_doc  # padded start position the view
        # right_bound = doc.elements[0].views[0].scale.end + padding_doc

        full_svg = base_svg + '<div class="custom-svg-container" style="overflow-y: hidden; width: 100%;">\n'

        for interval in exons_list:
            start = interval.begin
            end = interval.end
            # chrom = interval.chrom
            # strand = interval.strand
            
            if normalize_interval_width:
                display_width = (doc.width - reserved_width)/len(exons_list)
                padding = math.ceil((end - start) * padding_perc)
            else:
                display_width = math.floor((end - start + 2*padding) * per_base_size)

            slice_x_start = doc.elements[0].views[0].scale.topixels(start - padding) + doc.margin_x
            slice_x_end = doc.elements[0].views[0].scale.topixels(end + padding) + doc.margin_x

            slice_width = slice_x_end - slice_x_start

            interval_offset = math.floor((left_bound - start + padding) * per_base_size_doc) + doc.margin_x

            full_svg += f"""
            <svg width="{display_width}" height="{doc_actual_height}" viewBox="{slice_x_start} 0 {slice_width} {doc_actual_height}" preserveAspectRatio="none" style="display: inline-block">
                <use href="#{doc.id}" x="0" y="0" width="{slice_x_end}" height="{doc_actual_height}"/>
            </svg>
        """

        full_svg += f"""
        </div>
        <style>
            .custom-svg-container {{
                white-space: nowrap;
                overflow-x: auto;
            }}
            .custom-svg-container svg {{
                max-width: none;
                height: auto;
                display: inline-block;
            }}
        </style>
        """

        return full_svg



    def plot_exons(self, 
                   feature,
                   merge_exons = True,
                   N_per_row = 99999,
                   view_width = 1600,
                   as_widget = False,
                   **kwargs):
        """
        Parameters
        ----------
        feature : str
            the feature of interest
        merge_exons : bool, optional
            control whether to merge overlapping exons from different isoforms or keep them separate. Merging prevents redundence of parts of the view. (Default: True)
        shared_max_coverage : bool, optional
            control if the y-axis of the coverage tracks should be consistant across a row. Alternatively can provide an arbitrary value to use as the max for the y-axis. (Default: True)
        normalize_interval_width : bool, optional
            control if all exons should have the same view width rather than being proportional to the sequence length. (Default: False)
        N_per_row : int, optional
            How many exons to plot per row at most. (Default: 99999, basically all exons on a single row)
        view_width : int, optional
            width of the view in "pixels". (Defautl: 1600)
        padding_perc : int, optional
            Padding around each exon in proportion of the size of the smallest exon. (Default: 0.1)
        as_widget : bool, optional
            Control if the result is returned as a dropdown selection menu ipywidget where each tab is an exon. (Default: False)
        **kwargs: 
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        A new Document where only the windows of the exons of the provided feature are displayed. : :py:class:`integrative_transcriptomics_viewer.Document`
        """

        (feature_id, feature_type) = self.get_feature_info(feature)
        exons_list = None

        if feature_type == "exon":
            return self.plot_feature(feature_id, **kwargs)
        elif feature_type == "transcript":
            exons_list = sorted(self.transcript_to_exons[feature_id])
        elif feature_type == "gene":
            if merge_exons:
                tmp_exons = self.gene_to_exons[feature_id].copy()
                tmp_exons.merge_overlaps(data_reducer = interval_data_reduce)
                exons_list = sorted(tmp_exons)
            else:
                exons_list = sorted(self.gene_to_exons[feature_id])
        else:
            print("Error, could not find exons for feature")
            print("feature_id = ", feature_id, " ; feature_type = ", feature_type)
            return

        if as_widget:
            all_views = []
            all_titles = []
            for exon in exons_list:
                doc = itv.Document(view_width)
                self.make_intervals_row_through_virtual(doc, [exon], **kwargs)
                all_views.append(widgets.HTML(doc._repr_svg__()))
                all_titles.append("Exon:: " + exon.data + " : " + str(exon.begin) + " - " + str(exon.end))

            stack = widgets.Stack(all_views, selected_index=0)
            dropdown = widgets.Dropdown(options=all_titles)
            widgets.jslink((dropdown, 'index'), (stack, 'selected_index'))
            return(widgets.VBox([dropdown, stack]))

        else:
            doc = itv.Document(view_width)
            for i in range(0, len(exons_list), N_per_row):
               self.make_intervals_row_through_virtual(doc, exons_list[i:i+N_per_row], **kwargs)
            return doc



    def plot_splice_junctions(self, 
                              feature,
                              view_width = 1600,
                              as_widget = True, 
                              **kwargs):
        """
        Returns a Document where each known splice junction (according to isoforms) is displayed as a side by side exons view, one per row.

        Parameters
        ----------
        feature : str
            The feature of interest
        shared_max_coverage : bool, optional
            Control if the y-axis of the coverage tracks should be consistant across a row. Alternatively can provide an arbitrary value to use as the max for the y-axis. (Default: True)
        view_width : int, optional
            Width of the view in "pixels". (Defautl: 1600)
        padding_perc : float, optional
            Padding around each exon in proportion of the size of the smallest exon. (Default: 0.1)
        as_widget : bool, optional
            Control if the result is returned as a dropdown selection menu ipywidget where each tab is a splice junction. (Default: True)
        **kwargs: 
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        A new Document where each known splice junction (according to isoforms) is displayed as a side by side exons view, one per row, or an ipywidgets dropdown menu with one splice junction view option. : :py:class:`integrative_transcriptomics_viewer.Document` or :py:class:`ipywidgets.widgets.widget_box.VBox`
        """


        (feature_id, feature_type) = self.get_feature_info(feature)
        
        if feature_type == "exon":
            print("Error, feature type provided is an exon, hence there is no splice junction within it.")
            return

        exons_pairs = []
        if feature_type == "transcript":
            all_exons = sorted(self.transcript_to_exons[feature_id])
            for i in range(1, len(all_exons)):
                exons_pairs.append((all_exons[i-1], all_exons[i]))

        elif feature_type == "gene":
            for transcript_id in self.gene_to_transcripts[feature_id]:
                all_exons = sorted(self.transcript_to_exons[transcript_id])
                for i in range(1, len(all_exons)):
                    exons_pairs.append((all_exons[i-1], all_exons[i]))

        if len(exons_pairs) == 0:
            print("No splice junctions for the requested feature, this probably means it only has 1 exon")
            return

        seen = set()
        seen_add = seen.add
        exons_pairs = [x for x in exons_pairs if not (x in seen or seen_add(x))]

        if as_widget:
            all_views = []
            all_titles = []
            for pair in exons_pairs:
                doc = itv.Document(view_width)
                self.make_intervals_row_through_virtual(doc, pair, **kwargs)
                all_views.append(widgets.HTML(doc._repr_svg__()))
                all_titles.append("Splice junction btw: exon:: " + pair[0].data + ":" + str(pair[0].begin) + "-" + str(pair[0].end) + " and exon::" + pair[1].data + " : " + str(pair[1].begin) + " - " + str(pair[1].end))

            stack = widgets.Stack(all_views, selected_index=0)
            dropdown = widgets.Dropdown(options=all_titles)
            widgets.jslink((dropdown, 'index'), (stack, 'selected_index'))
            return(widgets.VBox([dropdown, stack]))

        else:
            doc = itv.Document(view_width)
            for pair in exons_pairs:
               self.make_intervals_row_through_virtual(doc, pair, **kwargs)
            return doc
     

    def make_intervals_row_through_virtual(self,
                                           doc,
                                           intervals_list,
                                           bams_dict,
                                           padding_perc = 0.1, 
                                           add_track_label = "auto",
                                           add_reads_label = "auto",
                                           add_coverage_label = "auto",
                                           with_coverage = True,
                                           include_secondary = False,
                                           row = None, 
                                           normalize_interval_width = False,
                                           shared_max_coverage = True,
                                           with_bed_label = True,
                                           **kwargs):
        """
        Adds a row to the provided doc for the given intervals such that reads ordering and spacing is consistant across them.

        Parameters
        ----------
        doc : py:class:`integrative_transcriptomics_viewer.Document`
            The Document to add to.
        intervals_list : list
            List of Intervals to make the views for.
        bams_dict : dict
            Dict of (virtual) BAMs to display the reads from.
        normalize_interval_width : bool, optional
            Control if all intervals should have the same view width rather than being proportional to the sequence length. (Default: False)
        shared_max_coverage : bool, optional
            Control if the y-axis of the coverage tracks should be consistant across a row. Alternatively can provide an arbitrary value to use as the max for the y-axis. (Default: True)
        with_bed_label : bool, optional
            Control if a label is displayed for each BED source. (Default: False)
        **kwargs: 
            anything that can be passed to :meth:`build_view_row()`
        """

        if row is None:
            row = itv.ViewRow("row")

        total_interval_size = 0
        left_bound = math.inf
        right_bound = -math.inf

        smallest_interval_size = math.inf
        for interval in intervals_list:
            total_interval_size += interval.end - interval.begin
            smallest_interval_size = min(smallest_interval_size, interval.end - interval.begin)
            left_bound = min(left_bound, interval.begin)
            right_bound = max(right_bound, interval.end)

        reserved_width = (len(intervals_list) - 1) * row.space_between + doc.margin_x * 2
        padding = math.ceil(smallest_interval_size * padding_perc)
        left_bound -= padding
        right_bound += padding
        if not normalize_interval_width:
            total_interval_size = total_interval_size + (padding * len(intervals_list))
            per_base_size = (doc.width - reserved_width)/total_interval_size
            padding_perc = 0

        max_coverage_dict = {}
        virtual_bams_dict = {}
        bam_track_to_series = {}
        for key, value in bams_dict.items():

            all_reads_for_coverage = set()

            bam_refs = None
            opener_fn = get_bam_opener(value)
            with opener_fn(value) as bam:
                bam_refs = bam.references
                virtual_bam = itv.bamtrack.VirtualBAM([], bam_refs)
                virtual_bam.dumb_fetch = True
                # bam_track.quick_consensus = False
                
                for read in bam.fetch(intervals_list[0].chrom, left_bound, right_bound):
                    if not include_secondary and read.is_secondary:
                        continue
                    # add check that it does overlap with intervals_list at least somewhere
                    virtual_bam.reads.append(read)

                virtual_bams_dict[key] = virtual_bam

            tmp_view = itv.GenomeView(intervals_list[0].chrom, left_bound, right_bound, intervals_list[0].strand)
            coverage_track_series = itv.BAMCoverageTrack(value, opener_fn=opener_fn)
            if "priming_orientation" in kwargs:
                coverage_track_series.priming_orientation = kwargs["priming_orientation"]
            if "coverage_bin_size" in kwargs:
                coverage_track_series.bin_size = kwargs["coverage_bin_size"]
            if "coverage_peak_min_distance" in kwargs:
                coverage_track_series.min_dist = kwargs["coverage_peak_min_distance"]
            if "coverage_tag" in kwargs:
                coverage_track_series.tag = kwargs["coverage_tag"]
            if "coverage_tag_fn" in kwargs:
                coverage_track_series.tag_fn = kwargs["coverage_tag_fn"]
            if "coverage_by_strand" in kwargs:
                coverage_track_series.stranded_coverage = kwargs["coverage_by_strand"]
            coverage_track_series.layout(tmp_view.scale)
            bam_track_to_series[virtual_bam] = coverage_track_series.series
            max_coverage_dict[key] = coverage_track_series.max_y


        bed_config = self  # default, if not with_bed_labels, use the original since not modifying it
        new_beds = {}
        new_bed_labels = {}
        secondary_new_bed_labels = {}
        if with_bed_label:
            bed_config = self.shallow_copy()
            seen_bed_entries = set()
            if type(self.bed_annotation) is dict:
                for bed_name, bed_path in self.bed_annotation.items():
                    is_not_first = 0
                    # seen_bed_entries = set()
                    for interval in intervals_list:
                        for bed_entry in itv.bedtrack.bed_fetch(bed_path, interval.chrom, interval.begin, interval.end):
                            if bed_entry.name in seen_bed_entries:
                                continue
                            label = None
                            if is_not_first:
                                new_key = "__tmp_" + str(is_not_first) + bed_name
                            else:
                                new_key = bed_name
                                label = bed_name
                            is_not_first += 1
                            new_beds[new_key] = itv.VirtualBEDTrack(transcripts=[bed_entry], name=None)
                            new_bed_labels[new_key] = label
                            secondary_new_bed_labels[new_key] = "" if label else None
                            seen_bed_entries.add(bed_entry.name)
                bed_config.update_bed(new_beds)
            elif type(self.bed_annotation) is list:
                new_bed_labels = False
                secondary_new_bed_labels = False
                i = 0
                for bed_path in self.bed_annotation:
                    # seen_bed_entries = set()
                    for interval in intervals_list:
                        for bed_entry in itv.bedtrack.bed_fetch(bed_path, interval.chrom, interval.begin, interval.end):
                            if bed_entry.name in seen_bed_entries:
                                continue
                            new_key = "__tmp_" + str(i)
                            new_beds[new_key] = itv.VirtualBEDTrack(transcripts=[bed_entry], name=None)
                            seen_bed_entries.add(bed_entry.name)
                            i += 1
                bed_config.update_bed(new_beds)
            else:  # just a string
                new_bed_labels = False
                secondary_new_bed_labels = False
                i = 0
                
                for interval in intervals_list:
                    for bed_entry in itv.bedtrack.bed_fetch(self.bed_annotation, interval.chrom, left_bound, right_bound):
                        if bed_entry.name in seen_bed_entries:
                            continue
                        seen_bed_entries.add(bed_entry.name)
                        new_key = "__tmp_" + str(i)
                        new_beds[new_key] = itv.VirtualBEDTrack(transcripts=[bed_entry], name=None)
                        i += 1
                bed_config.update_bed(new_beds)
        else:
            new_bed_labels = False


        if bams_dict and shared_max_coverage:
            if type(shared_max_coverage) is bool:  # is True and actually a Boolean rather than an Integer
                max_coverage_dict = max(max_coverage_dict.values())
            else:
                max_coverage_dict = shared_max_coverage

        for interval in intervals_list:
            start = interval.begin
            end = interval.end
            chrom = interval.chrom
            strand = interval.strand
            
            if normalize_interval_width:
                # padding = math.ceil((end - start) * padding_perc)
                interval_width = (doc.width - reserved_width)/len(intervals_list)
            else:
                interval_width = math.floor((interval.end - interval.begin + padding) * per_base_size)
                start = start - padding
                end = end + padding

            row = bed_config.build_view_row(start = start, 
                                                 end = end,
                                                 chrom = chrom,
                                                 strand = strand,
                                                 bams_dict = virtual_bams_dict, 
                                                 vertical_layout_reads = True,
                                                 padding_perc = padding_perc, 
                                                 add_track_label = add_track_label,
                                                 add_reads_label = add_reads_label,
                                                 add_coverage_label = add_coverage_label,
                                                 with_bed_label = new_bed_labels,
                                                 row = row, 
                                                 view_width = interval_width, 
                                                 view_margin_y = 0,
                                                 coverage_track_max_y = max_coverage_dict,
                                                 **kwargs)


            if add_track_label:
                add_track_label = "\n"
            new_bed_labels = secondary_new_bed_labels
            add_reads_label = None
            add_coverage_label = None

            for track in row.views[-1].get_tracks():
                if isinstance(track, itv.bamtrack.BAMCoverageTrack):
                    # track.series = bam_path_to_series[track.bam_path]
                    track.series = bam_track_to_series[track.bam_path]
                    track.cached_series = True

        doc.elements.append(row)
        return doc


    def get_gene_tab_title(self, feature):
        """
        Parameters
        ----------
        feature : str
            Name or id of the feature.

        Returns
        -------
        The name of the gene the feature is part of, or returns back feature. : str
        """

        if isinstance(feature, Interval):
            return feature.data + ":" + str(feature.begin) + "-" + str(feature.end)

        (feature_id, feature_type) = self.get_feature_info(feature)
        gene_name = self.get_gene_name((feature_id, feature_type))
        if gene_name != feature:
            feature = gene_name  # + "_" + feature
        return feature
        

    # custom bed dict accepts a dict with the same keys are bams_dict only
    def organize_tab_section(self,
                             bams_dict,
                             interval,
                             tab_name,
                             tab_id,
                             custom_bed_dict = None,
                             with_bed = True,
                             with_coverage = True,
                             fill_coverage = True,
                             coverage_bin_size = 0,
                             with_reads = True,
                             expended = False,
                             **kwargs):
        """
        Generate views and subtabs for a given tab and organize them so they can be used with a template.

        Parameters
        ----------
        bams_dict : dict
            Dict of (virtual) BAMs to display the reads from.
        interval : :class:`intervaltree.Interval`
            Interval to draw.
        tab_name : str
            Unique name of this tab.
        custom_bed_dict : dict, optional
            Dict of (virtual) BED entries to use instead of the internally reference BED for the shared annotations of this tab. For keys that are the same as the bams_dict, the specific BED entries will be used for the matching BAM subtab.
        expended : bool, optional
            Control if the subtabs should be expended by default within the tab. (Default: False)
        **kwargs:
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        A dict that stores the tab contents so it can be used with a template. : dict
        """

        bed_config = self.shallow_copy()

        if custom_bed_dict is not None:
            all_bed_entries = itv.VirtualBEDTrack()
            for virtual_bed in custom_bed_dict.values():
                for transcript in virtual_bed.transcripts:
                    if transcript not in all_bed_entries.transcripts:
                        all_bed_entries.transcripts.append(transcript)
            bed_config.update_bed(all_bed_entries)
        else:
            bed_config.update_bed(self.bed_annotation)

        shared_static_svg = bed_config.plot_interval(bams_dict={},
                                                    interval = interval,
                                                    with_bed = with_bed,
                                                    with_reads = False,
                                                    with_coverage = False,
                                                    add_track_label = False,
                                                    **kwargs
                                                    )._repr_svg__()
        
        tab_sections = []
        for key, bam in bams_dict.items():
            if bam is None:
                continue
            # unique_id = f"{tab_id}_{key}"

            static_svg = ""
            resizable_svg = ""
            if custom_bed_dict is not None and key in custom_bed_dict:
                bed_config.update_bed(custom_bed_dict[key])
                static_svg += bed_config.plot_interval(bams_dict={},
                                                       interval = interval,
                                                       with_bed = with_bed,
                                                       with_reads = False,
                                                       with_coverage = False,
                                                       add_track_label = False,
                                                       **kwargs
                                                      )._repr_svg__() + "</br>"

            if with_coverage:
                static_svg += self.plot_interval(bams_dict = {key: bam},
                                                 interval = interval,
                                                 with_reads = False,
                                                 with_coverage = True,
                                                 with_axis = len(bams_dict) > 1,
                                                 with_bed = False,
                                                 add_track_label = False,
                                                 fill_coverage = fill_coverage,
                                                 coverage_bin_size = coverage_bin_size,
                                                 **kwargs
                                                )._repr_svg__() + "</br>"

            if with_reads:
                resizable_svg += self.plot_interval(bams_dict = {key: bam},
                                                    interval = interval,
                                                    with_reads = with_reads,
                                                    with_coverage = False,
                                                    with_axis = False,
                                                    with_bed = False,
                                                    add_track_label = False,
                                                    add_reads_label = False,
                                                    vertical_layout_reads = True,
                                                    **kwargs
                                                )._repr_svg__() + "</br>"

            tab_sections.append({
                'unique_id': f"{tab_id}_{key}",
                'name': f"{tab_name}_{key}",
                'static_svg': static_svg,
                'resizable_svg': resizable_svg,
                'expended': expended
            })

        return {'tab_name': tab_name,
                'tab_id': tab_id,
                'shared_static_svg': shared_static_svg,
                'tab_sections': tab_sections}



    def organize_tabs_by_feature(self,
                                 bams_dict_dict,
                                 custom_bed_dict_dict = None,
                                 tab_title_fn = None,  # function that takes a feature_name/id or tuple(feature_id, feature_type) as input
                                 **kwargs):
        """
        Generate views and subtabs for all tabs and organize them so they can be used with a template.

        Parameters
        ----------
        bams_dict_dict : dict of dict
            Dict of Dict of (virtual) BAMs to display the reads from. Primary keys are the feature keys, secondary keys are the BAM names.
        interval : :class:`intervaltree.Interval`
            Interval to draw.
        custom_bed_dict_dict : dict of dict, optional
            Dict of Dict of (virtual) BED entries to use instead of the internally reference BED for the shared annotations for each tab and each BAM. For subkeys that are the same as the bams_dict_dict subkeys, the specific BED entries will be used for the matching BAM subtab.
        tab_title_fn : typing.Callable[[str], str], optional
            Method to use to generate unique tab names when provided a feature_id. (Default: :meth:`get_gene_tab_title()`)
        expended : bool, optional
            Control if the subtabs should be expended by default within the tab. (Default: False)
        **kwargs:
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        A dict that stores all the tab contents so it can be used with a template. : dict
        """

        if tab_title_fn is None:  # workaround because can't set a self.method as default parameter
            tab_title_fn = self.get_gene_tab_title

        tabs = []
        for feature, bams_dict in bams_dict_dict.items():
            #feature_name = str(feature)
            # print(feature)
            # (feature_id, feature_type) = self.get_feature_info(feature_name)
            # interval = self.id_to_coordinates[feature_id]
            interval = self.get_interval_from_feature(feature)

            custom_bed_dict = None
            if custom_bed_dict_dict is not None and feature in custom_bed_dict_dict:
                custom_bed_dict = custom_bed_dict_dict[feature]

            tab_name = tab_title_fn(feature)
            tab_id = tab_name + "_" + str(time.time())

            tabs.append(self.organize_tab_section(bams_dict = bams_dict,
                                                  interval = interval, 
                                                  tab_name = tab_name,
                                                  tab_id = tab_id,
                                                  custom_bed_dict = custom_bed_dict,
                                                  **kwargs))

        return tabs


    def organize_tabs_by_classification(self,
                                        bams_dict_dict,
                                        interval,
                                        custom_bed_dict_dict = None,
                                        tab_title_fn = None,  # function that takes a feature_name/id or tuple(feature_id, feature_type) as input
                                        **kwargs):
        """
        Generate views and subtabs for all tabs and organize them so they can be used with a template.

        Parameters
        ----------
        bams_dict_dict : dict of dict
            Dict of Dict of (virtual) BAMs to display the reads from. Primary keys are the feature keys, secondary keys are the BAM names.
        interval : :class:`intervaltree.Interval`
            Interval to draw.
        custom_bed_dict_dict : dict of dict, optional
            Dict of Dict of (virtual) BED entries to use instead of the internally reference BED for the shared annotations for each tab and each BAM. For subkeys that are the same as the bams_dict_dict subkeys, the specific BED entries will be used for the matching BAM subtab.
        tab_title_fn : typing.Callable[[str], str], optional
            Method to use to generate unique tab names when provided a classification. (Default: lambda x: x ; returns itself)
        expended : bool, optional
            Control if the subtabs should be expended by default within the tab. (Default: False)
        **kwargs:
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        A dict that stores all the tab contents so it can be used with a template. : dict
        """

        if tab_title_fn is None:  # workaround because can't set a self.method as default parameter
            tab_title_fn = lambda x: x


        tabs = []
        for classification, bams_dict in bams_dict_dict.items():
            custom_bed_dict = None
            if custom_bed_dict_dict is not None and classification in custom_bed_dict_dict:
                custom_bed_dict = custom_bed_dict_dict[classification]

            tab_name = tab_title_fn(classification)
            tab_id = tab_name + "_" + str(time.time())
            tabs.append(self.organize_tab_section(bams_dict = bams_dict,
                                                  interval = interval, 
                                                  tab_name = tab_name,
                                                  tab_id = tab_id,
                                                  custom_bed_dict = custom_bed_dict,
                                                  **kwargs))

        return tabs



    def plot_by_features_as_tab(self,
                                features_list,
                                bams_dict,
                                page_title = "Plot by feature", # and split by barcode whitelist category if provided
                                cellbarcode_whitelist = None,
                                cellbarcode_from = itv.StandardCellBarcode(),
                                **kwargs):
        """
        Returns an HTML object to output or display, with views separated over tabs by the list of features provided. Within tabs, split in subtabs over BAMs.
        This is mostly intended to use with non overlapping features since there is no separation by isoforms within a gene for example.
        Supports providing a cell barcodes whitelist as a list or dict of lists to only keep the indexed barcodes, and split barcodes according to dict keys into different groupings.

        Parameters
        ----------
        features_list : list
            List of features to generate views for, one per tab.
        bams_dict : dict
            Dict of (virtual) BAMs to display the reads from.
        page_title : str
            HTML page title
        cellbarcode_whitelist : dict of lists, or list, optional
            Dict of lists or list with whitelisted cell barcodes. If a dict, reads will be split according to which key they are listed under. If a list, they will be considered under an "all" key. (Default: None)
        cellbarcode_from : :py:class:`integrative_transcriptomics_viewer.CellBarcode`
            An instance of an implemented :py:class:`integrative_transcriptomics_viewer.CellBarcode` class with the :meth:integrative_transcriptomics_viewer.CellBarcode.get_barcode(self, read) method where read is a pysam.AlignmentSegment object and returns the cell barcode of the given read. (Default: :py:class:`integrative_transcriptomics_viewer.StandardCellBarcode()`)
        **kwargs:
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        An HTML document : str
        """

        virtual_bams_dict = {}
        for feature in features_list:
            # feature_name = str(feature)
            virtual_bams_dict[feature] = {}

            if cellbarcode_from is not None and cellbarcode_whitelist is not None:

                interval = self.get_interval_from_feature(feature)
                for bam_name, bam_file in bams_dict.items():

                    #for bam_name, bam_file in bams_dict.items():
                    virtual_bams_dict[feature].update(itv.split_bam_by_cellbarcode_whitelist(bam_name,
                                                                                                                            bam_file,
                                                                                                                            interval,
                                                                                                                            cellbarcode_whitelist = cellbarcode_whitelist,
                                                                                                                            cellbarcode_from = cellbarcode_from))
            else:
                for bam_name, bam_file in bams_dict.items():
                    virtual_bams_dict[feature][bam_name] = bam_file


        tabs = self.organize_tabs_by_feature(virtual_bams_dict, **kwargs)

        return itv.templates.render_tab_titles(tabs, page_title)


    # classified_dict has the classification as keys
    def match_classification_to_bed_entries(self, classifications, interval, annotation_matching):
        """
        Returns a virtual BED dict where keys are the classifications and entries are a list of all annotations that have been found to match according to the annotation_matching.match check.
        Classifications "ambiguous" and "unclassifed" are special classifications that are not checked.

        Parameters
        ----------
        classifications : list
            A list of the classifications to find the annotations for, "ambiguous" and "unclassified" are not checked against and will not exist as keys in the returned dict.
        interval : :class:`intervaltree.Interval`
            Interval to fetch the internal BED entries within to then try to match.
        annotation_matching : py:class:`integrative_transcriptomics_viewer.AnnotationMatching`
            An instance of an implemented AnnotationMatching class with the match(self, query_annotation, target_annotation) method where the query_annotation will be an entry of classifications and target_annotation is the name of a BED entry.
        
        Returns
        -------
        A dict of virtual BEDs where keys are the classifications and entries are a list of all annotations that have been found to match according to the annotation_matching.match check. : dict
        """

        # parse known annotations
        all_known_annotations = self.get_bed_entries(interval)

        virtual_bed_dict = {}
        for known_annotation, virtual_bed in all_known_annotations.items():
            for classification in classifications:
                if "ambiguous" in classification or "unclassified" in classification:
                    continue
                else:
                    if annotation_matching.match(classification, known_annotation):
                        if classification not in virtual_bed_dict:
                            virtual_bed_dict[classification] = virtual_bed.copy()  # copy so that the appends don't affect other entries that started from this same bed
                        else:
                            for transcript in virtual_bed.transcripts:
                                if transcript not in virtual_bed_dict[classification].transcripts:
                                    virtual_bed_dict[classification].transcripts.append(transcript)

        return virtual_bed_dict


    # returns two dict with same keys, one with dicts of split BAMs, and one with associated BED entries
    # here whitelist is just that, a whitelist, there is no split based on which dict key a barcode is found in
    def split_bams_dict_by_classification(self,
                                          bams_dict,
                                          feature,
                                          classification_from,
                                          annotation_matching,
                                          **kwargs):
        """
        Parameters
        ----------
        bams_dict : dict
            Dict of (virtual) BAMs to display the reads from.
        gene : str
            The gene id or name to use the coordinates of for the window to fetch reads over.
        classification_from : :py:class:`integrative_transcriptomics_viewer.Classification`
            An instance of an implemented py:class:`integrative_transcriptomics_viewer.Classification` class with the :meth:integrative_transcriptomics_viewer.Classification.get_classification(self, read, gene_id) method where read is a pysam.AlignmentSegment object, gene_id the id of the provided gene, and which returns the classification of the given read.
        annotation_matching : :py:class:`integrative_transcriptomics_viewer.AnnotationMatching`
            An instance of an implemented py:class:`integrative_transcriptomics_viewer.AnnotationMatching` class with the :meth:`integrative_transcriptomics_viewer.AnnotationMatching.match(self, query_annotation, target_annotation)` method where the query_annotation will be an entry of classifications and target_annotation is the name of a BED entry.
        cellbarcode_whitelist : dict of list, or list, optional
            A dict of lists or list with whitelisted cell barcodes. Here the whitelist is not used to split according to keys in the dict, it is purely a whitelist. If this is provided, a cellbarcode_from needs to also be provided. (Default: None)
        cellbarcode_from : :py:class:`integrative_transcriptomics_viewer.CellBarcode`, optional
            An instance of an implemented :py:class:`integrative_transcriptomics_viewer.CellBarcode` class with the :meth:`integrative_transcriptomics_viewer.CellBarcode.get_barcode(self, read)` method where read is a :py:class:`pysam.AlignmentSegment` object and returns the cell barcode of the given read. (Default: :py:class:`integrative_transcriptomics_viewer.StandardCellBarcode()`)
        **kwargs:
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        A tuple with a dict of virtual BAMs and a dict of virtual BEDs with matching keys, where provided input (virtual) BAMs are split according to classification. : tuple(dict, dict)
        """
        
        (feature_id, feature_type) = self.get_feature_info(feature)
        # if feature_type != "gene":
        #     print("feature provided is not a gene")
        #     return -1

        interval = self.get_interval_from_feature((feature_id, feature_type))

        virtual_bams_dict = {}
        for bam_name, bam_file in bams_dict.items():
            virtual_bams_dict.update(itv.split_bam_by_classification(bam_file = bam_file,
                                                                                                    name_prefix = bam_name,
                                                                                                    feature_id = feature_id,
                                                                                                    interval = interval,
                                                                                                    classification_from = classification_from,
                                                                                                    **kwargs))

        # parse known annotations
        virtual_bed_dict = self.match_classification_to_bed_entries(virtual_bams_dict.keys(), interval, annotation_matching)

        return (virtual_bams_dict, virtual_bed_dict)
        


    def plot_by_classification_over_features(self,
                                             bams_dict,
                                             features_list,
                                             classification_from,
                                             annotation_matching,
                                             page_title = "Split by Classification, Plot by feature",
                                             **kwargs):
        """
        Returns an HTML object to output or display, with views separated over tabs by the list of features provided. Within tabs, split in subtabs by classifications and BAMs.
        This is mostly intended to use with non overlapping features since there is no separation by isoforms within a gene for example.
        Supports providing a cell barcodes whitelist as a list or dict of lists to only keep the indexed barcodes, and split barcodes according to dict keys into different groupings.

        Parameters
        ----------
        bams_dict : dict
            A dict of (virtual) BAMs to display the reads from.
        features_list : list
            A list of features to generate views for, one per tab.
        classification_from : :py:class:`integrative_transcriptomics_viewer.Classification`
            An instance of an implemented py:class:`integrative_transcriptomics_viewer.Classification` class with the :meth:integrative_transcriptomics_viewer.Classification.get_classification(self, read, gene_id) method where read is a pysam.AlignmentSegment object, gene_id the id of the provided gene, and which returns the classification of the given read.
        annotation_matching : py:class:`integrative_transcriptomics_viewer.AnnotationMatching`
            An instance of an implemented AnnotationMatching class with the :meth:`integrative_transcriptomics_viewer.AnnotationMatching.match(self, query_annotation, target_annotation)` method where the query_annotation will be an entry of classifications and target_annotation is the name of a BED entry.
        page_title : str
            HTML page title
        **kwargs:
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        An HTML document : str
        """

        virtual_bams_dict_dict = {}
        virtual_beds_dict_dict = {}
        for feature in features_list:
            (virtual_bams_dict_dict[feature], virtual_beds_dict_dict[feature]) = \
                            self.split_bams_dict_by_classification(bams_dict = bams_dict,
                                                                   feature = feature,
                                                                   classification_from = classification_from,
                                                                   annotation_matching = annotation_matching,
                                                                   **kwargs
                                                                   )

        tabs = self.organize_tabs_by_feature(bams_dict_dict = virtual_bams_dict_dict,
                                             custom_bed_dict_dict = virtual_beds_dict_dict,
                                             **kwargs)

        return itv.templates.render_tab_titles(tabs, page_title)



    def plot_by_classification_as_tabs(self,
                                       bams_dict,
                                       feature,
                                       classification_from,
                                       annotation_matching,
                                       page_title = "Split by Classification",
                                       add_all_tab = True,
                                       **kwargs):
        """
        Returns an HTML object to output or display, with views separated over tabs by the list of features provided. Within tabs, split in subtabs by BAMs.
        This is mostly intended to use with non overlapping features since there is no separation by isoforms within a gene for example.
        Supports providing a cell barcodes whitelist as a list or dict of lists to only keep the indexed barcodes, and split barcodes according to dict keys into different groupings.

        Parameters
        ----------
        bams_dict : dict
            A dict of (virtual) BAMs to display the reads from.
        gene : str
            The gene id or name to use the coordinates of for the window to fetch reads over.
        classification_from : :py:class:`integrative_transcriptomics_viewer.Classification`
            An instance of an implemented py:class:`integrative_transcriptomics_viewer.Classification` class with the :meth:integrative_transcriptomics_viewer.Classification.get_classification(self, read, gene_id) method where read is a pysam.AlignmentSegment object, gene_id the id of the provided gene, and which returns the classification of the given read.
        annotation_matching : py:class:`integrative_transcriptomics_viewer.AnnotationMatching`
            An instance of an implemented py:class:`integrative_transcriptomics_viewer.AnnotationMatching` class with the :meth:`integrative_transcriptomics_viewer.AnnotationMatching.match(self, query_annotation, target_annotation)` method where the query_annotation will be an entry of classifications and target_annotation is the name of a BED entry.
        page_title : str
            HTML page title
        **kwargs:
            anything that can be passed to :meth:`build_view_row()`

        Returns
        -------
        An HTML document : str
        """

        (feature_id, feature_type) = self.get_feature_info(feature)
        # if feature_type != "gene":
        #     print("feature provided is not a gene")
        #     return -1

        interval = self.get_interval_from_feature((feature_id, feature_type))

        virtual_bams_dict_dict = {}
        if add_all_tab:
            virtual_bams_dict_dict["all"] = {} 
        for bam_name, bam_file in bams_dict.items():
            for classification, virtual_bam in (itv.split_bam_by_classification(bam_file = bam_file,
                                                                                                               name_prefix = "",
                                                                                                               feature_id = feature_id,
                                                                                                               interval = interval,
                                                                                                               classification_from = classification_from,
                                                                                                               **kwargs)).items():
                if classification not in virtual_bams_dict_dict:
                    virtual_bams_dict_dict[classification] = {}
                virtual_bams_dict_dict[classification][bam_name] = virtual_bam     

            if add_all_tab:
                virtual_bams_dict_dict["all"][bam_name] = bam_file 

        # parse known annotations
        virtual_bed_dict = self.match_classification_to_bed_entries(virtual_bams_dict_dict.keys(), interval, annotation_matching)

        virtual_beds_dict_dict = {}
        for classification, virtual_bed in virtual_bed_dict.items():
            if classification in virtual_bams_dict_dict:
                if classification not in virtual_beds_dict_dict:
                    virtual_beds_dict_dict[classification] = {}
                virtual_beds_dict_dict[classification]['all'] = virtual_bed


        tabs = self.organize_tabs_by_classification(bams_dict_dict = virtual_bams_dict_dict,
                                                    interval = interval,
                                                    custom_bed_dict_dict = virtual_beds_dict_dict,
                                                    **kwargs)

        return itv.templates.render_tab_titles(tabs, page_title)



























