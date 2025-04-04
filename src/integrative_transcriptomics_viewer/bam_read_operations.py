from integrative_transcriptomics_viewer.bamtrack import VirtualBAM, CombinedVirtualBAM
from intervaltree import Interval

import pysam

def get_bam_opener(bam):
    if isinstance(bam, VirtualBAM):
        opener_fn = lambda x: x
    elif isinstance(bam, CombinedVirtualBAM):
        opener_fn = lambda x: x
    else:
        opener_fn = pysam.AlignmentFile
    return opener_fn


def get_virtualbam_max_coverage(coverage_bam):
    all_starts = []
    all_ends = []
    for inter in sorted(coverage_bam.reads_interval_tree):
        all_starts.append(inter.begin)
        all_ends.append(inter.end)
    all_ends.sort()
    current_max = current_coverage = i = j = 0
    while i < len(all_starts) or j < len(all_ends):
        if i < len(all_starts) and all_starts[i] < all_ends[j]:
            current_coverage += 1
            i += 1
            current_max = max(current_max, current_coverage)
        elif i < len(all_starts):
            current_coverage -= 1
            j += 1
        elif j < len(all_ends):
            current_coverage -= 1
            j += 1
    return current_max


def get_regions_by_read_id(bam_file, read_id, opener_fn=pysam.AlignmentFile):
    regions = []

    with opener_fn(bam) as bam_in:
        for read in bam_in.fetch():
            if read.query_name != read_id:
                continue

            regions.append(Interval(read.reference_start, read.reference_end, bam_in.get_reference_name(read.reference_id) + "+" if read.is_forward else "-"))

    return(regions)


def find_read_in_bam(read_id, bams_dict, silence_error=False):
    regions = []
    virtual_bams = []

    for key, value in bams_dict.items():
        if isinstance(value, VirtualBAM):
            opener_fn = lambda x: x
        else:
            opener_fn = pysam.AlignmentFile

        with opener_fn(value) as bam_in:
            bam_refs = bam_in.references
            for read in bam_in.fetch():
                if read.query_name != read_id:
                    continue
                regions.append(Interval(read.reference_start, read.reference_end, read.reference_name + ("+" if read.is_forward else "-")))
                virtual_bams.append(VirtualBAM([read], bam_refs))
    
    if len(regions) == 0 and not silence_error:
        print("Error: read either not found or not aligned")
        return (-1, -1)

    return (regions, virtual_bams)


# similar to regular pysam.AlignmentSegment.get_tag() but instead of raising a KeyError when tag is missing, just return None
def get_read_tag(read, tag):
    if not read.has_tag(tag):
        return None

    return read.get_tag(tag)


def is_in_whitelist(ident, whitelist):
    if isinstance(whitelist, dict):
        for key, whitelist in whitelist.items():
            if ident in whitelist:
                return key
    elif isinstance(whitelist, list) or isinstance(whitelist, set):
        if ident in whitelist:
            return "all"
    return None


# whitelist_dict in the form of {condition: whitelisted_barcodes} or a plain list/set that gives the name "all"
def split_bam_by_cellbarcode_whitelist(bam_name,
                                       bam_file,
                                       interval,
                                       cellbarcode_from,
                                       cellbarcode_whitelist=None):
    if cellbarcode_whitelist is None:
        return {bam_name + "_all": bam_file}

    opener_fn = get_bam_opener(bam_file)

    tmp_reads = {}
    with opener_fn(bam_file) as bam:
        refs = bam.references
        for read in bam.fetch(interval.chrom, interval.begin, interval.end):
            if read.is_unmapped or read.is_secondary:
                continue

            cell_barcode = cellbarcode_from.get_barcode(read)
            if not cell_barcode:
                continue

            key = is_in_whitelist(cell_barcode, cellbarcode_whitelist)
            if key:
                if key not in tmp_reads:
                    tmp_reads[key] = []
                tmp_reads[key].append(read)
    
    virtual_bams_dict = {}
    for key, reads in tmp_reads.items():
        virtual_bams_dict[bam_name + "_" + key] = VirtualBAM(reads, refs)

    return virtual_bams_dict


def split_bam_by_classification(bam_file,
                                name_prefix,
                                feature_id,
                                interval,
                                classification_from,
                                cellbarcode_from = None,
                                cellbarcode_whitelist = None,
                                **kwargs):
    if classification_from is None:
        print("No way of getting classification provided")
        return -1

    opener_fn = get_bam_opener(bam_file)

    if name_prefix and name_prefix != "":
        name_prefix += "_"

    tmp_reads = {}
    with opener_fn(bam_file) as bam:
        refs = bam.references
        for read in bam.fetch(interval.chrom, interval.begin, interval.end):
            if read.is_unmapped or read.is_secondary:
                continue

            ## maybe still keep key name of barcode to use in name
            whitelist = ""
            if cellbarcode_whitelist is not None:
                cell_barcode = cellbarcode_from.get_barcode(read)
                if not cell_barcode:
                    continue
                whitelist = is_in_whitelist(cell_barcode, cellbarcode_whitelist)
                if whitelist is None:
                    continue
                whitelist += "_"

            classifications = classification_from.get_classification(read, feature_id)
            if classifications is None:
                if name_prefix + whitelist + "unclassified" not in tmp_reads:
                    tmp_reads[name_prefix + whitelist + "unclassified"] = []
                tmp_reads[name_prefix + whitelist + "unclassified"].append(read)

            else:
                for classification in classifications:
                    if name_prefix + whitelist + classification not in tmp_reads:
                        tmp_reads[name_prefix + whitelist + classification] = []
                    tmp_reads[name_prefix + whitelist + classification].append(read)

    
    virtual_bams_dict = {}
    for key, reads in tmp_reads.items():
        virtual_bams_dict[key] = VirtualBAM(reads, refs)

    return virtual_bams_dict
