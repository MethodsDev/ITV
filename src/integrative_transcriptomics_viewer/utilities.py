import pysam
import os
import gzip
import bz2
from collections.abc import MutableMapping


def match_chrom_format(chrom, keys):
    if chrom in keys:
        return chrom
    if "chr" in chrom:
        chrom2 = chrom.replace("chr", "")
    else:
        chrom2 = "chr{}".format(chrom)
        
    if chrom2 in keys:
        return chrom2
    return chrom


def get_one_track(doc_or_view, name):
    """
    Convenience function to get a single track by name from a document 
    or a view. If more than one track is found matching the provided 
    track name, then the first one is returned. Raises IndexError 
    if no matching tracks are found.
    """
    tracks = doc_or_view.get_tracks(name)
    return tracks[0]
    

def is_paired_end(bam_path, n=100):
    bam = pysam.AlignmentFile(bam_path)

    for i, read in enumerate(bam.fetch()):
        if read.is_paired:
            return True
        if i >= n:
            break

    return False


def is_long_frag_dataset(bam_path, n=1000):
    bam = pysam.AlignmentFile(bam_path)

    for i, read in enumerate(bam.fetch()):
        if read.is_paired:
            return False

        if read.query_length > 1000:
            return True

        if i > n:
            break

    return False



def flatten(dictionary, separator='_'):
    items = []
    for key, value in dictionary.items():
        if isinstance(value, MutableMapping):
            for el in flatten(value, separator=separator):
                items.append(key + separator + el)
        elif isinstance(value, list):
            for el in value:
                items.append(key + separator + el)
        else:
            items.append(key + separator + value)
    return items


def my_hook_compressed(filename, mode):
    if 'b' not in mode:
        mode += 't'
    ext = os.path.splitext(filename)[1]
    if ext == '.gz':
        return gzip.open(filename, mode)
    elif ext == '.bz2':
        return bz2.open(filename, mode)
    else:
        return open(filename, mode)

