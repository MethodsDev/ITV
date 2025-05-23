import pysam
import os
import gzip
import bz2
import math
import random
from collections.abc import MutableMapping, Iterable, Iterator
from itertools import islice


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

        
def reservoir_sampling_from_iterable(iterable_to_sample, sample_size, iterable_length):
    if iterable_length > sample_size:
        iterable_iter = iterable_to_sample.__iter__()
        
        # initial fill of the reservoir
        reservoir = [(item, i) for i, item in enumerate(islice(iterable_iter, sample_size))]
        latest_index = sample_size

        i = sample_size
        n = iterable_length - 1
        W = math.exp(math.log(random.random()) / sample_size)
        while i < n:
            # jump to the next element that will replace another in the reservoir
            i += math.floor(math.log(random.random()) / math.log(1 - W)) + 1

            # if we didn't reach the end of the list of stuff to sample yet
            if i < n:
                # reservoir[random.randint(0, sample_size - 1)] = iterable_to_sample[i]  # random index between 1 and k, inclusive
                for _ in range(i - 1 - latest_index):
                    next(iterable_iter)
                reservoir[random.randint(0, sample_size - 1)] = (next(iterable_iter), i)  # random index between 1 and k, inclusive
                latest_index = i
                W = W * math.exp(math.log(random.random()) / sample_size)
        
        # sorting by original index to keep input sorting
        reservoir = [item for item, _ in sorted(reservoir, key=lambda x: x[1])]
        return reservoir
    else:
        return iterable_to_sample


def reservoir_sampling_from_iterator(iterator_to_sample, sample_size, iterator_length):
    # sampling_pool = [_ for _ in iterable_to_sample]

    if len(iterator_length) > sample_size:
        # initial fill of the reservoir
        reservoir = [(iterator_to_sample.next(), i) for i in range(sample_size)]
        latest_index = sample_size

        i = sample_size
        n = iterator_length - 1
        W = math.exp(math.log(random.random()) / sample_size)
        while i < n:
            # jump to the next element that will replace another in the reservoir
            i += math.floor(math.log(random.random()) / math.log(1 - W)) + 1

            # if we didn't reach the end of the list of stuff to sample yet
            if i < n:
                for j in range(latest_index, i-1):
                    iterator_to_sample.next()
                reservoir[random.randint(0, sample_size - 1)] = (iterator_to_sample.next(), i)  # random index between 1 and k, inclusive
                latest_index = i
                W = W * math.exp(math.log(random.random()) / sample_size)
                
        # sorting by original index to keep input sorting
        reservoir = [item for item, _ in sorted(reservoir, key=lambda x: x[1])]
        return reservoir
    else:
        return iterator_to_sample

    
def reservoir_sampling(iter_to_sample, sample_size, iter_len):
    # could add a check of whether it's a generator or not? to only run this first line if it is
    if isinstance(iter_to_sample, Iterator):
        return reservoir_sampling_from_iterator(iter_to_sample, sample_size, iter_len)
    else:
        return reservoir_sampling_from_iterable(iter_to_sample, sample_size, iter_len)
