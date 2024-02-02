"""This module contains helper functions for the pipeline."""
import sys

import pysam


def bam_mean_coverage(bam_file: pysam.AlignmentFile, chr_id: str) -> float:
    """Calculate mean coverage for a chromosome.

    For each position in a chromosome, we check how many reads start at this position using mpileup engine.
    We sum all these values and divide by the length of the chromosome.
    """
    coverage = 0
    for pileup_column in bam_file.pileup(chr_id):
        coverage += pileup_column.nsegments
    return coverage / bam_file.get_reference_length(chr_id)


def get_total_size(obj, seen=None):
    """Returns the memory footprint of an object and all of its contents."""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    seen.add(obj_id)
    if isinstance(obj, dict):
        size += sum([get_total_size(v, seen) for v in obj.values()])
        size += sum([get_total_size(k, seen) for k in obj.keys()])
    elif hasattr(obj, '__dict__'):
        size += get_total_size(obj.__dict__, seen)
    elif hasattr(obj, '__iter__') and not isinstance(obj, (str, bytes, bytearray)):
        size += sum([get_total_size(i, seen) for i in obj])
    return size


def get_list_of_chr(bam_file: pysam.AlignmentFile, chr_list: list) -> list[str]:
    """Get list of chromosomes to process.

    If chr_list is 'all' then all chromosomes in bam file are returned.
    Otherwise, only chromosomes from chr_list are returned.
    If any chromosome from chr_list is not found in bam file, ValueError is raised.
    """
    chrs_in_bam = [x['SN'] for x in bam_file.header['SQ']]
    if chr_list != 'all':
        for chr_id in chr_list:
            if chr_id not in chrs_in_bam:
                error_msg = f"Chromosome {chr_id} not found in bam file."
                raise ValueError(error_msg)
        return chr_list
    return chrs_in_bam
