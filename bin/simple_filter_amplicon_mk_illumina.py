#!/usr/bin/env python3

"""Script for down-sampling bam file trying to keep the coverage uniform."""

# ruff: noqa: T201, Q000, TRY0003, FA100

import argparse
from collections import defaultdict
from typing import Optional, Union, List, Tuple

import pysam
from bitarray import bitarray

from nf_pipeline_utils import bam_mean_coverage, get_list_of_chr


def read_pair_generator(bam: pysam.AlignmentFile, region_string: Optional[str] = None) -> (pysam.AlignedSegment, pysam.AlignedSegment):
    """Generate read pairs in a BAM file or within a region string.

    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def organize_reads(bam_file: pysam.AlignmentFile, chr_id: str) -> Tuple[List[Tuple[List[int], List[int]]], List[pysam.AlignedSegment]]:
    """Organize reads in a convenient data structure.

    In result, we got a list of the length of chromosome.
    For each position we have a pair of lists.
    First list contains indexes of reads starting at this position.
    Second list contains indexes of reads ending at this position.

    A second return value is a list of reads.
    """
    genome_length = bam_file.get_reference_length(chr_id)
    organized_reads = [([], []) for _ in range(genome_length)]
    reads = list(bam_file.fetch(region=chr_id))
    for i, read in enumerate(reads):
        organized_reads[read.reference_start][0].append(i)
        organized_reads[read.reference_end - 1][1].append(i)
    return organized_reads, reads


def down_sample_bam(input_bam_path: str,
                    output_bam_path: str,
                    cycles: int = 30,
                    chromosomes: Union[list[str], str] = 'all',
                    mode: str = 'single') -> None:
    """Down-sample bam file."""
    input_bam_file = pysam.AlignmentFile(input_bam_path, "rb")
    output_bam_file = pysam.AlignmentFile(output_bam_path, "wb", template=input_bam_file)
    list_of_chr = get_list_of_chr(input_bam_file, chromosomes)
    for chr_id in list_of_chr:
        print(f"Processing chromosome {chr_id}...")
        if mode == 'single':
            run_mode_single(input_bam_file, chr_id, cycles, output_bam_file)
        elif mode == 'paired':
            raise NotImplementedError("Paired mode not implemented yet.")
            run_mode_paired(input_bam_file, chr_id, cycles, output_bam_file)
        else:
            msg = f"Unknown mode: {mode}"
            raise ValueError(msg)
    output_bam_file.close()
    print(f"Output bam file: {output_bam_path} saved.")
    sorting_and_indexing(output_bam_path)
    for chr_id in list_of_chr:
        print(f"Mean coverage for {chr_id}", bam_mean_coverage(pysam.AlignmentFile(output_bam_path, 'rb'), chr_id))
    print("Done.")


def sorting_and_indexing(bam_path):
    print(f"Sorting and indexing {bam_path} file...")
    pysam.sort("-o", bam_path, bam_path)
    pysam.index(bam_path)


def run_mode_single(bam_file: pysam.AlignmentFile, chr_id: str, cycles: int, output_bam: pysam.AlignmentFile) -> None:
    """Run down-sampling in single reads mode.

    This is two-way mode for unpaired reads.
    Memory O(n) where n is the number of reads.
    Pretty fast.
    """
    number_of_reads = bam_file.count(chr_id)
    genome_length = bam_file.get_reference_length(chr_id)
    used_reads = bitarray(number_of_reads)
    organized_reads, reads = organize_reads(bam_file, chr_id)
    direction = '>'
    for c in range(cycles):
        if direction == '>':
            last_covered = 0
            while last_covered < genome_length:
                for i in organized_reads[last_covered][0]:
                    r = reads[i]
                    if used_reads[i]:
                        continue
                    else:
                        output_bam.write(r)
                        last_covered = r.reference_end - 1
                        used_reads[i] = 1
                        break
                last_covered += 1
        else:
            last_covered = genome_length
            while last_covered > 0:
                for i in organized_reads[last_covered - 1][1]:
                    r = reads[i]
                    if used_reads[i]:
                        continue
                    else:
                        output_bam.write(r)
                        last_covered = r.reference_start
                        used_reads[i] = 1
                        break
                last_covered -= 1
        direction = '>' if direction == '<' else '<'


def run_mode_paired(bam_file: pysam.AlignmentFile, chr_id: str, cycles: int, output_bam: pysam.AlignmentFile) -> None:
    """Run down-sampling in paired mode (Paired reads)."""
    # TODO memory inefficient
    last_covered = -1
    reads = []
    for r1, r2 in read_pair_generator(bam_file, region_string=chr_id):
        reads.append((r1, r2))
    number_of_reads = len(reads)
    direction = '>'
    for c in range(cycles):
        print(f"Cycle {c + 1} {direction}")
        if direction == '>':
            for i in range(number_of_reads):
                if reads[i] is None:
                    continue
                r1, r2 = reads[i]
                if r1.reference_start > last_covered:
                    reads[i] = None
                    output_bam.write(r1)
                    output_bam.write(r2)
                    last_covered = r1.reference_end
        else:
            for i in range(number_of_reads - 1, -1, -1):
                if reads[i] is None:
                    continue
                r1, r2 = reads[i]
                if r2.reference_end < last_covered:
                    reads[i] = None
                    output_bam.write(r1)
                    output_bam.write(r2)
                    last_covered = r2.reference_start
        direction = '>' if direction == '<' else '<'
        last_covered = -1 if direction == '>' else float('inf')
    output_bam.close()


def main() -> None:
    """Uruchomienie skryptu."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", "--cycles", help="Number of cycles", type=int, default=30)
    parser.add_argument("-m", "--mode", help="Mode", choices=["paired", "single"], default="single")
    parser.add_argument("--chr", help="List of chromosomes to process", nargs="+", default='all')
    parser.add_argument("input_bam", help="Input bam file")
    parser.add_argument("output_bam", help="Output bam file")
    args = parser.parse_args()
    down_sample_bam(input_bam_path=args.input_bam,
                    cycles=args.cycles,
                    chromosomes=args.chr,
                    mode=args.mode,
                    output_bam_path=args.output_bam)


if __name__ == "__main__":
    main()
