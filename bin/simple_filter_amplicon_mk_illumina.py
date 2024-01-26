#!/usr/bin/env python3

"""Script for down-sampling bam file trying to keep the coverage uniform."""

# ruff: noqa: T201, Q000, TRY0003, FA100

import argparse
from collections import defaultdict
from typing import Optional, Union

import pysam
from bitarray import bitarray


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
            run_mode_paired(input_bam_file, chr_id, cycles, output_bam_file)
        else:
            msg = f"Unknown mode: {mode}"
            raise ValueError(msg)
    print(f"Output bam file: {output_bam_path} saved.")
    print(f"Sorting and indexing {output_bam_path} file...")
    pysam.sort("-o", output_bam_path, output_bam_path)
    pysam.index(output_bam_path)
    print("Done.")


def run_mode_single(bam_file: pysam.AlignmentFile, chr_id: str, cycles: int, output_bam: pysam.AlignmentFile) -> None:
    """Run down-sampling in single reads mode.

    This mode is one-way mode. In each cycle, reads
    are processed from the beginning to the end."""
    number_of_reads = bam_file.count(chr_id)
    last_covered = -1
    used_reads = bitarray(number_of_reads)
    direction = '>'
    for c in range(cycles):
        print(f"Cycle {c + 1} {direction}")
        for i, r1 in enumerate(bam_file.fetch(region=chr_id)):
            if used_reads[i]:
                continue
            if r1.reference_start > last_covered:
                output_bam.write(r1)
                last_covered = r1.reference_end - 1
                used_reads[i] = 1
        last_covered = -1
    output_bam.close()


def run_mode_single_two_way(bam_file: pysam.AlignmentFile, chr_id: str, cycles: int, output_bam: pysam.AlignmentFile) -> None:
    """Run down-sampling in single reads mode."""
    # TODO memory inefficient
    # TODO Two way mode is used only for experiments.
    last_covered = -1
    reads = list(bam_file.fetch(region=chr_id))
    number_of_reads = len(reads)
    direction = '<'
    for c in range(cycles):
        print(f"Cycle {c + 1} {direction}")
        if direction == '>':
            for i in range(number_of_reads):
                r1 = reads[i]
                if r1 is None:
                    continue
                if r1.reference_start > last_covered:
                    output_bam.write(r1)
                    last_covered = r1.reference_end - 1
                    reads[i] = None
        else:
            for i in range(number_of_reads - 1, -1, -1):
                r1 = reads[i]
                if r1 is None:
                    continue
                if r1.reference_end < last_covered:
                    output_bam.write(r1)
                    last_covered = r1.reference_start + 1
                    reads[i] = None
        direction = '>' if direction == '<' else '<'
        last_covered = -1 if direction == '>' else float('inf')
    output_bam.close()


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
