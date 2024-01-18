#!/usr/bin/env python3

import argparse
import gzip
import json
import os

import pysam


def analyze_bam_file(bam_file: str, sample: str) -> dict:
    f = pysam.AlignmentFile(bam_file)
    stats = {}
    meta = {
        "filename": os.path.basename(bam_file),
        "sample": sample,
        "file_md5sum": os.popen(f"md5sum {bam_file}").read().split()[0],
        "file_size_bytes": os.path.getsize(bam_file),
    }
    stats["meta"] = meta

    # This situation should never happen, but we don't want to crash if it does.
    if f.mapped + f.unmapped != 0:
        mapped_fraction = f.mapped / (f.mapped + f.unmapped)
        unmapped_fraction = f.unmapped / (f.mapped + f.unmapped)
    else:
        mapped_fraction = "NaN"
        unmapped_fraction = "NaN"

    bam_stats = {
        "mapped": f.mapped,
        "unmapped": f.unmapped,
        "total": f.mapped + f.unmapped,
        "mapped_fraction": mapped_fraction,
        "unmapped_fraction": unmapped_fraction
    }
    stats["bam_stats"] = bam_stats
    return stats


def main():
    parser = argparse.ArgumentParser(description='Generate some simple stats about a bam file in json format')
    parser.add_argument('bam_file', help='alignment bam file')
    parser.add_argument('--output', '-o', default="bamStats.json.gz", help='output json file')
    parser.add_argument('-s', '--sample', default="unspecified", help='sample name')
    args = parser.parse_args()
    stats = analyze_bam_file(args.bam_file, args.sample)
    # with open(args.output, "w") as f:
    #     json.dump(stats, f, indent=4)
    with gzip.open(args.output, 'w') as f:
        f.write(json.dumps(stats, indent=4).encode('utf-8'))
    print(f"File {args.output} saved...")


if __name__ == '__main__':
    main()
