#!/usr/bin/env python3

import argparse
import gzip
import json
from typing import List


def combine(json_files: List[str], output: str):
    combined = []
    for json_file in json_files:
        with gzip.open(json_file, 'r') as f:
            combined.append(json.load(f))
    with gzip.open(output, 'w') as f:
        f.write(json.dumps(combined, indent=4).encode('utf-8'))
    print(f"File {output} saved...")


def main():
    parser = argparse.ArgumentParser(description='Combine two or more json.gz files')
    parser.add_argument('json_files', nargs='+', help='json.gz files')
    parser.add_argument('--output', '-o', default="combined.json.gz", help='output json file')
    args = parser.parse_args()
    combine(args.json_files, args.output)


if __name__ == '__main__':
    main()
