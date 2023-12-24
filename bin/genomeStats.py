#!/usr/bin/env python3

import argparse
import gzip
import json
import os
import magic

from Bio import SeqIO


def analyze_genome(genome: str) -> dict:
    """Analyze a genome and return a dictionary of stats"""
    chromosomes = []
    stats = {}
    k = 0

    m = magic.open(magic.MAGIC_MIME_ENCODING)
    m.load()
    with open(genome, 'rb') as f:
        encoding = m.buffer(f.read())

    for k, r in enumerate(SeqIO.parse(genome, format='fasta')):
        chromosome = {"number": k + 1,
                      "id": r.id,
                      "header": r.description,
                      "length": len(r.seq),
                      'A': r.seq.count('A'),
                      'C': r.seq.count('C'),
                      'T': r.seq.count('T'),
                      'G': r.seq.count('G'),
                      'N': r.seq.count('N'),
                      }
        chromosome['Other_bases'] = chromosome['length'] - (chromosome['A'] +
                                                            chromosome['C'] +
                                                            chromosome['T'] +
                                                            chromosome['G'] +
                                                            chromosome['N'])
        chromosome['GC_content'] = (chromosome['G'] + chromosome['C']) / chromosome['length']
        chromosomes.append(chromosome)
    meta = {
        "filename": os.path.basename(genome),
        "file_md5sum": os.popen(f"md5sum {genome}").read().split()[0],
        "file_size_bytes": os.path.getsize(genome),
        "encoding": encoding,
    }
    stats["meta"] = meta
    summary = {"number_of_chromosomes": k + 1}
    for c in range(len(chromosomes)):
        for k in ['length', 'A', 'C', 'T', 'G', 'N', 'Other_bases']:
            if k in summary:
                summary[k] += chromosomes[c][k]
            else:
                summary[k] = chromosomes[c][k]
    list_of_ids = [c['id'] for c in chromosomes]
    summary['chromosome_ids'] = list_of_ids
    stats["summary"] = summary
    stats["chromosomes"] = chromosomes
    return stats


def main():
    parser = argparse.ArgumentParser(description='Generate some simple stats about a genome in json format')
    parser.add_argument('genome', help='genome fasta file')
    parser.add_argument('--output', '-o', default="genomeStats.json.gz", help='output json file')
    args = parser.parse_args()
    stats = analyze_genome(args.genome)
    # with open(args.output, "w") as f:
    #     json.dump(stats, f, indent=4)
    with gzip.open(args.output, 'w') as f:
        f.write(json.dumps(stats, indent=4).encode('utf-8'))
    print(f"File {args.output} saved...")


if __name__ == '__main__':
    main()
