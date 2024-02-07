#!/usr/bin/env python

"""Parse alignment with masked reference sequence and output a fasta.

This script takes an alignment file with two seqences:
- a reference sequence with masked regions (N's)
- a sample sequence from SV/INDEL/SNP caller.
It outputs a sample sequence with masked regions where the reference
sequence is masked.

Example:
    Ref:    NNNNNCTGNNNACTGACTGA---CTGNNN
    Sample: ACTGACTGACTACT---TGAACTCTGACT
    Output: NNNNNCTGNNNACTTGAACTCTGNNN

Deletions takes priority over masking.

Biopython.SeqIO is used to handle fasta files.
"""

import argparse

import Bio.SeqIO


def parse_alignment(alignment_path, output_path):
    with open(alignment_path, 'r') as f:
        alignment = Bio.SeqIO.parse(f, 'fasta')
        ref = next(alignment)
        sample = next(alignment)

    output = ''
    # ref = ref.seq.lower()
    # sample = sample.seq.lower()
    for r, s in zip(ref.seq, sample.seq):
        if s == '-':
            continue
        else:
            if r == 'n':
                output += 'n'
            else:
                output += s
    record = Bio.SeqIO.SeqRecord(Bio.Seq.Seq(output), id=sample.id)

    with open(output_path, "w") as output_handle:
        Bio.SeqIO.write([record], output_handle, "fasta")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('alignment', help='alignment file')
    parser.add_argument('output', help='output fasta file')
    args = parser.parse_args()
    parse_alignment(args.alignment, args.output)


if __name__ == '__main__':
    main()
