#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO
import csv

def parse_args():
    parser = argparse.ArgumentParser(description="""Add columns to metadata for specific AAs and dels""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-fasta', dest = 'in_fasta', required=True, help='Aligned FASTA')
    parser.add_argument('--in-metadata', dest = 'in_metadata', required=True, help='CSV of metadata to add to')
    parser.add_argument('--out-metadata', dest = 'out_metadata', required=True, help='CSV to write out')
    parser.add_argument('--reference-fasta', dest = 'reference_fasta', required=True, help='Reference FASTA')
    parser.add_argument('--aas', dest = 'aas', required=False, help='CSV of AAs')
    parser.add_argument('--dels', dest = 'dels', required=False, help='CSV of deletions')
    parser.add_argument('--index-column', dest = 'index_column', required=False, default='sequence_name')

    args = parser.parse_args()
    return args


def parse_AA_file(file):
    """
    input is in the format:
    start (1-based)
    e.g.:
    D614G,1605

    ls is a list of length-2 tuples with the format (name, position)
    position is the 1-based starting position of the codon in Wuhan-Hu-1 coordinates
    It has the same number of entries as lines in file
    """
    ls = []
    if not file:
        return ls

    with open(file, 'r') as f:
        for line in f:
            l = line.rstrip().split(",")
            name, pos = l
            ls = ls + [(name, int(pos))]
    return(ls)

def parse_del_file(file, ref_fasta):
    """
    input is in the format:
    start (1-based), length of deletion
    e.g.:
    1605,3

    ls is a list of length-3 tuples with the format (position, length, ref_allele)
    It has the same number of entries as lines in file
    """
    ls = []
    if not file:
        return ls
    WuhanHu1 = SeqIO.read(ref_fasta, 'fasta')

    with open(file, 'r') as f:
        for line in f:
            l = line.rstrip().split(',')
            pos, length = l
            ref_allele = str(WuhanHu1.seq).upper()[int(pos) - 1: int(pos) - 1 + int(length)]
            ls = ls + [(int(pos), int(length), ref_allele)]

    return(ls)

def type_aas_and_dels(in_fasta, in_aa_file, in_del_file, reference_fasta, in_metadata, out_metadata, index_column):
    alignment = SeqIO.index(in_fasta, "fasta")
    AAs = parse_AA_file(in_aa_file)
    dels = parse_del_file(in_del_file, reference_fasta)

    new_aa_columns = [x[0] for x in AAs]
    new_del_columns = ["del_" + str(x[0]) + "_" + str(x[1]) for x in dels]

    with open(in_metadata, 'r', newline = '') as csv_in, \
         open(out_metadata, 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames + new_aa_columns + new_del_columns, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            id = row[index_column]
            seq = alignment[id].seq

            for entry in AAs:
                pos = entry[1]
                try:
                    QUERY_allele = seq[pos - 1: pos + 2].translate()
                except:
                    QUERY_allele = 'X'
                row[entry[0]] = QUERY_allele

            for entry in dels:
                pos = entry[0]
                length = entry[1]
                ref_allele = entry[2]
                column_name = "del_" + str(pos) + "_" + str(length)

                if seq[pos - 1: pos - 1 + length] == '-' * length:
                    genotype = 'del'
                elif seq[pos - 1: pos - 1 + length] == ref_allele:
                    genotype = 'ref'
                else:
                    genotype = 'X'

                row[column_name] = genotype

            writer.writerow(row)

def main():
    args = parse_args()
    type_aas_and_dels(args.in_fasta, args.aas, args.dels, args.reference_fasta, args.in_metadata, args.out_metadata, args.index_column)

if __name__ == '__main__':
    main()