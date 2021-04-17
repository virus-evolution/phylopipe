#!/usr/bin/env python3

from Bio import SeqIO
import csv
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="""Filter UK sequences based on metadata""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-csv', dest = 'in_csv', required=True, help='CSV of taxon,lineage')
    parser.add_argument('--out_tsv', dest = 'out_csv', required=True, help='TSV for usher pangolin')

    args = parser.parse_args()
    return args

def convert(in_csv, out_csv):
    """
    input is CSV, last column being the representative outgroups:
    """
    with open(in_csv, "r") as csv_in, open(out_tsv, "r") as tsv_out:
        line = csv_in.readline()
        while line:
            line = csv_in.readline()
            try:
                taxon,lineage = line.strip().split(",")
                tsv_out.write("%s\t%s\n" % (lineage, taxon))
            except:
                continue

def main():
    args = parse_args()
    convert(args.in_csv, args.out_tsv)

if __name__ == '__main__':
    main()
