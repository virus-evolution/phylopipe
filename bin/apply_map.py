#!/usr/bin/env python3

from Bio import SeqIO
import csv
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="""Filter UK sequences based on metadata""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-metadata', dest = 'in_metadata', required=True, help='CSV: if provided hash keeps most recent sequence as representative')
    parser.add_argument('--in-map', dest = 'in_map', required=True, help='Map old sequence_name to new')
    parser.add_argument('--out-metadata', dest = 'out_metadata', required=True, help='CSV: output')

    args = parser.parse_args()
    return args

def parse_map(map_file):
    """
    input is CSV, with columns original and new
    """
    map = {}
    with open(map_file, "r") as map_handle:
        line = map_handle.readline()
        while line:
            try:
                original,new = line.strip().split(",")
                map[original]=new
            except:
                continue
            line = map_handle.readline()
    return(map)


def correct(in_metadata, in_map, out_metadata):
    map = parse_map(in_map)

    with open(in_metadata, 'r') as csv_in, \
        open(out_metadata, 'w') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(csv_out, fieldnames = fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        fasta_header = "sequence_name"
        for row in reader:
            if row[fasta_header] in map:
                row[fasta_header] = map[row[fasta_header]]
            writer.writerow(row)


def main():
    args = parse_args()
    correct(args.in_metadata, args.in_map, args.out_metadata)

if __name__ == '__main__':
    main()
