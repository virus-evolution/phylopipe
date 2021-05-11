#!/usr/bin/env python3

from Bio import SeqIO
import csv
import argparse
import unicodedata

def parse_args():
    parser = argparse.ArgumentParser(description="""Filter UK sequences based on metadata""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-metadata', dest = 'in_metadata', required=True, help='CSV: if provided hash keeps most recent sequence as representative')
    parser.add_argument('--in-map', dest = 'in_map', required=True, help='Map old sequence_name to new')
    parser.add_argument('--to-clean', dest = 'to_clean', required=False, nargs='+', default=[], help='List of metadata fields to clean up')
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

def strip_accents(s):
    if not s.isalnum():
        s = ''.join((c for c in unicodedata.normalize('NFD', s) if unicodedata.category(c) != 'Mn'))
        s = s.replace("'", "_")
    return s

def correct(in_metadata, in_map, to_clean, out_metadata):
    map = parse_map(in_map)

    with open(in_metadata, 'r') as csv_in, \
        open(out_metadata, 'w') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        fieldnames = reader.fieldnames
        clean_columns = []
        for column in to_clean:
            if column not in fieldnames:
                print("Column '%s' not in input" %column)
            else:
                clean_columns.append(column)
        writer = csv.DictWriter(csv_out, fieldnames = fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        fasta_header = "sequence_name"
        for row in reader:
            if row[fasta_header] in map:
                row[fasta_header] = map[row[fasta_header]]
            for column in clean_columns:
                row[column] = strip_accents(row[column])
            writer.writerow(row)


def main():
    args = parse_args()
    correct(args.in_metadata, args.in_map, args.to_clean, args.out_metadata)

if __name__ == '__main__':
    main()
