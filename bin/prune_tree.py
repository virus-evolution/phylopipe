#!/usr/bin/env python3

import argparse
import csv
import sys
import ete3
from ete3 import Tree

def parse_args():
    parser = argparse.ArgumentParser(description="""Apply a mask to some bases of alignment""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-tree', dest = 'in_tree', required=True, help='Newick tree')
    parser.add_argument('--metadata', dest = 'metadata', required=True, help='CSV of metadata')
    parser.add_argument('--index-column', dest = 'index_column', required=False, default='sequence_name',
                        help='CSV of metadata')
    parser.add_argument('--out-tree', dest = 'out_tree', required=True, help='Newick tree')

    args = parser.parse_args()
    return args

def get_keep_tips(metadata, index_column):
    keep = set()
    with open(metadata, 'r', newline = '') as csv_in:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        if index_column not in reader.fieldnames:
            sys.exit("Index column %s not in CSV" %index_column)
        for row in reader:
            keep.add(row[index_column])
    return keep

def prune_tips_with_metadata(in_tree, metadata, index_column, out_tree):
    keep = get_keep_tips(metadata, index_column)

    t = Tree(in_tree)
    for tip in keep:
        try:
            t.prune(keep)
        except:
            print("%s not in tree" %tip)
    t.write(outfile=out_tree)

def main():
    args = parse_args()
    prune_tips_with_metadata(args.in_tree, args.metadata, args.index_column, args.out_tree)

if __name__ == '__main__':
    main()
