#!/usr/bin/env python3

from Bio import SeqIO
import re
import argparse
import unicodedata

def parse_args():
    parser = argparse.ArgumentParser(description="""Cleans FASTA headers and updates tree file with these substitutions""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-fasta', dest = 'in_fasta', required=True, help='Aligned FASTA')
    parser.add_argument('--in-tree', dest = 'in_tree', required=False, help='Tree file')
    parser.add_argument('--out-fasta', dest = 'out_fasta', required=False, help='FASTA to write out')
    parser.add_argument('--out-metadata', dest = 'out_metadata', required=False, help='CSV to write out')
    parser.add_argument('--out-tree', dest = 'out_tree', required=False, help='If in-tree specified, must specify out tree too')

    args = parser.parse_args()
    return args

def strip_accents(s):
    s = ''.join((c for c in unicodedata.normalize('NFD', s) if unicodedata.category(c) != 'Mn'))
    return s

def clean(in_fasta, in_tree, out_fasta, out_metadata, out_tree):
    if not out_fasta:
        out_fasta = in_fasta.replace(".fa",".clean.fa")
    if not out_metadata:
        out_metadata = re.sub(".fa[sta]*",".map.csv", in_fasta)
    if in_tree and not out_tree:
        parts = in_tree.split(".")
        out_tree = ".".join(parts[:-1] + ["clean"] + parts[-1:])

    records = SeqIO.index(in_fasta, "fasta")
    new_names = {}

    with open(out_fasta, 'w') as fa_out:

        for record in records:
            id_string = record
            if not id_string.isalnum():
                id_string = strip_accents(id_string)
                p = re.compile(r'[^a-zA-Z0-9_/-]')
                id_string = p.sub('_', id_string)
                if id_string != record:
                    new_names[record] = id_string
            fa_out.write(">%s\n" %id_string)
            fa_out.write("%s\n" %str(records[record].seq))

    with open(out_metadata, 'w') as csv_out:
        csv_out.write('original,new\n')
        for item in new_names:
            csv_out.write("%s,%s\n" % (item, new_names[item]))

    if in_tree:
        with open(in_tree, 'r') as tree_in, open(out_tree, 'w') as tree_out:
            for line in tree_in:
                for item in new_names:
                    line.replace(item, new_names[item])
                tree_out.write(line)

def main():
    args = parse_args()
    clean(args.in_fasta, args.in_tree, args.out_fasta, args.out_metadata, args.out_tree)

if __name__ == '__main__':
    main()
