#!/usr/bin/env python3

import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="""Apply a mask to some bases of alignment""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-alignment', dest = 'in_alignment', required=True, help='Aligned FASTA')
    parser.add_argument('--reference', dest = 'reference', required=True, help='FASTA of reference sequence')
    parser.add_argument('--out-alignment', dest = 'out_alignment', required=True, help='Aligned FASTA')

    args = parser.parse_args()
    return args

def remove_reference(fasta_file, out_file, reference_file):
    """
    remove fasta entry from alignment that matches with ref id
    """
    with open(reference_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                reference = line

    found = False
    skip = False
    with open(fasta_file, 'r') as f, \
         open(out_file, 'w') as g:
        for line in f:
            if skip:
                skip = False
            elif found:
                g.write(line)
            elif line == reference:
                skip = True
                found = True
            else:
                g.write(line)

def main():
    args = parse_args()
    remove_reference(args.in_alignment, args.out_alignment, args.reference)

if __name__ == '__main__':
    main()
