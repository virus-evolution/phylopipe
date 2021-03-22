#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description="""Apply a mask to some bases of alignment""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-alignment', dest = 'in_alignment', required=True, help='Aligned FASTA')
    parser.add_argument('--mask', dest = 'mask', required=True, help='Mask CSV of pos, mask character, regex')
    parser.add_argument('--out-alignment', dest = 'out_alignment', required=True, help='FASTA to write out')

    args = parser.parse_args()
    return args

def parse_mask_file(file):
    """
    input is in the format:
    start (1-based), mask character, regex-format string to match record.id
    e.g.:
    13402,?,^Belgium/
    d is a dictionary with the regex strings as keys and position,
    mask character and compiled regular expression as values.
    it has the same number of entries as lines in file
    """
    d = {}

    with open(file, 'r') as f:
        for line in f:
            l = line.rstrip().split(',')
            pos, mask_char, regex = l

            d[regex] = {'pos': int(pos),
                        'mask_char': mask_char,
                        'regex': re.compile(regex)}

    return(d)

def apply_mask(in_fasta, out_fasta, mask):
    mask_info = parse_mask_file(mask)

    with open(in_fasta, "r") as fasta_in, \
         open(out_fasta, "w") as fasta_out:

        for record in SeqIO.parse(fasta_in, 'fasta'):
            ID = record.id
            seq = str(record.seq)

            for entry in mask_info:
                regex = mask_info[entry]['regex']

                if re.search(regex, ID):
                    pos = mask_info[entry]['pos']
                    mask_char = mask_info[entry]['mask_char']
                    seq = seq[:pos - 1] + mask_char + seq[pos:]

            fasta_out.write('>' + ID + '\n')
            fasta_out.write(seq + '\n')

def main():
    args = parse_args()
    apply_mask(args.in_alignment, args.out_alignment, args.mask)

if __name__ == '__main__':
    main()
