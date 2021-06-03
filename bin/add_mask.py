#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description="""Apply a mask to some bases of alignment""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-alignment', dest = 'in_alignment', required=True, help='Aligned FASTA')
    parser.add_argument('--mask', dest = 'mask', required=False, help='Mask CSV of pos, mask character, regex, note')
    parser.add_argument('--vcf', dest = 'vcf', required=False, help='Mask VCF')
    parser.add_argument('--filters', dest = 'filters', required=False, default=None, nargs='+', help='Filter strings to mask for')
    parser.add_argument('--out-alignment', dest = 'out_alignment', required=True, help='FASTA to write out')

    args = parser.parse_args()
    return args

def parse_mask_file(file, filters=None):
    """
    input is in the format:
    start (1-based), mask character, regex-format string to match record.id
    e.g.:
    13402,?,^Belgium/
    d is a dictionary with the regex strings as keys and position,
    mask character and compiled regular expression as values.
    it has the same number of entries as lines in file
    """
    d = []

    with open(file, 'r') as f:
        for line in f:
            l = line.rstrip().split(',')
            pos, mask_char, regex, note = l

            if filters:
                for filter in filters:
                    if filter in note:
                        d.append({'pos': int(pos),
                                  'mask_char': mask_char,
                                  'regex': re.compile(regex),
                                  'note':note})
                        break
            else:
                d.append({'pos': int(pos),
                          'mask_char': mask_char,
                          'regex': re.compile(regex),
                          'note':note})

    return(d)

def parse_vcf_file(file, filters=None, mask_char='N', regex='\w'):
    """
    input is VCF
    """
    d = []

    with open(file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            l = line.rstrip().split('\t')
            chrom, pos, id, ref, alt, qual, filter, info = l

            if filter != "mask":
                continue

            note = None
            if filters:
                for filter in filters:
                    if filter in info:
                        note = filter
                        break

            if note or not filters:
                d.append({'pos': int(pos),
                          'mask_char': mask_char,
                          'regex': re.compile(regex),
                          'note':note})
                print("Adding mask %s to position %s with regex %s and info %s" %(mask_char, pos, regex, info))
            else:
                print("Skipping position %s with info %s" %(pos, info))

    return(d)

def apply_mask(in_fasta, out_fasta, mask=None, vcf=None, filters=None):
    if mask:
        mask_info = parse_mask_file(mask)
    elif vcf:
        mask_info = parse_vcf_file(vcf, filters)
    else:
        print("No mask provided!!")
        return

    records = SeqIO.index(in_fasta, "fasta")

    with open(out_fasta, "w") as fasta_out:
        for record in records:
            ID = record
            seq = str(records[ID].seq)

            for entry in mask_info:
                regex = entry['regex']

                if re.search(regex, ID):
                    pos = entry['pos']
                    mask_char = entry['mask_char']
                    seq = seq[:pos - 1] + mask_char + seq[pos:]
                    print("updated sequence %s at pos %d with char %s" %(ID, pos, mask_char))

            fasta_out.write('>' + ID + '\n')
            fasta_out.write(seq + '\n')

def main():
    args = parse_args()
    apply_mask(args.in_alignment, args.out_alignment, args.mask, args.vcf, args.filters)

if __name__ == '__main__':
    main()
