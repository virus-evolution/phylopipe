#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import re

def parse_args():
    parser = argparse.ArgumentParser(description="""Apply a mask to some bases of alignment""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-alignment', dest = 'in_alignment', required=True, help='Aligned FASTA')
    parser.add_argument('--sites', dest = 'sites', required=True, help='List of sites where if ambiguous should exclude sequence')
    parser.add_argument('--out-alignment', dest = 'out_alignment', required=True, help='FASTA to write out')
    parser.add_argument('--outgroups', dest = 'outgroups', required=False, help='Lineage splits file containing representative outgroups to protect')

    args = parser.parse_args()
    return args

def parse_outgroups(outgroup_file):
    """
    input is CSV, last column being the representative outgroups:
    """
    outgroups = []
    if not outgroup_file:
        return outgroups
    with open(outgroup_file, "r") as outgroup_handle:
        line = outgroup_handle.readline()
        while line:
            try:
                outgroup = line.strip().split(",")[-1]
                outgroups.append(outgroup)
            except:
                continue
            line = outgroup_handle.readline()
    return(outgroups)

def parse_sites(file):
    """
    get list of site positions
    """
    d = []

    with open(file, 'r') as f:
        for line in f:
            l = line.rstrip()
            d.append(int(l))
    return(d)

def apply_filter(in_fasta, out_fasta, sites_file, outgroup_file=None):
    outgroups = parse_outgroups(outgroup_file)
    sites = parse_sites(sites_file)

    records = SeqIO.index(in_fasta, "fasta")

    with open(out_fasta, "w") as fasta_out:
        for record in records:
            ID = record
            seq = str(records[ID].seq)
            keep = True

            if ID not in outgroups:
                for pos in sites:
                    if seq[pos] not in ['A','G', 'C', 'T']:
                        keep = False

            if keep:
                fasta_out.write('>' + ID + '\n')
                fasta_out.write(seq + '\n')
            else:
                print(ID)

def main():
    args = parse_args()
    apply_filter(args.in_alignment, args.out_alignment, args.sites, args.outgroups)

if __name__ == '__main__':
    main()
