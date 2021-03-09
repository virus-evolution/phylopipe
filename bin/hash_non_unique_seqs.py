#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import csv

def parse_args():
    parser = argparse.ArgumentParser(description="""Pick a representative sample for each unique sequence""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-fasta', dest = 'in_fasta', required=True, help='Aligned FASTA')
    parser.add_argument('--in-metadata', dest = 'in_metadata', required=False, help='CSV: if provided hash keeps most recent sequence as representative')
    parser.add_argument('--outgroups', dest = 'outgroups', required=False, help='Lineage splits file containing representative outgroups')
    parser.add_argument('--out-fasta', dest = 'out_fasta', required=True, help='FASTA to write out')
    parser.add_argument('--out-metadata', dest = 'out_metadata', required=True, help='CSV to write out')

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

def write_hash_dict(in_fasta, out_fasta, out_metadata, outgroup_file, in_metadata=None):
    outgroups = parse_outgroups(outgroup_file)
    records = SeqIO.index(in_fasta, "fasta")
    hash_dict = {}

    for record_id in records:
        seq = str(records[record_id].seq)

        if seq in hash_dict:
            hash_dict[seq] = hash_dict[seq] + [record_id]
        else:
            hash_dict[seq] = [record_id]

    date_dict = {}
    if in_metadata is not None:
        with open(in_metadata, 'r', newline = '') as csv_in:
            reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
            for row in reader:
                if row["epi_week"] not in [None, "None", ""]:
                    date_dict[row["fasta_header"]] = row["epi_week"]
        print("Found", len(date_dict), "epi_week dates")

    with open(out_fasta, "w") as fasta, open(out_metadata, "w") as metadata:
        metadata.write("tip,redundant\n")

        for key, value in hash_dict.items():
            if len(value) == 1:
                r = records[value[0]]
                fasta.write(">" + r.id + "\n")
                fasta.write(str(r.seq) + "\n")

            elif len(value) > 1:
                r = None
                for id in value:
                    if id in outgroups:
                        r = records[id]
                        fasta.write(">" + r.id + "\n")
                        fasta.write(str(r.seq) + "\n")
                        value.remove(id)

                if not r:
                    index = 0
                    date = 0
                    if in_metadata is not None:
                        for i,id in enumerate(value):
                            if id in date_dict and date_dict[id] > date:
                                index = i
                                date = date_dict[id]
                    r = records[value[index]]
                    fasta.write(">" + r.id + "\n")
                    fasta.write(str(r.seq) + "\n")
                    value.remove(value[index])

                metadata.write(r.id + ",")
                metadata.write("|".join(value) + "\n")

def main():
    args = parse_args()
    write_hash_dict(args.in_fasta, args.out_fasta, args.out_metadata, args.outgroups)

if __name__ == '__main__':
    main()
