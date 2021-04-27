#!/usr/bin/env python3

from Bio import SeqIO
import csv
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="""Filter UK sequences based on metadata""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-fasta', dest = 'in_fasta', required=True, help='Aligned FASTA')
    parser.add_argument('--in-metadata', dest = 'in_metadata', required=True, help='CSV: if provided hash keeps most recent sequence as representative')
    parser.add_argument('--outgroups', dest = 'outgroups', required=False, help='Lineage splits file containing representative outgroups to protect')
    parser.add_argument('--out-fasta', dest = 'out_fasta', required=True, help='FASTA to write out')
    parser.add_argument('--out-metadata', dest = 'out_metadata', required=True, help='CSV to write out')
    parser.add_argument('--size', required=False, default=100, help='Number to keep from each lineage')

    args = parser.parse_args()
    return args

def parse_outgroups(outgroup_file):
    """
    input is CSV, last column being the representative outgroups:
    """
    outgroups = {}
    if not outgroup_file:
        return outgroups
    with open(outgroup_file, "r") as outgroup_handle:
        line = outgroup_handle.readline()
        while line:
            try:
                lineage,outgroup = line.strip().split(",")
                outgroups[outgroup]=lineage
            except:
                continue
            line = outgroup_handle.readline()
    return(outgroups)

def get_cog_uk_status(row):
    if "country" in row and row["country"] in ["UK", "United_Kingdom", "United Kingdom"]:
        return "cog_uk"
    elif "adm1" in row and row["adm1"] in ["England", "Wales", "Scotland", "Northern_Ireland",
                                           "Falkland_Islands", "Gibraltar", "Jersey", "The_Isle_of_Man", "Guernsey"]:
        return "cog_uk"
    return "global"

def get_split(row, sorted_list_outgroups):
    if "lineage" not in row:
        return None
    for outgroup in sorted_list_outgroups:
        if outgroup == row["lineage"] or outgroup.startswith("%s." % row["lineage"]):
            return outgroup
    return None


def make_test(in_fasta, in_metadata, outgroup_file, out_fasta, out_metadata, size):
    outgroups = parse_outgroups(outgroup_file)
    sorted_list_outgroups = outgroups.keys().sort(key=lambda x: re.sub("[^A-Z0-9]", "", x), reverse=True)
    print(sorted_list_outgroups)
    records = SeqIO.index(in_fasta, "fasta")

    counts = {"cog_uk":{}, "global": {}}
    all_counts = 0
    for lineage in sorted_list_outgroups:
        counts["cog_uk"][lineage] = 0
        counts["global"][lineage] = 0

    with open(in_metadata, 'r') as csv_in, \
        open(out_metadata, 'w') as csv_out, \
        open(out_fasta, 'w') as fa_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            fasta_header = row["sequence_name"]
            if fasta_header not in records:
                continue

            status = get_cog_uk_status(row)
            split = get_split(row, sorted_list_outgroups)

            if fasta_header in outgroups:
                writer.writerow(row)
                fa_out.write(">%s\n" %row["sequence_name"])
                fa_out.write("%s\n" %str(records[fasta_header].seq))
                counts[status][split] += 1
                continue

            if "why_excluded" in reader.fieldnames and row["why_excluded"] not in [None, "None", ""]:
                continue

            if counts[status][split] < size:
                writer.writerow(row)
                fa_out.write(">%s\n" %row["sequence_name"])
                fa_out.write("%s\n" %str(records[fasta_header].seq))
                counts[status][split] += 1
                all_counts += 1

            if all_counts > len(outgroups) * size:
                break


def main():
    args = parse_args()
    make_test(args.in_fasta, args.in_metadata, args.outgroups, args.out_fasta, args.out_metadata, args.size)

if __name__ == '__main__':
    main()