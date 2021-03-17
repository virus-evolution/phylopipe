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
    parser.add_argument('--include_true', required=False, nargs='+', default=[], help='List of CSV columns, include rows only if True')
    parser.add_argument('--exclude_true', required=False, nargs='+', default=[], help='List of CSV columns, exclude if True')

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

def is_uk(row):
    if "country" in row and row["country"] in ["UK", "United_Kingdom", "United Kingdom"]:
        return True
    elif "adm1" in row and row["adm1"] in ["England", "Wales", "Scotland", "Northern_Ireland"]:
        return True
    return False

def filter(in_fasta, in_metadata, outgroup_file, out_fasta, out_metadata, include_true, exclude_true):
    outgroups = parse_outgroups(outgroup_file)
    records = SeqIO.index(in_fasta, "fasta")

    with open(in_metadata, 'r') as csv_in, \
        open(out_metadata, 'w') as csv_out, \
        open(out_fasta, 'w') as out_fasta:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        fieldnames = reader.fieldnames
        if "why_excluded" not in reader.fieldnames:
            fieldnames.append("why_excluded")
        writer = csv.DictWriter(csv_out, fieldnames = fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            fasta_header = row["sequence_name"]

            if fasta_header in outgroups:
                writer.writerow(row)
                out_fasta.write(">%s\n" %fasta_header)
                out_fasta.write("%s\n" %str(records[fasta_header].seq))
                continue

            if "why_excluded" in reader.fieldnames and row["why_excluded"] not in [None, "None", ""]:
                writer.writerow(row)
                continue

            if "(" in fasta_header:
                row["why_excluded"] = "dodgy name"
                writer.writerow(row)
                continue

            filtered = False
            if is_uk(row):
                for column in exclude_true:
                    if column in reader.fieldnames and row[column] in [True, "True", "Yes", "yes", "Y"]:
                        row["why_excluded"] = "Filtered because %s is True" %column
                        filtered=True
                        break
                else:
                    for column in include_true:
                        if column in reader.fieldnames and row[column] not in [True, "True", "Yes", "yes", "Y"]:
                            row["why_excluded"] = "Filtered because %s is not True" %column
                            filtered=True
                            break
            if filtered:
                writer.writerow(row)
                continue

            if fasta_header in records:
                out_fasta.write(">%s\n" %fasta_header)
                out_fasta.write("%s\n" %str(records[fasta_header].seq))
                row["why_excluded"] = ""
                writer.writerow(row)
            else:
                row["why_excluded"] = "No matching metadata"
                writer.writerow(row)


def main():
    args = parse_args()
    filter(args.in_fasta, args.in_metadata, args.outgroups, args.out_fasta, args.out_metadata, args.include_true, args.exclude_true)

if __name__ == '__main__':
    main()
