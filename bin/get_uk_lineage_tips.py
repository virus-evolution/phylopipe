#!/usr/bin/env python3

import csv
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="""Filter UK sequences based on metadata""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-metadata', dest = 'in_metadata', required=True, help='CSV: if provided hash keeps most recent sequence as representative')

    args = parser.parse_args()
    return args

def split(in_metadata):
    tips = {}
    sequence_name = "taxon"

    with open(in_metadata, 'r') as csv_in:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        for row in reader:
            if row["uk_lineage"] in ["None", None, ""]:
                continue
            elif row["uk_lineage"] in tips:
                tips[row["uk_lineage"]].append(row[sequence_name])
            else:
                tips[row["uk_lineage"]] = [row[sequence_name]]

    for uk_lineage in tips:
        if len(tips[uk_lineage]) > 2:
            with open("%s.txt" %uk_lineage, "w") as out:
                for it in tips[uk_lineage]:
                    out.write("%s\n" %it)


def main():
    args = parse_args()
    split(args.in_metadata)

if __name__ == '__main__':
    main()
