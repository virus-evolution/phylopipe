#!/usr/bin/env python3

import sys
import random
import copy
import argparse
import csv
import re

def get_random_name():
    r = random.randint(0, 999999)
    num_str = "{:06d}".format(r)
    anonymous_name = 'COG' + str(num_str)
    return(anonymous_name)


def anonymize_microreact(metadata_in, tree_in, metadata_out, tree_out, seed):
    if seed:
        random.seed(seed)

    admn2_counts = {}
    with open(metadata_in, 'r', newline = '') as csv_in:
        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        if not all(x in reader.fieldnames for x in ['adm2', 'is_cog_uk', 'sequence_name']):
            sys.exit('required columns not found in metadata')
        for row in reader:
            location = row['adm2']
            if location in admn2_counts:
                admn2_counts[location] += 1
            else:
                admn2_counts[location] = 1

    anonymous_locations = set()
    for location in admn2_counts:
        if admn2_counts[location] < 5:
            anonymous_locations.add(location)
    print("Found %i anonymous_locations" %len(anonymous_locations))
    del admn2_counts

    anonymous_names = {}
    with open(metadata_in, 'r', newline = '') as csv_in, \
         open(metadata_out, 'w', newline = '') as csv_out:

        reader = csv.DictReader(csv_in, delimiter=",", quotechar='\"', dialect = "unix")
        writer = csv.DictWriter(csv_out, fieldnames = reader.fieldnames, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()

        for row in reader:
            if row['adm2'] in anonymous_locations:
                row['adm2'] = ''

            if row['is_cog_uk'] in ["True", True]:
                anonymous_name = get_random_name()
                while anonymous_name in anonymous_names:
                    anonymous_name = get_random_name()
                anonymous_names[row['sequence_name']] = anonymous_name
                row['sequence_name'] = anonymous_name

            writer.writerow(row)
    print("Found %i anonymous_names" %len(anonymous_names))
    tree = open(tree_in, 'r').read()
    tip = re.compile("[A-Za-z0-9_-]+/[A-Za-z0-9_-]+/202[01]", re.U)
    start_pos = 0
    new_strings = []
    m = tip.search(tree, start_pos)
    count = 0
    while m:
        #if count % 10 == 0:
        #    print(count, start_pos, m.start(), m.end())
        new_strings.append(tree[start_pos:m.start()])
        if m.group(0) in anonymous_names:
            count += 1
            new_strings.append(anonymous_names[m.group(0)])
            del anonymous_names[m.group(0)]
        else:
            new_strings.append(m.group(0))
        start_pos = m.end()
        m = tip.search(tree, start_pos)

    new_strings.append(tree[start_pos:])
    print("Replaced", count, "sequence_names in tree")
    tree = "".join(new_strings)

    for old, new in anonymous_names.items():
        tree = re.sub(old, new, tree)

    tree_out = open(tree_out, 'w')
    tree_out.write(tree)
    tree_out.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""anonymize_microreact""")

    parser.add_argument('--input-tree',
                        help='newick format tree whose UK labels will be anonymized',
                        required=True,
                        dest='tree_in',
                        metavar='input.tree')
    parser.add_argument('--input-metadata',
                        help='metadata whose UK sequence names and rare admn2 values will be anonymized',
                        required=True,
                        dest='metadata_in',
                        metavar='input.csv')
    parser.add_argument('--output-tree',
                        help='anonymized newick format tree to write',
                        required=True,
                        dest='tree_out',
                        metavar='output.tree')
    parser.add_argument('--output-metadata',
                        help='anonymized metadata file to write',
                        required=True,
                        dest='metadata_out',
                        metavar='output.csv')
    parser.add_argument('--seed',
                         help='random seed to use',
                         required=False,
                         type=int,
                         metavar='N')
    args = parser.parse_args()

    anonymize_microreact(metadata_in = args.metadata_in,
                         tree_in = args.tree_in,
                         metadata_out = args.metadata_out,
                         tree_out = args.tree_out,
                         seed = args.seed)
