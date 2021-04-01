#!/usr/bin/env python3

import sys
import random
import copy
import argparse
import csv

def get_random_name():
    r = random.randint(0, 999999)
    num_str = "{:06d}".format(r)
    anonymous_name = 'COG' + str(num_str)
    return(anonymous_name)


def anonymize_microreact(metadata_in, tree_in, metadata_out, tree_out, seed):
    if seed:
        random.seed(seed)

    metadata_list = []
    with open(metadata_in, 'r', newline = '') as csvfile:
        reader = csv.reader(csvfile, delimiter=",", quotechar='\"', dialect = "unix")
        First = True
        for line in reader:
            if First:
                header = line
                if not all(x in line for x in ['adm2', 'sequence_name']):
                    sys.exit('required columns not found in metadata')
                First = False
                continue
            d = {x:y for x,y in zip(header, line)}

            metadata_list = metadata_list + [d]


    admn2_counts = {}
    for entry in metadata_list:
        location = entry['adm2']
        if location in admn2_counts:
            admn2_counts[location] = admn2_counts[location] + 1
        else:
            admn2_counts[location] = 1

    anonymous_locations = []
    for location in admn2_counts:
        if admn2_counts[location] < 5:
            anonymous_locations.append(location)

    anonymous_locations = set(anonymous_locations)

    for entry in metadata_list:
        if entry['adm2'] in anonymous_locations:
            entry['anonymized_admn2'] = ''
        else:
            entry['anonymized_admn2'] = entry['adm2']


    random_names = []
    for entry in metadata_list:
        if entry['sequence_name'].split('/')[0] in ['Wales', 'Scotland', 'Northern_Ireland', 'England']:
            anonymous_name = get_random_name()
            while anonymous_name in random_names:
                anonymous_name = get_random_name()
            random_names.append(anonymous_name)

            entry['anonymized_sequence_name'] = anonymous_name
        else:
            entry['anonymized_sequence_name'] = entry['sequence_name']

    new_columns = copy.copy(header)
    for i, item in enumerate(new_columns):
        if item == 'sequence_name':
            new_columns[i] = 'anonymized_sequence_name'
        if item == 'adm2':
            new_columns[i] = 'anonymized_admn2'

    with open(metadata_out, 'w', newline = '') as csvfile:
        writer = csv.writer(csvfile, delimiter=",", quotechar='\"',
                            quoting=csv.QUOTE_MINIMAL, dialect = "unix")

        # writer.write(','.join(header) + '\n')
        writer.writerow(header)

        for entry in metadata_list:
            # file_out.write(','.join([entry[x] for x in new_columns]) + '\n')
            writer.writerow([entry[x] for x in new_columns])

    label_mappings = {}
    for entry in metadata_list:
        if entry['sequence_name'] != entry['anonymized_sequence_name']:
            label_mappings[entry['sequence_name']] = entry['anonymized_sequence_name']

    tree = open(tree_in, 'r').read()

    for old, new in label_mappings.items():
        tree = tree.replace("\'" + old + "\'", "\'" + new + "\'")

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
