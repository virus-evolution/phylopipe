#!/usr/bin/env python3

import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="""Convert usher stdout to metadata CSV""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in', dest = 'in_file', required=True, help='Usher log file from STDOUT')
    parser.add_argument('--out', dest = 'out_file', required=True, help='Out CSV')

    args = parser.parse_args()
    return args

def parse_log(in_file, out_file):
    """
    Usher STDOUT > CSV
    """
    with open(in_file, 'r') as in_handle, \
         open(out_file, 'w') as out_handle:
        out_handle.write("sequence_name,parsimony_score,num_parsimony_optimal_placements,is_unreliable_in_tree\n")
        for line in in_handle:
            if "Number of parsimony-optimal placements: " not in line:
                continue
            print(line)
            fields = line.strip().split('\t')
            print(fields)
            assert len(fields) > 3
            row = {"sequence_name": fields[1].split(": ")[1], "parsimony_score": fields[2].split(": ")[1], "num_parsimony_optimal_placements": fields[3].split(": ")[1]}
            if int(row["num_parsimony_optimal_placements"]) > 1:
                row["is_unreliable_in_tree"] = "Y"
            else:
                row["is_unreliable_in_tree"] = "N"
            out_handle.write("%s,%s,%s,%s\n" %(row["sequence_name"], row["parsimony_score"], row["num_parsimony_optimal_placements"],row["is_unreliable_in_tree"]))

def main():
    args = parse_args()
    parse_log(args.in_file, args.out_file)

if __name__ == '__main__':
    main()
