#!/usr/bin/env python3

import argparse
import json
import subprocess
import os
import sys
import glob
import random
from anonymize_microreact import anonymize_microreact

class Error (Exception): pass

def parse_args():
    parser = argparse.ArgumentParser(description="""Create published files from config file""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--in-fasta', dest = 'in_fasta', required=False, help='FASTA')
    parser.add_argument('--min-metadata', dest = 'min_metadata', required=False, help='Metadata columns for all outputs')
    parser.add_argument('--full-metadata', dest = 'full_metadata', required=False, help='MASSIVE CSV')

    parser.add_argument('--variants', dest = 'variants', required=False, help='Mutations CSV')
    parser.add_argument('--seed', dest = 'seed', required=False, help='If anonymize, use this random seed')

    parser.add_argument('--newick-tree', dest = 'newick_tree', required=False, help='Newick tree')
    parser.add_argument('--nexus-tree', dest = 'nexus_tree', required=False, help='Nexus tree')

    parser.add_argument('--recipes', dest = 'recipes', required=True, help='JSON of recipes')
    parser.add_argument('--date', dest = 'date', required=True, help='Datestamp for published files')

    args = parser.parse_args()
    return args

#"metadata_fields": []
#"mutations": True or False to add columns from mutations
#"where": free text to be passed to fastafunk fetch --where-column
#"suffix": something to append to file names
#"exclude_uk": True or False to exclude samples from UK
#"uk_only": True or False to include only samples from UK
#"tree": "newick" or "nexus"
#"anonymize": True or False to anonymize COG sequences e.g. for microreact
#"seed": int, seed to use for anonymizing

def get_info_from_config(config_dict, outdir, date, in_fasta, min_csv, full_csv, var_csv, tree_dict):
    info_dict = {"suffix":None, "metadata_fields":None, "where": None,
                 "mutations":False, "exclude_uk":False, "tree":None, "fasta": None,
                 "anonymize":False, "date": date, "data": "cog_global",
                 "in_fa":in_fasta, "min_csv":min_csv, "full_csv":full_csv, "in_var":None, "in_tree":None,
                 "out_fasta":None, "out_csv":"tmp.csv", "out_var":"tmp_var.csv", "out_anon":None, "out_tree":None}
    info_dict.update(config_dict)

    if info_dict["tree"] in tree_dict.keys():
        info_dict["in_tree"] = tree_dict[info_dict["tree"]]

    start = "%s/%s_%s" %(outdir, info_dict["data"], info_dict["date"])
    tree_start = "%s_tree" %start
    fasta_start = "%s_tree" %start

    if info_dict["fasta"] or info_dict["tree"]:
        start += "_metadata"

    if info_dict["suffix"]:
        start += "_%s" %info_dict["suffix"]
        tree_start += "_%s" %info_dict["suffix"]
        fasta_start += "_%s" %info_dict["suffix"]


    if info_dict["tree"] in tree_dict.keys():
        info_dict["out_tree"] = "%s.%s" %(tree_start, info_dict["tree"])

    if info_dict["fasta"]:
        info_dict["out_fasta"] = "%s.fasta" %fasta_start

    csv_end = ".csv"
    if info_dict["anonymize"]:
        info_dict["out_anon"] = "%s%s" %(start, csv_end)
    elif info_dict["mutations"]:
        info_dict["out_var"] = "%s%s" %(start, csv_end)
    else:
        info_dict["out_csv"] = "%s%s" %(start, csv_end)

    if info_dict["mutations"] and info_dict["in_var"] is None:
        sys.exit("Please provide the appropriate mutations file")
    print(info_dict)
    return info_dict

def syscall(cmd_list, allow_fail=False):
    if None in cmd_list:
        print('None in list', cmd_list, file=sys.stderr)
        raise Error('Error in command. Cannot continue')
    command = ' '.join(cmd_list)
    print(command)
    completed_process = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    if (not allow_fail) and completed_process.returncode != 0:
        print('Error running this command:', command, file=sys.stderr)
        print('Return code:', completed_process.returncode, file=sys.stderr)
        print('\nOutput from stdout:', completed_process.stdout, sep='\n', file=sys.stderr)
        print('\nOutput from stderr:', completed_process.stderr, sep='\n', file=sys.stderr)
        raise Error('Error in system call. Cannot continue')
    print(completed_process.stdout)
    return completed_process

def publish_file(outdir, info_dict, seed):
    if info_dict["tree"] is not None and not info_dict["anonymize"]:
        cmd_list = ["cp", info_dict["in_tree"], info_dict["out_tree"]]
        syscall(cmd_list)
        return

    if info_dict["exclude_uk"]:
        cmd_list = ["head -n1", info_dict["min_csv"], "> no_uk.csv"]
        syscall(cmd_list)
        cmd_list = ["tail -n+2", info_dict["min_csv"], "| grep -v -E \"^England|^Northern_Ireland|^Wales|^Scotland\"", ">> no_uk.csv"]
        syscall(cmd_list)
        info_dict["min_csv"] = "no_uk.csv"

    if info_dict["uk_only"]:
        cmd_list = ["head -n1", info_dict["min_csv"], "> uk_only.csv"]
        syscall(cmd_list)
        cmd_list = ["tail -n+2", info_dict["min_csv"], "| grep -E \"^England|^Northern_Ireland|^Wales|^Scotland\"", ">> uk_only.csv"]
        syscall(cmd_list)
        info_dict["min_csv"] = "uk_only.csv"

    if info_dict["out_fasta"] is not None:
        cmd_list = ["fastafunk fetch --in-fasta", info_dict["in_fa"], "--in-metadata", info_dict["full_csv"],
                  "--index-column sequence_name --out-fasta", info_dict["out_fa"],
                  "--out-metadata", info_dict["out_csv"], "--restrict --low-memory --keep-omit-rows"]
        if info_dict["metadata_fields"]:
                cmd_list.append("--filter-column")
                cmd_list.extend(info_dict["metadata_fields"])
        if info_dict["where"]:
            cmd_list.append("--where-column %s" %info_dict["where"])
        syscall(cmd_list)
    else:
        cmd_list = ["fastafunk add_columns --in-metadata", info_dict["min_csv"],
            "--in-data", info_dict["full_csv"], "--index-column sequence_name",
            "--join-on sequence_name --out-metadata", info_dict["out_csv"]]
        syscall(cmd_list)
        if info_dict["metadata_fields"]:
                cmd_list.append("--new-columns")
                cmd_list.extend(info_dict["metadata_fields"])
        if info_dict["where"]:
            cmd_list.append("--where-column %s" %info_dict["where"])
        syscall(cmd_list)

    if info_dict["mutations"]:
        cmd_list = ["fastafunk add_columns --in-metadata", info_dict["out_csv"],
        "--in-data", info_dict["in_var"], "--index-column sequence_name",
        "--join-on query --out-metadata", info_dict["out_var"]]
        syscall(cmd_list)

    if info_dict["anonymize"]:
        if info_dict["mutations"]:
            in_csv = info_dict["out_var"]
        else:
            in_csv = info_dict["out_csv"]
        anonymize_microreact(metadata_in = in_csv,
                             tree_in = info_dict["in_tree"],
                             metadata_out = info_dict["out_anon"],
                             tree_out = info_dict["out_tree"],
                             seed = seed)

    #tmp = glob.glob("tmp.*")
    #if len(tmp) > 0:
    #    cmd_list = ["rm tmp.*"]
    #    syscall(cmd_list)

def main():
    args = parse_args()
    print(args)
    if args.seed:
        random.seed(args.seed)

    tree_dict = {"newick":args.newick_tree, "nexus":args.nexus_tree}
    print(tree_dict)

    recipes = {}
    with open(args.recipes, 'r') as f:
        recipes = json.load(f)

    for outdir in recipes.keys():
        os.makedirs(outdir,exist_ok=True)
        for recipe in recipes[outdir]:
            info_dict = get_info_from_config(recipe, outdir, args.date, args.fasta, args.min_metadata, args.full_metadata, args.variants, tree_dict)
            publish_file(outdir, info_dict, args.seed)

if __name__ == '__main__':
    main()
