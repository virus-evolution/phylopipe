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
    parser.add_argument('--in-metadata', dest = 'full_metadata', required=False, help='MASSIVE CSV')

    parser.add_argument('--mutations', dest = 'mutations', required=False, help='Mutations CSV')
    parser.add_argument('--constellations', dest = 'constellations', required=False, help='Constellations CSV')
    parser.add_argument('--seed', dest = 'seed', required=False, help='If anonymize, use this random seed')

    parser.add_argument('--newick-tree', dest = 'newick_tree', required=False, help='Newick tree')
    parser.add_argument('--pruned-newick-tree', dest = 'pruned_newick_tree', required=False, help='Newick tree with fewer tips')
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
#"tree": "newick" or "nexus" or "pruned_newick"
#"anonymize": True or False to anonymize COG sequences e.g. for microreact
#"seed": int, seed to use for anonymizing

def get_info_from_config(config_dict, outdir, date, in_fasta, full_csv, muts_csv, con_csv, tree_dict):
    info_dict = {"suffix":None, "metadata_fields":None, "where": None,
                 "mutations":False, "constellations":False,
                 "exclude_uk":False, "uk_only": False, "exclude_cog":False, "cog_only": False,
                 "tree":None, "fasta": None,
                 "anonymize":False, "drop_index": False, "date": date, "data": "cog_global",
                 "in_fa":in_fasta, "full_csv":full_csv, "in_muts":muts_csv, "in_con":con_csv, "in_tree":None,
                 "out_fa":"tmp.fa", "intermediate_csv":"tmp.csv", "out_csv":"tmp.csv", "out_tree":None}
    info_dict.update(config_dict)

    if info_dict["tree"] in tree_dict.keys():
        info_dict["in_tree"] = tree_dict[info_dict["tree"]]

    start = "%s/%s_%s" %(outdir, info_dict["data"], info_dict["date"])
    tree_start = "%s_tree" %start
    fasta_start = start

    if info_dict["fasta"] or info_dict["tree"]:
        start += "_metadata"

    if info_dict["suffix"]:
        start += "_%s" %info_dict["suffix"]
        tree_start += "_%s" %info_dict["suffix"]
        fasta_start += "_%s" %info_dict["suffix"]


    if info_dict["tree"] in tree_dict.keys():
        info_dict["out_tree"] = "%s.%s" %(tree_start, info_dict["tree"])
        if "pruned" in info_dict["tree"]:
            info_dict["out_tree"] = info_dict["out_tree"].replace("pruned_", "pruned.")


    if info_dict["fasta"]:
        info_dict["out_fa"] = "%s.fasta" %fasta_start

    csv_end = ".csv"
    info_dict["out_csv"] = "%s%s" %(start, csv_end)

    if info_dict["mutations"] and info_dict["in_muts"] is None:
        sys.exit("Please provide the appropriate mutations file")
    if info_dict["constellations"] and info_dict["in_con"] is None:
            sys.exit("Please provide the appropriate constellations file")
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

    if info_dict["out_fa"] is not "tmp.fa" and info_dict["metadata_fields"] is None:
        cmd_list = ["cp", info_dict["in_fa"], info_dict["out_fa"]]
        syscall(cmd_list)

    if info_dict["metadata_fields"] is None:
        return

    if info_dict["exclude_uk"]:
        cmd_list = ["fastafunk filter_column --in-metadata", info_dict["in_csv"],
        "--out-metadata tmp.no_uk.csv --column is_uk --is_true"]
        syscall(cmd_list)
        info_dict["in_csv"] = "tmp.no_uk.csv"

    if info_dict["exclude_cog"]:
        cmd_list = ["fastafunk filter_column --in-metadata", info_dict["in_csv"],
        "--out-metadata tmp.no_cog.csv --column is_cog_uk --is_true"]
        syscall(cmd_list)
        info_dict["in_csv"] = "tmp.no_cog.csv"

    if info_dict["uk_only"]:
        cmd_list = ["fastafunk filter_column --in-metadata", info_dict["in_csv"],
        "--out-metadata tmp.uk_only.csv --column is_uk --is_false"]
        syscall(cmd_list)
        info_dict["in_csv"] = "tmp.uk_only.csv"

    if info_dict["cog_only"]:
        cmd_list = ["fastafunk filter_column --in-metadata", info_dict["in_csv"],
        "--out-metadata tmp.cog_only.csv --column is_cog_uk --is_false"]
        syscall(cmd_list)
        info_dict["in_csv"] = "tmp.cog_only.csv"

    cmd_list = ["fastafunk fetch --in-fasta", info_dict["in_fa"], "--in-metadata", info_dict["full_csv"],
                "--index-column sequence_name --out-fasta", info_dict["out_fa"],
                "--out-metadata", info_dict["intermediate_csv"], "--low-memory --keep-omit-rows"]
    if info_dict["metadata_fields"]:
        cmd_list.append("--filter-column")
        cmd_list.extend(info_dict["metadata_fields"])
    if info_dict["where"]:
        cmd_list.append("--where-column %s" %info_dict["where"])
    syscall(cmd_list)

    if info_dict["mutations"]:
            cmd_list = ["fastafunk add_columns --in-metadata", info_dict["intermediate_csv"],
            "--in-data", info_dict["in_muts"], "--index-column sequence_name",
            "--join-on sequence_name --out-metadata tmp.muts.csv"]
            info_dict["intermediate_csv"] = "tmp.muts.csv"
            syscall(cmd_list)

    if info_dict["constellations"]:
            cmd_list = ["fastafunk add_columns --in-metadata", info_dict["intermediate_csv"],
            "--in-data", info_dict["in_con"], "--index-column sequence_name",
            "--join-on sequence_name --out-metadata tmp.constellations.csv"]
            info_dict["intermediate_csv"] = "tmp.constellations.csv"
            syscall(cmd_list)

    if info_dict["anonymize"]:
        anonymize_microreact(metadata_in = info_dict["intermediate_csv"],
                             tree_in = info_dict["in_tree"],
                             metadata_out = "tmp.anon.csv",
                             tree_out = info_dict["out_tree"],
                             seed = seed)
        info_dict["intermediate_csv"] = "tmp.anon.csv"

    if info_dict["drop_index"]:
        cmd_list = ["fastafunk drop_columns --in-metadata", info_dict["intermediate_csv"],
        "--columns", info_dict["drop_index"],
        "--out-metadata tmp.drop.csv"]
        info_dict["intermediate_csv"] = "tmp.drop.csv"
        syscall(cmd_list)

    cmd_list = ["mv", info_dict["intermediate_csv"], info_dict["out_csv"]]
    syscall(cmd_list)

    #tmp = glob.glob("tmp.*")
    #if len(tmp) > 0:
    #    cmd_list = ["rm tmp.*"]
    #    syscall(cmd_list)

def main():
    args = parse_args()
    print(args)
    if args.seed:
        random.seed(args.seed)

    tree_dict = {"newick":args.newick_tree, "nexus":args.nexus_tree, "pruned_newick": args.pruned_newick_tree}
    print(tree_dict)

    recipes = {}
    with open(args.recipes, 'r') as f:
        recipes = json.load(f)

    for outdir in recipes.keys():
        os.makedirs(outdir,exist_ok=True)
        for recipe in recipes[outdir]:
            info_dict = get_info_from_config(recipe, outdir, args.date, args.in_fasta, args.full_metadata, args.mutations, args.constellations, tree_dict)
            publish_file(outdir, info_dict, args.seed)

if __name__ == '__main__':
    main()
