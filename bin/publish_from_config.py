#!/usr/bin/env python3

import argparse
import json
import subprocess
import os
import sys
import glob

class Error (Exception): pass

def parse_args():
    parser = argparse.ArgumentParser(description="""Create published files from config file""",
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--unaligned_fasta', dest = 'unaligned_fasta', required=False, help='Raw FASTA')
    parser.add_argument('--aligned_fasta', dest = 'aligned_fasta', required=False, help='Aligned, masked, untrimmed FASTA')
    parser.add_argument('--trimmed_fasta', dest = 'trimmed_fasta', required=False, help='Aligned, masked, trimmed and filtered FASTA')
    parser.add_argument('--gisaid_fasta', dest = 'global_fasta', required=False, help='GISAID aligned FASTA')
    parser.add_argument('--cog_global_fasta', dest = 'cog_global_fasta', required=False, help='COG GISAID aligned FASTA')

    parser.add_argument('--cog_metadata', dest = 'cog_metadata', required=False, help='MASSIVE CSV')
    parser.add_argument('--gisaid_metadata', dest = 'global_metadata', required=False, help='MASSIVE CSV')
    parser.add_argument('--cog_global_metadata', dest = 'cog_global_metadata', required=False, help='MASSIVE CSV')

    parser.add_argument('--cog_variants', dest = 'cog_variants', required=False, help='Mutations CSV')
    parser.add_argument('--gisaid_variants', dest = 'global_variants', required=False, help='Mutations CSV')
    parser.add_argument('--cog_global_variants', dest = 'cog_global_variants', required=False, help='Mutations CSV')

    parser.add_argument('--recipes', dest = 'recipes', required=True, help='JSON of recipes')
    parser.add_argument('--date', dest = 'date', required=True, help='Datestamp for published files')

    args = parser.parse_args()
    return args

#"data": "cog", "gisaid" or "cog_global"
#"fasta": "unaligned", "aligned", "trimmed", "cog_global" or "gisaid"
#"metadata_fields": []
#"mutations": True or False to add columns from mutations
#"where": free text to be passed to fastafunk fetch --where-column
#"suffix": something to append to file names
#"exclude_uk": True or False to exclude samples from UK

def get_info_from_config(config_dict, outdir, date, fasta_dict, csv_dict, var_dict):
    info_dict = {"suffix":None, "data":None, "fasta":None, "metadata_fields":None,
                 "where": None, "mutations":False, "exclude_uk":False, "date": date,
                 "in_fa":None, "in_csv":None, "in_var":None,
                 "out_fa":"tmp.fa", "out_csv":"tmp.csv", "out_var":None}
    info_dict.update(config_dict)

    if info_dict["fasta"] in fasta_dict.keys():
        info_dict["in_fa"] = fasta_dict[info_dict["fasta"]]
    elif info_dict["data"] == "cog_global":
        info_dict["in_fa"] = fasta_dict["cog_global"]
    elif info_dict["data"] == "gisaid":
        info_dict["in_fa"] = fasta_dict["gisaid"]
    elif info_dict["data"] == "cog":
            info_dict["in_fa"] = fasta_dict["trimmed"]
    else:
        sys.exit("Config entries need to specify either fasta in ['unaligned', 'aligned', 'trimmed', 'cog_global', 'gisaid'] or data \
        in ['cog', 'cog_global', 'gisaid']")

    if info_dict["data"] is None:
        if info_dict["fasta"] == "cog_global":
            info_dict["data"] = "cog_global"
        elif info_dict["fasta"] == "gisaid":
            info_dict["data"] = "gisaid"
        else:
            info_dict["data"] = "cog"

    if info_dict["data"] == "cog_global":
            info_dict["in_csv"] = csv_dict["cog_global"]
            info_dict["in_var"] = var_dict["cog_global"]
    elif info_dict["data"] == "cog":
            info_dict["in_csv"] = csv_dict["cog"]
            info_dict["in_var"] = var_dict["cog"]
    elif info_dict["data"] == "gisaid":
            info_dict["in_csv"] = csv_dict["gisaid"]
            info_dict["in_var"] = var_dict["gisaid"]

    start = "%s/%s_%s" %(outdir, info_dict["data"], info_dict["date"])
    if info_dict["suffix"]:
        start += "_%s" %info_dict["suffix"]
    csv_end = ".csv"

    if info_dict["fasta"]:
        csv_end = "_metadata.csv"
        if info_dict["fasta"]=="aligned" or (info_dict["metadata_fields"] and info_dict["fasta"]!="unaligned"):
            info_dict["out_fa"] = "%s_alignment.fa" %start
        else:
            info_dict["out_fa"] = "%s.fa" %start

    if info_dict["mutations"]:
        info_dict["out_var"] = "%s%s" %(start, csv_end)
    else:
        info_dict["out_csv"] = "%s%s" %(start, csv_end)

    if info_dict["out_fa"] != "tmp.fa" and info_dict["in_fa"] is None:
        sys.exit("Please provide the appropriate FASTA file")
    if (info_dict["out_csv"] != "tmp.csv" or info_dict["out_var"] is not None) and info_dict["in_csv"] is None:
        sys.exit("Please provide the appropriate CSV file")
    if info_dict["out_var"] is not None and info_dict["in_var"] is None:
        sys.exit("Please provide the appropriate mutations file")
    print(info_dict)
    return info_dict

def syscall(cmd_list, allow_fail=False):
    if None in cmd_list:
        print('None in list', cmd_list, file=sys.stderr)
        raise Error('Error in command. Cannot continue')
    command = ' '.join(cmd_list)
    completed_process = subprocess.run(command, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
    if (not allow_fail) and completed_process.returncode != 0:
        print('Error running this command:', command, file=sys.stderr)
        print('Return code:', completed_process.returncode, file=sys.stderr)
        print('\nOutput from stdout:', completed_process.stdout, sep='\n', file=sys.stderr)
        print('\nOutput from stderr:', completed_process.stderr, sep='\n', file=sys.stderr)
        raise Error('Error in system call. Cannot continue')
    print(completed_process.stdout)
    return completed_process

def publish_file(outdir, info_dict):
    if info_dict["metadata_fields"] is None:
        cmd_list = ["cp", info_dict["in_fa"], info_dict["out_fa"]]
        syscall(cmd_list)
        return

    if info_dict["exclude_uk"]:
        cmd_list = ["head -n1", info_dict["in_csv"], "> no_uk.csv"]
        syscall(cmd_list)
        cmd_list = ["tail -n+2", info_dict["in_csv"], "| grep -v -E \"^England|^Northern_Ireland|^Wales|^Scotland\"", ">> no_uk.csv"]
        syscall(cmd_list)
        info_dict["in_csv"] = "no_uk.csv"

    cmd_list = ["fastafunk fetch --in-fasta", info_dict["in_fa"], "--in-metadata", info_dict["in_csv"],
              "--index-column sequence_name --out-fasta", info_dict["out_fa"],
              "--out-metadata", info_dict["out_csv"], "--restrict"]
    if info_dict["metadata_fields"]:
            cmd_list.append("--filter-column")
            cmd_list.extend(info_dict["metadata_fields"])
    if info_dict["where"]:
        cmd_list.append("--where-column %s" %info_dict["where"])
    syscall(cmd_list)

    if info_dict["mutations"]:
        cmd_list = ["fastafunk add_columns --in-metadata", info_dict["out_csv"],
        "--in-data", info_dict["in_var"], "--index-column sequence_name",
        "--join-on query --out-metadata", info_dict["out_var"]]
        syscall(cmd_list)

    tmp = glob.glob("tmp.*")
    if len(tmp) > 0:
        cmd_list = ["rm tmp.*"]
        syscall(cmd_list)

def main():
    args = parse_args()
    print(args)
    fasta_dict = {"unaligned":args.unaligned_fasta, "aligned":args.aligned_fasta, "trimmed":args.trimmed_fasta, "cog_global": args.cog_global_fasta, "gisaid": args.global_fasta}
    print(fasta_dict)
    csv_dict = {"cog":args.cog_metadata, "cog_global":args.cog_global_metadata, "gisaid": args.global_metadata}
    print(csv_dict)
    var_dict = {"cog":args.cog_variants, "cog_global":args.cog_global_variants, "gisaid": args.global_variants}
    print(var_dict)

    recipes = {}
    with open(args.recipes, 'r') as f:
        recipes = json.load(f)

    for outdir in recipes.keys():
        os.makedirs(outdir,exist_ok=True)
        for recipe in recipes[outdir]:
            info_dict = get_info_from_config(recipe, outdir, args.date, fasta_dict, csv_dict, var_dict)
            publish_file(outdir, info_dict)

if __name__ == '__main__':
    main()
