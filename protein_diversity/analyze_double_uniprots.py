from pathlib import Path
import pandas as pd
from argparse import ArgumentParser
import json
import numpy as np


arg_parser = ArgumentParser(
    description="Identify double uniprot IDs for the LOBSTER dataset"
)
arg_parser.add_argument(
    "--config",
    type=str,
    required=True,
    help="Path to config in JSON format including all file paths.",
)
args = arg_parser.parse_args()

config = json.load(open(args.config, "r"))

lobster_data = pd.read_csv(Path(config["broad_cluster_csv"]))
n_total = 0
n_double_uniprot = 0
n_double_famids = 0
reverse_pfam = {}
n_failed_work = 0
filtered_data = []
failed_data = []
columns = lobster_data.columns.values
for index, row in lobster_data.iterrows():
    n_total += 1
    Uniprotids = row["Uniprotids"].split(";")
    pfamids = row["PfamID"].split(";")
    pfamnames = row["PfamName"].split(";")
    if len(pfamnames) == len(pfamids):
        for i in range(len(pfamnames)):
            if pfamnames[i] in reverse_pfam.keys():
                reverse_pfam[pfamnames[i]].add(pfamids[i])
            else:
                reverse_pfam[pfamnames[i]] = {pfamids[i]}
    else:
        raise RuntimeError
    if len(Uniprotids) == 1 and Uniprotids[0] == "failed":
        n_failed_work += 1
        failed_data.append([x for x in row])
    if len([x for x in Uniprotids if x != "failed"]) > 1:
        n_double_uniprot += 1
        filtered_data.append([x for x in row])
    if len([x for x in pfamids if x != "failed"]) > 1:
        n_double_famids += 1

filtered_data = pd.DataFrame(filtered_data, columns=columns)
filtered_data.to_csv(config["double_uniprot_csv"])
failed_data = pd.DataFrame(failed_data, columns=columns)
failed_data.to_csv(config["failing_uniprot_csv"])
print(n_total, n_double_uniprot, n_double_famids, n_failed_work)


reverse_pfam = {}
lobster_data = pd.read_csv(Path(config["exact_cluster_csv"]))
for index, row in lobster_data.iterrows():
    pfamid = row["PfamID"]
    pfamname = row["Pfam Name"]
    if pfamname in reverse_pfam.keys():
        reverse_pfam[pfamname].add(pfamid)
    else:
        reverse_pfam[pfamname] = {pfamid}

for entity in reverse_pfam:
    if len(reverse_pfam[entity]) > 1:
        print(entity, reverse_pfam[entity])
