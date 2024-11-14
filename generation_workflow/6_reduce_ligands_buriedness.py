import pandas as pd
import progressbar, json
from collections import defaultdict
import argparse
from pathlib import Path


arg_parser = argparse.ArgumentParser(description="Filtering for buriedness.")
arg_parser.add_argument(
    "--config",
    type=str,
    required=True,
    help="Path to config in JSON format including all file paths.",
)
args = arg_parser.parse_args()

config = json.load(open(args.config, "r"))


if not Path(config["3_activities_csv"]).is_file():
    print("Input file initial_csv specified in config does not exist")
    raise FileNotFoundError(config["3_activities_csv"])

# create result dir if necessary
if not Path(config["6_filtered_buriedness_csv"]).parent.is_dir():
    Path(config["6_filtered_buriedness_csv"]).parent.mkdir()

if not Path(config["6_finally_filtered_pdbs_ligands_json"]).parent.is_dir():
    Path(config["6_finally_filtered_pdbs_ligands_json"]).parent.mkdir()

if not Path(config["log_file"]).parent.is_dir():
    Path(config["log_file"]).parent.mkdir()


# files to handle
master_table_path = Path(config["3_activities_csv"])
buriedness_file_path = Path(config["5_buriednes_csv"])
new_master_table_path = Path(config["6_filtered_buriedness_csv"])
pdb_file_path = Path(config["6_finally_filtered_pdbs_ligands_json"])
master_table_log_path = Path(config["log_file"])


# build data frames
df = pd.read_csv(master_table_path, delimiter=",")
buriedness_df = pd.read_csv(buriedness_file_path, delimiter=";")

# variables
drop_idx = []
buriedness_too_low = 0
nof_ligands = 0
pdbs = defaultdict(list)

for idx, row in progressbar.progressbar(df.iterrows()):
    ligand = row["name"]
    pdb = row["pdb"]
    buriedness_row = buriedness_df.where(
        (buriedness_df["pdb"] == pdb) & (buriedness_df["name"] == ligand)
    ).dropna(subset=["pdb"])
    buriedness = float(buriedness_row["buriedness"])
    if buriedness < 0.5:
        buriedness_too_low += 1
        drop_idx.append(idx)
    else:
        nof_ligands += 1
        pdbs[pdb].append(ligand)
df = df.drop(drop_idx)
df.to_csv(new_master_table_path, mode="w")
log = open(master_table_log_path, "a")  # Important: Append existing log!
log.write(f"Removed by burriedness: {buriedness_too_low}\n")
log.write(f"Number of remaining ligands: {nof_ligands}\n")
log.write(f"Number of remaining pdbs: {len(pdbs)}\n")

pdb_file = open(pdb_file_path, "w")
json.dump(pdbs, pdb_file)
pdb_file.close()
