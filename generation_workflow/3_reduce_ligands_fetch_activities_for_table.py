import math, sys, argparse, json
from pathlib import Path
import pandas as pd
from progressbar import progressbar
from collections import defaultdict

# This script was implemented to gather activity data from PDB for a "master table"
arg_parser = argparse.ArgumentParser(description="Filtering for ligand efficiency")
arg_parser.add_argument(
    "--config",
    type=str,
    required=True,
    help="Path to config in JSON format including all file paths.",
)
args = arg_parser.parse_args()

config = json.load(open(args.config, "r"))

if not Path(config["2_no_duplicates_csv"]).is_file():
    print("Input file initial_csv specified in config does not exist")
    raise FileNotFoundError(config["2_no_duplicates_csv"])

# create result dir if necessary
if not Path(config["3_activities_csv"]).parent.is_dir():
    Path(config["3_activities_csv"]).parent.mkdir()

if not Path(config["log_file"]).parent.is_dir():
    Path(config["log_file"]).parent.mkdir()

pdbBind = pd.read_csv(
    config["PDBbind_file"],
    dtype=str,
    comment="#",
    delimiter=r"\s+",
    names=[
        "PDB code",
        "resolution",
        "release year",
        "-logKd/Ki",
        "Kd/Ki",
        "?",
        "reference",
        "ligand name",
    ],
)

bindingDB = pd.read_csv(
    config["BindingDB_file"],
    sep="\t",
    dtype=str,
    usecols=[
        "Ki (nM)",
        "IC50 (nM)",
        "Kd (nM)",
        "EC50 (nM)",
        "Ligand HET ID in PDB",
        "PDB ID(s) for Ligand-Target Complex",
    ],
)
print("Finished loading BindingDB")

drop_idx = []  # ligands without activity or too low LE
no_activity = 0
too_low_le = 0

# Update data frame for master table
df = pd.read_csv(config["2_no_duplicates_csv"], delimiter=",", dtype=str)
df.insert(len(df.columns), "PDBbind activity", "nan")
df.insert(len(df.columns), "BindingDB values (nM)", "nan")
df.insert(len(df.columns), "BindingDB units", "nan")

# Prepare BindingDB in a dictionary for fast access
bindingDB_dict = defaultdict(
    lambda: {"values": [], "unaltered_values": [], "units": []}
)
bindingDB.dropna(
    subset=["PDB ID(s) for Ligand-Target Complex", "Ligand HET ID in PDB"], inplace=True
)
print("Rows with PDB ID(s) and HET codes remaining:", len(bindingDB))
print("Preparing BindingDB")

for idx, row in progressbar(bindingDB.iterrows()):
    pdb_codes = row["PDB ID(s) for Ligand-Target Complex"].split(",")
    ligands = row["Ligand HET ID in PDB"].split(";")
    unit = []
    value = []
    unaltered_value = []

    # This should not be the case...
    if len(ligands) > 1:
        print("More than one HET code detected:", pdb_code, ligand, "in row", idx)
        sys.exit(1)
    ligand = ligands[0]

    ki = str(row["Ki (nM)"]).strip().replace("<", "").replace("<=", "")
    if ki != "nan" and not ">" in ki and not "~" in ki:
        value.append(ki)
        unit.append("Ki")
        unaltered_value.append(str(row["Ki (nM)"]).strip())

    ic50 = str(row["IC50 (nM)"]).strip().replace("<", "").replace("<=", "")
    if ic50 != "nan" and not ">" in ic50 and not "~" in ic50:
        value.append(ic50)
        unit.append("IC50")
        unaltered_value.append(str(row["IC50 (nM)"]).strip())

    kd = str(row["Kd (nM)"]).strip().replace("<", "").replace("<=", "")
    if kd != "nan" and not ">" in kd and not "~" in kd:
        value.append(kd)
        unit.append("Kd")
        unaltered_value.append(str(row["Kd (nM)"]).strip())

    ec50 = str(row["EC50 (nM)"]).strip().replace("<", "").replace("<=", "")
    if ec50 != "nan" and not ">" in ec50 and not "~" in ec50:
        value.append(ec50)
        unit.append("EC50")
        unaltered_value.append(str(row["EC50 (nM)"]).strip())

    if len(value) != 0:
        for pdb_code in pdb_codes:
            try:
                bindingDB_dict[(pdb_code, ligand)]["values"].extend(value)
                bindingDB_dict[(pdb_code, ligand)]["unaltered_values"].extend(
                    unaltered_value
                )
                bindingDB_dict[(pdb_code, ligand)]["units"].extend(unit)
            except:
                raise ValueError(
                    f"Error for value {value}, unit {unit}, pdb {pdb_code},"
                    f" ligand {ligand}."
                )
# Compute min norm value
for key in bindingDB_dict.keys():
    # math.pow(10,9) bc of nM conversion
    bindingDB_dict[key]["min_norm"] = min(
        [
            -math.log10(float(x) / math.pow(10, 9)) if float(x) != 0 else 1000000
            for x in bindingDB_dict[key]["values"]
        ]
    )

# prepare PDBbind as a dictionary
pdb_bind_dict = {}
for idx, row in pdbBind.iterrows():
    pdb_code = row["PDB code"].upper()
    ligands = row["ligand name"].strip().replace("(", "").replace(")", "").split("&")
    for ligand in ligands:
        pdb_bind_dict[(pdb_code, ligand)] = {
            "normalized": row["-logKd/Ki"],
            "everything": row["Kd/Ki"],
        }

# iterate, filter, and annotate master table with activity values
for idx, row in progressbar(df.iterrows()):
    pdb_code = row["pdb"]
    wanted_ligand = row["hetcodes"]
    if (
        ((pdb_code, wanted_ligand) in pdb_bind_dict.keys())
        and not ">" in pdb_bind_dict[(pdb_code, wanted_ligand)]["everything"]
        and not "~" in pdb_bind_dict[(pdb_code, wanted_ligand)]["everything"]
    ):
        pdbBind_values = pdb_bind_dict[(pdb_code, wanted_ligand)]
        normalized_pdb = float(pdbBind_values["normalized"])
        if (pdb_code, wanted_ligand) in bindingDB_dict.keys():
            bindingDB_values = bindingDB_dict[(pdb_code, wanted_ligand)]
            normalized_min = min([normalized_pdb, bindingDB_values["min_norm"]])
            # why is this *1.37 magic value here?
            # See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7610969/
            if normalized_min * 1.37 / float(row["heavy_atoms"]) > 0.3:
                df.at[idx, "PDBbind activity"] = pdbBind_values["everything"]
                df.at[idx, "BindingDB values (nM)"] = bindingDB_values[
                    "unaltered_values"
                ]
                df.at[idx, "BindingDB units"] = bindingDB_values["units"]
            else:
                too_low_le += 1
                drop_idx.append(idx)
                continue
        elif normalized_pdb * 1.37 / float(row["heavy_atoms"]) > 0.3:
            df.at[idx, "PDBbind activity"] = pdbBind_values["everything"]
        else:
            too_low_le += 1
            drop_idx.append(idx)
            continue
    elif (pdb_code, wanted_ligand) in bindingDB_dict.keys():
        bindingDB_values = bindingDB_dict[(pdb_code, wanted_ligand)]
        bindingDB_min = bindingDB_values["min_norm"]
        if bindingDB_min * 1.37 / float(row["heavy_atoms"]) > 0.3:
            df.at[idx, "BindingDB values (nM)"] = bindingDB_values["unaltered_values"]
            df.at[idx, "BindingDB units"] = bindingDB_values["units"]
        else:
            too_low_le += 1
            drop_idx.append(idx)
            continue
    else:
        no_activity += 1
        drop_idx.append(idx)
        continue

df = df.drop(drop_idx)
df.to_csv(config["3_activities_csv"], mode="w")
with Path(config["log_file"]).open("a") as f:
    f.write(f"No activity: {no_activity}\n")
    f.write(f"Too low LE: {too_low_le}\n")
    f.write(f"Remaining: {len(df)}\n")
