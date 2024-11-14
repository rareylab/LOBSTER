# NOTE: Each ligand only is counted once in the statistics,
# even if it occurs in multiple ensembles

import pandas as pd
import json
from matplotlib import pyplot as plt
from rdkit import Chem
import progressbar
from argparse import ArgumentParser
from pathlib import Path


FONT_SIZE = 14


arg_parser = ArgumentParser(
    description="Create statistics on physchem properties of the molecules"
    " within the LOBSTER dataset"
)
arg_parser.add_argument(
    "--config",
    type=str,
    required=True,
    help="Path to config in JSON format including all file paths.",
)
args = arg_parser.parse_args()

config = json.load(open(args.config, "r"))

final_dataset_dir = Path(config["10_final_dataset_dir"])
if not final_dataset_dir.is_dir():
    print("Input directory 10_final_dataset_dir as specified in config does not exist")
    raise IOError(str(final_dataset_dir))

representatives_path = final_dataset_dir / "representatives"
clustering_log = final_dataset_dir / "clusters_log.txt"

table_path = Path(config["6_filtered_buriedness_csv"])

if not table_path.is_file():
    print("Input file 6_filtered_buriedness_csv as specified in config does not exist")
    raise FileNotFoundError(str(table_path))


filtered_master_table = pd.read_csv(table_path, dtype=str)
plot_dir = Path(config["plots_dir"])
if not plot_dir.is_dir():
    plot_dir.mkdir()

nof_ligands = 0
unique_ligands = set()
rotatable_bonds = []
mw = []
logP = []
hbd = []
hba = []
unique = []
all_pdbs = set()

pdbs_log_file = Path(config["stats_all_pdbs_in_LOBSTER"])


for file in progressbar.progressbar(representatives_path.iterdir()):
    molecules = Chem.SDMolSupplier(str(file))
    for mol in molecules:
        nof_ligands += 1
        ligand = mol.GetProp("_Name")
        pdb = mol.GetProp("pdb_id")
        all_pdbs.add(pdb)
        table_row = filtered_master_table.where(
            (filtered_master_table["pdb"] == pdb)
            & (filtered_master_table["name"] == ligand)
        ).dropna(subset=["pdb"])
        usmiles = table_row["usmiles"].values[0]
        if not (ligand, pdb) in unique:
            unique.append((ligand, pdb))
        else:
            print("WARNING: Ensembles are not disjoint!")
            print(ligand, pdb)
        if not usmiles in unique_ligands:
            unique_ligands.add(usmiles)
            rotatable_bonds.append(float(table_row["nrot_simple"]))
            mw.append(float(table_row["mw"]))
            logP.append(float(table_row["logP"]))
            hba.append(float(table_row["lipinskiAcceptors"]))
            hbd.append(float(table_row["lipinskiDonors"]))

plt.rc("legend", fontsize=FONT_SIZE)
plt.rc("axes", labelsize=FONT_SIZE)

rot_bonds_bars = []
rotatable_bonds.sort()
for i in range(0, int(rotatable_bonds[-1]) + 1):
    rot_bonds_bars.append(rotatable_bonds.count(i))
plt.bar(
    range(0, int(rotatable_bonds[-1]) + 1),
    rot_bonds_bars,
    color="b",
    edgecolor="black",
    linewidth=1.2,
    alpha=0.5,
)
plt.xlabel("Number of Rotatable Bonds")
plt.ylabel("Number of Ligands")

plt.savefig(str(plot_dir / "rotBonds.png"), dpi=300)
plt.clf()

plt.hist(mw, bins=15, color="b", edgecolor="black", linewidth=1.2, alpha=0.5)
plt.xlabel("Molecular Weight [Da]")
plt.ylabel("Number of Ligands")
plt.savefig(str(plot_dir / "mw.png"), dpi=300)
plt.clf()

plt.hist(logP, bins=15, color="b", edgecolor="black", linewidth=1.2, alpha=0.5)
plt.xlabel("cLogP")
plt.ylabel("Number of Ligands")
plt.savefig(str(plot_dir / "clogP.png"), dpi=300)
plt.clf()

hba_bars = []
hbd_bars = []

hba.sort()
hbd.sort()

for i in range(0, int(hba[-1]) + 1):
    hba_bars.append(hba.count(i))

for i in range(0, int(hbd[-1]) + 1):
    hbd_bars.append(hbd.count(i))

plt.bar(
    range(0, int(hba[-1]) + 1),
    hba_bars,
    color="b",
    edgecolor="black",
    linewidth=1.2,
    label="HBA",
    alpha=0.5,
)
plt.bar(
    range(0, int(hbd[-1]) + 1),
    hbd_bars,
    color="r",
    edgecolor="black",
    linewidth=1.2,
    label="HBD",
    alpha=0.5,
)
plt.xticks(range(0, int(max(hbd[-1] + 1, hba[-1] + 1))))
plt.xlabel("Number of H-bonding Atoms")
plt.ylabel("Number of Ligands")
plt.legend()
plt.savefig(str(plot_dir / "hba_hbd.png"), dpi=300)

print("Number of Ligands:", nof_ligands)
print("Number of Unique Ligands:", len(unique_ligands))
print("Number of PDBs:", len(all_pdbs))

with pdbs_log_file.open("w") as f:
    for pdb in all_pdbs:
        f.write(f"{pdb}\n")
