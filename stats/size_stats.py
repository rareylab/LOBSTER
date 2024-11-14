# NOTE: Each ligand only is counted PER ENSEMBLE,
# no matter if it occurs in other ensembles as well.

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from rdkit import Chem
from statistics import stdev, median, mean
import progressbar, json
from pathlib import Path
from argparse import ArgumentParser


FONT_SIZE = 14

arg_parser = ArgumentParser(
    description="Create statistics on the size of the LOBSTER dataset"
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

nof_mols_per_ensemble = []
std_devs_mol_size_per_ensemble = []
ensembles_per_cluster = []

# Stats about stdev of mol size per ensemble and ensemble size
for file in progressbar.progressbar(representatives_path.iterdir()):
    molecules = Chem.SDMolSupplier(str(file))
    ensemble_mol_sizes = []
    nof_mols = 0
    for mol in molecules:
        nof_mols += 1
        ligand = mol.GetProp("_Name")
        pdb = mol.GetProp("pdb_id")
        table_row = filtered_master_table.where(
            (filtered_master_table["pdb"] == pdb)
            & (filtered_master_table["name"] == ligand)
        ).dropna(subset=["pdb"])
        usmiles = table_row["usmiles"].values[0]
        ensemble_mol_sizes.append(float(table_row["heavy_atoms"].values[0]))
    nof_mols_per_ensemble.append(nof_mols)
    std_devs_mol_size_per_ensemble.append(stdev(ensemble_mol_sizes))

# Stats for whole cluster
ensembles_in_cluster = []
clustering_content = (
    open(clustering_log, "r")
    .read()
    .split("--------------------------------------------")
)
clusters_elements_split = [x.split("\n") for x in clustering_content[:-1]]
# -2 because of trailing \n and leading Representative tag
clustering_elements_sizes = [len(x) - 2 for x in clusters_elements_split]

plt.rc("legend", fontsize=FONT_SIZE)
plt.rc("axes", labelsize=FONT_SIZE)

print(
    f"Mols per ensemble: mean {mean(nof_mols_per_ensemble)},"
    f" meadian {median(nof_mols_per_ensemble)}"
)
sns.violinplot(y=nof_mols_per_ensemble, color="#920fa3", cut=2)
plt.ylabel("Ensemble Size")
plt.xlabel("Set of Representative Ensembles")

plt.grid()
plt.savefig(str(plot_dir / "mols_per_ensemble.png"), dpi=300)
plt.grid()
plt.tight_layout()
plt.clf()

print(
    f"Standard deviation in HAC: mean {mean(std_devs_mol_size_per_ensemble)},"
    f" median {median(std_devs_mol_size_per_ensemble)}"
)
sns.violinplot(
    y=std_devs_mol_size_per_ensemble, color="#c23c81", cut=0
)  # , color="b", alpha=0.5, edgecolor="black")
plt.ylabel("Standard Deviation in HAC")
plt.xlabel("Set of Representative Ensembles")

plt.tight_layout()
plt.grid()
plt.savefig(str(plot_dir / "hac_stdev.png"), dpi=300)
plt.clf()

print(
    f"Ensembles per cluster: mean {mean(clustering_elements_sizes)},"
    f" median {median(clustering_elements_sizes)}, min {min(clustering_elements_sizes)},"
    f"max {max(clustering_elements_sizes)}"
)
sns.violinplot(
    y=clustering_elements_sizes, color="#fa9c3d", cut=1
)  # ,  color="b", alpha=0.5, edgecolor="black")
plt.ylabel("Ensembles in Cluster")
plt.xlabel("Set of All Clusters")
plt.tight_layout()
plt.grid()
plt.savefig(str(plot_dir / "ensembles_in_cluster.png"), dpi=300)
plt.clf()
