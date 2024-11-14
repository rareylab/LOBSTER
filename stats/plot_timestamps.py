from pathlib import Path
from matplotlib import pyplot as plt
import json
from argparse import ArgumentParser

arg_parser = ArgumentParser(
    description="Create statistics on timestamps of the molecules"
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
timestamps = Path(config["timestamps_json"])
if not timestamps.is_file():
    print("Input directory 10_final_dataset_dir as specified in config does not exist")
    raise FileNotFoundError(str(timestamps))
plots_dir = Path(config["plots_dir"])
if not plots_dir.is_dir():
    print("Input directory 10_final_dataset_dir as specified in config does not exist")
    raise IOError(str(plots_dir))

with timestamps.open("r") as f:
    new_pairs = json.load(f)

years = []
counts = []
accumulated = []
pdb_resolutions = []
edia_m_values = []
mw_values = []
hac_values = []
rot_bonds = []

for i, year in enumerate(sorted(new_pairs.keys())):
    current_count = new_pairs[year]["count"]
    if current_count != 0:
        years.append(year)
        counts.append(current_count)
        if len(accumulated) > 0:
            accumulated.append(current_count + accumulated[len(accumulated) - 1])
        else:
            accumulated.append(current_count)


plt.barh(years, counts, color="b", edgecolor="black")
plt.yticks(years)
plt.ylabel("Year")
plt.xticks(rotation=90)
plt.xlabel("Number of Ligand Pairs")
plt.tight_layout()
plt.savefig(plots_dir / "timestamps.png", dpi=300)
plt.clf()

plt.barh(years, counts, color="b", edgecolor="black")
plt.yticks(years)
plt.ylabel("Year")
plt.xticks(rotation=90)
plt.xlabel("Number of Ligand Pairs (Log Scaled)")
plt.xscale("log")
plt.tight_layout()
plt.savefig(plots_dir / "timestamps_log.png", dpi=300)
plt.clf()

plt.barh(years, accumulated, color="b", edgecolor="black")
plt.yticks(years)
plt.ylabel("Year")
plt.xticks(rotation=90)
plt.xlabel("Number of Accumulated Ligand Pairs (Log Scaled)")
plt.xscale("log")
plt.tight_layout()
plt.savefig(plots_dir / "accumulated_timestamps_log.png", dpi=300)
plt.clf()

years = []
for year in sorted(new_pairs.keys()):
    years.append(year)
    pdb_resolutions.append(new_pairs[year]["pdb_resolutions"])
    edia_m_values.append(new_pairs[year]["edia_m_values"])
    mw_values.append(new_pairs[year]["mw_values"])
    hac_values.append(new_pairs[year]["hac_values"])
    rot_bonds.append(new_pairs[year]["rot_bonds"])

plt.boxplot(pdb_resolutions, vert=False, showmeans=True)
plt.yticks(range(0, len(years)), labels=years)
plt.ylabel("Year")
plt.xticks(rotation=90)
plt.xlabel("PDB Resolution (cutoff > 2.5 Å)")
plt.tight_layout()
plt.savefig(plots_dir / "timestamps_pdb_resolution.png", dpi=300)
plt.clf()

plt.boxplot(edia_m_values, vert=False, showmeans=True)
plt.yticks(range(0, len(years)), labels=years)
plt.ylabel("Year")
plt.xticks(rotation=90)
plt.xlabel("EDIAm (cutoff > 0.8)")
plt.tight_layout()
plt.savefig(plots_dir / "timestamps_edia_m.png", dpi=300)
plt.clf()

plt.boxplot(mw_values, vert=False, showmeans=True)
plt.yticks(range(0, len(years)), labels=years)
plt.ylabel("Year")
plt.xticks(rotation=90)
plt.xlabel("Molecular Weight (cutoff > 975 Da)")
plt.tight_layout()
plt.savefig(plots_dir / "timestamps_mw.png", dpi=300)
plt.clf()

plt.boxplot(hac_values, vert=False, showmeans=True)
plt.yticks(range(0, len(years)), labels=years)
plt.ylabel("Year")
plt.xticks(rotation=90)
plt.xlabel("Heavy Atom Count (cutoff < 10)")
plt.tight_layout()
plt.savefig(plots_dir / "timestamps_hac.png", dpi=300)
plt.clf()

plt.boxplot(rot_bonds, vert=False, showmeans=True)
plt.yticks(range(0, len(years)), labels=years)
plt.ylabel("Year")
plt.xticks(rotation=90)
plt.xlabel("Number of Rotatable Bonds (cutoff > 10)")
plt.tight_layout()
plt.savefig(plots_dir / "timestamps_rot_bonds.png", dpi=300)
plt.clf()
