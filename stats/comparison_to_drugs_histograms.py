import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import json
from pathlib import Path


arg_parser = ArgumentParser(
    description="Plot statistics on the druglikeness of the molecules from the"
    " LOBSTER dataset"
)
arg_parser.add_argument(
    "--config",
    type=str,
    required=True,
    help="Path to config in JSON format including all file paths.",
)
args = arg_parser.parse_args()

config = json.load(open(args.config, "r"))

df = pd.read_csv(config["known_drugs_csv"], sep=";", header=0)
print(df)
oral_df = df[df["source"] == "FDA_ORAL"]
our_df = df[df["source"] == "LOBSTER"]

patterns = ["", "/"]
colors1 = ["deepskyblue", "lightskyblue"]
colors2 = ["deepskyblue", "lightskyblue", "r", "coral"]

plt.rcParams["figure.figsize"] = (10, 7)
plt.rcParams.update({"font.size": 11})
fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(nrows=2, ncols=2)

ax0.hist(
    [oral_df["SlogP"], our_df["SlogP"]],
    bins=np.arange(math.floor(min(df["SlogP"])), max(df["SlogP"]) + 1, 1),
    label=("orally available drugs", "benchmark ligands"),
    density=True,
    color=colors1,
)
print(np.arange(math.floor(min(df["SlogP"])), max(df["SlogP"]) + 1, 1))
ax1.hist(
    [oral_df["NumRotatableBonds"], our_df["NumRotatableBonds"]],
    bins=np.arange(
        math.floor(min(df["NumRotatableBonds"])), max(df["NumRotatableBonds"]) + 1, 1
    ),
    label=("orally available drugs", "benchmark ligands"),
    density=True,
    color=colors1,
)
ax2.hist(
    [oral_df["ExactMW"], our_df["ExactMW"]],
    bins=np.arange(0, max(df["ExactMW"]) + 50, 50),
    label=("orally available drugs", "benchmark ligands"),
    density=True,
    color=colors1,
)
ax3.hist(
    [oral_df["NumHBD"], our_df["NumHBD"], oral_df["NumHBA"], our_df["NumHBA"]],
    bins=np.arange(math.floor(min(df["NumHBA"])), max(df["NumHBA"]) + 1, 1),
    label=(
        "orally available drugs (HBD)",
        "benchmark ligands (HBD)",
        "orally available drugs (HBA)",
        "benchmark ligands (HBA)",
    ),
    density=True,
    color=colors2,
)
ax0.legend(loc=2, fontsize="small")
ax1.legend(fontsize="small")
ax2.legend(fontsize="small")
ax3.legend(fontsize="small")
ax0.set_ylabel("Relative Frequency")
ax0.set_xlabel("SlogP")
ax1.set_ylabel("Relative Frequency")
ax1.set_xlabel("Number of Rotatable Bonds")
ax2.set_ylabel("Relative Frequency")
ax2.set_xlabel("Molecular Weight [g/mol]")
ax3.set_ylabel("Relative Frequency")
ax3.set_xlabel("Number of Hydrogen Bond Donors and Acceptors")

plt.tight_layout()
plt.savefig(Path(config["plots_dir"]) / "binned_descriptors_drugbank.png", dpi=600)
