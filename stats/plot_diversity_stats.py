import plotnine as p9
import pandas as pd
from matplotlib import pyplot as plt
import json
from argparse import ArgumentParser
from pathlib import Path


arg_parser = ArgumentParser(
    description="Plot statistics on diversity of the molecules"
    " within the LOBSTER and Astra Zeneca dataset"
)
arg_parser.add_argument(
    "--config",
    type=str,
    required=True,
    help="Path to config in JSON format including all file paths.",
)
args = arg_parser.parse_args()

config = json.load(open(args.config, "r"))
plots_path = Path(config["plots_dir"])

subset_data = pd.read_csv(config["diversity_per_subset"], index_col=0)

graph = (
    p9.ggplot(
        subset_data,
        p9.aes(
            x="Median MACCS FP-based Tanimoto Similarity",
            y="Scaffold AUC",
            color="Set",
            size="Size",
        ),
    )
    + p9.geom_point(alpha=0.5)
    + p9.xlim([0.3, 0.4])
    + p9.ylim([0.6, 0.7])
)
graph.save(plots_path / "diversity_subsets.png", dpi=800)

ax = subset_data.plot(
    x="Set",
    y="Size",
    kind="bar",
    color="blue",
    edgecolor="black",
    linewidth=1.2,
    alpha=0.5,
)
ax.set_xticklabels([f"≥{x/100}" for x in range(0, 100, 10)], rotation="horizontal")
ax.get_legend().remove()
ax.set_ylabel("Number of Ligand Pairs")
ax.set_xlabel("Shape Tversky Index")
plt.tight_layout()
plt.savefig(plots_path / "subset_size.png", dpi=300)

data_ensembles = pd.read_csv(config["diversity_per_ensemble"], index_col=0)
data_az = pd.read_csv(config["diversity_az"], index_col=0)

graph_1 = (
    p9.ggplot(
        data_ensembles,
        p9.aes(
            x="Median MACCS FP-based Tanimoto Similarity", y="Scaffold AUC", size="Size"
        ),
    )
    + p9.geom_point(alpha=0.5, color="red")
    + p9.xlim([0, 1])
    + p9.ylim([0.49, 0.76])
)
graph_1.save(plots_path / "diversity_ensemble.png", dpi=800)

graph_2 = (
    p9.ggplot(
        data_az,
        p9.aes(
            x="Median MACCS FP-based Tanimoto Similarity", y="Scaffold AUC", size="Size"
        ),
    )
    + p9.geom_point(alpha=0.5, color="blue")
    + p9.xlim([0, 1])
    + p9.ylim([0.49, 0.76])
)
graph_2.save(plots_path / "diversity_az.png", dpi=800)

data_ensembles["lobster"] = pd.Series(
    ["LOBSTER" for x in range(len(data_ensembles.index))], index=data_ensembles.index
)
data_az["az"] = pd.Series(
    ["Astra Zeneca" for x in range(len(data_az.index))], index=data_az.index
)

graph_2 = (
    p9.ggplot()
    + p9.geom_point(
        data_ensembles,
        p9.aes(
            x="Median MACCS FP-based Tanimoto Similarity",
            y="Scaffold AUC",
            size="Size",
            color="lobster",
        ),
        alpha=0.5,
    )
    + p9.geom_point(
        data_az,
        p9.aes(
            x="Median MACCS FP-based Tanimoto Similarity",
            y="Scaffold AUC",
            size="Size",
            color="az",
        ),
        alpha=0.5,
    )
    + p9.xlim([0, 1])
    + p9.ylim([0.49, 0.76])
)
graph_2.save(plots_path / "diversity_merged.png", dpi=800)

print("Finished")
