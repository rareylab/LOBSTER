from pathlib import Path
import pandas as pd
from matplotlib import pyplot as plt
import json
from argparse import ArgumentParser


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


results_csv = Path(config["stats_pairs"])

plot_dir = Path(config["plots_dir"])
if not plot_dir.is_dir():
    plot_dir.mkdir()

result_df = pd.read_csv(str(results_csv), delimiter=";")
all_morgan = result_df["morgan_fp_tanimoto"].tolist()
all_gobbi = result_df["gobbi_2D_pharmacophore_fp_tanimoto"].tolist()
bp = plt.boxplot([all_morgan, all_gobbi], showmeans=True)
plt.title("Tanimoto Fingerprint Similarity")
plt.xticks(ticks=[1, 2], labels=["ECFP 4", "Gobbi 2D"])
plt.savefig(str(plot_dir / "tanimoto_boxplots.png"), dpi=300)
