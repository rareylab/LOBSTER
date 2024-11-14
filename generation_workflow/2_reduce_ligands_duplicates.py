import pandas as pd
import progressbar
from collections import defaultdict
import argparse, json
from pathlib import Path


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="Deduplication by usmiles")
    arg_parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to config in JSON format including all file paths.",
    )
    args = arg_parser.parse_args()

    config = json.load(open(args.config, "r"))

    if not Path(config["1_main_filter_csv"]).is_file():
        print("Input file initial_csv specified in config does not exist")
        raise FileNotFoundError(config["1_main_filter_csv"])

    # create result dir if necessary
    if not Path(config["2_no_duplicates_csv"]).parent.is_dir():
        Path(config["2_no_duplicates_csv"]).parent.mkdir()

    if not Path(config["log_file"]).parent.is_dir():
        Path(config["log_file"]).parent.mkdir()

    pdb_ligand_dict = defaultdict(list)
    idx_to_keep = []

    df = pd.read_csv(config["1_main_filter_csv"], delimiter=",", dtype=str)
    for idx, row in progressbar.progressbar(df.iterrows()):
        pdb_ligand_dict[(row["pdb"], row["usmiles"])].append(
            (idx, float(row["EDIAm"]), row["name"])
        )

    for key in pdb_ligand_dict:
        idx_edia_list = pdb_ligand_dict[key]
        idx_edia_list.sort(key=lambda x: x[1], reverse=True)
        keep_idx = idx_edia_list[0][0]
        idx_to_keep.append(keep_idx)

    keeps_df = df.loc[idx_to_keep]
    keeps_df.to_csv(config["2_no_duplicates_csv"])

    with Path(config["log_file"]).open("a") as f:
        f.write(f"Removed deuplicates: {len(df) - len(keeps_df)}\n")
        f.write(f"Remaining: {len(keeps_df)}\n")
