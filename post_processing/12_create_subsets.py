from pathlib import Path
from argparse import ArgumentParser
import json
from progressbar import progressbar
import pandas as pd


def write_header(file):
    file.write(
        "query_file;query;pdb;template_file;template;template_pdb;ensemble;"
        "morgan_fp_tanimoto_r4;gobbi_2D_pharmacophore_fp_tanimoto;hac_difference;"
        "shape_tversky_index;shape_tanimoto_distance;shape_protrude_distance;\n"
    )


# ----------------------------------------------------------------------------------------


def split_dataset(subsets_path, pair_stats_csv):
    subset_90 = subsets_path / "subset_90.csv"
    subset_80 = subsets_path / "subset_80.csv"
    subset_70 = subsets_path / "subset_70.csv"
    subset_60 = subsets_path / "subset_60.csv"
    subset_50 = subsets_path / "subset_50.csv"
    subset_40 = subsets_path / "subset_40.csv"
    subset_30 = subsets_path / "subset_30.csv"
    subset_20 = subsets_path / "subset_20.csv"
    subset_10 = subsets_path / "subset_10.csv"
    subset_0 = subsets_path / "subset_0.csv"

    df = pd.read_csv(pair_stats_csv, delimiter=";")
    df = df.sort_values(by="shape_tversky_index", ascending=False)

    b90 = subset_90.open("w")
    write_header(b90)
    b80 = subset_80.open("w")
    write_header(b80)
    b70 = subset_70.open("w")
    write_header(b70)
    b60 = subset_60.open("w")
    write_header(b60)
    b50 = subset_50.open("w")
    write_header(b50)
    b40 = subset_40.open("w")
    write_header(b40)
    b30 = subset_30.open("w")
    write_header(b30)
    b20 = subset_20.open("w")
    write_header(b20)
    b10 = subset_10.open("w")
    write_header(b10)
    b0 = subset_0.open("w")
    write_header(b0)

    for idx, row in progressbar(df.iterrows()):
        if row["shape_tversky_index"] >= 0.9:
            df.loc[[idx]].to_csv(
                b90, sep=";", index=False, header=False, line_terminator="\n"
            )
        elif row["shape_tversky_index"] >= 0.8:
            df.loc[[idx]].to_csv(
                b80, sep=";", index=False, header=False, line_terminator="\n"
            )
        elif row["shape_tversky_index"] >= 0.7:
            df.loc[[idx]].to_csv(
                b70, sep=";", index=False, header=False, line_terminator="\n"
            )
        elif row["shape_tversky_index"] >= 0.6:
            df.loc[[idx]].to_csv(
                b60, sep=";", index=False, header=False, line_terminator="\n"
            )
        elif row["shape_tversky_index"] >= 0.5:
            df.loc[[idx]].to_csv(
                b50, sep=";", index=False, header=False, line_terminator="\n"
            )
        elif row["shape_tversky_index"] >= 0.4:
            df.loc[[idx]].to_csv(
                b40, sep=";", index=False, header=False, line_terminator="\n"
            )
        elif row["shape_tversky_index"] >= 0.3:
            df.loc[[idx]].to_csv(
                b30, sep=";", index=False, header=False, line_terminator="\n"
            )
        elif row["shape_tversky_index"] >= 0.2:
            df.loc[[idx]].to_csv(
                b20, sep=";", index=False, header=False, line_terminator="\n"
            )
        elif row["shape_tversky_index"] >= 0.1:
            df.loc[[idx]].to_csv(
                b10, sep=";", index=False, header=False, line_terminator="\n"
            )
        else:
            df.loc[[idx]].to_csv(
                b0, sep=";", index=False, header=False, line_terminator="\n"
            )


# ----------------------------------------------------------------------------------------


if __name__ == "__main__":
    arg_parser = ArgumentParser(
        description="Split the dataset of all pairs to csvs for the subsets"
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
    # Dataset should have already been built by the time this script is executed
    if not final_dataset_dir.is_dir():
        print(
            "Input directory 10_final_dataset_dir as specified in config does not exist"
        )
        raise IOError(str(final_dataset_dir))

    pair_stats = Path(config["stats_pairs"])
    if not pair_stats.is_file():
        print("Input file stats_pairs as specified in config does not exist")
        print("Please execute stats/all_csv_stats.py first!")
        raise FileNotFoundError(str(pair_stats))

    single_files_dir = final_dataset_dir / "all_ligands"
    if not single_files_dir.is_dir():
        single_files_dir.mkdir()

    subsets_dir = final_dataset_dir / "subsets"
    split_dataset(subsets_dir, pair_stats)
