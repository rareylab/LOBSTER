from pathlib import Path
from argparse import ArgumentParser
import json
from progressbar import progressbar
import pandas as pd


def write_pairs(all_pairs_csv_path, overlay_representatives, single_files):
    lig_pdb_combinations = set()
    with all_pairs_csv_path.open("w") as f:
        for file in overlay_representatives.iterdir():
            ligands_content = [
                x for x in file.open("r").read().split("$$$$") if x != ""
            ]
            ligand_pdb_list = []
            for entry in ligands_content:
                entry_lines = [x for x in entry.split("\n") if x != ""]
                if entry_lines != []:
                    pdb_idx = entry_lines.index("> <pdb_id>") + 1
                    pdb = entry_lines[pdb_idx]
                    ligand = entry_lines[0]
                    ligand_file_path = single_files / f"{ligand}-{pdb}.sdf"
                    ligand_file = ligand_file_path.open("w")
                    ligand_file.write(entry.strip())
                    ligand_file.write("\n\n$$$$\n")
                    ligand_file.close()
                    ligand_pdb_list.append(f"{ligand}-{pdb}.sdf")
                    if not (ligand, pdb) in lig_pdb_combinations:
                        lig_pdb_combinations.add((ligand, pdb))
                    else:
                        print("This should never happen.")
                        raise RuntimeError("Same ligand-pdb occurring twice!")
            # All pairs, including a with b AND b with a!
            for template_sdf in ligand_pdb_list:
                for ligand_sdf in ligand_pdb_list:
                    if template_sdf != ligand_sdf:
                        query = file.stem.replace("-ligands", "")
                        f.write(f"{template_sdf};{ligand_sdf};{query}\n")


# ----------------------------------------------------------------------------------------


if __name__ == "__main__":
    arg_parser = ArgumentParser(
        description="Split the dataset to single files and a csvs for all pairs"
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

    single_files_dir = final_dataset_dir / "all_ligands"
    if not single_files_dir.is_dir():
        single_files_dir.mkdir()

    subsets_dir = final_dataset_dir / "subsets"
    if not subsets_dir.is_dir():
        subsets_dir.mkdir()

    all_pairs_csv = subsets_dir / "all_pairs.csv"

    write_pairs(all_pairs_csv, final_dataset_dir / "representatives", single_files_dir)
