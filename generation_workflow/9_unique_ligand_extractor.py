import pandas as pd
import sys, json
import progressbar, argparse
from statistics import mean
from pathlib import Path


def write_ligand_file(ligand_path, unique_ligands_dict):
    ligand_file = open(ligand_path, "w")
    for key in unique_ligands_dict:
        ligand_file.write(unique_ligands_dict[key]["entry"].strip())
        ligand_file.write("\n\n$$$$\n")
    ligand_file.close()


if __name__ == "__main__":

    arg_parser = argparse.ArgumentParser(
        description="Filter for unique ligands within the ensembles."
    )
    arg_parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to config in JSON format including all file paths.",
    )
    args = arg_parser.parse_args()

    config = json.load(open(args.config, "r"))

    input_dir = Path(config["8_siena_results_dir"])
    unique_ligands_dir = Path(config["9_unique_ligands_dir"])

    if not input_dir.is_dir():
        print(
            "Input directory 6_finally_filtered_pdbs_ligands_json as specified in config"
            "does not exist"
        )
        raise IOError(str(input_dir))

    # create result dir if necessary
    if not unique_ligands_dir.is_dir():
        unique_ligands_dir.mkdir()

    filtered_master_table = pd.read_csv(config["6_filtered_buriedness_csv"])
    unique_queries = {}

    print("Uniqueness is checked...", flush=True)
    for file in progressbar.progressbar(input_dir.iterdir()):
        if file.suffix == ".sdf":
            query_ligand, query_pdb = file.name.split("-")[0:2]
            avg_bb_rmsd = sys.float_info.max
            bb_rmsds = []
            ligands_content = [
                x
                for x in file.open("r", encoding="utf8").read().split("$$$$")
                if x != ""
            ]
            unique_ligands = {}
            file_base_name = f"{query_ligand}-{query_pdb}"
            statistics_path = file.parent / f"{file_base_name}-resultStatistic.csv"
            stats_df = pd.read_csv(str(statistics_path), delimiter=";", dtype=str)
            for entry in ligands_content:
                entry_lines = [x for x in entry.split("\n") if x != ""]
                if entry_lines != []:
                    pdb_idx = entry_lines.index("> <pdb_id>") + 1
                    pdb = entry_lines[pdb_idx]
                    ligand = entry_lines[0]
                    table_entry = filtered_master_table.where(
                        (filtered_master_table["pdb"] == pdb)
                        & (filtered_master_table["name"] == ligand)
                    ).dropna(subset=["pdb"])
                    usmiles = table_entry["usmiles"].values[0]
                    uniprot = table_entry["uniprotID"].values[0]

                    # extend entry by uniprot-ID and usmiles
                    entry_lines.insert(pdb_idx + 1, "\n> <uniprot_id>")
                    entry_lines.insert(pdb_idx + 2, f"{uniprot}")
                    entry_lines.insert(pdb_idx + 3, "\n> <usmiles>")
                    entry_lines.insert(pdb_idx + 4, f"{usmiles}")
                    final_entry = "\n".join(entry_lines)

                    bbrmsd = float(
                        stats_df.where(
                            (stats_df["PDB code"] == pdb)
                            & (stats_df["Ligand PDB code"] == ligand)
                        )
                        .dropna(subset=["PDB code"])["Backbone RMSD"]
                        .values[0]
                    )
                    if (not usmiles in unique_ligands.keys()) or (
                        bbrmsd < unique_ligands[usmiles]["bbrmsd"]
                    ):
                        unique_ligands[usmiles] = {
                            "entry": final_entry,
                            "bbrmsd": bbrmsd,
                        }
                    # ignore template with bbrmsd = 0
                    if ligand != query_ligand:
                        bb_rmsds.append(bbrmsd)
                    # ensure query ligand from query pdb is in dictionary
                    elif pdb.lower() == query_pdb.lower():
                        unique_ligands[usmiles] = {
                            "entry": final_entry,
                            "bbrmsd": bbrmsd,
                        }

            if len(bb_rmsds) > 0:
                avg_bb_rmsd = mean(bb_rmsds)
            ligand_path = unique_ligands_dir / file.name
            if (not query_ligand in unique_queries.keys()) or (
                unique_queries[query_ligand]["avg_bbrmsd"] > avg_bb_rmsd
            ):
                unique_queries[query_ligand] = {
                    "avg_bbrmsd": avg_bb_rmsd,
                    "ligand_path": ligand_path,
                    "content": unique_ligands,
                }

    print(f"Rewriting {len(unique_queries.keys())} ensembles...")
    for key in progressbar.progressbar(unique_queries.keys()):
        ligand_dict = unique_queries[key]
        write_ligand_file(ligand_dict["ligand_path"], ligand_dict["content"])
