from collections import defaultdict
import os, argparse
import pandas as pd
import progressbar
import shutil, json
from statistics import mean
from pathlib import Path


# Filter out ensembles containing less than min_nof_ligands ligands
# from directory siena_dir
def clean_up(
    min_nof_ligands,
    sdf_dir,
):
    keep_files = []
    no_reference = 0
    too_few_ligands = 0
    cluster_index = defaultdict(lambda: [])

    print(f"Sorting out ensembles with less than {min_nof_ligands} ligands")
    for ligand_path in progressbar.progressbar(sdf_dir.iterdir()):
        # Filter for number of ligands
        content_read = ligand_path.open("r")
        content = content_read.read()
        nof_ligands = content.count("$$$$")
        query_name, pdb_code = ligand_path.name.split("-")[:2]
        nof_references = content.count(f"{query_name}")
        content_read.close()
        if nof_ligands < min_nof_ligands:
            too_few_ligands += 1
            continue
        elif nof_references < 1:
            no_reference += 1
            continue
        keep_files.append((ligand_path, nof_ligands))
        # Create Index if ensemble is accepted
        ligands_content = [
            x for x in ligand_path.open("r").read().split("$$$$") if x != ""
        ]
        file_base_name = f"{query_name}-{pdb_code}"

        for entry in ligands_content:
            entry_lines = [x for x in entry.split("\n") if x != ""]
            if entry_lines != []:
                pdb = entry_lines[entry_lines.index("> <pdb_id>") + 1]
                ligand = entry_lines[0]
                # cluster by ligand and pdb
                cluster_index[file_base_name].append((ligand, pdb))
    keep_files.sort(key=lambda x: (x[1], x[0]), reverse=True)
    keep_files = [x for x, _ in keep_files]
    print(f"{no_reference} ensemble without reference ligand filtered out", flush=True)
    print(
        f"{too_few_ligands} ensembles contained less than {min_nof_ligands} ligands",
        flush=True,
    )
    print(f"Processing {len(keep_files)} remaining ensembles", flush=True)
    return keep_files, cluster_index


# ----------------------------------------------------------------------------------------


def pick_representative_by_mean_bbrmsd(
    cluster,
    stats_dir,
    sdfs_dir,
):
    representative_queue = []
    for file_base_name in cluster:
        statistics_path = f"{stats_dir}/{file_base_name}-resultStatistic.csv"
        stats_df = pd.read_csv(str(statistics_path), delimiter=";", dtype=str)
        rmsds = []
        ligands_content = [
            x
            for x in open(f"{sdfs_dir}/{file_base_name}-ligands.sdf", "r")
            .read()
            .split("$$$$")
            if x != ""
        ]

        for entry in ligands_content:
            entry_lines = [x for x in entry.split("\n") if x != ""]
            if entry_lines != []:
                pdb = entry_lines[entry_lines.index("> <pdb_id>") + 1]
                ligand = entry_lines[0]
                # leave out reference with bb rmsd 0
                if file_base_name != "-".join([ligand, pdb.lower()]):
                    bb_rmsd = float(
                        stats_df.where(
                            (stats_df["PDB code"] == pdb)
                            & (stats_df["Ligand PDB code"] == ligand)
                        )
                        .dropna(subset=["PDB code"])["Backbone RMSD"]
                        .values[0]
                    )
                    rmsds.append(bb_rmsd)
        representative_queue.append((file_base_name, mean(rmsds)))
    representative_queue.sort(key=lambda x: x[1])
    print(representative_queue)
    return representative_queue[0][0]


# ----------------------------------------------------------------------------------------


def pick_representative_by_highest_query_HAC(cluster, filtered_table_path):
    filtered_master_table = pd.read_csv(str(filtered_table_path), dtype=str)
    representative_queue = []

    for file_base_name in cluster:
        query, pdb = file_base_name.split("-")
        table_row = filtered_master_table.where(
            (filtered_master_table["pdb"] == pdb.upper())
            & (filtered_master_table["name"] == query)
        ).dropna(subset=["pdb"])
        representative_queue.append(
            (file_base_name, float(table_row["heavyAtoms"].values[0]))
        )
    representative_queue.sort(key=lambda x: x[1], reverse=True)
    return representative_queue[0][0]


# ------------------------------------------------------------------------------


# Note: filtered_table_path is only needed for clustering by query HAC
def cluster(dir, sdf_dir, res_dir, keep_files, cluster_index, filtered_table_path):
    cluster_dir = res_dir / "clusters"
    if not cluster_dir.is_dir():
        cluster_dir.mkdir()

    representative_dir = res_dir / "representatives"
    if not representative_dir.is_dir():
        representative_dir.mkdir()

    cluster_set = []

    for ligand_path in progressbar.progressbar(keep_files):
        query_name, pdb_code = ligand_path.name.split("-")[:2]
        file_base_name = f"{query_name}-{pdb_code}"

        # Cluster Alignments
        ensemble_added = []
        for elem_list in cluster_set:
            for elem in elem_list:
                if bool(
                    set(cluster_index[elem]) & (set(cluster_index[file_base_name]))
                ):
                    elem_list.append(file_base_name)
                    ensemble_added.append(elem_list)
                    break
        if ensemble_added == []:
            cluster_set.append([file_base_name])
        if len(ensemble_added) > 1:
            all_elems = set()
            for elem in ensemble_added:
                all_elems = all_elems.union(set(elem))
                cluster_set.remove(elem)
            cluster_set.append(list(all_elems))

    print(f"Writing {len(cluster_set)} clusters to {res_dir}", flush=True)
    cluster_log = res_dir / "clusters_log.txt"
    cluster_file = cluster_log.open("w")

    for elem_list in progressbar.progressbar(cluster_set):
        # Choose representative
        representative = pick_representative_by_highest_query_HAC(
            cluster=elem_list, filtered_table_path=filtered_table_path
        )
        # Copy representative into extra directory
        # and move clustered files to a directory wrt their cluster.
        shutil.copy(
            f"{sdf_dir}/{representative}-ligands.sdf",
            f"{representative_dir}/{representative}-ligands.sdf",
        )
        single_cluster_dir = f"{cluster_dir}/{representative}"
        cluster_file.write(f"Representative: {representative}\n")

        if not os.path.exists(single_cluster_dir):
            os.makedirs(single_cluster_dir)

        for elem in elem_list:
            cluster_file.write(f"{elem}\n")
            shutil.copy(
                f"{sdf_dir}/{elem}-ligands.sdf",
                f"{single_cluster_dir}/{elem}-ligands.sdf",
            )
            shutil.copy(
                f"{dir}/{elem}-alignment.txt",
                f"{single_cluster_dir}/{elem}-alignment.txt",
            )
            shutil.copy(
                f"{dir}/{elem}-resultStatistic.csv",
                f"{single_cluster_dir}/{elem}-resultStatistic.csv",
            )
        cluster_file.write("--------------------------------------------\n")
    cluster_file.close()


# ----------------------------------------------------------------------------------------


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(
        description="Clustering by PDB ID and ligand hetcode."
    )
    arg_parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to config in JSON format including all file paths.",
    )
    args = arg_parser.parse_args()

    config = json.load(open(args.config, "r"))

    siena_dir = Path(config["8_siena_results_dir"])
    unique_sdfs_dir = Path(config["9_unique_ligands_dir"])
    final_dataset_dir = Path(config["10_final_dataset_dir"])
    table_path = Path(config["6_filtered_buriedness_csv"])

    if not table_path.is_file():
        print(
            "Input file 6_filtered_buriedness_csv as specified in config does not exist"
        )
        raise FileNotFoundError(str(table_path))

    if not siena_dir.is_dir():
        print(
            "Input directory 6_finally_filtered_pdbs_ligands_json"
            " as specified in config does not exist"
        )
        raise IOError(str(siena_dir))

    if not unique_sdfs_dir.is_dir():
        print(
            "Input directory 6_finally_filtered_pdbs_ligands_json"
            " as specified in config does not exist"
        )
        raise IOError(str(unique_sdfs_dir))

    # create result dir if necessary
    if not final_dataset_dir.is_dir():
        final_dataset_dir.mkdir()

    files, index = clean_up(2, unique_sdfs_dir)
    cluster(siena_dir, unique_sdfs_dir, final_dataset_dir, files, index, table_path)
