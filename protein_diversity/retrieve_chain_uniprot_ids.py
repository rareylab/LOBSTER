import requests
from pathlib import Path
import pandas as pd
import json
from argparse import ArgumentParser
from progressbar import progressbar


def get_uniprot_for_chain(pdb_id, chain_id):
    new_results = requests.get(
        f"https://data.rcsb.org/rest/v1/core/polymer_entity_instance/"
        f"{pdb_id.lower()}/{chain_id}"
    )
    new_entity_id = new_results.json()[
        "rcsb_polymer_entity_instance_container_identifiers"
    ]["entity_id"]
    uniprot_id_raw = requests.get(
        f"https://data.rcsb.org/rest/v1/core/uniprot/{pdb_id.lower()}/{new_entity_id}"
    )
    uniprot_id = uniprot_id_raw.json()[0]["rcsb_uniprot_container_identifiers"][
        "uniprot_id"
    ]
    return uniprot_id


# ----------------------------------------------------------------------------------------


def get_identifier_from_uniprot(uniprot_id):
    result = requests.get(f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}")
    pfam_id = None
    for reference in result.json()["dbReferences"]:
        if reference["type"] == "Pfam":
            pfam_id = reference["id"]
    return pfam_id


# ----------------------------------------------------------------------------------------


def get_name_from_pfam(pfam_id):
    results = requests.get(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam_id}")
    full_name = results.json()["metadata"]["name"]["name"]
    return full_name


# ----------------------------------------------------------------------------------------


def get_pfam_name_for_pdb_chain(pdb_id, chain_id):
    try:
        uniprot_id = get_uniprot_for_chain(pdb_id, chain_id)
        pfam_id = get_identifier_from_uniprot(uniprot_id)
        full_pfam_name = get_name_from_pfam(pfam_id)
    except:
        full_pfam_name = "failed"
    return full_pfam_name


# ----------------------------------------------------------------------------------------


if __name__ == "__main__":
    arg_parser = ArgumentParser(
        description="Retrieve chain uniprot IDs for the LOBSTER dataset"
    )
    arg_parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to config in JSON format including all file paths.",
    )
    args = arg_parser.parse_args()

    config = json.load(open(args.config, "r"))

    asym_chain_ids = {}
    asym_data = pd.read_csv(config["asym_id_csv"], dtype=str)
    for index, row in progressbar(asym_data.iterrows()):
        auth_asym = row["auth_asym"]
        asym = row["asym"]
        pdb = row["entry"].upper()
        if auth_asym != asym:
            if not pdb in asym_chain_ids:
                asym_chain_ids[pdb] = {auth_asym: asym}
            else:
                asym_chain_ids[pdb][auth_asym] = asym

    print("Finished processing asym data, continuing with chain ID assignment")

    master_table_path = config["6_filtered_buriedness_csv"]
    lobster_path = Path(config["10_final_dataset_dir"])

    lobster_molecules = {}
    all_pdbs = set()
    all_ligands = set()
    for cluster in (lobster_path / "clusters").glob("*"):
        lobster_molecules[cluster.name] = []
        for sdf in cluster.glob("*.sdf"):
            name = sdf.name.split("-")[0]
            pdb = sdf.name.split("-")[1].upper()
            all_pdbs.add(pdb)
            all_ligands.add(name)
            lobster_molecules[cluster.name].append([pdb, name])

    for chunk in pd.read_csv(master_table_path, chunksize=1000, dtype=str):
        for index, row in chunk.iterrows():
            pdb = row["pdb"].upper()
            chains = row["chains"]
            chain = row["chain"]
            if pdb not in all_pdbs:
                continue
            if pd.isna(chains):
                if pd.isna(chain):
                    continue
                chains = row["chain"]
            full_name = (
                row["HETCode"] + "_" + row["chain"] + "_" + row["ID"].split(".")[0]
            )
            if not full_name in all_ligands:
                continue
            for cluster in lobster_molecules:
                for molecule in lobster_molecules[cluster]:
                    if pdb == molecule[0].upper() and full_name == molecule[1]:
                        molecule.append(chains)
    n_skipped = 0
    cluster_uniprot_dict = {}
    cluster_pfam_id_dict = {}
    cluster_pfam_name_dict = {}
    cluster_data = []

    for cluster in progressbar(lobster_molecules):
        cluster_uniprot_dict[cluster] = set()
        cluster_pfam_id_dict[cluster] = set()
        cluster_pfam_name_dict[cluster] = set()

        for molecule in lobster_molecules[cluster]:
            pdb = molecule[0].upper()
            ligand = molecule[1]
            chains = molecule[2]
            if len(molecule) == 2:
                n_skipped += 1
                continue
            for current_chain in molecule[2].split(","):
                if (
                    pdb in asym_chain_ids.keys()
                    and current_chain in asym_chain_ids[pdb].keys()
                ):
                    correct_chain = asym_chain_ids[pdb][current_chain]
                else:
                    correct_chain = current_chain
                current_uniprot_id = None
                pfam_id = None
                full_pfam_name = None
                try:
                    current_uniprot_id = get_uniprot_for_chain(pdb, correct_chain)
                    pfam_id = get_identifier_from_uniprot(current_uniprot_id)
                    full_pfam_name = get_name_from_pfam(pfam_id)
                except:
                    if not current_uniprot_id:
                        current_uniprot_id = "failed"
                    if not pfam_id:
                        pfam_id = "failed"
                    if not full_pfam_name:
                        full_pfam_name = "failed"
                cluster_uniprot_dict[cluster].add(current_uniprot_id)
                cluster_pfam_id_dict[cluster].add(pfam_id)
                cluster_pfam_name_dict[cluster].add(full_pfam_name)
                cluster_data.append(
                    [
                        cluster,
                        molecule[0],
                        molecule[1],
                        correct_chain,
                        current_uniprot_id,
                        pfam_id,
                        full_pfam_name,
                    ]
                )

    final_data = []
    for cluster in cluster_uniprot_dict:
        final_data.append(
            [
                cluster,
                ";".join(cluster_uniprot_dict[cluster]),
                ";".join(cluster_pfam_id_dict[cluster]),
                ";".join(cluster_pfam_name_dict[cluster]),
            ]
        )
    final_data = pd.DataFrame(
        final_data, columns=["Cluster", "Uniprotids", "PfamID", "PfamName"]
    )
    final_data.to_csv(config["broad_cluster_csv"])
    cluster_data = pd.DataFrame(
        cluster_data,
        columns=[
            "Cluster",
            "PDB",
            "Ligand",
            "chain",
            "Uniprotid",
            "PfamID",
            "Pfam Name",
        ],
    )
    cluster_data.to_csv(config["exact_cluster_csv"])
    print("Number of skipped entries:", n_skipped)
