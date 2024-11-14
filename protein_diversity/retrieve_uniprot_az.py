from pathlib import Path
import pandas as pd
import re, json
from retrieve_chain_uniprot_ids import (
    get_uniprot_for_chain,
    get_pfam_name_for_pdb_chain,
    get_name_from_pfam,
    get_identifier_from_uniprot,
)
from argparse import ArgumentParser
from progressbar import progressbar


def get_exact_name(pdb, hetcode, az_data):
    for index, row in az_data.iterrows():
        current_hetcode = row["hetcode"]
        current_pdb = row["pdb"].upper()
        ligand = row["ligand"]
        if current_pdb == pdb and current_hetcode == hetcode:
            return ligand
    for index, row in az_data.iterrows():
        current_pdb = row["pdb"]
        ligand = row["ligand"]
        if current_pdb == pdb:
            return ligand
    return None


# ----------------------------------------------------------------------------------------


if __name__ == "__main__":
    arg_parser = ArgumentParser(
        description="Retrieve chain uniprot IDs for the AstraZeneca dataset"
    )
    arg_parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to config in JSON format including all file paths.",
    )
    args = arg_parser.parse_args()

    config = json.load(open(args.config, "r"))

    obsolete_pdb_ids = {"2ZCP": "3W7F", "3P25": "4DZ7", "3P29": "4DZ9", "1QY5": "6D28"}
    inverse_obsolete_pdbs = {obsolete_pdb_ids[x]: x for x in obsolete_pdb_ids}

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

initial_data = Path(config["initial_csv"])
print("Finished processing asym data, continuing with chain ID assignment")

az_overlay_path = Path(config["path_to_az_data"])
az_ligand_identifier_path = Path(config["az_ligand_identifiers"])
az_identifier_data = pd.read_csv(az_ligand_identifier_path, dtype=str)
az_clusters = {}
for sdf in az_overlay_path.glob("*"):
    az_clusters[sdf.name] = []
    with open(str(sdf), "r") as f:
        for line in f.readlines():
            line = line.strip()
            if re.search("_lig_", line) != None:
                pdb = line.split("_")[0].upper()
                hetcode = line.split("_")[2]
                if "." in hetcode:
                    continue
                az_clusters[sdf.name].append([pdb, hetcode])

cluster_uniprot_dict = {}
cluster_pfam_id_dict = {}
cluster_pfam_name_dict = {}
cluster_data = []
not_found = []
all_pdbs = set()
all_ligands = set()
pdb_ligand_data = []
az_molecules = {}
for cluster in az_clusters.keys():
    cluster_uniprot_dict[cluster] = set()
    cluster_pfam_id_dict[cluster] = set()
    cluster_pfam_name_dict[cluster] = set()
    az_molecules[cluster] = []

    for pdb, hetcode in az_clusters[cluster]:
        ligand = get_exact_name(pdb, hetcode, az_identifier_data)
        az_molecules[cluster].append([pdb, ligand])
        all_pdbs.add(pdb)
        all_ligands.add(ligand)

for chunk in pd.read_csv(initial_data, chunksize=10000, dtype=str):
    for index, row in chunk.iterrows():
        pdb = row["pdb"].upper()
        if pdb in inverse_obsolete_pdbs:
            pdb = inverse_obsolete_pdbs[pdb]
            print("inverse pdb found")
        chains = row["chains"]
        chain = row["chain"]
        if pdb not in all_pdbs:
            continue
        if pd.isna(chains):
            if pd.isna(chain):
                continue
            chains = row["chain"]
        full_name = row["HETCode"] + "_" + row["chain"] + "_" + row["ID"]
        if not full_name in all_ligands:
            continue
        for cluster in az_molecules:
            for molecule in az_molecules[cluster]:
                if pdb == molecule[0] and full_name == molecule[1]:
                    molecule.append(chains)
n_skipped = 0
for cluster in progressbar(az_molecules):
    cluster_uniprot_dict[cluster] = set()
    cluster_pfam_id_dict[cluster] = set()
    cluster_pfam_name_dict[cluster] = set()
    for molecule in az_molecules[cluster]:
        if len(molecule) == 2:
            n_skipped += 1
            molecule.append(str(molecule[1].split("_")[1]))
        pdb = molecule[0].upper()
        old_pdb = pdb
        ligand = molecule[1]
        chains = molecule[2]
        for current_chain in chains.split(","):
            if (
                pdb in asym_chain_ids.keys()
                and current_chain in asym_chain_ids[pdb].keys()
            ):
                correct_chain = asym_chain_ids[pdb][current_chain]
                print(f"chain correction happening: {correct_chain} {current_chain}")
            else:
                correct_chain = current_chain
            if pdb in obsolete_pdb_ids:
                pdb = obsolete_pdb_ids[old_pdb]
                print("pdb correction happening")
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
                    old_pdb,
                    ligand,
                    correct_chain,
                    current_chain,
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
final_data.to_csv(config["az_broad_cluster_csv"])
cluster_data = pd.DataFrame(
    cluster_data,
    columns=[
        "Cluster",
        "PDB",
        "Ligand",
        "correct chain",
        "Author chain",
        "Uniprotid",
        "PfamID",
        "Pfam Name",
    ],
)
cluster_data.to_csv(config["az_exact_cluster_csv"])
print("Number of skipped entries:", n_skipped)
