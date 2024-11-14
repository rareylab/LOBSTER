import pandas as pd
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys
import statistics
from collections import Counter
from sklearn import metrics
from pathlib import Path
from progressbar import progressbar
from argparse import ArgumentParser
import json


def calc_median_maccs(mol_list):
    maccs_list = []
    maccs_keys = []
    for mol in mol_list:
        maccs_keys.append(MACCSkeys.GenMACCSKeys(mol))

    for i, key1 in enumerate(maccs_keys):
        for j, key2 in enumerate(maccs_keys):
            if i == j:
                continue
            maccs_list.append(DataStructs.FingerprintSimilarity(key1, key2))
    return statistics.median(maccs_list)


# ----------------------------------------------------------------------------------------


def calc_scaffold_diversity(mol):
    all_scaffolds = []
    n_compounds = len(mol)
    for mol in mol:
        scaffold = Chem.MolToSmiles(Chem.MurckoDecompose(mol))
        all_scaffolds.append(scaffold)
    scaffold_counter = Counter(all_scaffolds)
    number_scaffolds = len(scaffold_counter)
    values = list(scaffold_counter.values())
    values.sort(reverse=True)
    first_rate = [0.0]
    second_rate = [0.0]
    total = 0
    i = 0
    for value in values:
        total += value
        i += 1
        first_rate.append(i / number_scaffolds)
        second_rate.append(total / n_compounds)
    auc = round(metrics.auc(x=first_rate, y=second_rate), 3)
    return auc


# ----------------------------------------------------------------------------------------


def diversity_subset(subset_path, result_path, dataset_dir):
    total_data = []
    for subset in progressbar(subset_path.iterdir()):
        if subset.is_file() and subset.name != "all_pairs.csv":
            df = pd.read_csv(subset, delimiter=";")
            threshold = subset.stem
            mols = []
            smiles = set()
            for idx, row in df.iterrows():
                file_1 = dataset_dir / row["query_file"]
                file_2 = dataset_dir / row["template_file"]
                mol_1 = [x for x in Chem.SDMolSupplier(str(file_1))][0]
                mol_2 = [x for x in Chem.SDMolSupplier(str(file_2))][0]
                if not mol_1.GetProp("usmiles") in smiles:
                    smiles.add(mol_1.GetProp("usmiles"))
                    mols.append(mol_1)
                if not mol_2.GetProp("usmiles") in smiles:
                    smiles.add(mol_2.GetProp("usmiles"))
                    mols.append(mol_2)
            scaffold_stuff = calc_scaffold_diversity(mols)
            median_maccs = calc_median_maccs(mols)
            local_data = [
                row["query_file"],
                row["template_file"],
                threshold,
                len(df.index),
                median_maccs,
                scaffold_stuff,
            ]
            total_data.append(local_data)
    total_data.sort(key=lambda x: x[2])
    total_data = pd.DataFrame(
        total_data,
        columns=[
            "Ligand",
            "Template",
            "Set",
            "Size",
            "Median MACCS FP-based Tanimoto Similarity",
            "Scaffold AUC",
        ],
    )
    total_data.to_csv(result_path)


# ----------------------------------------------------------------------------------------


def diversity_ensembles(ensembles_path, results_path):
    total_data = []
    all_mols = []
    for ensemble in progressbar(ensembles_path.iterdir()):
        ensemble_name = ensemble.stem.replace("-ligands", "")
        mols = [x for x in Chem.SDMolSupplier(str(ensemble))]
        all_mols.extend(mols)
        scaffold_stuff = calc_scaffold_diversity(mols)
        median_maccs = calc_median_maccs(mols)
        local_data = [
            ensemble_name,
            len(mols),
            median_maccs,
            scaffold_stuff,
        ]
        total_data.append(local_data)
    total_data.sort(key=lambda x: x[2])
    total_data = pd.DataFrame(
        total_data,
        columns=[
            "Ensemble",
            "Size",
            "Median MACCS FP-based Tanimoto Similarity",
            "Scaffold AUC",
        ],
    )
    total_data.to_csv(results_path)
    print("Overall Ensemble Scaffold Diversity:", calc_scaffold_diversity(all_mols))
    print("Overall Ensemble Median MACCS Similarity:", calc_median_maccs(all_mols))


# ----------------------------------------------------------------------------------------


def diversity_for_az(path_az, results_path):
    total_data = []
    all_mols = []
    none_molecules = 0
    for ensemble in progressbar(path_az.iterdir()):
        ensemble_name = ensemble.stem.replace("-ligands", "")
        mols = []
        for mol in Chem.SDMolSupplier(str(ensemble)):
            if mol != None:
                mols.append(mol)
            else:
                none_molecules += 1
        all_mols.extend(mols)
        scaffold_stuff = calc_scaffold_diversity(mols)
        median_maccs = calc_median_maccs(mols)
        local_data = [
            ensemble_name,
            len(mols),
            median_maccs,
            scaffold_stuff,
        ]
        total_data.append(local_data)
    print("Molecules not processed by RDKit:", none_molecules)
    total_data.sort(key=lambda x: x[2])
    total_data = pd.DataFrame(
        total_data,
        columns=[
            "Ensemble",
            "Size",
            "Median MACCS FP-based Tanimoto Similarity",
            "Scaffold AUC",
        ],
    )
    total_data.to_csv(results_path)
    print("Overall AZ Scaffold Diversity:", calc_scaffold_diversity(all_mols))
    print("Overall AZ Median MACCS Similarity:", calc_median_maccs(all_mols))


# ----------------------------------------------------------------------------------------


if __name__ == "__main__":

    arg_parser = ArgumentParser(
        description="Create statistics on diversity of the molecules"
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

    final_dataset_dir = Path(config["10_final_dataset_dir"])
    if not final_dataset_dir.is_dir():
        print(
            "Input directory 10_final_dataset_dir as specified in config does not exist"
        )
        raise IOError(str(final_dataset_dir))

    # diversity for each subset
    subsets = final_dataset_dir / "subsets"
    result_subsets = Path(config["diversity_per_subset"])
    diversity_subset(subsets, result_subsets, final_dataset_dir)
    print(f"Finished processing subsets, diversity csv: {result_subsets}")

    # diversity for each ensemble
    ensembles_path = final_dataset_dir / "representatives"
    result_ensembles = Path(config["diversity_per_ensemble"])
    diversity_ensembles(ensembles_path, result_ensembles)
    print(f"Finished processing ensembles, diversity csv: {result_ensembles}")

    # diversity for Astra Zeneca Overlays dataset
    az_path = Path(config["path_to_az_data"])
    if not az_path.is_dir():
        print("Input directory path_to_az_data as specified in config does not exist")
        raise IOError(str(az_path))

    result_az = Path(config["diversity_az"])
    diversity_for_az(az_path, result_az)
    print(f"Finished processing AZ dataset, diversity csv: {result_az}")
