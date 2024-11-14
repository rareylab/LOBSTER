import pandas as pd
from rdkit import Chem
import progressbar
from pathlib import Path
import argparse
import json


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="Initial filtering steps")
    arg_parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to config in JSON format including all file paths.",
    )
    args = arg_parser.parse_args()

    config = json.load(open(args.config, "r"))

    if not Path(config["initial_csv"]).is_file():
        print("Input file initial_csv specified in config does not exist")
        raise FileNotFoundError(config["initial_csv"])

    # create result dir if necessary due to unaltered config paths
    if not Path(config["1_main_filter_csv"]).parent.is_dir():
        Path(config["1_main_filter_csv"]).parent.mkdir()

    if not Path(config["log_file"]).parent.is_dir():
        Path(config["log_file"]).parent.mkdir()

    # variable initialization
    unprocessable = 0
    skipped_by_flachsi = 0
    resolution_removed = 0
    ediaM_removed = 0
    mw_removed = 0
    hac_removed = 0
    rot_bond_removed = 0
    atom_type_removed = 0
    occurring_unallowed = set()
    drop_idx = []
    smiles_removed = set()
    allowed_elements = {"C", "N", "O", "S", "P", "Br", "Cl", "F", "I", "B"}

    zipped = set()
    total_lig_pdbs = 0

    df = pd.read_csv(config["initial_csv"], delimiter=",", dtype=str)
    df.dropna(subset=["name"], inplace=True)  # remove rows without ligands

    for idx, row in progressbar.progressbar(df.iterrows()):
        ligand = row["name"].split("_")[0]

        # Filter for skipreason
        if row["skip_reason"] != "NotSkipped" or row["ligand_structure"] != "both":
            skipped_by_flachsi += 1
            drop_idx.append(idx)
            continue

        # Filter for PDB Resolution
        if float(row["resolution"]) > 2.5:
            resolution_removed += 1
            drop_idx.append(idx)
            continue

        # Filter for EdiaM
        # not >= 0.8, to filter out examples, where EDIAm is nan
        if not float(row["EDIAm"]) >= 0.8:
            ediaM_removed += 1
            drop_idx.append(idx)
            continue

        # Filter for molecular weight >975
        if float(row["mw"]) > 975:
            mw_removed += 1
            drop_idx.append(idx)
            continue

        # Filter for rotatable bonds > 10
        if float(row["nrot_simple"]) > 10:
            rot_bond_removed += 1
            drop_idx.append(idx)
            continue

        try:
            molecule = Chem.MolFromSmiles(row["usmiles"])

            # Filter for atom types, only include C,N,O,S,P,Br,Cl,F,I,B
            unallowed_elems = [
                atom.GetSymbol()
                for atom in molecule.GetAtoms()
                if atom.GetSymbol() not in allowed_elements
            ]
            if len(unallowed_elems) > 0:
                occurring_unallowed.update(unallowed_elems)
                atom_type_removed += 1
                drop_idx.append(idx)
                continue

        except:
            unprocessable += 1
            drop_idx.append(idx)
            continue

        # Filter for too small molecules
        if float(row["heavyAtoms"]) < 10:
            hac_removed += 1
            drop_idx.append(idx)
            continue

    df = df.drop(drop_idx)
    df.to_csv(config["1_main_filter_csv"], mode="w")
    log = open(config["log_file"], "w")
    log.write(f"Skipped by construction: {skipped_by_flachsi}\n")
    log.write(f"Molecules removed with resolution > 2.5 A: {resolution_removed}\n")
    log.write(f"Molecules removed with EDIAm < 0.8: {ediaM_removed}\n")
    log.write(f"Molecules removed with MW > 975: {mw_removed}\n")
    log.write(f"Molecules removed by nof rot bonds > 10: {rot_bond_removed}\n")
    log.write(f"Unprocessable molecules removed: {unprocessable}\n")
    log.write(f"Molecules removed by atom type: {atom_type_removed}\n")
    log.write(f"Unallowed elements occurring: {occurring_unallowed}\n")
    log.write(f"Molecules removed HAC < 10: {hac_removed}\n")
