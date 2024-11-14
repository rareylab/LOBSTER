import pandas as pd
from pathlib import Path
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
import progressbar, json
from argparse import ArgumentParser
from statistics import mean, stdev


arg_parser = ArgumentParser(
    description="Create CSV files with statistics on pairs / molecules of the"
    " LOBSTER dataset"
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
representatives_path = final_dataset_dir / "representatives"
clusters_path = final_dataset_dir / "clusters"

if not final_dataset_dir.is_dir():
    print("Input directory 10_final_dataset_dir as specified in config does not exist")
    raise IOError(str(final_dataset_dir))

table_path = Path(config["6_filtered_buriedness_csv"])


if not table_path.is_file():
    print("Input file 6_filtered_buriedness_csv as specified in config does not exist")
    raise FileNotFoundError(str(table_path))

filtered_master_table = pd.read_csv(table_path, dtype=str)

molecule_stats_csv_path = Path(config["stats_molecules"])
pairs_stats_csv_path = Path(config["stats_pairs"])
ensemble_stats_csv_path = Path(config["stats_ensembles"])

# check for parent directory due to default configuration to store files in results/stats
if not molecule_stats_csv_path.parent.is_dir():
    molecule_stats_csv_path.parent.mkdir()

if not pairs_stats_csv_path.parent.is_dir():
    pairs_stats_csv_path.parent.mkdir()

# write headers
with molecule_stats_csv_path.open("w") as f:
    f.write(
        "ligand;pdb;ensemble;mw;hac;rot_bonds;hba;hbd;logP;usmiles;uniprot_id;edia_m;"
        "pdb_resolution\n"
    )

with pairs_stats_csv_path.open("w") as f:
    f.write(
        "query_file;query;pdb;template_file;template;template_pdb;ensemble;"
        "morgan_fp_tanimoto;gobbi_2D_pharmacophore_fp_tanimoto;hac_difference;"
        "shape_tversky_index;shape_tanimoto_distance;shape_protrude_distance;\n"
    )

# initialize columns of ensemble result csv
information = {
    "Query Ligand": [],
    "Query Protein": [],
    "Average Morgan Fingerprint Tanimoto": [],
    "Number of Molecules": [],
    "HAC Standard Deviation": [],
    "Average Backbone RMSD": [],
}

# single molecule stats
for file in progressbar.progressbar(representatives_path.iterdir()):
    ensemble = file.stem.replace("-ligands", "")
    split_file_name = ensemble.split("-")
    query_ligand = split_file_name[0]
    query_pdb = split_file_name[1]

    # stats for ensembles
    ensemble_mol_sizes = []
    bb_rmsds = []
    ref_fp = None
    mol_fps = []

    for mol in Chem.SDMolSupplier(str(file)):
        ligand = mol.GetProp("_Name")
        pdb = mol.GetProp("pdb_id")
        table_row = filtered_master_table.where(
            (filtered_master_table["pdb"] == pdb)
            & (filtered_master_table["name"] == ligand)
        ).dropna(subset=["pdb"])
        hac = float(table_row["heavyAtoms"].values[0])
        # write the table
        with molecule_stats_csv_path.open("a") as f:
            f.write(f"{ligand};")
            f.write(f"{pdb.upper()};")
            f.write(f"{ensemble};")
            f.write(f"{table_row['mw'].values[0]};")
            f.write(f"{hac};")
            f.write(f"{table_row['nrot_simple'].values[0]};")
            f.write(f"{table_row['lipinskiAcceptors'].values[0]};")
            f.write(f"{table_row['lipinskiDonors'].values[0]};")
            f.write(f"{table_row['logP'].values[0]};")
            f.write(f"{mol.GetProp('usmiles')};")
            f.write(f"{mol.GetProp('uniprot_id')};")
            f.write(f"{table_row['EDIAm'].values[0]};")
            f.write(f"{table_row['resolution'].values[0]}\n")

        ensemble_mol_sizes.append(hac)
        # for average bb rmsd per ensemble
        statistics_path = clusters_path / ensemble / f"{ensemble}-resultStatistic.csv"
        stats_df = pd.read_csv(statistics_path, delimiter=";", dtype=str)
        if ligand == query_ligand and pdb.lower() == query_pdb.lower():
            ref_fp = AllChem.GetMorganFingerprint(mol, 2)
        else:
            bb_rmsd = float(
                stats_df.where(
                    (stats_df["PDB code"] == pdb)
                    & (stats_df["Ligand PDB code"] == ligand)
                )
                .dropna(subset=["PDB code"])["Backbone RMSD"]
                .values[0]
            )
            bb_rmsds.append(bb_rmsd)
            mol_fps.append(AllChem.GetMorganFingerprint(mol, 2))

        mol_morgan_fp = AllChem.GetMorganFingerprint(mol, 2)
        mol_gobbi_pharm_fp = fp = Generate.Gen2DFingerprint(mol, Gobbi_Pharm2D.factory)
        # multi molecule stats
        for paired_mol in Chem.SDMolSupplier(str(file)):
            paired_ligand = paired_mol.GetProp("_Name")
            paired_pdb = paired_mol.GetProp("pdb_id")
            paired_table_row = filtered_master_table.where(
                (filtered_master_table["pdb"] == paired_pdb)
                & (filtered_master_table["name"] == paired_ligand)
            ).dropna(subset=["pdb"])
            paired_hac = float(paired_table_row["heavyAtoms"].values[0])

            if paired_ligand != ligand and paired_pdb != pdb:
                # Calculate FP Tanimoto
                paired_mol_morgan_fp = AllChem.GetMorganFingerprint(paired_mol, 2)
                paired_mol_gobbi_pharm_fp = Generate.Gen2DFingerprint(
                    paired_mol, Gobbi_Pharm2D.factory
                )

                # calculate shape distances
                shape_tversky = rdShapeHelpers.ShapeTverskyIndex(
                    mol, paired_mol, alpha=1, beta=0
                )
                shape_tani = rdShapeHelpers.ShapeTanimotoDist(mol, paired_mol)
                shape_protrude = rdShapeHelpers.ShapeProtrudeDist(mol, paired_mol)

                with pairs_stats_csv_path.open("a") as f:
                    f.write(f"all_ligands/{paired_ligand}-{paired_pdb.upper()}.sdf;")
                    f.write(f"{paired_ligand};")
                    f.write(f"{paired_pdb.lower()};")
                    f.write(f"all_ligands/{ligand}-{pdb.upper()}.sdf;")
                    f.write(f"{ligand};")
                    f.write(f"{pdb.lower()};")
                    f.write(f"{ensemble};")
                    tanimoto_sim = round(
                        DataStructs.TanimotoSimilarity(
                            mol_morgan_fp, paired_mol_morgan_fp
                        ),
                        3,
                    )
                    f.write(f"{tanimoto_sim};")
                    gobbi_sim = round(
                        DataStructs.TanimotoSimilarity(
                            mol_gobbi_pharm_fp, paired_mol_gobbi_pharm_fp
                        ),
                        3,
                    )
                    f.write(f"{gobbi_sim};")
                    f.write(f"{abs(hac - paired_hac)};")
                    f.write(f"{round(shape_tversky, 3)};")
                    f.write(f"{round(shape_tani, 3)};")
                    f.write(f"{round(shape_protrude, 3)};\n")
    avg_similarity = 0
    if ref_fp != None:
        for non_ref_fp in mol_fps:
            avg_similarity += DataStructs.TanimotoSimilarity(ref_fp, non_ref_fp)
        avg_similarity /= len(mol_fps)
    else:
        print("No reference found!", ensemble)
        print("This should never happen.")

    information["Query Ligand"].append(query_ligand)
    information["Query Protein"].append(query_pdb)
    information["Average Morgan Fingerprint Tanimoto"].append(avg_similarity)
    information["Number of Molecules"].append(len(ensemble_mol_sizes))
    information["HAC Standard Deviation"].append(
        "{0}".format(round(stdev(ensemble_mol_sizes), 3))
    )
    information["Average Backbone RMSD"].append("{0}".format(round(mean(bb_rmsds), 3)))

df = pd.DataFrame.from_dict(information)
df.to_csv(ensemble_stats_csv_path, sep=";", index=False)
